#!/usr/bin/env python3
"""Tracker multi-fit closure and joint beta/R upper-bound study for C12/C13.

This script never changes the production estimator.  It joins tracker-internal
fit variants to the already frozen C12/C13 controlled train/validation split,
then asks two questions:

1. Can tracker-only redundancy predict the event-level rigidity error on files
   not used for fitting?
2. If an oracle tracker correction exists, how much additional isotope
   separation is available when it is combined with the frozen beta results?

MC truth is used only for the explicit oracle regression and final metrics.
The feature transform itself uses tracker fit disagreements, fit errors and
track geometry, not isotope labels, beta, dE/dx, or MC momentum.
"""

import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import uproot

CHARGE = 6.0
MASSES = np.array([11.17793, 12.11221], dtype=float)
ISOTOPE_NAMES = np.array(["C12", "C13"])
SEED = 131512

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "c_isotope_proportion_closure"
CACHE_PATH = BASE / "controlled_cache.npz"
MANIFEST_PATH = BASE / "manifest.tsv"
PREDICTIONS_PATH = BASE / "predictions.npz"

TRACKER_BRANCHES = [
    "run",
    "event",
    "tk_rigidity1",
    "tk_rigidityk",
    "tk_erigidity1",
    "tk_rigidityi",
    "tk_chisi",
    "tk_bcor",
    "tk_cdif",
    "tk_fzm",
    "tk_qrms",
    "tk_qhit",
]


def log(message):
    print(message, flush=True)


def w68(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan")
    lo, hi = np.percentile(values, [16.0, 84.0])
    return float((hi - lo) / 2.0)


def rms(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.sqrt(np.mean(values * values))) if len(values) else float("nan")


def auc(score, label):
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=np.int8)
    finite = np.isfinite(score) & ((label == 0) | (label == 1))
    score = score[finite]
    label = label[finite]
    n1 = int(np.sum(label == 1))
    n0 = int(np.sum(label == 0))
    if n0 == 0 or n1 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    sorted_label = label[order]
    starts = np.r_[0, np.flatnonzero(sorted_score[1:] != sorted_score[:-1]) + 1]
    ends = np.r_[starts[1:], len(sorted_score)]
    ranks = np.repeat((starts + 1.0 + ends) / 2.0, ends - starts)
    return float(
        (np.sum(ranks[sorted_label == 1]) - n1 * (n1 + 1.0) / 2.0)
        / (float(n0) * n1)
    )


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as source:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_controlled_cache(path):
    with np.load(str(path), allow_pickle=False) as source:
        train, valid = {}, {}
        for name in source.files:
            prefix, key = name.split("__", 1)
            (train if prefix == "train" else valid)[key] = source[name]
    return train, valid


def load_manifest(path):
    groups = {(role, label): [] for role in (0, 1) for label in (0, 1)}
    label_map = {"C12": 0, "C13": 1}
    with open(str(path), newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            groups[(int(row["role"]), label_map[row["isotope"]])].append(
                row["path"]
            )
    return groups


def allocate_tracker_arrays(n):
    return {
        "r1": np.full((n, 3, 3, 7), np.nan, dtype=np.float64),
        "rk": np.full((n, 4, 7), np.nan, dtype=np.float64),
        "er1": np.full((n, 3, 7), np.nan, dtype=np.float64),
        "ri": np.full((n, 9), np.nan, dtype=np.float64),
        "chisi": np.full((n, 9), np.nan, dtype=np.float64),
        "bcor": np.full(n, np.nan, dtype=np.float64),
        "cdif": np.full((n, 2, 3), np.nan, dtype=np.float64),
        "fzm": np.full((n, 2), np.nan, dtype=np.float64),
        "qrms": np.full((n, 2), np.nan, dtype=np.float64),
        "qhit": np.full((n, 2), np.nan, dtype=np.float64),
    }


def extract_tracker_for_split(data, files_by_group, split_name):
    n = len(data["run"])
    output = allocate_tracker_arrays(n)
    found = np.zeros(n, dtype=bool)
    audit = []

    for label in (0, 1):
        rows_label = np.flatnonzero(data["label"] == label)
        role_values = np.unique(data["role"][rows_label])
        if len(role_values) != 1:
            raise RuntimeError("Mixed roles inside %s/%s" % (split_name, label))
        role = int(role_values[0])
        files = files_by_group[(role, label)]
        for source_index in sorted(int(v) for v in np.unique(data["source"][rows_label])):
            target = np.flatnonzero(
                (data["label"] == label) & (data["source"] == source_index)
            )
            if source_index >= len(files):
                raise RuntimeError("Source index outside manifest")
            path = files[source_index]
            started = time.time()
            with uproot.open(path) as root_file:
                tree = root_file["amstreea"]
                arrays = tree.arrays(TRACKER_BRANCHES, library="np")

            target_lookup = {
                (int(data["run"][row]), int(data["event"][row])): int(row)
                for row in target
            }
            matched_entries = []
            matched_rows = []
            for entry, (run, event) in enumerate(
                zip(arrays["run"], arrays["event"])
            ):
                row = target_lookup.get((int(run), int(event)))
                if row is not None:
                    matched_entries.append(entry)
                    matched_rows.append(row)
            if len(matched_rows) != len(target):
                missing = len(target) - len(matched_rows)
                raise RuntimeError(
                    "%s %s source=%d missing %d/%d selected events"
                    % (split_name, ISOTOPE_NAMES[label], source_index, missing, len(target))
                )
            entries = np.asarray(matched_entries, dtype=np.int64)
            rows = np.asarray(matched_rows, dtype=np.int64)
            output["r1"][rows] = arrays["tk_rigidity1"][entries]
            output["rk"][rows] = arrays["tk_rigidityk"][entries]
            output["er1"][rows] = arrays["tk_erigidity1"][entries]
            output["ri"][rows] = arrays["tk_rigidityi"][entries]
            output["chisi"][rows] = arrays["tk_chisi"][entries]
            output["bcor"][rows] = arrays["tk_bcor"][entries]
            output["cdif"][rows] = arrays["tk_cdif"][entries]
            output["fzm"][rows] = arrays["tk_fzm"][entries]
            output["qrms"][rows] = arrays["tk_qrms"][entries]
            output["qhit"][rows] = arrays["tk_qhit"][entries]
            found[rows] = True
            audit.append(
                {
                    "split": split_name,
                    "isotope": ISOTOPE_NAMES[label],
                    "source": source_index,
                    "path": path,
                    "tree_entries": int(len(arrays["run"])),
                    "matched": int(len(rows)),
                    "seconds": time.time() - started,
                }
            )
            log(
                "extract %s %s file=%d/%d matched=%d elapsed=%.1fs"
                % (
                    split_name,
                    ISOTOPE_NAMES[label],
                    source_index + 1,
                    len(files),
                    len(rows),
                    time.time() - started,
                )
            )

    if not np.all(found):
        raise RuntimeError("Unfilled tracker rows in " + split_name)
    baseline = output["rk"][:, 1, 1]
    relative = np.abs(baseline - data["R"]) / np.maximum(np.abs(data["R"]), 1e-12)
    maximum = float(np.nanmax(relative))
    if maximum > 2e-6:
        raise RuntimeError("Baseline rigidity mismatch %.6g" % maximum)
    return output, audit


def candidate_matrix(tracker):
    columns = []
    names = []
    for algorithm in range(3):
        for alignment in range(3):
            for span in range(7):
                columns.append(tracker["r1"][:, algorithm, alignment, span])
                names.append("r1_a%d_v%d_s%d" % (algorithm, alignment, span))
    for plane in range(4):
        for span in range(7):
            columns.append(tracker["rk"][:, plane, span])
            names.append("rk_z%d_s%d" % (plane, span))
    for index in range(9):
        columns.append(tracker["ri"][:, index])
        names.append("ri_%d" % index)
    return np.column_stack(columns), names


def raw_quality_features(data, tracker, candidates, candidate_names):
    baseline = data["R"].astype(np.float64)
    values = []
    names = []

    with np.errstate(divide="ignore", invalid="ignore"):
        candidate_ratio = np.log(np.abs(candidates / baseline[:, None]))
    candidate_ratio[~np.isfinite(candidate_ratio)] = np.nan
    candidate_ratio[np.abs(candidate_ratio) > 1.5] = np.nan
    values.append(candidate_ratio)
    names.extend("logratio_" + name for name in candidate_names)

    with np.errstate(divide="ignore", invalid="ignore"):
        relative_error = np.log(np.abs(tracker["er1"] * baseline[:, None, None]))
    relative_error[~np.isfinite(relative_error)] = np.nan
    relative_error[np.abs(relative_error) > 12.0] = np.nan
    values.append(relative_error.reshape(len(baseline), -1))
    names.extend(
        "log_relerr_a%d_s%d" % (algorithm, span)
        for algorithm in range(3)
        for span in range(7)
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        chisi = np.log(np.abs(tracker["chisi"]))
    chisi[~np.isfinite(chisi)] = np.nan
    chisi[np.abs(chisi) > 12.0] = np.nan
    values.append(chisi)
    names.extend("log_chisi_%d" % index for index in range(9))

    values.extend(
        [
            tracker["bcor"][:, None],
            tracker["cdif"].reshape(len(baseline), -1),
            tracker["fzm"],
            tracker["qrms"],
            tracker["qhit"].astype(np.float64),
            data["continuous"][:, :14].astype(np.float64),
            data["relR"][:, None].astype(np.float64),
        ]
    )
    names.extend(["bcor"])
    names.extend("cdif_%d_%d" % (side, axis) for side in range(2) for axis in range(3))
    names.extend("fzm_%d" % index for index in range(2))
    names.extend("qrms_%d" % index for index in range(2))
    names.extend("qhit_%d" % index for index in range(2))
    names.extend("geometry_%d" % index for index in range(14))
    names.append("baseline_relR")
    return np.column_stack(values), names


def fit_transform_state(raw):
    raw = np.asarray(raw, dtype=np.float64)
    finite_fraction = np.mean(np.isfinite(raw), axis=0)
    keep = finite_fraction >= 0.60
    raw_kept = raw[:, keep]
    medians = np.nanmedian(raw_kept, axis=0)
    filled = np.where(np.isfinite(raw_kept), raw_kept, medians)
    low, high = np.percentile(filled, [1.0, 99.0], axis=0)
    clipped = np.minimum(np.maximum(filled, low), high)
    center = np.median(clipped, axis=0)
    q25, q75 = np.percentile(clipped, [25.0, 75.0], axis=0)
    scale = 0.741301109 * (q75 - q25)
    keep_scale = np.isfinite(scale) & (scale > 1e-7)
    kept_indices = np.flatnonzero(keep)[keep_scale]
    return {
        "indices": kept_indices,
        "median": medians[keep_scale],
        "low": low[keep_scale],
        "high": high[keep_scale],
        "center": center[keep_scale],
        "scale": scale[keep_scale],
    }


def apply_transform_state(raw, state):
    selected = np.asarray(raw, dtype=np.float64)[:, state["indices"]]
    filled = np.where(np.isfinite(selected), selected, state["median"])
    clipped = np.minimum(np.maximum(filled, state["low"]), state["high"])
    standardized = np.clip(
        (clipped - state["center"]) / state["scale"], -6.0, 6.0
    )
    return np.column_stack([np.ones(len(standardized)), standardized])


def huber_ridge(x, y, ridge, iterations=20):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    penalty = np.full(x.shape[1], float(ridge), dtype=np.float64)
    penalty[0] = 1e-8
    theta = np.linalg.solve(x.T.dot(x) + np.diag(penalty), x.T.dot(y))
    for _ in range(iterations):
        residual = y - x.dot(theta)
        delta = max(0.015, 1.8 * w68(residual))
        weight = np.minimum(1.0, delta / np.maximum(np.abs(residual), 1e-12))
        weighted = x * weight[:, None]
        updated = np.linalg.solve(
            x.T.dot(weighted) + np.diag(penalty), x.T.dot(weight * y)
        )
        if np.max(np.abs(updated - theta)) < 1e-9:
            theta = updated
            break
        theta = updated
    return theta


def group_folds(label, source, n_folds=5):
    groups = sorted(set(zip(label.tolist(), source.tolist())))
    rng = np.random.RandomState(SEED)
    permutation = rng.permutation(len(groups))
    assignment = {}
    for order, group_index in enumerate(permutation):
        assignment[groups[group_index]] = order % n_folds
    return np.asarray(
        [assignment[(int(a), int(b))] for a, b in zip(label, source)],
        dtype=np.int8,
    )


def select_ridge(x, y, label, source):
    candidates = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]
    folds = group_folds(label, source)
    records = []
    for ridge in candidates:
        widths = []
        medians = []
        for fold in range(5):
            fit = folds != fold
            hold = folds == fold
            theta = huber_ridge(x[fit], y[fit], ridge)
            residual = y[hold] - x[hold].dot(theta)
            widths.append(w68(residual))
            medians.append(abs(float(np.median(residual))))
        record = {
            "ridge": ridge,
            "mean_w68": float(np.mean(widths)),
            "mean_abs_median": float(np.mean(medians)),
            "score": float(np.mean(widths) + 0.25 * np.mean(medians)),
        }
        records.append(record)
        log("ridge CV %s" % record)
    best = min(records, key=lambda item: item["score"])
    return float(best["ridge"]), records


def fit_logistic(x, label):
    label = label.astype(np.float64)

    def objective(theta):
        score = x.dot(theta)
        value = float(np.mean(np.logaddexp(0.0, score) - label * score))
        probability = 1.0 / (1.0 + np.exp(-np.clip(score, -40.0, 40.0)))
        gradient = x.T.dot(probability - label) / len(label)
        value += 0.5e-3 * float(np.sum(theta[1:] * theta[1:]))
        gradient[1:] += 1e-3 * theta[1:]
        return value, gradient

    result = minimize(
        objective,
        np.zeros(x.shape[1], dtype=float),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 800, "ftol": 1e-11, "gtol": 1e-6},
    )
    if not result.success:
        log("warning: label diagnostic optimizer: %s" % result.message)
    return result.x


def rigidity_metrics(name, corrected_r, data):
    truth = data["truth_valid"] & np.isfinite(data["R_true"]) & (data["R_true"] > 0)
    rows = []
    for group_name, group_mask in [
        ("all", np.ones(len(corrected_r), dtype=bool)),
        ("C12", data["label"] == 0),
        ("C13", data["label"] == 1),
    ]:
        select = truth & group_mask & np.isfinite(corrected_r) & (corrected_r > 0)
        error = np.log(corrected_r[select] / data["R_true"][select])
        rows.append(
            {
                "method": name,
                "group": group_name,
                "n": int(np.sum(select)),
                "median_dlogR": float(np.median(error)),
                "w68_dlogR": w68(error),
                "rms_dlogR": rms(error),
                "improvement_vs_raw_percent": float("nan"),
            }
        )
    return rows


def mass_metrics(name, u, rigidity, data):
    truth = (
        data["truth_valid"]
        & np.isfinite(u)
        & (u > 0)
        & np.isfinite(rigidity)
        & (rigidity > 0)
    )
    mass2 = (CHARGE * rigidity) ** 2 * (u * u - 1.0)
    true_mass2 = MASSES[data["label"]] ** 2
    light = truth & (data["label"] == 0)
    heavy = truth & (data["label"] == 1)
    denominator = math.sqrt(w68(mass2[light]) ** 2 + w68(mass2[heavy]) ** 2)
    return {
        "method": name,
        "n": int(np.sum(truth)),
        "mass2_auc": auc(mass2[truth], data["label"][truth]),
        "mass2_separation": (
            (float(np.median(mass2[heavy])) - float(np.median(mass2[light])))
            / denominator
            if denominator > 0
            else float("nan")
        ),
        "median_mass2_error": float(np.median(mass2[truth] - true_mass2[truth])),
        "w68_mass2_error": w68(mass2[truth] - true_mass2[truth]),
        "rms_mass2_error": rms(mass2[truth] - true_mass2[truth]),
    }


def write_tsv(path, rows):
    if not rows:
        return
    with open(str(path), "w", newline="") as output:
        writer = csv.DictWriter(
            output, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def single_candidate_metrics(candidates_train, candidates_valid, names, train, valid):
    truth_train = train["truth_valid"] & (train["R_true"] > 0)
    truth_valid = valid["truth_valid"] & (valid["R_true"] > 0)
    rows = []
    for column, name in enumerate(names):
        train_value = candidates_train[:, column]
        valid_value = candidates_valid[:, column]
        train_select = truth_train & np.isfinite(train_value) & (train_value > 0)
        valid_select = truth_valid & np.isfinite(valid_value) & (valid_value > 0)
        if np.mean(train_select[truth_train]) < 0.80 or np.mean(valid_select[truth_valid]) < 0.80:
            continue
        train_error = np.log(train_value[train_select] / train["R_true"][train_select])
        calibration = float(np.median(train_error))
        valid_error = (
            np.log(valid_value[valid_select] / valid["R_true"][valid_select])
            - calibration
        )
        rows.append(
            {
                "candidate": name,
                "train_coverage": float(np.mean(train_select[truth_train])),
                "valid_coverage": float(np.mean(valid_select[truth_valid])),
                "train_median_dlogR": calibration,
                "train_w68_dlogR": w68(train_error),
                "valid_median_after_train_calibration": float(np.median(valid_error)),
                "valid_w68_after_train_calibration": w68(valid_error),
                "valid_rms_after_train_calibration": rms(valid_error),
            }
        )
    rows.sort(key=lambda item: item["train_w68_dlogR"])
    return rows


def main():
    started = time.time()
    train, valid = load_controlled_cache(CACHE_PATH)
    files_by_group = load_manifest(MANIFEST_PATH)

    tracker_train, audit_train = extract_tracker_for_split(train, files_by_group, "train")
    tracker_valid, audit_valid = extract_tracker_for_split(valid, files_by_group, "validation")
    write_tsv(HERE / "tracker_extraction_audit.tsv", audit_train + audit_valid)

    np.savez_compressed(
        str(HERE / "tracker_multifit_cache.npz"),
        **{"train__" + key: value for key, value in tracker_train.items()},
        **{"valid__" + key: value for key, value in tracker_valid.items()},
    )

    candidates_train, candidate_names = candidate_matrix(tracker_train)
    candidates_valid, candidate_names_valid = candidate_matrix(tracker_valid)
    if candidate_names != candidate_names_valid:
        raise RuntimeError("Candidate schema mismatch")
    single_rows = single_candidate_metrics(
        candidates_train, candidates_valid, candidate_names, train, valid
    )
    write_tsv(HERE / "single_fit_metrics.tsv", single_rows)
    log("top single fits: %s" % single_rows[:5])

    raw_train, raw_names = raw_quality_features(
        train, tracker_train, candidates_train, candidate_names
    )
    raw_valid, raw_names_valid = raw_quality_features(
        valid, tracker_valid, candidates_valid, candidate_names
    )
    if raw_names != raw_names_valid:
        raise RuntimeError("Raw feature schema mismatch")
    state = fit_transform_state(raw_train)
    x_train_all = apply_transform_state(raw_train, state)
    x_valid_all = apply_transform_state(raw_valid, state)
    selected_names = [raw_names[int(index)] for index in state["indices"]]

    truth_train = train["truth_valid"] & np.isfinite(train["R_true"]) & (train["R_true"] > 0)
    truth_valid = valid["truth_valid"] & np.isfinite(valid["R_true"]) & (valid["R_true"] > 0)
    target = np.log(train["R"][truth_train] / train["R_true"][truth_train])
    x_train = x_train_all[truth_train]
    best_ridge, cv_records = select_ridge(
        x_train,
        target,
        train["label"][truth_train],
        train["source"][truth_train],
    )
    theta = huber_ridge(x_train, target, best_ridge)
    predicted_error_valid = x_valid_all.dot(theta)
    corrected_r = valid["R"] * np.exp(-predicted_error_valid)

    rigidity_rows = []
    rigidity_rows.extend(rigidity_metrics("raw", valid["R"], valid))
    rigidity_rows.extend(rigidity_metrics("tracker_multifit_oracle", corrected_r, valid))
    raw_width = {
        row["group"]: row["w68_dlogR"]
        for row in rigidity_rows
        if row["method"] == "raw"
    }
    for row in rigidity_rows:
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_dlogR"] / raw_width[row["group"]]
        )
    write_tsv(HERE / "rigidity_metrics.tsv", rigidity_rows)

    classifier = fit_logistic(x_train_all, train["label"])
    label_auc = auc(x_valid_all.dot(classifier), valid["label"])

    with np.load(str(PREDICTIONS_PATH), allow_pickle=False) as predictions:
        for key in ("run", "event", "label", "bin"):
            if not np.array_equal(predictions[key], valid[key]):
                raise RuntimeError("Prediction ordering mismatch: " + key)
        u_methods = {
            "raw_beta": predictions["u_raw"],
            "weak_no_edep_beta": predictions["u_weak_no_edep"],
            "oracle_full_beta": predictions["u_oracle_full"],
        }

    mass_rows = []
    for beta_name, u in u_methods.items():
        mass_rows.append(mass_metrics(beta_name + "__raw_R", u, valid["R"], valid))
        mass_rows.append(
            mass_metrics(beta_name + "__tracker_multifit_oracle_R", u, corrected_r, valid)
        )
    write_tsv(HERE / "joint_mass_metrics.tsv", mass_rows)

    model = {
        "truth_policy": (
            "Tracker feature preprocessing is label/truth free. MC truth is used only "
            "for the explicitly named tracker_multifit_oracle regression and final metrics."
        ),
        "feature_names": selected_names,
        "transform_state": state,
        "theta": theta,
        "ridge": best_ridge,
        "cv_records": cv_records,
        "candidate_names": candidate_names,
    }
    (HERE / "tracker_oracle_model.json").write_text(
        json.dumps(jsonable(model), indent=2, sort_keys=True)
    )

    summary = {
        "train_events": int(len(train["R"])),
        "validation_events": int(len(valid["R"])),
        "train_truth_events": int(np.sum(truth_train)),
        "validation_truth_events": int(np.sum(truth_valid)),
        "tracker_feature_dimensions": int(x_train_all.shape[1]),
        "tracker_feature_label_auc_posthoc": label_auc,
        "best_ridge": best_ridge,
        "rigidity_metrics": {
            row["method"] + ":" + row["group"]: row for row in rigidity_rows
        },
        "joint_mass_metrics": {row["method"]: row for row in mass_rows},
        "top_single_fits": single_rows[:20],
        "elapsed_seconds": time.time() - started,
        "script_sha256": sha256_file(Path(__file__)),
        "input_sha256": {
            "controlled_cache": sha256_file(CACHE_PATH),
            "manifest": sha256_file(MANIFEST_PATH),
            "predictions": sha256_file(PREDICTIONS_PATH),
        },
        "limitation": (
            "C12 is B1308 and C13 is B1315. The post-hoc label AUC measures "
            "production-version separability of tracker features; no flight transfer "
            "claim is permitted without a same-production isotope pair."
        ),
    }
    (HERE / "tracker_probe_results.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )
    log("DONE elapsed=%.1fs" % (time.time() - started))


if __name__ == "__main__":
    main()
