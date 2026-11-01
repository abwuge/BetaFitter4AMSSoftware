#!/usr/bin/env python3
"""Train and evaluate an unconstrained supervised C12/C13 mass toy.

The toy deliberately asks for the strongest event-level isotope inference that
can be obtained from current reconstructed detector information.  It is an MC
truth upper-bound study, not a flight-trainable calibration.  Train/validation/
test separation is by whole source file, isotope priors are balanced separately
in every reconstructed-rigidity bin, and all metrics are computed on untouched
source files.
"""
from __future__ import print_function

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

import mass_toy_common as mt

HERE = Path(__file__).resolve().parent
JOINT = HERE.parent / "c_joint_beta_rigidity_closure"
EXPANDED = JOINT / "expanded_validation"
MANIFEST = EXPANDED / "manifest.tsv"
EXPANDED_PREDICTIONS = EXPANDED / "predictions.npz"
FEATURE_CACHE = HERE / "feature_cache.npz"
FEATURE_NAMES = HERE / "feature_names.json"

METHODS = [
    "raw_linear",
    "raw_hybrid",
    "mass_proxy",
    "mass_energy",
    "all_physics",
    "all_detector",
    "nuisance_control",
    "shuffle_control",
]
DISPLAY_NAMES = {
    "raw_linear": "Linear TOF + raw R",
    "raw_hybrid": "Hybrid TOF + raw R",
    "mass_proxy": "All beta/R mass proxies",
    "mass_energy": "+ TOF/tracker/TRD/ECAL energy loss",
    "all_physics": "+ RICH beta",
    "all_detector": "+ geometry and fit-quality variables",
    "nuisance_control": "Geometry/quality only",
    "shuffle_control": "Shuffled-label control",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-features", action="store_true")
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--bootstrap", type=int, default=500)
    return parser.parse_args()


def load_truth_valid_data():
    with np.load(str(EXPANDED_PREDICTIONS), allow_pickle=False) as source:
        keep = (
            source["truth_valid"]
            & np.isfinite(source["u_true"])
            & np.isfinite(source["R_true"])
            & (source["R_true"] > 0.0)
        )
        keys = ["run", "event", "source", "label", "bin", "u", "R", "u_true", "R_true"]
        return {key: source[key][keep] for key in keys}


def assign_roles(data, role_by_source):
    role = np.full(len(data["u"]), "", dtype="U5")
    for label in (0, 1):
        for source, split_role in role_by_source[label].items():
            rows = (data["label"] == label) & (data["source"] == source)
            role[rows] = split_role
    if np.any(role == ""):
        raise RuntimeError("Some events did not receive a source-file split role")
    return role


def build_or_load_features(data, manifest, rebuild, processes):
    if FEATURE_CACHE.exists() and FEATURE_NAMES.exists() and not rebuild:
        mt.log("loading cached detector features")
        with np.load(str(FEATURE_CACHE), allow_pickle=False) as source:
            groups = {group: source[group] for group in mt.GROUP_ORDER}
            hybrid_beta = source["hybrid_beta"]
            baseline_relR = source["baseline_relR"]
            cache_run = source["run"]
            cache_event = source["event"]
            cache_label = source["label"]
        if not (
            np.array_equal(cache_run, data["run"])
            and np.array_equal(cache_event, data["event"])
            and np.array_equal(cache_label, data["label"])
        ):
            raise RuntimeError("Feature cache event identity mismatch")
        group_names = json.loads(FEATURE_NAMES.read_text())
        audits = []
        return groups, group_names, hybrid_beta, baseline_relR, audits

    mt.log("extracting common-branch detector features from %d files" % sum(len(x) for x in manifest.values()))
    groups, group_names, hybrid_beta, baseline_relR, audits = mt.extract_all_features(
        data, manifest, processes=processes
    )
    payload = {
        "run": data["run"],
        "event": data["event"],
        "label": data["label"],
        "hybrid_beta": hybrid_beta,
        "baseline_relR": baseline_relR,
    }
    payload.update(groups)
    np.savez_compressed(str(FEATURE_CACHE), **payload)
    FEATURE_NAMES.write_text(json.dumps(group_names, indent=2, sort_keys=True) + "\n")
    mt.write_tsv(HERE / "extraction_audit.tsv", sorted(audits, key=lambda row: (row["label"], row["source"])))
    return groups, group_names, hybrid_beta, baseline_relR, audits


def metric_from_rows(method, probability, labels, masses_true, rows):
    p = np.asarray(probability[rows], dtype=np.float64)
    y = np.asarray(labels[rows], dtype=np.int8)
    mass_true = np.asarray(masses_true[rows], dtype=np.float64)
    soft_mass = mt.MASSES[0] + p * (mt.MASSES[1] - mt.MASSES[0])
    hard = p >= 0.5
    confident = (p <= 0.10) | (p >= 0.90)
    return {
        "method": method,
        "n": int(len(rows)),
        "n_c12": int(np.sum(y == 0)),
        "n_c13": int(np.sum(y == 1)),
        "auc": mt.auc(p, y),
        "balanced_accuracy": mt.balanced_accuracy(p, y),
        "accuracy": float(np.mean(hard == (y == 1))),
        "rejection50": mt.rejection_at_efficiency(p, y, 0.50),
        "rejection70": mt.rejection_at_efficiency(p, y, 0.70),
        "brier": float(np.mean((p - y) ** 2)),
        "balanced_logloss": mt.stable_logloss(p, y),
        "soft_mass_mae": float(np.mean(np.abs(soft_mass - mass_true))),
        "soft_mass_rmse": float(np.sqrt(np.mean((soft_mass - mass_true) ** 2))),
        "soft_mass_error_w68": mt.w68(soft_mass - mass_true),
        "confident_fraction": float(np.mean(confident)),
        "confident_accuracy": float(np.mean(hard[confident] == (y[confident] == 1))) if np.any(confident) else float("nan"),
        "c12_to_c13": int(np.sum(hard & (y == 0))),
        "c13_to_c12": int(np.sum((~hard) & (y == 1))),
    }


def event_metric_rows(probabilities, data, test_mask):
    labels = data["label"].astype(np.int8)
    masses_true = mt.MASSES[labels]
    bins = data["bin"].astype(np.int8)
    rows_out = []
    for bin_id in [-1, 0, 1, 2, 3]:
        selected = np.flatnonzero(test_mask & ((bins == bin_id) if bin_id >= 0 else True))
        for method in METHODS:
            row = metric_from_rows(method, probabilities[method], labels, masses_true, selected)
            row["bin"] = int(bin_id)
            rows_out.append(row)
    return rows_out


def physical_metric_rows(physical_mass2, data, test_mask):
    labels = data["label"].astype(np.int8)
    bins = data["bin"].astype(np.int8)
    true_mass2 = mt.MASSES[labels] ** 2
    rows_out = []
    for bin_id in [-1, 0, 1, 2, 3]:
        mask = test_mask & ((bins == bin_id) if bin_id >= 0 else True)
        for method, mass2 in physical_mass2.items():
            values = np.asarray(mass2[mask], dtype=np.float64)
            y = labels[mask]
            error = values - true_mass2[mask]
            positive = values >= 0.0
            mass = np.full(len(values), np.nan, dtype=np.float64)
            mass[positive] = np.sqrt(values[positive])
            mass_error = mass - mt.MASSES[y]
            rows_out.append({
                "method": method,
                "bin": int(bin_id),
                "n": int(np.sum(mask)),
                "auc": mt.auc(values, y),
                "rejection50": mt.rejection_at_efficiency(values, y, 0.50),
                "mass2_error_median": float(np.median(error)),
                "mass2_error_w68": mt.w68(error),
                "mass_error_median": float(np.nanmedian(mass_error)),
                "mass_error_w68": mt.w68(mass_error),
                "negative_mass2_fraction": float(np.mean(~positive)),
            })
    return rows_out


def calibration_rows(probabilities, labels, test_mask, bins_count=10):
    rows_out = []
    for method in ["raw_linear", "mass_proxy", "mass_energy", "all_physics", "all_detector", "nuisance_control"]:
        p = probabilities[method][test_mask]
        y = labels[test_mask]
        edges = np.linspace(0.0, 1.0, bins_count + 1)
        for index in range(bins_count):
            selected = (p >= edges[index]) & (
                p < edges[index + 1] if index + 1 < bins_count else p <= edges[index + 1]
            )
            if not np.any(selected):
                continue
            rows_out.append({
                "method": method,
                "probability_low": float(edges[index]),
                "probability_high": float(edges[index + 1]),
                "n": int(np.sum(selected)),
                "mean_prediction": float(np.mean(p[selected])),
                "observed_c13_fraction": float(np.mean(y[selected])),
            })
    return rows_out


def source_bootstrap(probabilities, data, test_mask, replicates):
    labels = data["label"].astype(np.int8)
    sources = data["source"].astype(np.int32)
    masses_true = mt.MASSES[labels]
    groups = {}
    keys = {0: [], 1: []}
    for label in (0, 1):
        for source in sorted(int(value) for value in np.unique(sources[test_mask & (labels == label)])):
            key = (label, source)
            groups[key] = np.flatnonzero(test_mask & (labels == label) & (sources == source))
            keys[label].append(key)
    rng = np.random.RandomState(290929)
    fields = ["auc", "balanced_accuracy", "soft_mass_mae", "rejection50"]
    samples = {method: {field: [] for field in fields} for method in probabilities}
    improvements = {field: [] for field in fields}
    for replicate in range(int(replicates)):
        chunks = []
        for label in (0, 1):
            draw = rng.randint(0, len(keys[label]), size=len(keys[label]))
            chunks.extend(groups[keys[label][index]] for index in draw)
        rows = np.concatenate(chunks)
        local = {}
        for method, probability in probabilities.items():
            metric = metric_from_rows(method, probability, labels, masses_true, rows)
            local[method] = metric
            for field in fields:
                samples[method][field].append(metric[field])
        for field in fields:
            sign = -1.0 if field == "soft_mass_mae" else 1.0
            improvements[field].append(
                sign * (local["all_physics"][field] - local["raw_linear"][field])
            )
    summary = {
        "replicates": int(replicates),
        "methods": {},
        "all_physics_minus_raw_improvement": {},
    }
    for method in samples:
        summary["methods"][method] = {}
        for field in fields:
            values = np.asarray(samples[method][field], dtype=np.float64)
            summary["methods"][method][field] = {
                "median": float(np.median(values)),
                "low95": float(np.percentile(values, 2.5)),
                "high95": float(np.percentile(values, 97.5)),
            }
    for field, values in improvements.items():
        values = np.asarray(values, dtype=np.float64)
        summary["all_physics_minus_raw_improvement"][field] = {
            "median": float(np.median(values)),
            "low95": float(np.percentile(values, 2.5)),
            "high95": float(np.percentile(values, 97.5)),
            "positive_fraction": float(np.mean(values > 0.0)),
        }
    return summary


def coefficient_audit(model_states, transformed_names):
    rows_out = []
    for method, state in model_states.items():
        if method == "shuffle_control":
            continue
        names = []
        groups = []
        for group in state["groups"]:
            names.extend(transformed_names[group])
            groups.extend([group] * len(transformed_names[group]))
        coefficient_matrix = np.vstack([
            np.asarray(state["models"][str(bin_id)]["parameters"], dtype=np.float64)[1:]
            for bin_id in range(len(mt.R_EDGES) - 1)
        ])
        rms = np.sqrt(np.mean(coefficient_matrix ** 2, axis=0))
        maximum = np.max(np.abs(coefficient_matrix), axis=0)
        order = np.argsort(-rms)
        for rank, index in enumerate(order[:50], 1):
            rows_out.append({
                "method": method,
                "rank": int(rank),
                "group": groups[index],
                "feature": names[index],
                "coefficient_rms": float(rms[index]),
                "coefficient_max_abs": float(maximum[index]),
            })
    return rows_out


def main():
    args = parse_args()
    started = time.time()
    manifest = mt.load_manifest(MANIFEST)
    role_by_source, split_rows = mt.split_sources(manifest)
    mt.write_tsv(HERE / "split_manifest.tsv", split_rows)

    data = load_truth_valid_data()
    role = assign_roles(data, role_by_source)
    train_mask = role == "train"
    valid_mask = role == "valid"
    test_mask = role == "test"
    mt.log(
        "truth-valid events train=%d valid=%d test=%d total=%d"
        % (np.sum(train_mask), np.sum(valid_mask), np.sum(test_mask), len(role))
    )

    groups, group_names, hybrid_beta, baseline_relR, audits = build_or_load_features(
        data, manifest, args.rebuild_features, args.processes
    )
    if audits:
        max_u_replay = max(row["u_replay_max"] for row in audits)
        max_r_replay = max(row["R_replay_max"] for row in audits)
    else:
        max_u_replay = 0.0
        max_r_replay = 0.0
    if max_u_replay > 2e-6 or max_r_replay > 2e-6:
        raise RuntimeError("Measurement replay mismatch u=%g R=%g" % (max_u_replay, max_r_replay))

    preprocessors = {}
    transformed = {}
    transformed_names = {}
    feature_counts = {}
    for group in mt.GROUP_ORDER:
        state = mt.fit_preprocessor(groups[group], group_names[group], train_mask)
        preprocessors[group] = state
        transformed[group] = mt.apply_preprocessor(groups[group], state)
        transformed_names[group] = state["output_names"]
        feature_counts[group] = {
            "raw": int(groups[group].shape[1]),
            "retained": int(len(state["feature_names"])),
            "transformed": int(transformed[group].shape[1]),
        }
        mt.log(
            "group %-10s raw=%d retained=%d transformed=%d"
            % (
                group,
                feature_counts[group]["raw"],
                feature_counts[group]["retained"],
                feature_counts[group]["transformed"],
            )
        )
    del groups

    labels = data["label"].astype(np.int8)
    bins = data["bin"].astype(np.int8)
    model_states = {}
    raw_logits = {}
    scans = {}
    regularization_candidates = [3e-4, 1e-3, 3e-3]

    for method in ["mass_proxy", "mass_energy", "all_physics", "all_detector", "nuisance_control"]:
        x = np.column_stack([transformed[group] for group in mt.MODEL_GROUPS[method]]).astype(np.float32)
        mt.log("fit %s with %d transformed features" % (method, x.shape[1]))
        best, scan = mt.choose_regularization(
            x, labels, bins, train_mask, valid_mask, regularization_candidates
        )
        model_states[method] = {
            "groups": mt.MODEL_GROUPS[method],
            "feature_count": int(x.shape[1]),
            "regularization": float(best["regularization"]),
            "models": best["models"],
        }
        scans[method] = [
            {
                "regularization": item["regularization"],
                "train_macro_auc": item["train_macro_auc"],
                "valid_macro_auc": item["valid_macro_auc"],
            }
            for item in scan
        ]
        raw_logits[method] = mt.predict_binned(x, bins, best["models"])
        del x

    # Negative control: retain the full feature map but shuffle isotope labels
    # independently within each reconstructed-rigidity bin.
    x_all = np.column_stack([
        transformed[group] for group in mt.MODEL_GROUPS["all_detector"]
    ]).astype(np.float32)
    shuffled = labels.copy()
    rng = np.random.RandomState(130929)
    for bin_id in range(len(mt.R_EDGES) - 1):
        rows = np.flatnonzero(train_mask & (bins == bin_id))
        shuffled[rows] = shuffled[rng.permutation(rows)]
    shuffle_models = {}
    for bin_id in range(len(mt.R_EDGES) - 1):
        rows = train_mask & (bins == bin_id)
        shuffle_models[str(bin_id)] = mt.fit_logistic(
            x_all[rows], shuffled[rows],
            model_states["all_detector"]["regularization"], maxiter=120,
        )
    model_states["shuffle_control"] = {
        "groups": mt.MODEL_GROUPS["all_detector"],
        "feature_count": int(x_all.shape[1]),
        "regularization": model_states["all_detector"]["regularization"],
        "models": shuffle_models,
        "training_labels": "shuffled independently inside each R bin",
    }
    raw_logits["shuffle_control"] = mt.predict_binned(x_all, bins, shuffle_models)
    del x_all

    raw_m2 = (mt.CHARGE * data["R"]) ** 2 * (data["u"] ** 2 - 1.0)
    hybrid_valid = np.isfinite(hybrid_beta) & (hybrid_beta > 0.2) & (hybrid_beta < 1.5)
    hybrid_u = np.where(hybrid_valid, 1.0 / hybrid_beta, data["u"])
    hybrid_m2 = (mt.CHARGE * data["R"]) ** 2 * (hybrid_u ** 2 - 1.0)
    raw_logits["raw_linear"] = raw_m2
    raw_logits["raw_hybrid"] = hybrid_m2

    calibrations = {}
    probabilities = {}
    for method in METHODS:
        score_mode = method in ("raw_linear", "raw_hybrid")
        calibration, probability = mt.fit_calibration(
            raw_logits[method], labels, bins, valid_mask, score_mode=score_mode
        )
        calibrations[method] = calibration
        probabilities[method] = probability

    metric_rows = event_metric_rows(probabilities, data, test_mask)
    mt.write_tsv(HERE / "event_metrics.tsv", metric_rows)

    physical_mass2 = {
        "linear_raw": raw_m2,
        "hybrid_raw": hybrid_m2,
        "perfect_beta_raw_R": (mt.CHARGE * data["R"]) ** 2 * (data["u_true"] ** 2 - 1.0),
        "raw_beta_perfect_R": (mt.CHARGE * data["R_true"]) ** 2 * (data["u"] ** 2 - 1.0),
        "perfect_beta_perfect_R": (mt.CHARGE * data["R_true"]) ** 2 * (data["u_true"] ** 2 - 1.0),
    }
    physical_rows = physical_metric_rows(physical_mass2, data, test_mask)
    mt.write_tsv(HERE / "physical_mass_metrics.tsv", physical_rows)
    mt.write_tsv(HERE / "calibration.tsv", calibration_rows(probabilities, labels, test_mask))
    mt.write_tsv(HERE / "coefficient_audit.tsv", coefficient_audit(model_states, transformed_names))

    bootstrap_methods = {
        method: probabilities[method]
        for method in METHODS
        if method != "shuffle_control"
    }
    bootstrap = source_bootstrap(bootstrap_methods, data, test_mask, args.bootstrap)
    (HERE / "bootstrap.json").write_text(
        json.dumps(mt.jsonable(bootstrap), indent=2, sort_keys=True) + "\n"
    )

    model_document = {
        "seed_text": mt.SEED_TEXT,
        "masses_GeV": mt.MASSES,
        "rigidity_edges_GV": mt.R_EDGES,
        "common_input_branches": mt.BRANCHES,
        "feature_counts": feature_counts,
        "preprocessors": preprocessors,
        "models": model_states,
        "calibrations": calibrations,
        "regularization_scans": scans,
        "input_exclusions": [
            "run", "event", "source index", "production version",
            "mpar/mch isotope codes", "MC momentum", "MC beta",
            "tof_betahmc", "all production-exclusive branches",
        ],
    }
    (HERE / "frozen_model.json").write_text(
        json.dumps(mt.jsonable(model_document), indent=2, sort_keys=True) + "\n"
    )

    prediction_payload = {
        "run": data["run"],
        "event": data["event"],
        "source": data["source"],
        "label": labels,
        "bin": bins,
        "role": role.astype("S5"),
        "u": data["u"],
        "R": data["R"],
        "u_true": data["u_true"],
        "R_true": data["R_true"],
        "hybrid_beta": hybrid_beta,
        "baseline_relR": baseline_relR,
        "mass2_raw_linear": raw_m2,
        "mass2_raw_hybrid": hybrid_m2,
        "mass2_perfect_beta_raw_R": physical_mass2["perfect_beta_raw_R"],
        "mass2_raw_beta_perfect_R": physical_mass2["raw_beta_perfect_R"],
    }
    for method in METHODS:
        prediction_payload["score_" + method] = raw_logits[method]
        prediction_payload["p_c13_" + method] = probabilities[method]
        prediction_payload["mass_soft_" + method] = (
            mt.MASSES[0] + probabilities[method] * (mt.MASSES[1] - mt.MASSES[0])
        )
    np.savez_compressed(str(HERE / "predictions.npz"), **prediction_payload)

    overall = {row["method"]: row for row in metric_rows if row["bin"] == -1}
    physical_overall = {row["method"]: row for row in physical_rows if row["bin"] == -1}
    source_files = {
        split_role: {
            isotope: int(sum(row["role"] == split_role and row["isotope"] == isotope for row in split_rows))
            for isotope in ("C12", "C13")
        }
        for split_role in ("train", "valid", "test")
    }
    summary = {
        "events": {
            "total_truth_valid": int(len(role)),
            "train": int(np.sum(train_mask)),
            "valid": int(np.sum(valid_mask)),
            "test": int(np.sum(test_mask)),
            "test_c12": int(np.sum(test_mask & (labels == 0))),
            "test_c13": int(np.sum(test_mask & (labels == 1))),
        },
        "source_files": source_files,
        "feature_counts": feature_counts,
        "measurement_replay": {
            "max_abs_u_difference": float(max_u_replay),
            "max_abs_R_difference_GV": float(max_r_replay),
            "hybrid_beta_fallback_fraction": float(np.mean(~hybrid_valid)),
        },
        "overall_event_metrics": overall,
        "overall_physical_mass_metrics": physical_overall,
        "bootstrap": bootstrap,
        "interpretation_guardrail": (
            "Soft mass is the equal-prior posterior mean m12 + P(C13|x)(m13-m12). "
            "Its visual peak width is not a physical mass resolution; event AUC, confusion, "
            "calibration and held-out-source bootstrap are the microscopic tests."
        ),
        "known_limitation": (
            "C12 B1308 and C13 B1315 are different productions. Common branches and whole-file "
            "holdout remove direct schema and event leakage, but cannot remove production/isotope "
            "confounding. The nuisance-only model quantifies the remaining shortcut."
        ),
        "elapsed_seconds": float(time.time() - started),
    }
    (HERE / "results.json").write_text(
        json.dumps(mt.jsonable(summary), indent=2, sort_keys=True) + "\n"
    )

    freeze = {
        "common_script_sha256": mt.sha256_file(HERE / "mass_toy_common.py"),
        "run_script_sha256": mt.sha256_file(Path(__file__)),
        "expanded_predictions_sha256": mt.sha256_file(EXPANDED_PREDICTIONS),
        "expanded_manifest_sha256": mt.sha256_file(MANIFEST),
        "split_manifest_sha256": mt.sha256_file(HERE / "split_manifest.tsv"),
        "feature_cache_sha256": mt.sha256_file(FEATURE_CACHE),
        "frozen_model_sha256": mt.sha256_file(HERE / "frozen_model.json"),
        "predictions_sha256": mt.sha256_file(HERE / "predictions.npz"),
    }
    (HERE / "freeze.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")

    mt.log("COMPLETED overall metrics")
    for method in METHODS:
        row = overall[method]
        mt.log(
            "%-18s AUC %.5f  bal.acc %.5f  MAE %.4f GeV  rejection50 %.2f"
            % (method, row["auc"], row["balanced_accuracy"], row["soft_mass_mae"], row["rejection50"])
        )
    mt.log("elapsed %.1fs" % (time.time() - started))


if __name__ == "__main__":
    main()
