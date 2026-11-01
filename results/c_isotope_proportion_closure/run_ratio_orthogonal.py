#!/usr/bin/env python3
"""Label-free ratio-orthogonal energy-loss closure for C12/C13."""
import json
import time
from pathlib import Path

import numpy as np

import run_closure as rc


def save_cache(path, train, valid):
    arrays = {}
    for prefix, data in (("train", train), ("valid", valid)):
        for key, value in data.items():
            arrays[prefix + "__" + key] = value
    np.savez_compressed(str(path), **arrays)


def load_cache(path):
    source = np.load(str(path))
    train, valid = {}, {}
    for name in source.files:
        prefix, key = name.split("__", 1)
        (train if prefix == "train" else valid)[key] = source[name]
    return train, valid


def ratio_null_projection(x_train, train):
    """Constrain mean correction to be equal across ratio bags in every R bin."""
    z = x_train[:, 1:]
    rows = []
    for b in range(rc.N_BINS):
        low = (train["bin"] == b) & (train["bag"] == 0)
        high = (train["bin"] == b) & (train["bag"] == 1)
        rows.append(np.mean(z[low], axis=0) - np.mean(z[high], axis=0))
    difference = np.vstack(rows)
    _, singular, vt = np.linalg.svd(difference, full_matrices=True)
    tolerance = max(difference.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > tolerance))
    null_basis = vt[rank:].T
    return null_basis, difference, singular, rank


def projected_matrix(x, null_basis):
    return np.column_stack([np.ones(len(x)), x[:, 1:].dot(null_basis)])


def main():
    destination = Path(__file__).resolve().parent
    started = time.time()
    cache = destination / "controlled_cache.npz"

    if cache.exists():
        train, valid = load_cache(cache)
        rc.log("loaded controlled cache train=%d valid=%d" % (len(train["u"]), len(valid["u"])))
    else:
        c12_files = rc.choose_files(rc.C12_PATTERN, rc.FILES_PER_ROLE)
        c13_files = rc.choose_files(rc.C13_PATTERN, rc.C13_FILES_PER_ROLE)
        pools = {0: {}, 1: {}}
        audit = []
        for role in (0, 1):
            pools[role][0] = rc.extract_pool(
                c12_files[role], 0, 67, rc.MASSES[0], role, audit
            )
            pools[role][1] = rc.extract_pool(
                c13_files[role], 1, 117, rc.MASSES[1], role, audit
            )
        train, _ = rc.controlled_mix(
            pools[0], rc.TRAIN_FRACTIONS, rc.TRAIN_TARGET_PER_BAG_BIN, rc.SEED
        )
        valid, _ = rc.controlled_mix(
            pools[1], rc.VALID_FRACTIONS, rc.VALID_TARGET_PER_BIN, rc.SEED + 1
        )
        save_cache(cache, train, valid)
        rc.log("saved controlled cache %s" % cache)

    # Verify that the cache is exactly the validation sample from the first closure.
    first_predictions = np.load(str(destination / "predictions.npz"))
    for key in ("run", "event", "label", "bin"):
        if not np.array_equal(valid[key], first_predictions[key]):
            raise RuntimeError("Controlled-cache validation ordering mismatch for " + key)

    full_state = rc.make_feature_state(train, "full")
    full_train = rc.apply_feature_state(train, full_state)
    full_valid = rc.apply_feature_state(valid, full_state)
    null_basis, bag_difference, singular, rank = ratio_null_projection(full_train, train)
    x_train = projected_matrix(full_train, null_basis)
    x_valid = projected_matrix(full_valid, null_basis)

    transformed_differences = []
    for b in range(rc.N_BINS):
        low = (train["bin"] == b) & (train["bag"] == 0)
        high = (train["bin"] == b) & (train["bag"] == 1)
        transformed_differences.append(
            np.mean(x_train[low, 1:], axis=0) - np.mean(x_train[high, 1:], axis=0)
        )
    transformed_differences = np.vstack(transformed_differences)
    max_feature_mean_difference = float(np.max(np.abs(transformed_differences)))
    rc.log(
        "ratio-null projection original=%d rank=%d projected=%d max_bag_feature_mean_diff=%.3g"
        % (full_train.shape[1], rank, x_train.shape[1], max_feature_mean_difference)
    )

    weak_train = {key: train[key] for key in ["u", "R", "relR", "bin"]}
    weak_valid = {key: valid[key] for key in ["u", "R", "relR", "bin"]}
    frozen = json.load(open(str(destination / "weak_frozen_models.json")))
    initial_logsigma = np.asarray(frozen["models"]["raw"]["logsigma"], dtype=float)

    model = rc.fit_density(
        "weak_full_ratio_orthogonal",
        weak_train,
        x_train,
        train["prior_c13"],
        True,
        initial_logsigma,
        starts=2,
    )
    wrong_prior = np.clip(train["prior_c13"] + 0.10, 0.01, 0.99)
    wrong_model = rc.fit_density(
        "weak_full_ratio_orthogonal_wrong_fraction",
        weak_train,
        x_train,
        wrong_prior,
        True,
        model["logsigma"],
        starts=1,
    )
    oracle = rc.fit_oracle(
        "oracle_full_ratio_orthogonal", train, x_train, regularization=20.0
    )

    models = {
        "weak_full_ratio_orthogonal": model,
        "weak_full_ratio_orthogonal_wrong_fraction": wrong_model,
        "oracle_full_ratio_orthogonal": oracle,
        "null_basis": null_basis,
        "bag_difference_singular_values": singular,
        "bag_difference_rank": rank,
        "max_transformed_bag_feature_mean_difference": max_feature_mean_difference,
        "script_sha256": rc.sha256_file(Path(__file__)),
    }
    (destination / "ratio_orthogonal_models.json").write_text(
        json.dumps(rc.jsonable(models), indent=2, sort_keys=True)
    )

    corrected = {
        "raw": valid["u"].copy(),
        "weak_full_ratio_orthogonal": valid["u"] - x_valid.dot(model["theta"]),
        "weak_full_ratio_orthogonal_wrong_fraction": (
            valid["u"] - x_valid.dot(wrong_model["theta"])
        ),
        "oracle_full_ratio_orthogonal": valid["u"] - x_valid.dot(oracle["theta"]),
    }
    nominal_density = rc.MixtureDensity(
        weak_valid, x_valid, valid["prior_c13"]
    )
    heldout_nll = {
        "raw": float("nan"),
        "weak_full_ratio_orthogonal": float(
            np.mean(nominal_density.event_nll(model["parameters"]))
        ),
        "weak_full_ratio_orthogonal_wrong_fraction": float(
            np.mean(nominal_density.event_nll(wrong_model["parameters"]))
        ),
        "oracle_full_ratio_orthogonal": float(
            np.mean(
                nominal_density.event_nll(
                    model["parameters"],
                    corrected=corrected["oracle_full_ratio_orthogonal"],
                )
            )
        ),
    }

    rows = []
    bin_rows = []
    mass2 = {}
    for name in (
        "raw",
        "weak_full_ratio_orthogonal",
        "weak_full_ratio_orthogonal_wrong_fraction",
        "oracle_full_ratio_orthogonal",
    ):
        method_rows, values = rc.method_metrics(
            name, corrected[name], valid, heldout_nll[name]
        )
        rows.extend(method_rows)
        bin_rows.extend(rc.bin_metrics(name, corrected[name], values, valid))
        mass2[name] = values
    raw_width = {
        row["group"]: row["w68_du"]
        for row in rows
        if row["method"] == "raw"
    }
    for row in rows:
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_du"] / raw_width[row["group"]]
        )
    rc.write_tsv(destination / "ratio_orthogonal_metrics.tsv", rows)
    rc.write_tsv(destination / "ratio_orthogonal_bin_metrics.tsv", bin_rows)

    classifier = rc.fit_logistic(x_train, train["label"])
    correction = valid["u"] - corrected["weak_full_ratio_orthogonal"]
    correction_wrong = (
        valid["u"] - corrected["weak_full_ratio_orthogonal_wrong_fraction"]
    )
    correction_oracle = (
        valid["u"] - corrected["oracle_full_ratio_orthogonal"]
    )
    diagnostics = {
        "projected_feature_label_auc": rc.auc(
            x_valid.dot(classifier), valid["label"]
        ),
        "weak_correction_label_auc": rc.auc(correction, valid["label"]),
        "wrong_fraction_correction_label_auc": rc.auc(
            correction_wrong, valid["label"]
        ),
        "oracle_correction_label_auc": rc.auc(
            correction_oracle, valid["label"]
        ),
        "weak_correction_rms": rc.rms(correction),
        "wrong_fraction_correction_rms": rc.rms(correction_wrong),
        "oracle_correction_rms": rc.rms(correction_oracle),
    }
    overall = {
        row["method"]: row for row in rows if row["group"] == "all"
    }
    summary = {
        "projection": {
            "original_dimensions": int(full_train.shape[1]),
            "constraint_rank": rank,
            "projected_dimensions": int(x_train.shape[1]),
            "singular_values": singular,
            "max_transformed_bag_feature_mean_difference": max_feature_mean_difference,
        },
        "diagnostics": diagnostics,
        "overall_metrics": overall,
        "elapsed_seconds": time.time() - started,
        "script_sha256": rc.sha256_file(Path(__file__)),
    }
    (destination / "ratio_orthogonal_results.json").write_text(
        json.dumps(rc.jsonable(summary), indent=2, sort_keys=True)
    )
    np.savez_compressed(
        str(destination / "ratio_orthogonal_predictions.npz"),
        run=valid["run"],
        event=valid["event"],
        label=valid["label"],
        bin=valid["bin"],
        u_true=valid["u_true"],
        truth_valid=valid["truth_valid"],
        **{"u_" + key: value for key, value in corrected.items()}
    )
    rc.log(
        "DONE ratio-orthogonal elapsed=%.1fs results=%s"
        % (time.time() - started, destination / "ratio_orthogonal_results.json")
    )


if __name__ == "__main__":
    main()
