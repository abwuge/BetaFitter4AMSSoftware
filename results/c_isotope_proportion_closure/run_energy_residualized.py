#!/usr/bin/env python3
"""Energy-loss residualization ablation for C12/C13 weak beta closure."""
import json
import time
from pathlib import Path

import numpy as np

import run_closure as rc


def load_cache(path):
    source = np.load(str(path))
    train, valid = {}, {}
    for name in source.files:
        prefix, key = name.split("__", 1)
        (train if prefix == "train" else valid)[key] = source[name]
    return train, valid


def fit_bin_response(train, value):
    """Robust quadratic response versus log(raw inverse beta), independently per R bin."""
    coefficients = np.zeros((rc.N_BINS, 3, value.shape[1]), dtype=float)
    for b in range(rc.N_BINS):
        select = train["bin"] == b
        x = np.log(train["u"][select])
        center = np.median(x)
        scale = max(1e-6, 0.741301109 * np.diff(np.percentile(x, [25, 75]))[0])
        z = np.clip((x - center) / scale, -4.0, 4.0)
        design = np.column_stack([np.ones(len(z)), z, z * z])
        target = value[select]
        penalty = np.diag([1e-6, 1e-3, 1e-2])
        beta = np.linalg.solve(design.T.dot(design) + penalty, design.T.dot(target))
        for _ in range(15):
            residual = target - design.dot(beta)
            width = np.maximum(
                1e-4,
                1.4826 * np.median(np.abs(residual - np.median(residual, axis=0)), axis=0),
            )
            scaled = np.abs(residual) / (2.5 * width[None, :])
            weight = np.minimum(1.0, 1.0 / np.maximum(scaled, 1e-12))
            updated = np.empty_like(beta)
            for j in range(target.shape[1]):
                weighted = design * weight[:, j, None]
                updated[:, j] = np.linalg.solve(
                    design.T.dot(weighted) + penalty,
                    design.T.dot(weight[:, j] * target[:, j]),
                )
            if np.max(np.abs(updated - beta)) < 1e-9:
                beta = updated
                break
            beta = updated
        coefficients[b] = beta
        coefficients[b, 0, :] = coefficients[b, 0, :]
        # Store center/scale in an extra side table returned separately.
    return coefficients


def response_state(train):
    energy = train["continuous"][:, 18:29]
    common = np.column_stack([np.mean(energy[:, :4], axis=1), np.mean(energy[:, 4:], axis=1)])
    coefficients = np.zeros((rc.N_BINS, 3, 2), dtype=float)
    centers = np.zeros(rc.N_BINS, dtype=float)
    scales = np.zeros(rc.N_BINS, dtype=float)
    for b in range(rc.N_BINS):
        select = train["bin"] == b
        x = np.log(train["u"][select])
        centers[b] = np.median(x)
        scales[b] = max(1e-6, 0.741301109 * np.diff(np.percentile(x, [25, 75]))[0])
        z = np.clip((x - centers[b]) / scales[b], -4.0, 4.0)
        design = np.column_stack([np.ones(len(z)), z, z * z])
        target = common[select]
        penalty = np.diag([1e-6, 1e-3, 1e-2])
        beta = np.linalg.solve(design.T.dot(design) + penalty, design.T.dot(target))
        for _ in range(15):
            residual = target - design.dot(beta)
            width = np.maximum(
                1e-4,
                1.4826 * np.median(
                    np.abs(residual - np.median(residual, axis=0)), axis=0
                ),
            )
            scaled = np.abs(residual) / (2.5 * width[None, :])
            weight = np.minimum(1.0, 1.0 / np.maximum(scaled, 1e-12))
            updated = np.empty_like(beta)
            for j in range(2):
                weighted = design * weight[:, j, None]
                updated[:, j] = np.linalg.solve(
                    design.T.dot(weighted) + penalty,
                    design.T.dot(weight[:, j] * target[:, j]),
                )
            if np.max(np.abs(updated - beta)) < 1e-9:
                beta = updated
                break
            beta = updated
        coefficients[b] = beta
    state = {"coefficients": coefficients, "centers": centers, "scales": scales}
    residual = apply_response(train, state)
    bin_median = np.zeros((rc.N_BINS, 2), dtype=float)
    for b in range(rc.N_BINS):
        bin_median[b] = np.median(residual[train["bin"] == b], axis=0)
    state["bin_median"] = bin_median
    return state


def apply_response(data, state):
    energy = data["continuous"][:, 18:29]
    common = np.column_stack([np.mean(energy[:, :4], axis=1), np.mean(energy[:, 4:], axis=1)])
    residual = np.zeros_like(common)
    for b in range(rc.N_BINS):
        select = data["bin"] == b
        z = np.clip(
            (np.log(data["u"][select]) - state["centers"][b]) / state["scales"][b],
            -4.0,
            4.0,
        )
        design = np.column_stack([np.ones(len(z)), z, z * z])
        residual[select] = common[select] - design.dot(state["coefficients"][b])
    if "bin_median" in state:
        residual -= state["bin_median"][data["bin"]]
    return residual


def energy_raw(data, kind, response=None):
    energy = data["continuous"][:, 18:29]
    tof = energy[:, :4]
    tracker = energy[:, 4:]
    tof_shape = tof[:, :3] - np.mean(tof, axis=1, keepdims=True)
    tracker_shape = tracker[:, :6] - np.mean(tracker, axis=1, keepdims=True)
    shape = np.column_stack([tof_shape, tracker_shape])
    if kind == "shape":
        return shape
    if kind == "fixed":
        physical = 2.0 * np.log(data["u"])
        common = np.column_stack(
            [np.mean(tof, axis=1) - physical, np.mean(tracker, axis=1) - physical]
        )
        return common
    if kind == "fitted":
        return apply_response(data, response)
    if kind == "shape_fitted":
        return np.column_stack([shape, apply_response(data, response)])
    raise ValueError(kind)


def standardizer(raw):
    low, high = np.percentile(raw, [1.0, 99.0], axis=0)
    clipped = np.minimum(np.maximum(raw, low), high)
    center = np.median(clipped, axis=0)
    quartiles = np.percentile(clipped, [25.0, 75.0], axis=0)
    scale = np.maximum(1e-6, 0.741301109 * (quartiles[1] - quartiles[0]))
    return {"low": low, "high": high, "center": center, "scale": scale}


def apply_standardizer(raw, state):
    clipped = np.minimum(np.maximum(raw, state["low"]), state["high"])
    return np.clip((clipped - state["center"]) / state["scale"], -4.0, 4.0)


def append_energy(base_train, base_valid, raw_train, raw_valid):
    state = standardizer(raw_train)
    z_train = apply_standardizer(raw_train, state)
    z_valid = apply_standardizer(raw_valid, state)
    derived_train = np.column_stack([z_train, z_train * z_train])
    derived_valid = np.column_stack([z_valid, z_valid * z_valid])
    center = np.mean(derived_train, axis=0)
    scale = np.std(derived_train, axis=0)
    scale[scale < 1e-6] = 1.0
    derived_train = np.clip((derived_train - center) / scale, -5.0, 5.0)
    derived_valid = np.clip((derived_valid - center) / scale, -5.0, 5.0)
    return (
        np.column_stack([base_train, derived_train]),
        np.column_stack([base_valid, derived_valid]),
    )


def ratio_project(x_train, x_valid, train):
    z = x_train[:, 1:]
    difference = []
    for b in range(rc.N_BINS):
        low = (train["bin"] == b) & (train["bag"] == 0)
        high = (train["bin"] == b) & (train["bag"] == 1)
        difference.append(np.mean(z[low], axis=0) - np.mean(z[high], axis=0))
    difference = np.vstack(difference)
    _, singular, vt = np.linalg.svd(difference, full_matrices=True)
    tolerance = max(difference.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > tolerance))
    null = vt[rank:].T
    return (
        np.column_stack([np.ones(len(x_train)), z.dot(null)]),
        np.column_stack([np.ones(len(x_valid)), x_valid[:, 1:].dot(null)]),
        singular,
        rank,
    )


def main():
    destination = Path(__file__).resolve().parent
    started = time.time()
    train, valid = load_cache(destination / "controlled_cache.npz")
    first = np.load(str(destination / "predictions.npz"))
    for key in ("run", "event", "label", "bin"):
        if not np.array_equal(valid[key], first[key]):
            raise RuntimeError("Validation ordering mismatch: " + key)

    base_state = rc.make_feature_state(train, "no_edep")
    base_train = rc.apply_feature_state(train, base_state)
    base_valid = rc.apply_feature_state(valid, base_state)
    response = response_state(train)

    matrices_train = {}
    matrices_valid = {}
    for kind in ("shape", "fixed", "fitted", "shape_fitted"):
        matrices_train[kind], matrices_valid[kind] = append_energy(
            base_train,
            base_valid,
            energy_raw(train, kind, response),
            energy_raw(valid, kind, response),
        )
        rc.log("%s dimensions=%d" % (kind, matrices_train[kind].shape[1]))
    projected_train, projected_valid, singular, rank = ratio_project(
        matrices_train["shape_fitted"], matrices_valid["shape_fitted"], train
    )
    matrices_train["shape_fitted_ratio"] = projected_train
    matrices_valid["shape_fitted_ratio"] = projected_valid
    rc.log(
        "shape_fitted_ratio dimensions=%d constraint_rank=%d"
        % (projected_train.shape[1], rank)
    )

    weak_train = {key: train[key] for key in ["u", "R", "relR", "bin"]}
    weak_valid = {key: valid[key] for key in ["u", "R", "relR", "bin"]}
    frozen = json.load(open(str(destination / "weak_frozen_models.json")))
    initial_logsigma = np.asarray(frozen["models"]["raw"]["logsigma"], dtype=float)

    weak_models = {}
    for kind in ("shape", "fixed", "fitted", "shape_fitted", "shape_fitted_ratio"):
        weak_models[kind] = rc.fit_density(
            "weak_energy_" + kind,
            weak_train,
            matrices_train[kind],
            train["prior_c13"],
            True,
            initial_logsigma,
            starts=1,
        )

    oracle_models = {}
    for kind in ("shape", "fitted", "shape_fitted", "shape_fitted_ratio"):
        oracle_models[kind] = rc.fit_oracle(
            "oracle_energy_" + kind,
            train,
            matrices_train[kind],
            regularization=20.0,
        )

    corrected = {
        "raw": valid["u"].copy(),
        "weak_no_edep_existing": first["u_weak_no_edep"],
        "oracle_full_existing": first["u_oracle_full"],
    }
    heldout = {
        "raw": float("nan"),
        "weak_no_edep_existing": float("nan"),
        "oracle_full_existing": float("nan"),
    }
    diagnostics = {}
    for kind, model in weak_models.items():
        name = "weak_energy_" + kind
        x_train = matrices_train[kind]
        x_valid = matrices_valid[kind]
        corrected[name] = valid["u"] - x_valid.dot(model["theta"])
        density = rc.MixtureDensity(weak_valid, x_valid, valid["prior_c13"])
        heldout[name] = float(np.mean(density.event_nll(model["parameters"])))
        classifier = rc.fit_logistic(x_train, train["label"])
        correction = valid["u"] - corrected[name]
        truth = valid["truth_valid"]
        diagnostics[name] = {
            "feature_label_auc": rc.auc(x_valid.dot(classifier), valid["label"]),
            "correction_label_auc": rc.auc(correction, valid["label"]),
            "correction_rms": rc.rms(correction),
            "correlation_with_raw_truth_residual": float(
                np.corrcoef(
                    correction[truth],
                    valid["u"][truth] - valid["u_true"][truth],
                )[0, 1]
            ),
        }
    for kind, model in oracle_models.items():
        name = "oracle_energy_" + kind
        x_valid = matrices_valid[kind]
        corrected[name] = valid["u"] - x_valid.dot(model["theta"])
        density = rc.MixtureDensity(weak_valid, x_valid, valid["prior_c13"])
        heldout[name] = float(
            np.mean(
                density.event_nll(
                    weak_models[kind]["parameters"], corrected=corrected[name]
                )
            )
        )

    rows = []
    bin_rows = []
    method_order = list(corrected)
    for name in method_order:
        method_rows, mass2 = rc.method_metrics(name, corrected[name], valid, heldout[name])
        rows.extend(method_rows)
        bin_rows.extend(rc.bin_metrics(name, corrected[name], mass2, valid))
    raw_width = {
        row["group"]: row["w68_du"] for row in rows if row["method"] == "raw"
    }
    for row in rows:
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_du"] / raw_width[row["group"]]
        )
    rc.write_tsv(destination / "energy_residualized_metrics.tsv", rows)
    rc.write_tsv(destination / "energy_residualized_bin_metrics.tsv", bin_rows)

    overall = {row["method"]: row for row in rows if row["group"] == "all"}
    output = {
        "response_state": response,
        "ratio_projection": {"singular_values": singular, "rank": rank},
        "diagnostics": diagnostics,
        "overall_metrics": overall,
        "weak_models": weak_models,
        "oracle_models": oracle_models,
        "elapsed_seconds": time.time() - started,
        "script_sha256": rc.sha256_file(Path(__file__)),
    }
    (destination / "energy_residualized_results.json").write_text(
        json.dumps(rc.jsonable(output), indent=2, sort_keys=True)
    )
    rc.log(
        "DONE energy residualization elapsed=%.1fs results=%s"
        % (time.time() - started, destination / "energy_residualized_results.json")
    )


if __name__ == "__main__":
    main()
