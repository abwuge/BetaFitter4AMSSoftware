#!/usr/bin/env python3
"""Anchored weak beta/R scan for C12/C13.

Both detector heads are forced to have zero mean in every reconstructed-
rigidity bin on the training sample, and linear directions whose bag means
change with the synthetic isotope proportion are projected out.  Thus the
existing beta and rigidity scales are retained; the mass likelihood may only
learn event-level resolution corrections from within-bin detector redundancy.
Event labels and MC truth are not consulted until all models are frozen.
"""

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "c_isotope_proportion_closure"
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(HERE))

import run_closure as rc
import run_joint_identifiability_scan as js


def log(message):
    print(message, flush=True)


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


def anchored_projection(x_train, x_valid, train, valid, name):
    """Per-bin zero mean plus null projection of bag-mean differences."""
    x_train = np.asarray(x_train, dtype=np.float64)
    x_valid = np.asarray(x_valid, dtype=np.float64)
    centers = np.zeros((rc.N_BINS, x_train.shape[1]), dtype=np.float64)
    centered_train = np.empty_like(x_train)
    centered_valid = np.empty_like(x_valid)
    for b in range(rc.N_BINS):
        selected = train["bin"] == b
        centers[b] = np.mean(x_train[selected], axis=0)
        centered_train[selected] = x_train[selected] - centers[b]
        selected_valid = valid["bin"] == b
        centered_valid[selected_valid] = x_valid[selected_valid] - centers[b]

    constraints = []
    for b in range(rc.N_BINS):
        bag0 = (train["bin"] == b) & (train["bag"] == 0)
        bag1 = (train["bin"] == b) & (train["bag"] == 1)
        constraints.append(
            np.mean(centered_train[bag0], axis=0)
            - np.mean(centered_train[bag1], axis=0)
        )
    constraints = np.asarray(constraints)
    _, singular, vt = np.linalg.svd(constraints, full_matrices=True)
    tolerance = max(constraints.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > tolerance))
    null = vt[rank:].T
    projected_train = centered_train.dot(null)
    projected_valid = centered_valid.dot(null)

    scale = np.std(projected_train, axis=0)
    keep = np.isfinite(scale) & (scale > 1e-8)
    projected_train = projected_train[:, keep] / scale[keep]
    projected_valid = projected_valid[:, keep] / scale[keep]
    null = null[:, keep]
    scale = scale[keep]

    # A leading all-zero column occupies the intercept slot expected by the
    # shared regularizer without permitting any scale shift.
    transformed_train = np.column_stack(
        [np.zeros(len(projected_train)), projected_train]
    )
    transformed_valid = np.column_stack(
        [np.zeros(len(projected_valid)), projected_valid]
    )

    max_bin_mean = 0.0
    max_bag_difference = 0.0
    for b in range(rc.N_BINS):
        in_bin = train["bin"] == b
        max_bin_mean = max(
            max_bin_mean,
            float(np.max(np.abs(np.mean(projected_train[in_bin], axis=0)))),
        )
        bag0 = in_bin & (train["bag"] == 0)
        bag1 = in_bin & (train["bag"] == 1)
        max_bag_difference = max(
            max_bag_difference,
            float(
                np.max(
                    np.abs(
                        np.mean(projected_train[bag0], axis=0)
                        - np.mean(projected_train[bag1], axis=0)
                    )
                )
            ),
        )
    state = {
        "name": name,
        "bin_centers": centers,
        "constraint_singular_values": singular,
        "constraint_rank": rank,
        "null_projection": null,
        "scale": scale,
        "input_dimensions": int(x_train.shape[1]),
        "output_dimensions": int(transformed_train.shape[1]),
        "max_training_bin_feature_mean": max_bin_mean,
        "max_training_bag_feature_mean_difference": max_bag_difference,
    }
    return transformed_train, transformed_valid, state


def write_tsv(path, rows):
    if not rows:
        return
    with open(str(path), "w", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    started = time.time()
    train, valid = js.load_controlled_cache(BASE / "controlled_cache.npz")
    tracker_train, tracker_valid = js.load_tracker_cache(
        HERE / "tracker_multifit_cache.npz"
    )

    beta_state_base = rc.make_feature_state(train, "no_edep")
    beta_train_base = rc.apply_feature_state(train, beta_state_base)[:, 1:]
    beta_valid_base = rc.apply_feature_state(valid, beta_state_base)[:, 1:]
    rigidity_train_base, rigidity_valid_base, rigidity_state_base = (
        js.tracker_feature_matrices(train, valid, tracker_train, tracker_valid)
    )
    rigidity_train_base = rigidity_train_base[:, 1:]
    rigidity_valid_base = rigidity_valid_base[:, 1:]

    x_beta_train, x_beta_valid, beta_anchor_state = anchored_projection(
        beta_train_base, beta_valid_base, train, valid, "beta"
    )
    x_r_train, x_r_valid, rigidity_anchor_state = anchored_projection(
        rigidity_train_base, rigidity_valid_base, train, valid, "rigidity"
    )
    log(
        "anchored features beta=%d rigidity=%d beta_mean=%.3g beta_bag=%.3g "
        "R_mean=%.3g R_bag=%.3g"
        % (
            x_beta_train.shape[1],
            x_r_train.shape[1],
            beta_anchor_state["max_training_bin_feature_mean"],
            beta_anchor_state["max_training_bag_feature_mean_difference"],
            rigidity_anchor_state["max_training_bin_feature_mean"],
            rigidity_anchor_state["max_training_bag_feature_mean_difference"],
        )
    )

    frozen_beta = json.loads((BASE / "weak_frozen_models.json").read_text())
    raw_logsigma = np.asarray(
        frozen_beta["models"]["raw"]["logsigma"], dtype=np.float64
    )
    p_beta = x_beta_train.shape[1]
    p_r = x_r_train.shape[1]
    zero_beta = np.zeros(p_beta, dtype=np.float64)
    zero_r = np.zeros(p_r, dtype=np.float64)
    raw_parameters = np.r_[zero_beta, zero_r, raw_logsigma]

    small_rows = np.arange(min(400, len(train["u"])))
    small_data = {key: value[small_rows] for key, value in train.items()}
    check_density = js.JointDensity(
        small_data,
        x_beta_train[small_rows],
        x_r_train[small_rows],
        train["prior_c13"][small_rows],
    )
    check_parameters = raw_parameters.copy()
    if p_beta > 1:
        check_parameters[1] = 0.001
    if p_r > 1:
        check_parameters[p_beta + 1] = 0.002
    gradient_error = js.gradient_check(check_density, check_parameters)

    scenarios = [
        ("anchored_raw", False, False, 1.0, 1.0, train["prior_c13"]),
        ("anchored_beta_only", True, False, 1.0, 1.0, train["prior_c13"]),
        ("anchored_rigidity_only", False, True, 1.0, 1.0, train["prior_c13"]),
        ("anchored_joint_R_stiff", True, True, 1.0, 10.0, train["prior_c13"]),
        ("anchored_joint_balanced", True, True, 1.0, 1.0, train["prior_c13"]),
        ("anchored_joint_R_loose", True, True, 1.0, 0.1, train["prior_c13"]),
        (
            "anchored_joint_wrong_fraction",
            True,
            True,
            1.0,
            1.0,
            np.clip(train["prior_c13"] + 0.10, 0.01, 0.99),
        ),
    ]

    models = {}
    beta_parameters = None
    r_parameters = None
    for name, fit_beta, fit_r, beta_multiplier, r_multiplier, prior in scenarios:
        density = js.JointDensity(
            {key: train[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_train,
            x_r_train,
            prior,
            beta_multiplier,
            r_multiplier,
        )
        starts = []
        if fit_beta and beta_parameters is not None:
            starts.append(beta_parameters)
        if fit_r and r_parameters is not None:
            starts.append(r_parameters)
        model = js.fit_model(
            name,
            density,
            fit_beta,
            fit_r,
            raw_parameters,
            additional_starts=starts,
        )
        model["beta_multiplier"] = beta_multiplier
        model["rigidity_multiplier"] = r_multiplier
        model["training_prior_shift"] = (
            0.10 if name == "anchored_joint_wrong_fraction" else 0.0
        )
        models[name] = model
        if name == "anchored_beta_only":
            beta_parameters = model["parameters"].copy()
        if name == "anchored_rigidity_only":
            r_parameters = model["parameters"].copy()

    frozen = {
        "truth_policy": (
            "All anchored models were frozen before isotope labels, u_true, or "
            "R_true were consulted. Per-bin zero means and bag-difference null "
            "projections were computed from unlabeled reconstructed features."
        ),
        "gradient_check_max_relative_error": gradient_error,
        "beta_base_state": beta_state_base,
        "rigidity_base_state": rigidity_state_base,
        "beta_anchor_state": beta_anchor_state,
        "rigidity_anchor_state": rigidity_anchor_state,
        "models": models,
        "script_sha256": sha256_file(Path(__file__)),
    }
    (HERE / "anchored_joint_frozen_models.json").write_text(
        json.dumps(jsonable(frozen), indent=2, sort_keys=True)
    )
    log("anchored models frozen before truth evaluation")

    rows = []
    predictions = {
        "run": valid["run"],
        "event": valid["event"],
        "label": valid["label"],
        "bin": valid["bin"],
        "u_true": valid["u_true"],
        "R_true": valid["R_true"],
        "truth_valid": valid["truth_valid"],
    }
    for name, model in models.items():
        beta_multiplier = model["beta_multiplier"]
        r_multiplier = model["rigidity_multiplier"]
        train_prior = (
            np.clip(train["prior_c13"] + 0.10, 0.01, 0.99)
            if model["training_prior_shift"]
            else train["prior_c13"]
        )
        train_density = js.JointDensity(
            {key: train[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_train,
            x_r_train,
            train_prior,
            beta_multiplier,
            r_multiplier,
        )
        valid_density = js.JointDensity(
            {key: valid[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_valid,
            x_r_valid,
            valid["prior_c13"],
            beta_multiplier,
            r_multiplier,
        )
        row = js.evaluate(name, model, train_density, valid_density, valid)
        rows.append(row)
        values = valid_density.components(model["parameters"])
        predictions["u_" + name] = values["corrected_u"]
        predictions["R_" + name] = values["corrected_r"]
        predictions["du_correction_" + name] = values["beta_correction"]
        predictions["dlogR_correction_" + name] = values["log_r_correction"]

    write_tsv(HERE / "anchored_joint_metrics.tsv", rows)
    np.savez_compressed(str(HERE / "anchored_joint_predictions.npz"), **predictions)

    joint_names = [
        "anchored_joint_R_stiff",
        "anchored_joint_balanced",
        "anchored_joint_R_loose",
    ]
    pairwise = []
    with np.load(str(HERE / "anchored_joint_predictions.npz"), allow_pickle=False) as source:
        for first_index, first in enumerate(joint_names):
            for second in joint_names[first_index + 1 :]:
                pairwise.append(
                    {
                        "first": first,
                        "second": second,
                        "rms_delta_beta_correction": rc.rms(
                            source["du_correction_" + first]
                            - source["du_correction_" + second]
                        ),
                        "rms_delta_logR_correction": rc.rms(
                            source["dlogR_correction_" + first]
                            - source["dlogR_correction_" + second]
                        ),
                    }
                )
    write_tsv(HERE / "anchored_joint_pairwise.tsv", pairwise)

    summary = {
        "metrics": {row["method"]: row for row in rows},
        "pairwise_joint_solution_differences": pairwise,
        "elapsed_seconds": time.time() - started,
        "script_sha256": sha256_file(Path(__file__)),
        "limitation": (
            "C12 is B1308 and C13 is B1315. Anchoring protects the existing "
            "scales but does not by itself establish flight transfer."
        ),
    }
    (HERE / "anchored_joint_results.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )
    log("DONE elapsed=%.1fs" % (time.time() - started))


if __name__ == "__main__":
    main()
