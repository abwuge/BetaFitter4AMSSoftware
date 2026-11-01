#!/usr/bin/env python3
"""Weak C12/C13 mass-likelihood scan with separate beta and rigidity heads.

The weak objective receives only reconstructed measurements, tracker-internal
redundancy and known bag-level isotope fractions.  Event labels and MC truth
are excluded until every weak model is frozen.  Several relative penalties on
the beta and rigidity heads are scanned to expose or reject the beta/R
identifiability degeneracy of the mass spectrum.
"""

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "c_isotope_proportion_closure"
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(HERE))

import run_closure as rc
import run_tracker_multifit_probe as rt

NU = 6.0
T_NORMALIZATION = (
    math.lgamma((NU + 1.0) / 2.0)
    - math.lgamma(NU / 2.0)
    - 0.5 * math.log(NU * math.pi)
)
TRACKER_PCS = 24


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


def load_controlled_cache(path):
    with np.load(str(path), allow_pickle=False) as source:
        train, valid = {}, {}
        for name in source.files:
            prefix, key = name.split("__", 1)
            (train if prefix == "train" else valid)[key] = source[name]
    return train, valid


def load_tracker_cache(path):
    with np.load(str(path), allow_pickle=False) as source:
        train, valid = {}, {}
        for name in source.files:
            prefix, key = name.split("__", 1)
            (train if prefix == "train" else valid)[key] = source[name]
    return train, valid


def tracker_feature_matrices(train, valid, tracker_train, tracker_valid):
    candidates_train, candidate_names = rt.candidate_matrix(tracker_train)
    candidates_valid, candidate_names_valid = rt.candidate_matrix(tracker_valid)
    if candidate_names != candidate_names_valid:
        raise RuntimeError("Tracker candidate schema mismatch")
    raw_train, names = rt.raw_quality_features(
        train, tracker_train, candidates_train, candidate_names
    )
    raw_valid, names_valid = rt.raw_quality_features(
        valid, tracker_valid, candidates_valid, candidate_names
    )
    if names != names_valid:
        raise RuntimeError("Tracker raw feature schema mismatch")

    # Keep this head detector-specific: remove TOF-position/path geometry that
    # the tracker probe used only as an oracle diagnostic.
    keep = np.asarray([not name.startswith("geometry_") for name in names])
    raw_train = raw_train[:, keep]
    raw_valid = raw_valid[:, keep]
    kept_names = [name for name, selected in zip(names, keep) if selected]

    state = rt.fit_transform_state(raw_train)
    standardized_train = rt.apply_transform_state(raw_train, state)[:, 1:]
    standardized_valid = rt.apply_transform_state(raw_valid, state)[:, 1:]
    selected_names = [kept_names[int(index)] for index in state["indices"]]

    _, singular, vt = np.linalg.svd(standardized_train, full_matrices=False)
    rank = int(np.sum(singular > singular[0] * 1e-9))
    dimensions = min(TRACKER_PCS, rank)
    components = vt[:dimensions]
    scales = singular[:dimensions] / math.sqrt(len(standardized_train))
    train_pc = standardized_train.dot(components.T) / scales
    valid_pc = standardized_valid.dot(components.T) / scales
    return (
        np.column_stack([np.ones(len(train_pc)), train_pc]),
        np.column_stack([np.ones(len(valid_pc)), valid_pc]),
        {
            "raw_state": state,
            "raw_feature_names": selected_names,
            "pca_components": components,
            "pca_scales": scales,
            "rank": rank,
            "dimensions": dimensions,
            "explained_variance_fraction": float(
                np.sum(singular[:dimensions] ** 2) / np.sum(singular ** 2)
            ),
        },
    )


class JointDensity(object):
    def __init__(
        self,
        data,
        x_beta,
        x_rigidity,
        prior_c13,
        beta_multiplier=1.0,
        rigidity_multiplier=1.0,
    ):
        self.u = np.asarray(data["u"], dtype=np.float64)
        self.r = np.asarray(data["R"], dtype=np.float64)
        self.rel_r = np.asarray(data["relR"], dtype=np.float64)
        self.bin = np.asarray(data["bin"], dtype=np.int32)
        self.x_beta = np.asarray(x_beta, dtype=np.float64)
        self.x_rigidity = np.asarray(x_rigidity, dtype=np.float64)
        self.prior = np.asarray(prior_c13, dtype=np.float64)
        self.p_beta = self.x_beta.shape[1]
        self.p_rigidity = self.x_rigidity.shape[1]
        self.beta_multiplier = float(beta_multiplier)
        self.rigidity_multiplier = float(rigidity_multiplier)

        raw_means = np.sqrt(
            1.0 + (rc.MASSES[None, :] / (rc.CHARGE * self.r[:, None])) ** 2
        )
        raw_sensitivity = (raw_means * raw_means - 1.0) / raw_means
        self.fixed_tracker_variance = (
            ((NU - 2.0) / NU)
            * raw_sensitivity
            * raw_sensitivity
            * self.rel_r[:, None] ** 2
        )

    def unpack(self, parameters):
        first = self.p_beta
        second = first + self.p_rigidity
        return (
            parameters[:first],
            parameters[first:second],
            parameters[second : second + rc.N_BINS],
        )

    def components(self, parameters):
        theta_beta, theta_rigidity, logsigma = self.unpack(parameters)
        beta_correction = self.x_beta.dot(theta_beta)
        log_r_correction = self.x_rigidity.dot(theta_rigidity)
        corrected_u = self.u - beta_correction
        corrected_r = self.r * np.exp(-log_r_correction)
        means = np.sqrt(
            1.0
            + (rc.MASSES[None, :] / (rc.CHARGE * corrected_r[:, None])) ** 2
        )
        sensitivity = (means * means - 1.0) / means
        residual = corrected_u[:, None] - means
        sigma2 = np.exp(2.0 * logsigma)[self.bin]
        variance = sigma2[:, None] + self.fixed_tracker_variance
        logpdf = (
            T_NORMALIZATION
            - 0.5 * np.log(variance)
            - ((NU + 1.0) / 2.0)
            * np.log1p(residual * residual / (NU * variance))
        )
        return {
            "logpdf": logpdf,
            "residual": residual,
            "variance": variance,
            "sigma2": sigma2,
            "sensitivity": sensitivity,
            "beta_correction": beta_correction,
            "log_r_correction": log_r_correction,
            "corrected_u": corrected_u,
            "corrected_r": corrected_r,
        }

    def value_gradient(self, parameters):
        values = self.components(parameters)
        logpdf = values["logpdf"]
        residual = values["residual"]
        variance = values["variance"]
        sigma2 = values["sigma2"]
        sensitivity = values["sensitivity"]
        with np.errstate(divide="ignore"):
            log_mixture = np.logaddexp(
                np.log1p(-self.prior) + logpdf[:, 0],
                np.log(self.prior) + logpdf[:, 1],
            )
        posterior = np.exp(logpdf - log_mixture[:, None]) * np.column_stack(
            [1.0 - self.prior, self.prior]
        )
        denominator = NU * variance + residual * residual
        component_gradient = (
            -posterior * (NU + 1.0) * residual / denominator
        )
        event_beta_gradient = np.sum(component_gradient, axis=1)
        event_rigidity_gradient = np.sum(
            component_gradient * sensitivity, axis=1
        )
        sigma_gradient = np.sum(
            posterior
            * (sigma2[:, None] / variance)
            * (1.0 - (NU + 1.0) * residual * residual / denominator),
            axis=1,
        )
        n = float(len(self.u))
        gradient_beta = self.x_beta.T.dot(event_beta_gradient) / n
        gradient_rigidity = self.x_rigidity.T.dot(event_rigidity_gradient) / n
        gradient_sigma = np.bincount(
            self.bin, weights=sigma_gradient, minlength=rc.N_BINS
        ) / n
        value = -float(np.mean(log_mixture))

        theta_beta, theta_rigidity, _ = self.unpack(parameters)
        if self.p_beta > 1:
            coefficient = theta_beta[1:] / 0.02
            value += (
                0.5
                * self.beta_multiplier
                * 5e-4
                * float(np.sum(coefficient * coefficient))
            )
            gradient_beta[1:] += (
                self.beta_multiplier
                * 5e-4
                * theta_beta[1:]
                / (0.02 ** 2)
            )
        if self.p_rigidity > 1:
            coefficient = theta_rigidity[1:] / 0.08
            value += (
                0.5
                * self.rigidity_multiplier
                * 5e-4
                * float(np.sum(coefficient * coefficient))
            )
            gradient_rigidity[1:] += (
                self.rigidity_multiplier
                * 5e-4
                * theta_rigidity[1:]
                / (0.08 ** 2)
            )

        beta_prediction = values["beta_correction"]
        rigidity_prediction = values["log_r_correction"]
        value += (
            0.5
            * self.beta_multiplier
            * 2e-3
            * float(np.mean((beta_prediction / 0.03) ** 2))
        )
        value += (
            0.5
            * self.rigidity_multiplier
            * 2e-3
            * float(np.mean((rigidity_prediction / 0.15) ** 2))
        )
        gradient_beta += (
            self.beta_multiplier
            * 2e-3
            / (0.03 ** 2)
            * self.x_beta.T.dot(beta_prediction)
            / n
        )
        gradient_rigidity += (
            self.rigidity_multiplier
            * 2e-3
            / (0.15 ** 2)
            * self.x_rigidity.T.dot(rigidity_prediction)
            / n
        )
        return value, np.concatenate(
            [gradient_beta, gradient_rigidity, gradient_sigma]
        )

    def event_nll(self, parameters):
        logpdf = self.components(parameters)["logpdf"]
        with np.errstate(divide="ignore"):
            return -np.logaddexp(
                np.log1p(-self.prior) + logpdf[:, 0],
                np.log(self.prior) + logpdf[:, 1],
            )


def gradient_check(density, parameters):
    value, analytic = density.value_gradient(parameters)
    indices = [
        0,
        min(3, density.p_beta - 1),
        density.p_beta,
        density.p_beta + min(3, density.p_rigidity - 1),
        density.p_beta + density.p_rigidity,
    ]
    errors = []
    for index in sorted(set(indices)):
        plus = parameters.copy()
        minus = parameters.copy()
        step = 1e-6
        plus[index] += step
        minus[index] -= step
        numeric = (
            density.value_gradient(plus)[0]
            - density.value_gradient(minus)[0]
        ) / (2.0 * step)
        errors.append(abs(analytic[index] - numeric) / (1.0 + abs(numeric)))
    maximum = float(max(errors))
    log("gradient check value=%.8g max_relative_error=%.3g" % (value, maximum))
    if maximum > 5e-5:
        raise RuntimeError("Joint gradient check failed")
    return maximum


def fit_model(
    name,
    density,
    fit_beta,
    fit_rigidity,
    initial_parameters,
    additional_starts=None,
):
    starts = [np.asarray(initial_parameters, dtype=float).copy()]
    for start in additional_starts or []:
        starts.append(np.asarray(start, dtype=float).copy())
    bounds = (
        [(-0.05, 0.05)] * density.p_beta
        + [(-0.25, 0.25)] * density.p_rigidity
        + [(math.log(0.0015), math.log(0.08))] * rc.N_BINS
    )
    if not fit_beta:
        for index in range(density.p_beta):
            bounds[index] = (0.0, 0.0)
    if not fit_rigidity:
        for index in range(
            density.p_beta, density.p_beta + density.p_rigidity
        ):
            bounds[index] = (0.0, 0.0)

    records = []
    solutions = []
    for start_index, initial in enumerate(starts):
        result = minimize(
            density.value_gradient,
            initial,
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": 1800,
                "maxls": 80,
                "ftol": 1e-12,
                "gtol": 2e-6,
            },
        )
        gradient_norm = float(np.max(np.abs(result.jac)))
        record = {
            "start": start_index,
            "success": bool(result.success),
            "status": int(result.status),
            "objective": float(result.fun),
            "iterations": int(result.nit),
            "gradient_inf_norm": gradient_norm,
            "message": str(result.message),
        }
        records.append(record)
        log("fit %s %s" % (name, record))
        if np.isfinite(result.fun) and np.all(np.isfinite(result.x)):
            solutions.append(result)
    if not solutions:
        raise RuntimeError("All starts failed for " + name)
    best = min(solutions, key=lambda item: item.fun)
    values = density.components(best.x)
    return {
        "name": name,
        "parameters": best.x,
        "theta_beta": density.unpack(best.x)[0],
        "theta_rigidity": density.unpack(best.x)[1],
        "logsigma": density.unpack(best.x)[2],
        "objective": float(best.fun),
        "optimizer": records,
        "beta_correction_rms_training": rc.rms(values["beta_correction"]),
        "log_r_correction_rms_training": rc.rms(values["log_r_correction"]),
    }


def roc_rejection(score, label, signal_efficiency=0.50):
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=np.int8)
    finite = np.isfinite(score)
    signal = score[finite & (label == 1)]
    background = score[finite & (label == 0)]
    if not len(signal) or not len(background):
        return float("nan")
    threshold = np.quantile(signal, 1.0 - signal_efficiency)
    acceptance = float(np.mean(background >= threshold))
    return float(1.0 / acceptance) if acceptance > 0.0 else float("inf")


def evaluate(name, model, train_density, valid_density, valid):
    train_values = train_density.components(model["parameters"])
    values = valid_density.components(model["parameters"])
    corrected_u = values["corrected_u"]
    corrected_r = values["corrected_r"]
    truth = (
        valid["truth_valid"]
        & np.isfinite(valid["u_true"])
        & np.isfinite(valid["R_true"])
        & (valid["R_true"] > 0.0)
        & np.isfinite(corrected_u)
        & (corrected_u > 0.0)
        & np.isfinite(corrected_r)
        & (corrected_r > 0.0)
    )
    du = corrected_u[truth] - valid["u_true"][truth]
    dlogr = np.log(corrected_r[truth] / valid["R_true"][truth])
    mass2 = (rc.CHARGE * corrected_r) ** 2 * (corrected_u * corrected_u - 1.0)
    true_mass2 = rc.MASSES[valid["label"]] ** 2
    light = truth & (valid["label"] == 0)
    heavy = truth & (valid["label"] == 1)
    denominator = math.sqrt(rc.w68(mass2[light]) ** 2 + rc.w68(mass2[heavy]) ** 2)
    return {
        "method": name,
        "n": int(np.sum(truth)),
        "training_objective": float(model["objective"]),
        "heldout_nll": float(np.mean(valid_density.event_nll(model["parameters"]))),
        "median_du": float(np.median(du)),
        "w68_du": rc.w68(du),
        "rms_du": rc.rms(du),
        "median_dlogR": float(np.median(dlogr)),
        "w68_dlogR": rc.w68(dlogr),
        "rms_dlogR": rc.rms(dlogr),
        "beta_correction_median": float(np.median(values["beta_correction"][truth])),
        "beta_correction_rms": rc.rms(values["beta_correction"][truth]),
        "logR_correction_median": float(np.median(values["log_r_correction"][truth])),
        "logR_correction_rms": rc.rms(values["log_r_correction"][truth]),
        "correction_correlation": float(
            np.corrcoef(
                values["beta_correction"][truth],
                values["log_r_correction"][truth],
            )[0, 1]
        ),
        "mass2_auc": rc.auc(mass2[truth], valid["label"][truth]),
        "mass2_separation": (
            (float(np.median(mass2[heavy])) - float(np.median(mass2[light])))
            / denominator
            if denominator > 0.0
            else float("nan")
        ),
        "c12_rejection_at_50pct_c13_eff": roc_rejection(
            mass2[truth], valid["label"][truth], 0.50
        ),
        "median_mass2_error": float(
            np.median(mass2[truth] - true_mass2[truth])
        ),
        "w68_mass2_error": rc.w68(mass2[truth] - true_mass2[truth]),
        "rms_mass2_error": rc.rms(mass2[truth] - true_mass2[truth]),
        "train_beta_correction_rms": rc.rms(train_values["beta_correction"]),
        "train_logR_correction_rms": rc.rms(train_values["log_r_correction"]),
    }


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
    train, valid = load_controlled_cache(BASE / "controlled_cache.npz")
    tracker_train, tracker_valid = load_tracker_cache(HERE / "tracker_multifit_cache.npz")

    beta_state = rc.make_feature_state(train, "no_edep")
    x_beta_train = rc.apply_feature_state(train, beta_state)
    x_beta_valid = rc.apply_feature_state(valid, beta_state)
    x_r_train, x_r_valid, tracker_state = tracker_feature_matrices(
        train, valid, tracker_train, tracker_valid
    )
    log(
        "feature dimensions beta=%d rigidity=%d rigidity_variance=%.5f"
        % (
            x_beta_train.shape[1],
            x_r_train.shape[1],
            tracker_state["explained_variance_fraction"],
        )
    )

    frozen_beta = json.loads((BASE / "weak_frozen_models.json").read_text())
    beta_model = frozen_beta["models"]["weak_no_edep"]
    raw_model = frozen_beta["models"]["raw"]
    beta_theta = np.asarray(beta_model["theta"], dtype=float)
    raw_logsigma = np.asarray(raw_model["logsigma"], dtype=float)
    beta_logsigma = np.asarray(beta_model["logsigma"], dtype=float)
    if len(beta_theta) != x_beta_train.shape[1]:
        raise RuntimeError("Frozen beta feature dimension mismatch")

    p_beta = x_beta_train.shape[1]
    p_r = x_r_train.shape[1]
    zero_r = np.zeros(p_r, dtype=float)
    zero_beta = np.zeros(p_beta, dtype=float)
    raw_parameters = np.r_[zero_beta, zero_r, raw_logsigma]
    beta_parameters = np.r_[beta_theta, zero_r, beta_logsigma]

    small_rows = np.arange(min(400, len(train["u"])))
    small_data = {key: value[small_rows] for key, value in train.items()}
    check_density = JointDensity(
        small_data,
        x_beta_train[small_rows],
        x_r_train[small_rows],
        train["prior_c13"][small_rows],
    )
    check_parameters = raw_parameters.copy()
    check_parameters[0] = 0.001
    check_parameters[p_beta] = 0.002
    gradient_error = gradient_check(check_density, check_parameters)

    scenarios = [
        ("raw", False, False, 1.0, 1.0, train["prior_c13"]),
        ("beta_only_frozen", True, False, 1.0, 1.0, train["prior_c13"]),
        ("rigidity_only", False, True, 1.0, 1.0, train["prior_c13"]),
        ("joint_R_stiff", True, True, 1.0, 10.0, train["prior_c13"]),
        ("joint_balanced", True, True, 1.0, 1.0, train["prior_c13"]),
        ("joint_R_loose", True, True, 1.0, 0.1, train["prior_c13"]),
        ("joint_beta_loose", True, True, 0.1, 1.0, train["prior_c13"]),
        (
            "joint_balanced_wrong_fraction",
            True,
            True,
            1.0,
            1.0,
            np.clip(train["prior_c13"] + 0.10, 0.01, 0.99),
        ),
    ]

    models = {}
    r_only_parameters = None
    for name, fit_beta, fit_r, beta_multiplier, r_multiplier, prior in scenarios:
        density = JointDensity(
            {key: train[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_train,
            x_r_train,
            prior,
            beta_multiplier,
            r_multiplier,
        )
        if name == "raw":
            initial = raw_parameters
            model = fit_model(name, density, False, False, initial)
        elif name == "beta_only_frozen":
            # Keep the pre-existing frozen weak-beta solution unchanged.  Its
            # objective is recomputed under this joint implementation.
            values = density.components(beta_parameters)
            model = {
                "name": name,
                "parameters": beta_parameters,
                "theta_beta": beta_theta,
                "theta_rigidity": zero_r,
                "logsigma": beta_logsigma,
                "objective": float(density.value_gradient(beta_parameters)[0]),
                "optimizer": [{"source": "pre-existing frozen weak_no_edep"}],
                "beta_correction_rms_training": rc.rms(values["beta_correction"]),
                "log_r_correction_rms_training": 0.0,
            }
        elif name == "rigidity_only":
            model = fit_model(name, density, False, True, raw_parameters)
            r_only_parameters = model["parameters"].copy()
        else:
            starts = []
            if r_only_parameters is not None:
                starts.append(r_only_parameters)
            model = fit_model(
                name,
                density,
                fit_beta,
                fit_r,
                beta_parameters,
                additional_starts=starts,
            )
        model["beta_multiplier"] = beta_multiplier
        model["rigidity_multiplier"] = r_multiplier
        model["training_prior_shift"] = (
            0.10 if name == "joint_balanced_wrong_fraction" else 0.0
        )
        models[name] = model

    frozen = {
        "truth_policy": (
            "All models in this file were fit before event isotope labels, u_true, "
            "or R_true were consulted. Labels only constructed controlled bags with "
            "known fractions. Truth is used after this file is written."
        ),
        "gradient_check_max_relative_error": gradient_error,
        "beta_feature_state": beta_state,
        "tracker_feature_state": tracker_state,
        "models": models,
        "script_sha256": sha256_file(Path(__file__)),
        "inputs": {
            "controlled_cache_sha256": sha256_file(BASE / "controlled_cache.npz"),
            "tracker_multifit_cache_sha256": sha256_file(
                HERE / "tracker_multifit_cache.npz"
            ),
            "weak_frozen_beta_sha256": sha256_file(
                BASE / "weak_frozen_models.json"
            ),
        },
    }
    (HERE / "joint_weak_frozen_models.json").write_text(
        json.dumps(jsonable(frozen), indent=2, sort_keys=True)
    )
    log("joint weak models frozen before truth evaluation")

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
        train_density = JointDensity(
            {key: train[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_train,
            x_r_train,
            train_prior,
            beta_multiplier,
            r_multiplier,
        )
        valid_density = JointDensity(
            {key: valid[key] for key in ["u", "R", "relR", "bin"]},
            x_beta_valid,
            x_r_valid,
            valid["prior_c13"],
            beta_multiplier,
            r_multiplier,
        )
        row = evaluate(name, model, train_density, valid_density, valid)
        rows.append(row)
        values = valid_density.components(model["parameters"])
        predictions["u_" + name] = values["corrected_u"]
        predictions["R_" + name] = values["corrected_r"]
        predictions["du_correction_" + name] = values["beta_correction"]
        predictions["dlogR_correction_" + name] = values["log_r_correction"]

    write_tsv(HERE / "joint_identifiability_metrics.tsv", rows)
    np.savez_compressed(
        str(HERE / "joint_identifiability_predictions.npz"), **predictions
    )

    joint_names = [
        "joint_R_stiff",
        "joint_balanced",
        "joint_R_loose",
        "joint_beta_loose",
    ]
    pairwise = []
    with np.load(
        str(HERE / "joint_identifiability_predictions.npz"), allow_pickle=False
    ) as source:
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
    write_tsv(HERE / "joint_solution_pairwise.tsv", pairwise)

    summary = {
        "metrics": {row["method"]: row for row in rows},
        "pairwise_joint_solution_differences": pairwise,
        "elapsed_seconds": time.time() - started,
        "script_sha256": sha256_file(Path(__file__)),
        "limitation": (
            "C12 is B1308 and C13 is B1315. This is an identifiability closure, "
            "not evidence that a model transfers to flight or a same-production pair."
        ),
    }
    (HERE / "joint_identifiability_results.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )
    log("DONE elapsed=%.1fs" % (time.time() - started))


if __name__ == "__main__":
    main()
