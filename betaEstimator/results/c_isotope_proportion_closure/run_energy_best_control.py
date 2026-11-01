#!/usr/bin/env python3
"""Long convergence and wrong-ratio controls for the best residualized-energy model."""
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import run_closure as rc
import run_energy_residualized as er


def optimize(name, data, x, prior, starts):
    density = rc.MixtureDensity(data, x, prior)
    p = x.shape[1]
    bounds = [(-0.05, 0.05)] * p + [
        (math.log(0.0015), math.log(0.08))
    ] * rc.N_BINS
    records = []
    solutions = []
    for index, initial in enumerate(starts):
        result = minimize(
            density.value_gradient,
            np.asarray(initial, dtype=float),
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": 5000,
                "maxls": 80,
                "ftol": 1e-13,
                "gtol": 5e-7,
            },
        )
        gradient = density.value_gradient(result.x)[1]
        record = {
            "start": index,
            "success": bool(result.success),
            "status": int(result.status),
            "objective": float(result.fun),
            "iterations": int(result.nit),
            "message": str(result.message),
            "gradient_inf_norm": float(np.max(np.abs(gradient))),
        }
        records.append(record)
        rc.log("%s %s" % (name, record))
        if np.isfinite(result.fun) and np.all(np.isfinite(result.x)):
            solutions.append(result)
    if not solutions:
        raise RuntimeError("No finite fit for " + name)
    best = min(solutions, key=lambda result: result.fun)
    return {
        "name": name,
        "parameters": best.x,
        "theta": best.x[:p],
        "logsigma": best.x[p:],
        "objective": float(best.fun),
        "optimizer": records,
    }


def main():
    destination = Path(__file__).resolve().parent
    started = time.time()
    train, valid = er.load_cache(destination / "controlled_cache.npz")
    base_state = rc.make_feature_state(train, "no_edep")
    base_train = rc.apply_feature_state(train, base_state)
    base_valid = rc.apply_feature_state(valid, base_state)
    response = er.response_state(train)
    x_train, x_valid = er.append_energy(
        base_train,
        base_valid,
        er.energy_raw(train, "fitted", response),
        er.energy_raw(valid, "fitted", response),
    )
    weak_train = {key: train[key] for key in ["u", "R", "relR", "bin"]}
    weak_valid = {key: valid[key] for key in ["u", "R", "relR", "bin"]}

    earlier = json.load(open(str(destination / "energy_residualized_results.json")))
    initial = np.asarray(
        earlier["weak_models"]["fitted"]["parameters"], dtype=float
    )
    rng = np.random.RandomState(4317)
    perturb = initial.copy()
    perturb[: x_train.shape[1]] += rng.normal(
        0.0, 2e-4, size=x_train.shape[1]
    )
    nominal = optimize(
        "weak_energy_fitted_converged",
        weak_train,
        x_train,
        train["prior_c13"],
        [initial, perturb],
    )
    plus_prior = np.clip(train["prior_c13"] + 0.10, 0.01, 0.99)
    minus_prior = np.clip(train["prior_c13"] - 0.10, 0.01, 0.99)
    wrong_plus = optimize(
        "weak_energy_fitted_wrong_plus_0p10",
        weak_train,
        x_train,
        plus_prior,
        [nominal["parameters"]],
    )
    wrong_minus = optimize(
        "weak_energy_fitted_wrong_minus_0p10",
        weak_train,
        x_train,
        minus_prior,
        [nominal["parameters"]],
    )

    models = {
        "nominal": nominal,
        "wrong_plus": wrong_plus,
        "wrong_minus": wrong_minus,
    }
    corrected = {"raw": valid["u"].copy()}
    heldout = {"raw": float("nan")}
    nominal_density = rc.MixtureDensity(
        weak_valid, x_valid, valid["prior_c13"]
    )
    diagnostics = {}
    for key, model in models.items():
        name = model["name"]
        correction = x_valid.dot(model["theta"])
        corrected[name] = valid["u"] - correction
        heldout[name] = float(
            np.mean(nominal_density.event_nll(model["parameters"]))
        )
        diagnostics[name] = {
            "correction_rms": rc.rms(correction),
            "correction_label_auc": rc.auc(correction, valid["label"]),
            "correction_median_c12": float(
                np.median(correction[valid["label"] == 0])
            ),
            "correction_median_c13": float(
                np.median(correction[valid["label"] == 1])
            ),
        }

    rows = []
    bin_rows = []
    for name, values in corrected.items():
        method_rows, mass2 = rc.method_metrics(
            name, values, valid, heldout[name]
        )
        rows.extend(method_rows)
        bin_rows.extend(rc.bin_metrics(name, values, mass2, valid))
    raw_width = {
        row["group"]: row["w68_du"]
        for row in rows
        if row["method"] == "raw"
    }
    for row in rows:
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_du"] / raw_width[row["group"]]
        )
    rc.write_tsv(destination / "energy_best_control_metrics.tsv", rows)
    rc.write_tsv(destination / "energy_best_control_bin_metrics.tsv", bin_rows)
    overall = {
        row["method"]: row for row in rows if row["group"] == "all"
    }
    output = {
        "models": models,
        "diagnostics": diagnostics,
        "overall_metrics": overall,
        "elapsed_seconds": time.time() - started,
        "script_sha256": rc.sha256_file(Path(__file__)),
    }
    (destination / "energy_best_control_results.json").write_text(
        json.dumps(rc.jsonable(output), indent=2, sort_keys=True)
    )
    rc.log(
        "DONE best-control elapsed=%.1fs results=%s"
        % (
            time.time() - started,
            destination / "energy_best_control_results.json",
        )
    )


if __name__ == "__main__":
    main()
