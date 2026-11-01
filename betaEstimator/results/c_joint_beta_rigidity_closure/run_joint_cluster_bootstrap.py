#!/usr/bin/env python3
"""Paired source-file bootstrap for the C12/C13 joint closure."""

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "c_isotope_proportion_closure"
sys.path.insert(0, str(BASE))
import run_closure as rc

N_BOOTSTRAP = 3000
SEED = 260929


def load_cache():
    with np.load(str(BASE / "controlled_cache.npz"), allow_pickle=False) as source:
        valid = {
            name.split("__", 1)[1]: source[name]
            for name in source.files
            if name.startswith("valid__")
        }
    return valid


def rejection(score, label):
    signal = score[label == 1]
    background = score[label == 0]
    threshold = np.quantile(signal, 0.5)
    acceptance = float(np.mean(background >= threshold))
    return 1.0 / acceptance if acceptance > 0 else float("inf")


def metrics(u, rigidity, data, rows):
    label = data["label"][rows]
    u = u[rows]
    rigidity = rigidity[rows]
    u_true = data["u_true"][rows]
    r_true = data["R_true"][rows]
    mass2 = (rc.CHARGE * rigidity) ** 2 * (u * u - 1.0)
    true_mass2 = rc.MASSES[label] ** 2
    return {
        "auc": rc.auc(mass2, label),
        "rejection50": rejection(mass2, label),
        "w68_mass2_error": rc.w68(mass2 - true_mass2),
        "w68_du": rc.w68(u - u_true),
        "w68_dlogR": rc.w68(np.log(rigidity / r_true)),
        "median_du": float(np.median(u - u_true)),
        "median_dlogR": float(np.median(np.log(rigidity / r_true))),
    }


def summarize(values):
    values = np.asarray(values, dtype=float)
    return {
        "median": float(np.median(values)),
        "p2p5": float(np.percentile(values, 2.5)),
        "p97p5": float(np.percentile(values, 97.5)),
        "probability_positive": float(np.mean(values > 0.0)),
    }


def main():
    data = load_cache()
    truth = (
        data["truth_valid"]
        & np.isfinite(data["u_true"])
        & np.isfinite(data["R_true"])
        & (data["R_true"] > 0)
    )

    with np.load(str(BASE / "predictions.npz"), allow_pickle=False) as source:
        methods = {
            "raw": (source["u_raw"].copy(), data["R"].copy()),
            "beta_only_frozen": (
                source["u_weak_no_edep"].copy(),
                data["R"].copy(),
            ),
        }
    with np.load(
        str(HERE / "joint_identifiability_predictions.npz"), allow_pickle=False
    ) as source:
        methods["joint_unanchored"] = (
            source["u_joint_balanced"].copy(),
            source["R_joint_balanced"].copy(),
        )
    with np.load(str(HERE / "anchored_joint_predictions.npz"), allow_pickle=False) as source:
        methods["joint_anchored"] = (
            source["u_anchored_joint_balanced"].copy(),
            source["R_anchored_joint_balanced"].copy(),
        )

    base_rows = np.flatnonzero(truth)
    observed = {
        name: metrics(u, rigidity, data, base_rows)
        for name, (u, rigidity) in methods.items()
    }

    clusters = {}
    for label in (0, 1):
        source_ids = sorted(
            int(item) for item in np.unique(data["source"][truth & (data["label"] == label)])
        )
        clusters[label] = [
            np.flatnonzero(
                truth & (data["label"] == label) & (data["source"] == source_id)
            )
            for source_id in source_ids
        ]

    rng = np.random.RandomState(SEED)
    samples = {
        name: {metric: [] for metric in next(iter(observed.values()))}
        for name in methods
    }
    differences = {
        comparison: {metric: [] for metric in next(iter(observed.values()))}
        for comparison in [
            "joint_anchored-minus-beta_only_frozen",
            "joint_unanchored-minus-beta_only_frozen",
            "joint_anchored-minus-raw",
            "beta_only_frozen-minus-raw",
        ]
    }

    for _ in range(N_BOOTSTRAP):
        pieces = []
        for label in (0, 1):
            choices = rng.randint(0, len(clusters[label]), size=len(clusters[label]))
            pieces.extend(clusters[label][choice] for choice in choices)
        rows = np.concatenate(pieces)
        current = {}
        for name, (u, rigidity) in methods.items():
            current[name] = metrics(u, rigidity, data, rows)
            for metric, value in current[name].items():
                samples[name][metric].append(value)

        comparisons = {
            "joint_anchored-minus-beta_only_frozen": (
                "joint_anchored",
                "beta_only_frozen",
            ),
            "joint_unanchored-minus-beta_only_frozen": (
                "joint_unanchored",
                "beta_only_frozen",
            ),
            "joint_anchored-minus-raw": ("joint_anchored", "raw"),
            "beta_only_frozen-minus-raw": ("beta_only_frozen", "raw"),
        }
        for comparison, (first, second) in comparisons.items():
            for metric in current[first]:
                # For width/error metrics, positive means improvement by using
                # second-first.  For AUC/rejection, positive means first-second.
                if metric in ("auc", "rejection50"):
                    delta = current[first][metric] - current[second][metric]
                elif metric in ("w68_mass2_error", "w68_du", "w68_dlogR"):
                    delta = current[second][metric] - current[first][metric]
                else:
                    delta = abs(current[second][metric]) - abs(current[first][metric])
                differences[comparison][metric].append(delta)

    result = {
        "bootstrap_note": (
            "Paired source-file cluster bootstrap, separately resampling 32 C12 and "
            "8 C13 validation files. Models are fixed; training and production-version "
            "uncertainty are not included. Positive comparison deltas mean the first "
            "method is better under the metric convention documented in the script."
        ),
        "replicates": N_BOOTSTRAP,
        "observed": observed,
        "method_intervals": {
            name: {metric: summarize(values) for metric, values in metric_values.items()}
            for name, metric_values in samples.items()
        },
        "paired_differences": {
            comparison: {
                metric: summarize(values) for metric, values in metric_values.items()
            }
            for comparison, metric_values in differences.items()
        },
    }

    per_bin = {}
    for b in range(rc.N_BINS):
        rows = np.flatnonzero(truth & (data["bin"] == b))
        per_bin[str(b)] = {
            name: metrics(u, rigidity, data, rows)
            for name, (u, rigidity) in methods.items()
        }
    result["per_reconstructed_rigidity_bin"] = per_bin

    (HERE / "joint_cluster_bootstrap_results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True)
    )
    print(json.dumps(result["paired_differences"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
