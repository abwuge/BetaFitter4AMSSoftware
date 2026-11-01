#!/usr/bin/env python3
"""One-shot confirmation on source files unused by all model development.

The selected weak model, detector response residualization, feature scaling and
coefficients are reconstructed entirely from the original training cache and
frozen JSON files. Blind corrections are written before isotope labels or MC
truth are used for metrics.
"""
import csv
import glob
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

import run_closure as rc
import run_energy_residualized as er

N_C12_FILES = 32
N_C13_FILES = 8
BLIND_TARGET_PER_BIN = 1600
BLIND_SEED = 20260929
BOOTSTRAP_REPLICATES = 5000


def file_sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value):
    value = np.asarray(value)
    digest = hashlib.sha256()
    digest.update(str(value.dtype).encode("ascii"))
    digest.update(str(value.shape).encode("ascii"))
    digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def choose_untouched(pattern, used, count):
    candidates = [path for path in sorted(glob.glob(pattern)) if path not in used]
    # The development scans only touched the lexicographic beginning. Taking the
    # tail additionally avoids those exploratory files, not merely the manifest.
    selected = candidates[-count:]
    if len(selected) != count:
        raise RuntimeError("Not enough untouched source files")
    if set(selected).intersection(used):
        raise RuntimeError("Blind file overlaps development manifest")
    return selected


def width68(values):
    q = np.percentile(values, [16.0, 84.0])
    return float((q[1] - q[0]) / 2.0)


def bootstrap(valid, corrected, replicates, seed):
    rng = np.random.RandomState(seed)
    truth = valid["truth_valid"]
    strata = [
        np.flatnonzero(truth & (valid["label"] == label) & (valid["bin"] == b))
        for label in (0, 1)
        for b in range(rc.N_BINS)
    ]
    if any(len(rows) == 0 for rows in strata):
        raise RuntimeError("Empty blind bootstrap stratum")
    names = list(corrected)
    values = {
        group: {
            "no_edep_vs_raw_percent": [],
            "best_vs_raw_percent": [],
            "best_vs_no_edep_percent": [],
        }
        for group in ("all", "C12", "C13")
    }
    for _ in range(replicates):
        sample = np.concatenate(
            [rows[rng.randint(0, len(rows), len(rows))] for rows in strata]
        )
        for group in ("all", "C12", "C13"):
            select = sample
            if group == "C12":
                select = sample[valid["label"][sample] == 0]
            elif group == "C13":
                select = sample[valid["label"][sample] == 1]
            widths = {
                name: width68(corrected[name][select] - valid["u_true"][select])
                for name in names
            }
            values[group]["no_edep_vs_raw_percent"].append(
                100.0 * (1.0 - widths["no_edep"] / widths["raw"])
            )
            values[group]["best_vs_raw_percent"].append(
                100.0 * (1.0 - widths["best"] / widths["raw"])
            )
            values[group]["best_vs_no_edep_percent"].append(
                100.0 * (1.0 - widths["best"] / widths["no_edep"])
            )
    summary = {}
    for group, metrics in values.items():
        summary[group] = {}
        for metric, sequence in metrics.items():
            sequence = np.asarray(sequence)
            summary[group][metric] = {
                "median": float(np.median(sequence)),
                "p2p5": float(np.percentile(sequence, 2.5)),
                "p97p5": float(np.percentile(sequence, 97.5)),
                "probability_positive": float(np.mean(sequence > 0.0)),
            }
    return summary


def main():
    destination = Path(__file__).resolve().parent
    started = time.time()

    with open(str(destination / "manifest.tsv")) as source:
        development_manifest = list(csv.DictReader(source, delimiter="\t"))
    used = {row["path"] for row in development_manifest}
    c12_files = choose_untouched(rc.C12_PATTERN, used, N_C12_FILES)
    c13_files = choose_untouched(rc.C13_PATTERN, used, N_C13_FILES)
    blind_manifest = [
        {"isotope": isotope, "path": path, "selection": "tail_of_unused_sorted_paths"}
        for isotope, paths in (("C12", c12_files), ("C13", c13_files))
        for path in paths
    ]
    rc.write_tsv(destination / "blind_manifest.tsv", blind_manifest)

    audit = []
    pools = {
        0: rc.extract_pool(c12_files, 0, 67, rc.MASSES[0], 2, audit),
        1: rc.extract_pool(c13_files, 1, 117, rc.MASSES[1], 2, audit),
    }
    rc.write_tsv(destination / "blind_extraction_audit.tsv", audit)
    blind, mixture_audit = rc.controlled_mix(
        pools,
        rc.VALID_FRACTIONS,
        BLIND_TARGET_PER_BIN,
        BLIND_SEED,
    )
    rc.write_tsv(destination / "blind_mixture_audit.tsv", mixture_audit)

    development_ids = {
        (int(run), int(event), int(label))
        for run, event, label in zip(
            np.load(str(destination / "predictions.npz"))["run"],
            np.load(str(destination / "predictions.npz"))["event"],
            np.load(str(destination / "predictions.npz"))["label"],
        )
    }
    blind_ids = {
        (int(run), int(event), int(label))
        for run, event, label in zip(blind["run"], blind["event"], blind["label"])
    }
    if development_ids.intersection(blind_ids):
        raise RuntimeError("Blind event overlaps development validation")
    if set(c12_files + c13_files).intersection(used):
        raise RuntimeError("Blind file overlap detected after extraction")

    train, _ = er.load_cache(destination / "controlled_cache.npz")
    base_state = rc.make_feature_state(train, "no_edep")
    base_train = rc.apply_feature_state(train, base_state)
    base_blind = rc.apply_feature_state(blind, base_state)
    response = er.response_state(train)
    x_train, x_blind = er.append_energy(
        base_train,
        base_blind,
        er.energy_raw(train, "fitted", response),
        er.energy_raw(blind, "fitted", response),
    )
    del x_train

    weak_frozen = json.load(open(str(destination / "weak_frozen_models.json")))
    best_control = json.load(open(str(destination / "energy_best_control_results.json")))
    theta_no_edep = np.asarray(
        weak_frozen["models"]["weak_no_edep"]["theta"], dtype=float
    )
    best_parameters = np.asarray(
        best_control["models"]["nominal"]["parameters"], dtype=float
    )
    theta_best = np.asarray(
        best_control["models"]["nominal"]["theta"], dtype=float
    )
    if base_blind.shape[1] != len(theta_no_edep) or x_blind.shape[1] != len(theta_best):
        raise RuntimeError("Frozen model feature dimension mismatch")

    corrected = {
        "raw": blind["u"].copy(),
        "no_edep": blind["u"] - base_blind.dot(theta_no_edep),
        "best": blind["u"] - x_blind.dot(theta_best),
    }
    correction_best = blind["u"] - corrected["best"]

    # Freeze predictions before opening labels/truth for performance metrics.
    frozen_prediction_path = destination / "blind_predictions_frozen.npz"
    np.savez_compressed(
        str(frozen_prediction_path),
        run=blind["run"],
        event=blind["event"],
        R=blind["R"],
        relR=blind["relR"],
        bin=blind["bin"],
        prior_c13=blind["prior_c13"],
        u_raw=corrected["raw"],
        u_no_edep=corrected["no_edep"],
        u_best=corrected["best"],
    )
    freeze = {
        "blind_events": int(len(blind["u"])),
        "blind_source_files": len(c12_files) + len(c13_files),
        "development_manifest_overlap": 0,
        "selection_rule": "last lexicographic paths among files absent from development manifest",
        "weak_frozen_models_sha256": file_sha256(destination / "weak_frozen_models.json"),
        "energy_best_control_results_sha256": file_sha256(
            destination / "energy_best_control_results.json"
        ),
        "theta_no_edep_sha256": array_sha256(theta_no_edep),
        "theta_best_sha256": array_sha256(theta_best),
        "predictions_sha256": file_sha256(frozen_prediction_path),
        "truth_policy": "Corrections and this freeze record were written before blind labels or MC truth were used for metrics.",
    }
    (destination / "blind_freeze.json").write_text(
        json.dumps(freeze, indent=2, sort_keys=True)
    )
    rc.log("BLIND PREDICTIONS FROZEN %s" % freeze)

    weak_blind = {key: blind[key] for key in ["u", "R", "relR", "bin"]}
    density_best = rc.MixtureDensity(
        weak_blind, x_blind, blind["prior_c13"]
    )
    heldout = {
        "raw": float("nan"),
        "no_edep": float("nan"),
        "best": float(np.mean(density_best.event_nll(best_parameters))),
    }

    metric_rows = []
    bin_rows = []
    mass2_by_method = {}
    for name in ("raw", "no_edep", "best"):
        rows, mass2 = rc.method_metrics(name, corrected[name], blind, heldout[name])
        metric_rows.extend(rows)
        bin_rows.extend(rc.bin_metrics(name, corrected[name], mass2, blind))
        mass2_by_method[name] = mass2
    raw_width = {
        row["group"]: row["w68_du"]
        for row in metric_rows
        if row["method"] == "raw"
    }
    for row in metric_rows:
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_du"] / raw_width[row["group"]]
        )
    rc.write_tsv(destination / "blind_metrics.tsv", metric_rows)
    rc.write_tsv(destination / "blind_bin_metrics.tsv", bin_rows)

    bootstrap_result = bootstrap(
        blind, corrected, BOOTSTRAP_REPLICATES, BLIND_SEED + 1
    )
    overall = {
        row["method"]: row for row in metric_rows if row["group"] == "all"
    }
    output = {
        "freeze": freeze,
        "truth_valid_events": int(np.sum(blind["truth_valid"])),
        "versions": {
            "C12": sorted(int(value) for value in np.unique(pools[0]["version"])),
            "C13": sorted(int(value) for value in np.unique(pools[1]["version"])),
        },
        "overall_metrics": overall,
        "correction_diagnostics_post_freeze": {
            "best_correction_rms": rc.rms(correction_best),
            "best_correction_label_auc": rc.auc(correction_best, blind["label"]),
            "best_correction_median_c12": float(
                np.median(correction_best[blind["label"] == 0])
            ),
            "best_correction_median_c13": float(
                np.median(correction_best[blind["label"] == 1])
            ),
        },
        "paired_stratified_event_bootstrap_5000": bootstrap_result,
        "bootstrap_note": "Fixed frozen models; validation events resampled within isotope x reconstructed-rigidity strata. Does not include training, file-cluster, or production-version uncertainty.",
        "elapsed_seconds": time.time() - started,
        "script_sha256": rc.sha256_file(Path(__file__)),
    }
    (destination / "blind_results.json").write_text(
        json.dumps(rc.jsonable(output), indent=2, sort_keys=True)
    )

    prediction_output = {
        "run": blind["run"],
        "event": blind["event"],
        "label": blind["label"],
        "version": blind["version"],
        "bin": blind["bin"],
        "R": blind["R"],
        "relR": blind["relR"],
        "truth_valid": blind["truth_valid"],
        "u_true": blind["u_true"],
        "R_true": blind["R_true"],
    }
    for name, values in corrected.items():
        prediction_output["u_" + name] = values
        prediction_output["mass2_" + name] = mass2_by_method[name]
    np.savez_compressed(
        str(destination / "blind_predictions_with_truth.npz"), **prediction_output
    )
    rc.log(
        "DONE blind confirmation events=%d truth=%d elapsed=%.1fs"
        % (len(blind["u"]), np.sum(blind["truth_valid"]), time.time() - started)
    )


if __name__ == "__main__":
    main()
