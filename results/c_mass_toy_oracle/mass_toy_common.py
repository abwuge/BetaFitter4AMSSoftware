#!/usr/bin/env python3
"""Shared extraction and fitting utilities for the C12/C13 mass toy."""
from __future__ import print_function

import csv
import hashlib
import json
import math
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import uproot

CHARGE = 6.0
MASSES = np.asarray([11.17793, 12.11221], dtype=np.float64)
R_EDGES = np.asarray([2.5, 3.0, 4.0, 5.0, 6.0], dtype=np.float64)
C_LIGHT = 29.9792458
SEED_TEXT = "c-mass-toy-oracle-20260929-v2-common-branches"

GROUP_ORDER = ["mass", "timing", "energy", "rich", "nuisance"]
MODEL_GROUPS = {
    "mass_proxy": ["mass", "timing"],
    "mass_energy": ["mass", "timing", "energy"],
    "all_physics": ["mass", "timing", "energy", "rich"],
    "all_detector": ["mass", "timing", "energy", "rich", "nuisance"],
    "nuisance_control": ["nuisance"],
}

# Every branch below was checked to exist in both B1308 C12 and B1315 C13.
# Production-exclusive TOF variants are deliberately excluded: their mere
# presence would reveal the isotope through the production version.
BRANCHES = [
    "run", "event",
    "tof_tl", "tof_leng", "tof_pos", "tof_barid", "tof_edep",
    "tof_betah", "tof_ebetah", "tof_ebetah_n", "tof_tres", "tof_cres",
    "tof_chist", "tof_chisc", "tof_chist_n", "tof_chisc_n",
    "tof_qrms", "tof_nhits", "tof_nclhl", "tof_oncl", "tof_oq",
    "tk_ql", "tk_rigidity1", "tk_rigidityk", "tk_erigidity1",
    "tk_rigidityi", "tk_chisi", "tk_bcor", "tk_cdif", "tk_fzm",
    "tk_qrms", "tk_qhit",
    "rich_beta", "rich_ebeta", "rich_ebeta1", "rich_good", "rich_clean",
    "rich_NaF", "rich_npe", "rich_q", "rich_width", "rich_hit",
    "rich_used", "rich_usedm", "rich_stat", "rich_BetaConsistency",
    "rich_tile", "rich_distb",
    "trd_qk", "trd_qnhk", "trd_nhitk", "trd_lik", "trd_like",
    "trd_amplk", "trd_pathlk", "trd_statk", "trd_onhitk", "trd_oampk",
    "ecal_q", "ecal_en", "ecal_bdt", "ecal_stat", "ecal_nhit",
    "ecal_el", "ecal_eh",
]


def log(message):
    print(message, flush=True)


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
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


def w68(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan")
    low, high = np.percentile(values, [16.0, 84.0])
    return float((high - low) / 2.0)


def auc(score, label):
    score = np.asarray(score, dtype=np.float64)
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


def sigmoid(values):
    values = np.asarray(values, dtype=np.float64)
    output = np.empty_like(values)
    positive = values >= 0.0
    output[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exponent = np.exp(values[~positive])
    output[~positive] = exponent / (1.0 + exponent)
    return output


def balanced_weights(label):
    label = np.asarray(label, dtype=np.int8)
    n = float(len(label))
    n0 = max(1, int(np.sum(label == 0)))
    n1 = max(1, int(np.sum(label == 1)))
    return np.where(label == 1, 0.5 * n / n1, 0.5 * n / n0).astype(np.float64)


def rejection_at_efficiency(score, label, efficiency=0.50):
    score = np.asarray(score, dtype=np.float64)
    label = np.asarray(label, dtype=np.int8)
    finite = np.isfinite(score)
    score = score[finite]
    label = label[finite]
    signal = score[label == 1]
    background = score[label == 0]
    if not len(signal) or not len(background):
        return float("nan")
    threshold = float(np.percentile(signal, 100.0 * (1.0 - efficiency)))
    acceptance = float(np.mean(background >= threshold))
    return float("inf") if acceptance <= 0.0 else 1.0 / acceptance


def balanced_accuracy(probability, label, threshold=0.5):
    probability = np.asarray(probability, dtype=np.float64)
    label = np.asarray(label, dtype=np.int8)
    prediction = probability >= threshold
    tpr = np.mean(prediction[label == 1]) if np.any(label == 1) else np.nan
    tnr = np.mean(~prediction[label == 0]) if np.any(label == 0) else np.nan
    return float(0.5 * (tpr + tnr))


def stable_logloss(probability, label):
    probability = np.clip(np.asarray(probability, dtype=np.float64), 1e-9, 1.0 - 1e-9)
    label = np.asarray(label, dtype=np.float64)
    weights = balanced_weights(label)
    loss = -(label * np.log(probability) + (1.0 - label) * np.log1p(-probability))
    return float(np.mean(weights * loss))


def safe_log_positive(values):
    values = np.asarray(values, dtype=np.float64)
    output = np.full(values.shape, np.nan, dtype=np.float64)
    good = np.isfinite(values) & (values > 0.0)
    output[good] = np.log(values[good])
    return output


def safe_valid(values, low=-np.inf, high=np.inf):
    values = np.asarray(values, dtype=np.float64)
    output = values.copy()
    good = np.isfinite(output) & (output > low) & (output < high)
    output[~good] = np.nan
    return output


def flatten_columns(array, prefix):
    array = np.asarray(array, dtype=np.float64)
    if array.ndim == 1:
        return array[:, None], [prefix]
    flat = array.reshape(len(array), -1)
    return flat, ["%s_%d" % (prefix, index) for index in range(flat.shape[1])]


def add_block(blocks, names, group, values, prefix, transform=None):
    values = np.asarray(values, dtype=np.float64)
    if transform is not None:
        values = transform(values)
    matrix, local_names = flatten_columns(values, prefix)
    blocks[group].append(matrix)
    names[group].extend(local_names)


def build_file_features(arrays, entries, target_u, target_r):
    """Build common-branch feature groups for selected entries in one file."""
    a = {
        key: np.asarray(value[entries])
        for key, value in arrays.items()
        if key not in ("run", "event")
    }
    n = len(entries)
    blocks = {name: [] for name in GROUP_ORDER}
    names = {name: [] for name in GROUP_ORDER}

    times = a["tof_tl"].astype(np.float64)
    path_ns = a["tof_leng"].astype(np.float64) / C_LIGHT
    centered_path = path_ns - path_ns.mean(axis=1, keepdims=True)
    centered_time = times - times.mean(axis=1, keepdims=True)
    denominator = np.sum(centered_path * centered_path, axis=1)
    u_recomputed = np.sum(centered_path * centered_time, axis=1) / denominator
    residual = centered_time - u_recomputed[:, None] * centered_path
    span = np.ptp(path_ns, axis=1)
    positions = a["tof_pos"].astype(np.float64)
    dz = positions[:, 3, 2] - positions[:, 0, 2]
    direction_x = (positions[:, 3, 0] - positions[:, 0, 0]) / dz
    direction_y = (positions[:, 3, 1] - positions[:, 0, 1]) / dz
    normalized_residual = residual / span[:, None]

    maximum_u_difference = float(np.max(np.abs(u_recomputed - target_u))) if n else 0.0
    baseline_r = a["tk_rigidityk"][:, 1, 1].astype(np.float64)
    maximum_r_difference = float(np.max(np.abs(baseline_r - target_r))) if n else 0.0

    # Detector geometry and topology.  This group intentionally carries no
    # beta, rigidity, energy-loss or RICH measurement and is used as a source/
    # production shortcut diagnostic.
    geometry = np.column_stack([
        positions[:, :, 0] / 120.0,
        positions[:, :, 1] / 120.0,
        direction_x,
        direction_y,
        centered_path,
        np.log(span / 4.0),
    ])
    add_block(blocks, names, "nuisance", geometry, "geometry")
    bars = a["tof_barid"].astype(np.int16)
    bar_onehot = np.column_stack([
        (bars[:, layer] == bar).astype(np.float64)
        for layer in range(4) for bar in range(10)
    ])
    add_block(blocks, names, "nuisance", bar_onehot, "bar")

    # Baseline and hybrid mass observables.
    raw_m2 = (CHARGE * target_r) ** 2 * (target_u * target_u - 1.0)
    raw_signed_mass = np.sign(raw_m2) * np.sqrt(np.abs(raw_m2))
    core = np.column_stack([
        target_u,
        1.0 / target_u,
        target_r,
        np.log(target_r),
        raw_m2,
        raw_signed_mass,
        target_u * target_r,
    ])
    add_block(blocks, names, "mass", core, "linear_core")

    hybrid_beta = safe_valid(a["tof_betah"], low=0.2, high=1.5)
    hybrid_u = 1.0 / hybrid_beta
    hybrid_m2 = (CHARGE * target_r) ** 2 * (hybrid_u * hybrid_u - 1.0)
    hybrid = np.column_stack([
        hybrid_u,
        hybrid_u - target_u,
        hybrid_m2,
        np.sign(hybrid_m2) * np.sqrt(np.abs(hybrid_m2)),
    ])
    add_block(blocks, names, "mass", hybrid, "hybrid_core")
    add_block(blocks, names, "mass", safe_log_positive(a["tof_ebetah"]), "log_tof_ebetah")
    add_block(blocks, names, "mass", safe_log_positive(a["tof_ebetah_n"]), "log_tof_ebetah_n")

    # Raw TOF residuals retain event-level timing information not captured by a
    # single beta estimate.
    add_block(blocks, names, "timing", normalized_residual, "linear_time_residual")
    add_block(blocks, names, "timing", a["tof_tres"], "tof_tres")
    add_block(blocks, names, "timing", a["tof_cres"], "tof_cres")

    # Tracker multi-fit redundancy.  The inputs are differences among measured
    # fits and candidate masses, never their differences to MC truth.
    tracker_candidates = []
    tracker_candidate_names = []
    r1 = np.asarray(a["tk_rigidity1"], dtype=np.float64)
    rk = np.asarray(a["tk_rigidityk"], dtype=np.float64)
    ri = np.asarray(a["tk_rigidityi"], dtype=np.float64)
    for algorithm in range(r1.shape[1]):
        for alignment in range(r1.shape[2]):
            for fit_span in range(r1.shape[3]):
                tracker_candidates.append(r1[:, algorithm, alignment, fit_span])
                tracker_candidate_names.append("r1_a%d_v%d_s%d" % (algorithm, alignment, fit_span))
    for plane in range(rk.shape[1]):
        for fit_span in range(rk.shape[2]):
            tracker_candidates.append(rk[:, plane, fit_span])
            tracker_candidate_names.append("rk_z%d_s%d" % (plane, fit_span))
    for index in range(ri.shape[1]):
        tracker_candidates.append(ri[:, index])
        tracker_candidate_names.append("ri_%d" % index)
    tracker_candidates = np.column_stack(tracker_candidates)
    ratio = tracker_candidates / target_r[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ratio = np.log(ratio)
    valid_ratio = np.isfinite(log_ratio) & (np.abs(log_ratio) < 1.5)
    log_ratio[~valid_ratio] = np.nan
    candidate_m2 = raw_m2[:, None] * ratio * ratio
    candidate_m2[~valid_ratio] = np.nan
    blocks["mass"].extend([log_ratio, candidate_m2])
    names["mass"].extend("logratio_%s" % item for item in tracker_candidate_names)
    names["mass"].extend("m2_%s" % item for item in tracker_candidate_names)
    candidate_summary = np.column_stack([
        np.nanmedian(log_ratio, axis=1),
        np.nanstd(log_ratio, axis=1),
        np.nanpercentile(log_ratio, 10.0, axis=1),
        np.nanpercentile(log_ratio, 90.0, axis=1),
        np.nanmedian(candidate_m2, axis=1),
        np.nanstd(candidate_m2, axis=1),
    ])
    add_block(blocks, names, "mass", candidate_summary, "tracker_candidate_summary")

    relative_error = np.abs(np.asarray(a["tk_erigidity1"], dtype=np.float64) * target_r[:, None, None])
    add_block(blocks, names, "mass", safe_log_positive(relative_error), "log_tracker_relerr")
    baseline_rel_error = relative_error[:, 2, 1]
    add_block(blocks, names, "mass", baseline_rel_error, "baseline_relR")

    # Measured TOF and tracker energy loss.
    tof_dedx = np.log(np.asarray(a["tof_edep"], dtype=np.float64) / (CHARGE * CHARGE))
    tracker_dedx = 2.0 * np.log(np.asarray(a["tk_ql"], dtype=np.float64)[:, 1:8] / CHARGE)
    add_block(blocks, names, "energy", tof_dedx, "log_tof_edep")
    add_block(blocks, names, "energy", tracker_dedx, "log_tracker_dedx")
    energy_summary = np.column_stack([
        np.mean(tof_dedx, axis=1),
        np.std(tof_dedx, axis=1),
        np.mean(tracker_dedx, axis=1),
        np.std(tracker_dedx, axis=1),
        np.median(tracker_dedx, axis=1),
    ])
    add_block(blocks, names, "energy", energy_summary, "dedx_summary")

    # RICH beta is kept as an independent physical mass proxy.  Sentinel values
    # become missing; the preprocessor learns event-validity indicators using
    # training events only.
    rich_beta = safe_valid(a["rich_beta"], low=0.2, high=1.2)
    rich_u = 1.0 / rich_beta
    rich_du = rich_u - target_u[:, None]
    rich_m2 = (CHARGE * target_r[:, None]) ** 2 * (rich_u * rich_u - 1.0)
    blocks["rich"].extend([rich_u, rich_du, rich_m2])
    rich_width = rich_beta.reshape(n, -1).shape[1]
    names["rich"].extend("rich_u_%d" % index for index in range(rich_width))
    names["rich"].extend("rich_du_%d" % index for index in range(rich_width))
    names["rich"].extend("rich_m2_%d" % index for index in range(rich_width))
    add_block(blocks, names, "rich", safe_log_positive(a["rich_ebeta"]), "log_rich_ebeta")
    add_block(blocks, names, "rich", safe_log_positive(a["rich_ebeta1"]), "log_rich_ebeta1")
    for branch in [
        "rich_good", "rich_clean", "rich_NaF", "rich_hit", "rich_used",
        "rich_usedm", "rich_stat", "rich_BetaConsistency", "rich_tile", "rich_distb",
    ]:
        add_block(blocks, names, "rich", a[branch], branch)
    add_block(blocks, names, "rich", safe_log_positive(a["rich_npe"]), "log_rich_npe")
    add_block(blocks, names, "rich", safe_valid(a["rich_q"], low=0.0, high=100.0), "rich_q")
    add_block(blocks, names, "rich", safe_valid(a["rich_width"], low=0.0, high=100.0), "rich_width")

    # TRD and ECAL energy/charge response.
    add_block(blocks, names, "energy", safe_log_positive(a["trd_qk"]), "log_trd_qk")
    add_block(blocks, names, "energy", safe_valid(a["trd_lik"], low=-2.5, high=100.0), "trd_lik")
    add_block(blocks, names, "energy", safe_valid(a["trd_like"], low=-0.5, high=100.0), "trd_like")
    trd_amp = np.asarray(a["trd_amplk"], dtype=np.float64)
    trd_path = np.asarray(a["trd_pathlk"], dtype=np.float64)
    trd_dedx = np.full(trd_amp.shape, np.nan, dtype=np.float64)
    good_trd = np.isfinite(trd_amp) & np.isfinite(trd_path) & (trd_amp > 0.0) & (trd_path > 0.0)
    trd_dedx[good_trd] = np.log(trd_amp[good_trd] / trd_path[good_trd])
    add_block(blocks, names, "energy", trd_dedx, "log_trd_dedx")
    add_block(blocks, names, "energy", safe_log_positive(a["ecal_q"]), "log_ecal_q")
    add_block(blocks, names, "energy", safe_log_positive(a["ecal_en"]), "log_ecal_en")
    add_block(blocks, names, "energy", safe_valid(a["ecal_bdt"], low=-2.5, high=2.5), "ecal_bdt")
    add_block(blocks, names, "energy", safe_log_positive(a["ecal_el"]), "log_ecal_el")
    add_block(blocks, names, "energy", safe_log_positive(a["ecal_eh"]), "log_ecal_eh")

    # Generic fit/cluster quality and topology.  This is deliberately isolated
    # so its standalone AUC measures residual production/source shortcuts.
    for branch in ["tof_chist", "tof_chisc", "tof_chist_n", "tof_chisc_n"]:
        add_block(blocks, names, "nuisance", safe_log_positive(a[branch]), "log_" + branch)
    add_block(blocks, names, "nuisance", safe_log_positive(a["tof_qrms"]), "log_tof_qrms")
    for branch in ["tof_nhits", "tof_nclhl", "tof_oncl"]:
        add_block(blocks, names, "nuisance", a[branch], branch)
    add_block(blocks, names, "nuisance", safe_log_positive(np.abs(a["tk_chisi"])), "log_tk_chisi")
    for branch in ["tk_bcor", "tk_cdif", "tk_fzm", "tk_qhit"]:
        add_block(blocks, names, "nuisance", a[branch], branch)
    add_block(blocks, names, "nuisance", safe_log_positive(a["tk_qrms"]), "log_tk_qrms")
    for branch in ["trd_qnhk", "trd_nhitk", "trd_statk", "trd_onhitk", "trd_oampk", "ecal_stat", "ecal_nhit"]:
        add_block(blocks, names, "nuisance", a[branch], branch)

    output = {}
    for group in GROUP_ORDER:
        output[group] = np.column_stack(blocks[group]).astype(np.float32)
        if output[group].shape[1] != len(names[group]):
            raise RuntimeError("Feature name mismatch in %s" % group)

    special = {
        "u_replay_max": maximum_u_difference,
        "R_replay_max": maximum_r_difference,
        "hybrid_beta": np.asarray(a["tof_betah"], dtype=np.float64),
        "baseline_relR": baseline_rel_error.astype(np.float64),
    }
    return output, names, special


def extract_one_file(task):
    label, source, path, rows, target_run, target_event, target_u, target_r = task
    started = time.time()
    with uproot.open(path) as root_file:
        tree = root_file["amstreea"]
        missing = [name for name in BRANCHES if name not in tree.keys()]
        if missing:
            raise RuntimeError("Missing common branches in %s: %s" % (path, ", ".join(missing)))
        arrays = tree.arrays(BRANCHES, library="np")
    ids = np.rec.fromarrays(
        [arrays["run"].astype(np.int64), arrays["event"].astype(np.int64)],
        names="run,event",
    )
    target_ids = np.rec.fromarrays(
        [target_run.astype(np.int64), target_event.astype(np.int64)],
        names="run,event",
    )
    common, entries, targets = np.intersect1d(
        ids, target_ids, assume_unique=False, return_indices=True
    )
    if len(common) != len(rows):
        raise RuntimeError(
            "Event match failure label=%d source=%d got=%d need=%d"
            % (label, source, len(common), len(rows))
        )
    order = np.argsort(targets)
    entries = entries[order]
    target_order = targets[order]
    rows_out = rows[target_order]
    features, names, special = build_file_features(
        arrays, entries, target_u[target_order], target_r[target_order]
    )
    audit = {
        "label": int(label),
        "source": int(source),
        "path": path,
        "events": int(len(rows_out)),
        "seconds": float(time.time() - started),
        "u_replay_max": special["u_replay_max"],
        "R_replay_max": special["R_replay_max"],
    }
    return rows_out, features, names, special, audit


def load_manifest(path):
    mapping = {0: {}, 1: {}}
    label_map = {"C12": 0, "C13": 1}
    with open(str(path), newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            mapping[label_map[row["isotope"]]][int(row["source"])] = row["path"]
    return mapping


def split_sources(manifest):
    target_counts = {
        0: {"train": 160, "valid": 48, "test": 48},
        1: {"train": 20, "valid": 6, "test": 6},
    }
    role_by_source = {0: {}, 1: {}}
    rows = []
    for label in (0, 1):
        items = list(manifest[label].items())
        items.sort(
            key=lambda item: hashlib.sha256(
                (SEED_TEXT + ":" + item[1]).encode("utf-8")
            ).hexdigest()
        )
        expected = sum(target_counts[label].values())
        if len(items) != expected:
            raise RuntimeError(
                "Expected %d files for label %d, found %d" % (expected, label, len(items))
            )
        cursor = 0
        for role in ("train", "valid", "test"):
            count = target_counts[label][role]
            for source, path in items[cursor : cursor + count]:
                role_by_source[label][source] = role
                rows.append({
                    "isotope": "C12" if label == 0 else "C13",
                    "source": int(source),
                    "role": role,
                    "path": path,
                    "hash": hashlib.sha256((SEED_TEXT + ":" + path).encode("utf-8")).hexdigest(),
                })
            cursor += count
    return role_by_source, rows


def write_tsv(path, rows):
    if not rows:
        return
    with open(str(path), "w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def extract_all_features(data, manifest, processes=4):
    n = len(data["u"])
    tasks = []
    for label in (0, 1):
        for source, path in sorted(manifest[label].items()):
            rows = np.flatnonzero((data["label"] == label) & (data["source"] == source))
            if not len(rows):
                continue
            tasks.append((
                label, source, path, rows,
                data["run"][rows], data["event"][rows],
                data["u"][rows], data["R"][rows],
            ))

    groups = None
    group_names = None
    hybrid_beta = np.full(n, np.nan, dtype=np.float64)
    baseline_relR = np.full(n, np.nan, dtype=np.float64)
    filled = np.zeros(n, dtype=np.bool_)
    audits = []
    started = time.time()
    with mp.Pool(processes=max(1, int(processes))) as pool:
        iterator = pool.imap_unordered(extract_one_file, tasks, chunksize=1)
        for completed, (rows, local, names, special, audit) in enumerate(iterator, 1):
            if groups is None:
                group_names = names
                groups = {
                    group: np.empty((n, local[group].shape[1]), dtype=np.float32)
                    for group in GROUP_ORDER
                }
            for group in GROUP_ORDER:
                if names[group] != group_names[group]:
                    raise RuntimeError("Feature schema changed for group %s" % group)
                groups[group][rows] = local[group]
            hybrid_beta[rows] = special["hybrid_beta"]
            baseline_relR[rows] = special["baseline_relR"]
            filled[rows] = True
            audits.append(audit)
            if completed % 24 == 0 or completed == len(tasks):
                log(
                    "feature extraction %d/%d files, elapsed %.1fs"
                    % (completed, len(tasks), time.time() - started)
                )
    if groups is None or not np.all(filled):
        raise RuntimeError("Feature extraction left unfilled events")
    return groups, group_names, hybrid_beta, baseline_relR, audits


def fit_preprocessor(raw, names, train_mask):
    train = np.asarray(raw[train_mask], dtype=np.float64)
    finite_fraction = np.mean(np.isfinite(train), axis=0)
    keep = finite_fraction >= 0.005
    kept_indices = np.flatnonzero(keep)
    train = train[:, keep]
    kept_names = [name for name, selected in zip(names, keep) if selected]
    finite_fraction = finite_fraction[keep]

    low = np.nanpercentile(train, 0.5, axis=0)
    high = np.nanpercentile(train, 99.5, axis=0)
    clipped = np.minimum(np.maximum(train, low), high)
    center = np.nanmedian(clipped, axis=0)
    q25, q75 = np.nanpercentile(clipped, [25.0, 75.0], axis=0)
    robust_scale = 0.741301109 * (q75 - q25)
    standard_scale = np.nanstd(clipped, axis=0)
    scale = np.maximum(robust_scale, 0.10 * standard_scale)
    informative = np.isfinite(center) & np.isfinite(scale) & (scale > 1e-7)

    kept_indices = kept_indices[informative]
    kept_names = [name for name, selected in zip(kept_names, informative) if selected]
    low = low[informative]
    high = high[informative]
    center = center[informative]
    scale = scale[informative]
    finite_fraction = finite_fraction[informative]
    missing_fraction = 1.0 - finite_fraction
    missing_indicator = (missing_fraction >= 0.005) & (missing_fraction <= 0.995)

    # Add z^2 only for genuinely continuous variables.  Binary/topological
    # columns would merely duplicate themselves after standardization.
    selected_train = train[:, informative]
    square_mask = np.zeros(selected_train.shape[1], dtype=np.bool_)
    for column in range(selected_train.shape[1]):
        values = selected_train[:, column]
        values = values[np.isfinite(values)]
        if len(values) and len(np.unique(values[: min(len(values), 20000)])) > 6:
            square_mask[column] = True

    output_names = list(kept_names)
    output_names += ["square_" + kept_names[index] for index in np.flatnonzero(square_mask)]
    output_names += ["missing_" + kept_names[index] for index in np.flatnonzero(missing_indicator)]
    return {
        "indices": kept_indices.astype(np.int32),
        "low": low,
        "high": high,
        "center": center,
        "scale": scale,
        "square_mask": square_mask,
        "missing_indicator": missing_indicator,
        "feature_names": kept_names,
        "output_names": output_names,
        "finite_fraction": finite_fraction,
    }


def apply_preprocessor(raw, state):
    values = np.asarray(raw[:, state["indices"]], dtype=np.float64)
    missing = ~np.isfinite(values)
    values = np.minimum(np.maximum(values, state["low"]), state["high"])
    values = np.where(missing, state["center"], values)
    standardized = np.clip((values - state["center"]) / state["scale"], -6.0, 6.0)
    columns = [standardized]
    square_mask = state["square_mask"]
    if np.any(square_mask):
        columns.append(standardized[:, square_mask] ** 2)
    indicator = state["missing_indicator"]
    if np.any(indicator):
        columns.append(missing[:, indicator].astype(np.float64))
    return np.column_stack(columns).astype(np.float32)


def fit_logistic(x, y, regularization, initial=None, maxiter=180):
    # scipy's optimizer owns float64 parameters.  Keeping x in float32 makes
    # numpy cast the complete matrix to float64 on every dot product, which is
    # far more expensive than one conversion per fit.
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n, dimensions = x.shape
    weights = balanced_weights(y)

    def objective(parameters):
        intercept = parameters[0]
        coefficients = parameters[1:]
        logits = intercept + x.dot(coefficients)
        losses = np.logaddexp(0.0, logits) - y * logits
        probability = sigmoid(logits)
        residual = weights * (probability - y)
        value = float(np.mean(weights * losses))
        value += 0.5 * regularization * float(np.dot(coefficients, coefficients))
        gradient = np.empty(dimensions + 1, dtype=np.float64)
        gradient[0] = np.mean(residual)
        gradient[1:] = x.T.dot(residual) / float(n) + regularization * coefficients
        return value, gradient

    if initial is None:
        initial = np.zeros(dimensions + 1, dtype=np.float64)
    fit = minimize(
        objective,
        initial,
        jac=True,
        method="L-BFGS-B",
        options={
            "maxiter": int(maxiter), "maxls": 40, "ftol": 1e-10,
            "gtol": 2e-6, "maxcor": 20,
        },
    )
    gradient_inf = float(np.max(np.abs(objective(fit.x)[1])))
    return {
        "parameters": fit.x,
        "success": bool(fit.success),
        "message": str(fit.message),
        "iterations": int(fit.nit),
        "objective": float(fit.fun),
        "gradient_inf": gradient_inf,
    }


def predict_binned(x, bins, models):
    output = np.full(len(x), np.nan, dtype=np.float64)
    for bin_id in range(len(R_EDGES) - 1):
        rows = bins == bin_id
        parameters = np.asarray(models[str(bin_id)]["parameters"], dtype=np.float64)
        output[rows] = parameters[0] + np.asarray(x[rows], dtype=np.float64).dot(parameters[1:])
    if not np.all(np.isfinite(output)):
        raise RuntimeError("Non-finite binned predictions")
    return output


def macro_auc(score, labels, bins, mask):
    values = []
    for bin_id in range(len(R_EDGES) - 1):
        rows = mask & (bins == bin_id)
        values.append(auc(score[rows], labels[rows]))
    return float(np.mean(values))


def choose_regularization(x, labels, bins, train_mask, valid_mask, candidates):
    scan = []
    best = None
    warm = {str(bin_id): None for bin_id in range(len(R_EDGES) - 1)}
    for regularization in candidates:
        models = {}
        for bin_id in range(len(R_EDGES) - 1):
            rows = train_mask & (bins == bin_id)
            fit = fit_logistic(
                x[rows], labels[rows], regularization,
                initial=warm[str(bin_id)], maxiter=180,
            )
            models[str(bin_id)] = fit
            warm[str(bin_id)] = fit["parameters"]
        logits = predict_binned(x, bins, models)
        train_auc = macro_auc(logits, labels, bins, train_mask)
        valid_auc = macro_auc(logits, labels, bins, valid_mask)
        item = {
            "regularization": float(regularization),
            "train_macro_auc": train_auc,
            "valid_macro_auc": valid_auc,
            "models": models,
        }
        scan.append(item)
        log(
            "lambda %.3g train macro AUC %.6f valid macro AUC %.6f"
            % (regularization, train_auc, valid_auc)
        )
        if best is None or valid_auc > best["valid_macro_auc"] + 1e-12 or (
            abs(valid_auc - best["valid_macro_auc"]) <= 1e-12
            and regularization > best["regularization"]
        ):
            best = item
    return best, scan


def fit_calibration(logits, labels, bins, valid_mask, score_mode=False):
    """Fit equal-prior one-dimensional logistic calibration in each R bin."""
    state = {}
    output = np.full(len(logits), np.nan, dtype=np.float64)
    for bin_id in range(len(R_EDGES) - 1):
        rows = valid_mask & (bins == bin_id) & np.isfinite(logits)
        x = np.asarray(logits[rows], dtype=np.float64)
        y = np.asarray(labels[rows], dtype=np.float64)
        center = float(np.median(x))
        scale = float(max(1e-9, 0.741301109 * (np.percentile(x, 75.0) - np.percentile(x, 25.0))))
        z = (x - center) / scale
        weights = balanced_weights(y)

        def objective(parameters):
            calibrated = parameters[0] + parameters[1] * z
            loss = np.logaddexp(0.0, calibrated) - y * calibrated
            probability = sigmoid(calibrated)
            value = float(np.mean(weights * loss)) + 5e-5 * float(parameters[1] ** 2)
            gradient = np.asarray([
                np.mean(weights * (probability - y)),
                np.mean(weights * (probability - y) * z) + 1e-4 * parameters[1],
            ])
            return value, gradient

        slope_bounds = (-20.0, 20.0) if score_mode else (0.02, 20.0)
        fit = minimize(
            objective,
            np.asarray([0.0, 1.0]),
            jac=True,
            method="L-BFGS-B",
            bounds=[(-20.0, 20.0), slope_bounds],
            options={"maxiter": 200, "ftol": 1e-12, "gtol": 1e-8},
        )
        state[str(bin_id)] = {
            "center": center,
            "scale": scale,
            "intercept": float(fit.x[0]),
            "slope": float(fit.x[1]),
            "success": bool(fit.success),
            "objective": float(fit.fun),
        }
    for bin_id in range(len(R_EDGES) - 1):
        rows = bins == bin_id
        item = state[str(bin_id)]
        z = (logits[rows] - item["center"]) / item["scale"]
        output[rows] = sigmoid(item["intercept"] + item["slope"] * z)
    return state, output
