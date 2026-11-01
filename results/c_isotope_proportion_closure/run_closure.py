#!/usr/bin/env python3
"""C12/C13 weak isotope-proportion closure for the TOF inverse-beta estimator.

Weak fits receive measured detector quantities, reconstructed rigidity, bag id,
and a known C13 fraction for each (bag, rigidity-bin). They never receive the
per-event isotope label or MC momentum. MC truth is used only after weak
models are frozen, for closure metrics and separate oracle benchmarks.
"""
import csv
import glob
import hashlib
import json
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import uproot

C_LIGHT = 29.9792458
CHARGE = 6.0
NU = 6.0
MASSES = np.array([11.17793, 12.11221], dtype=float)
ISOTOPE_NAMES = np.array(["C12", "C13"])
R_EDGES = np.array([2.5, 3.0, 4.0, 5.0, 6.0], dtype=float)
N_BINS = len(R_EDGES) - 1
T_NORMALIZATION = (
    math.lgamma((NU + 1.0) / 2.0)
    - math.lgamma(NU / 2.0)
    - 0.5 * math.log(NU * math.pi)
)
SEED = 1301315

C12_PATTERN = "/lustre02/user/yachen/public/hxwu/amsd69n_mc_z1_z2_z6_z8_20260903/Z6_C12/*.root"
C13_PATTERN = "/lustre03/data/eos/ams/group/mit/amsd69nMCC13l1_c13.B1315_612000_2022/*.root"
FILES_PER_ROLE = {0: 32, 1: 32}
C13_FILES_PER_ROLE = {0: 8, 1: 8}
TRAIN_FRACTIONS = np.array(
    [[0.20, 0.25, 0.30, 0.35], [0.50, 0.50, 0.50, 0.50]], dtype=float
)
VALID_FRACTIONS = np.array([[0.35, 0.35, 0.35, 0.35]], dtype=float)
TRAIN_TARGET_PER_BAG_BIN = 900
VALID_TARGET_PER_BIN = 1600

BRANCHES = [
    "run", "event", "version", "mpar", "mch", "mevmom1",
    "tof_betah", "tof_bpatt", "tof_tl", "tof_leng", "tof_pos",
    "tof_edep", "tof_barid", "tk_ql", "tk_q", "tk_qin", "tk_pos",
    "tk_rigidityk", "tk_erigidity1", "tk_chis1",
]

CONTINUOUS_NAMES = (
    ["tof_x%d" % i for i in range(4)]
    + ["tof_y%d" % i for i in range(4)]
    + ["tof_dx_dz", "tof_dy_dz"]
    + ["path_centered%d" % i for i in range(4)]
    + ["time_residual%d" % i for i in range(4)]
    + ["log_tof_edep%d" % i for i in range(4)]
    + ["log_tracker_dedx%d" % i for i in range(1, 8)]
    + ["log_path_span"]
)
FULL_INDICES = list(range(30))
NO_EDEP_INDICES = list(range(18)) + [29]
LEGACY_NAMES = [
    "mean_tof_x", "mean_tof_y", "mean_tof_x_sq", "mean_tof_y_sq",
    "log_path_span", "upper_time_contrast", "lower_time_contrast",
]
FULL_INTERACTIONS = (
    [(14 + i, 18 + i) for i in range(4)]
    + [(i, 18 + i) for i in range(4)]
    + [(4 + i, 18 + i) for i in range(4)]
    + [(10 + i, 14 + i) for i in range(4)]
)
NO_EDEP_INTERACTIONS = [(10 + i, 14 + i) for i in range(4)]


def log(message):
    print(message, flush=True)


def w68(values):
    values = np.asarray(values)
    if values.size == 0:
        return float("nan")
    q = np.percentile(values, [16.0, 84.0])
    return float((q[1] - q[0]) / 2.0)


def rms(values):
    values = np.asarray(values)
    return float(np.sqrt(np.mean(values * values))) if values.size else float("nan")


def auc(score, label):
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=np.int32)
    finite = np.isfinite(score) & ((label == 0) | (label == 1))
    score, label = score[finite], label[finite]
    n1, n0 = int(np.sum(label == 1)), int(np.sum(label == 0))
    if n0 == 0 or n1 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    s, y = score[order], label[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    ends = np.r_[starts[1:], len(s)]
    ranks = np.repeat((starts + 1.0 + ends) / 2.0, ends - starts)
    return float(
        (np.sum(ranks[y == 1]) - n1 * (n1 + 1.0) / 2.0) / (float(n0) * n1)
    )


def subset(data, rows):
    return {key: value[rows] for key, value in data.items()}


def concatenate(chunks):
    if not chunks:
        raise ValueError("No chunks to concatenate")
    return {
        key: np.concatenate([chunk[key] for chunk in chunks], axis=0)
        for key in chunks[0]
    }


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
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


def choose_files(pattern, per_role):
    roles = {0: [], 1: []}
    for path in glob.glob(pattern):
        name = os.path.basename(path)
        role = int(hashlib.sha256(name.encode("utf-8")).hexdigest(), 16) % 2
        roles[role].append(path)
    for role in roles:
        roles[role] = sorted(roles[role])[: per_role[role]]
        if len(roles[role]) != per_role[role]:
            raise RuntimeError("Not enough input files for role %d" % role)
    if set(roles[0]).intersection(roles[1]):
        raise RuntimeError("File split overlap")
    return roles


def extract_pool(paths, isotope_label, expected_mpar, mass, role, audit_rows):
    chunks = []
    for source_index, path in enumerate(paths):
        started = time.time()
        with uproot.open(path) as root_file:
            tree = root_file["amstreea"]
            arrays = tree.arrays(BRANCHES, library="np")
            raw_entries = int(tree.num_entries)

        q = np.where(
            arrays["tk_qin"][:, 0, 2] < 2.5,
            arrays["tk_q"][:, 1],
            arrays["tk_qin"][:, 0, 2],
        )
        reconstructed_charge = np.where(np.isfinite(q), q + 0.5, -99).astype(np.int16)
        times = arrays["tof_tl"].astype(np.float64)
        path_ns = arrays["tof_leng"].astype(np.float64) / C_LIGHT
        centered_path = path_ns - path_ns.mean(axis=1, keepdims=True)
        centered_time = times - times.mean(axis=1, keepdims=True)
        denominator = np.sum(centered_path * centered_path, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            inverse_beta = np.sum(centered_path * centered_time, axis=1) / denominator
        residual = centered_time - inverse_beta[:, None] * centered_path
        span = np.ptp(path_ns, axis=1)
        positions = arrays["tof_pos"].astype(np.float64)
        dz = positions[:, 3, 2] - positions[:, 0, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            direction_x = (positions[:, 3, 0] - positions[:, 0, 0]) / dz
            direction_y = (positions[:, 3, 1] - positions[:, 0, 1]) / dz

        rigidity = arrays["tk_rigidityk"][:, 1, 1].astype(np.float64)
        relative_rigidity_error = np.abs(
            arrays["tk_erigidity1"][:, 2, 1].astype(np.float64) * rigidity
        )
        chi_x = arrays["tk_chis1"][:, 2, 2, 1, 0].astype(np.float64)
        chi_y = arrays["tk_chis1"][:, 2, 2, 1, 1].astype(np.float64)
        momentum = arrays["mevmom1"][:, 10].astype(np.float64)

        good = (arrays["mch"] == CHARGE) & (arrays["mpar"] == expected_mpar)
        good &= (reconstructed_charge == int(CHARGE)) & (arrays["tof_bpatt"] == 4444)
        good &= np.all(np.isfinite(times) & np.isfinite(arrays["tof_leng"]), axis=1)
        good &= np.isfinite(inverse_beta) & (inverse_beta > 0.0) & (denominator > 1e-12)
        good &= np.isfinite(span) & (span > 0.0)
        good &= (
            np.all(np.isfinite(positions), axis=(1, 2))
            & np.isfinite(dz)
            & (np.abs(dz) > 1.0)
        )
        good &= np.isfinite(direction_x) & np.isfinite(direction_y)
        good &= np.all(
            np.isfinite(arrays["tof_edep"]) & (arrays["tof_edep"] > 0.0), axis=1
        )
        good &= np.all(
            (arrays["tof_barid"] >= 0) & (arrays["tof_barid"] < 10), axis=1
        )
        good &= np.all(
            np.isfinite(arrays["tk_ql"][:, 1:8])
            & (arrays["tk_ql"][:, 1:8] > 0.0),
            axis=1,
        )
        good &= np.all(np.isfinite(arrays["tk_pos"][:, 1:8, :]), axis=(1, 2))
        good &= (arrays["tk_pos"][:, 4, 2] > 0.0) & (
            arrays["tk_pos"][:, 5, 2] < 0.0
        )
        good &= (
            np.isfinite(rigidity)
            & (rigidity >= R_EDGES[0])
            & (rigidity < R_EDGES[-1])
        )
        good &= (
            np.isfinite(relative_rigidity_error)
            & (relative_rigidity_error > 0.0)
            & (relative_rigidity_error < 0.20)
        )
        good &= (
            np.isfinite(chi_x)
            & np.isfinite(chi_y)
            & (chi_x > 0.0)
            & (chi_y > 0.0)
            & (chi_x <= 10.0)
            & (chi_y <= 10.0)
        )
        index = np.flatnonzero(good)

        pos_x = positions[index, :, 0] / 120.0
        pos_y = positions[index, :, 1] / 120.0
        with np.errstate(divide="ignore", invalid="ignore"):
            normalized_residual = residual[index] / span[index, None]
            log_tof_edep = np.log(
                arrays["tof_edep"][index].astype(np.float64) / (CHARGE * CHARGE)
            )
            log_tracker_dedx = 2.0 * np.log(
                arrays["tk_ql"][index, 1:8].astype(np.float64) / CHARGE
            )
            log_span = np.log(span[index] / 4.0)
        continuous = np.column_stack(
            [
                pos_x,
                pos_y,
                direction_x[index],
                direction_y[index],
                centered_path[index],
                normalized_residual,
                log_tof_edep,
                log_tracker_dedx,
                log_span,
            ]
        ).astype(np.float64)
        if continuous.shape[1] != 30:
            raise RuntimeError("Continuous feature schema mismatch: %s" % (continuous.shape,))

        mean_x = np.sum(positions[index, :, 0] / 240.0, axis=1)
        mean_y = np.sum(positions[index, :, 1] / 240.0, axis=1)
        legacy = np.column_stack(
            [
                mean_x,
                mean_y,
                mean_x * mean_x,
                mean_y * mean_y,
                log_span,
                (residual[index, 0] - residual[index, 1]) / span[index],
                (residual[index, 2] - residual[index, 3]) / span[index],
            ]
        ).astype(np.float64)

        truth_valid = np.isfinite(momentum[index]) & (momentum[index] > 0.0)
        truth_inverse_beta = np.full(len(index), np.nan, dtype=np.float64)
        truth_rigidity = np.full(len(index), np.nan, dtype=np.float64)
        truth_inverse_beta[truth_valid] = np.sqrt(
            1.0 + (mass / momentum[index][truth_valid]) ** 2
        )
        truth_rigidity[truth_valid] = momentum[index][truth_valid] / CHARGE
        version = arrays["version"][index].astype(np.int32)
        reconstructed_bin = np.searchsorted(R_EDGES, rigidity[index], side="right") - 1

        chunk = {
            "run": arrays["run"][index].astype(np.int64),
            "event": arrays["event"][index].astype(np.int64),
            "source": np.full(len(index), source_index, dtype=np.int32),
            "role": np.full(len(index), role, dtype=np.int8),
            "label": np.full(len(index), isotope_label, dtype=np.int8),
            "version": version,
            "u": inverse_beta[index].astype(np.float64),
            "R": rigidity[index].astype(np.float64),
            "relR": relative_rigidity_error[index].astype(np.float64),
            "bin": reconstructed_bin.astype(np.int8),
            "continuous": continuous,
            "legacy": legacy,
            "bars": arrays["tof_barid"][index].astype(np.int8),
            "truth_valid": truth_valid.astype(np.bool_),
            "u_true": truth_inverse_beta,
            "R_true": truth_rigidity,
        }
        if not np.all(np.isfinite(continuous)) or not np.all(np.isfinite(legacy)):
            raise RuntimeError("Non-finite selected detector features")
        chunks.append(chunk)
        bincount = np.bincount(reconstructed_bin, minlength=N_BINS)
        audit_rows.append(
            {
                "role": role,
                "isotope": ISOTOPE_NAMES[isotope_label],
                "path": path,
                "raw_entries": raw_entries,
                "selected_entries": len(index),
                "truth_entries": int(np.sum(truth_valid)),
                "version_values": ",".join(str(int(v)) for v in np.unique(version)),
                "bin_counts": ",".join(str(int(v)) for v in bincount),
                "seconds": "%.3f" % (time.time() - started),
            }
        )
        log(
            "extract role=%d %s file=%d/%d raw=%d selected=%d"
            % (
                role,
                ISOTOPE_NAMES[isotope_label],
                source_index + 1,
                len(paths),
                raw_entries,
                len(index),
            )
        )
    return concatenate(chunks)


def controlled_mix(pools, fractions, target_per_bag_bin, seed):
    rng = np.random.RandomState(seed)
    selected = []
    consumption = {
        (label, b): 0 for label in (0, 1) for b in range(N_BINS)
    }
    permutations = {}
    for label in (0, 1):
        for b in range(N_BINS):
            available = np.flatnonzero(pools[label]["bin"] == b)
            permutations[(label, b)] = available[rng.permutation(len(available))]

    audit = []
    for b in range(N_BINS):
        demand_c12 = float(np.sum(1.0 - fractions[:, b]))
        demand_c13 = float(np.sum(fractions[:, b]))
        available_c12 = len(permutations[(0, b)])
        available_c13 = len(permutations[(1, b)])
        n_common = min(
            target_per_bag_bin,
            int(
                math.floor(
                    (available_c12 - len(fractions)) / max(demand_c12, 1e-9)
                )
            ),
            int(
                math.floor(
                    (available_c13 - len(fractions)) / max(demand_c13, 1e-9)
                )
            ),
        )
        if n_common < 300:
            raise RuntimeError(
                "Too few controlled events in rigidity bin %d: %d" % (b, n_common)
            )

        for bag in range(len(fractions)):
            desired = float(fractions[bag, b])
            n13 = int(round(desired * n_common))
            n12 = n_common - n13
            pieces = []
            for label, count in ((0, n12), (1, n13)):
                start = consumption[(label, b)]
                stop = start + count
                rows = permutations[(label, b)][start:stop]
                if len(rows) != count:
                    raise RuntimeError("Controlled sampling exhausted")
                consumption[(label, b)] = stop
                piece = subset(pools[label], rows)
                piece["bag"] = np.full(count, bag, dtype=np.int8)
                piece["prior_c13"] = np.full(
                    count, n13 / float(n_common), dtype=np.float64
                )
                pieces.append(piece)
            bag_data = concatenate(pieces)
            order = rng.permutation(len(bag_data["u"]))
            selected.append(subset(bag_data, order))
            audit.append(
                {
                    "bin": b,
                    "bag": bag,
                    "n_total": n_common,
                    "n_c12": n12,
                    "n_c13": n13,
                    "fraction_c13": n13 / float(n_common),
                }
            )

    data = concatenate(selected)
    order = rng.permutation(len(data["u"]))
    return subset(data, order), audit


def make_feature_state(data, mode):
    if mode == "geometry":
        raw = data["legacy"].astype(np.float64)
        names = list(LEGACY_NAMES)
        low, high = np.percentile(raw, [1.0, 99.0], axis=0)
        clipped = np.minimum(np.maximum(raw, low), high)
        center = np.median(clipped, axis=0)
        quartiles = np.percentile(clipped, [25.0, 75.0], axis=0)
        scale = np.maximum(1e-6, 0.741301109 * (quartiles[1] - quartiles[0]))
        return {
            "mode": mode,
            "low": low,
            "high": high,
            "center": center,
            "scale": scale,
            "feature_names": ["intercept"] + names,
        }

    selected_indices = FULL_INDICES if mode == "full" else NO_EDEP_INDICES
    interactions = FULL_INTERACTIONS if mode == "full" else NO_EDEP_INTERACTIONS
    raw = data["continuous"][:, selected_indices].astype(np.float64)
    low, high = np.percentile(raw, [1.0, 99.0], axis=0)
    clipped = np.minimum(np.maximum(raw, low), high)
    center = np.median(clipped, axis=0)
    quartiles = np.percentile(clipped, [25.0, 75.0], axis=0)
    scale = np.maximum(1e-6, 0.741301109 * (quartiles[1] - quartiles[0]))
    local = {original: position for position, original in enumerate(selected_indices)}
    local_pairs = [(local[a], local[b]) for a, b in interactions]
    z = np.clip((clipped - center) / scale, -4.0, 4.0)
    columns = [z, z * z]
    names = [CONTINUOUS_NAMES[i] for i in selected_indices] + [
        "square_" + CONTINUOUS_NAMES[i] for i in selected_indices
    ]
    if local_pairs:
        products = np.column_stack([z[:, a] * z[:, b] for a, b in local_pairs])
        columns.append(products)
        names += [
            "interaction_%s__%s"
            % (CONTINUOUS_NAMES[interactions[k][0]], CONTINUOUS_NAMES[interactions[k][1]])
            for k in range(len(interactions))
        ]
    one_hot_columns = []
    for layer in range(4):
        for bar in range(1, 10):
            one_hot_columns.append(
                (data["bars"][:, layer] == bar).astype(np.float64)
            )
            names.append("bar_layer%d_%d" % (layer, bar))
    columns.append(np.column_stack(one_hot_columns))
    derived = np.column_stack(columns)
    derived_center = np.mean(derived, axis=0)
    derived_scale = np.std(derived, axis=0)
    derived_scale[derived_scale < 1e-6] = 1.0
    return {
        "mode": mode,
        "selected_indices": np.asarray(selected_indices, dtype=np.int32),
        "interactions": np.asarray(interactions, dtype=np.int32),
        "low": low,
        "high": high,
        "center": center,
        "scale": scale,
        "derived_center": derived_center,
        "derived_scale": derived_scale,
        "feature_names": ["intercept"] + names,
    }


def apply_feature_state(data, state):
    mode = state["mode"]
    if mode == "geometry":
        raw = data["legacy"].astype(np.float64)
        clipped = np.minimum(np.maximum(raw, state["low"]), state["high"])
        z = np.clip((clipped - state["center"]) / state["scale"], -4.0, 4.0)
        return np.column_stack([np.ones(len(z)), z])

    selected_indices = state["selected_indices"].astype(np.int32)
    raw = data["continuous"][:, selected_indices].astype(np.float64)
    clipped = np.minimum(np.maximum(raw, state["low"]), state["high"])
    z = np.clip((clipped - state["center"]) / state["scale"], -4.0, 4.0)
    columns = [z, z * z]
    interactions = state["interactions"].astype(np.int32)
    local = {int(original): position for position, original in enumerate(selected_indices)}
    if len(interactions):
        columns.append(
            np.column_stack(
                [
                    z[:, local[int(a)]] * z[:, local[int(b)]]
                    for a, b in interactions
                ]
            )
        )
    one_hot_columns = []
    for layer in range(4):
        for bar in range(1, 10):
            one_hot_columns.append(
                (data["bars"][:, layer] == bar).astype(np.float64)
            )
    columns.append(np.column_stack(one_hot_columns))
    derived = np.column_stack(columns)
    standardized = np.clip(
        (derived - state["derived_center"]) / state["derived_scale"], -5.0, 5.0
    )
    return np.column_stack([np.ones(len(standardized)), standardized])



class MixtureDensity(object):
    def __init__(
        self,
        data,
        x,
        prior_c13,
        ridge_strength=5e-4,
        coefficient_scale=0.02,
        prediction_strength=2e-3,
        prediction_scale=0.03,
    ):
        self.x = np.asarray(x, dtype=np.float64)
        self.u = data["u"].astype(np.float64)
        self.b = data["bin"].astype(np.int32)
        self.prior = np.asarray(prior_c13, dtype=np.float64)
        self.means = np.sqrt(
            1.0 + (MASSES[None, :] / (CHARGE * data["R"][:, None])) ** 2
        )
        sensitivity = (self.means * self.means - 1.0) / self.means
        self.tracker = (
            ((NU - 2.0) / NU)
            * sensitivity
            * sensitivity
            * data["relR"][:, None] ** 2
        )
        self.p = self.x.shape[1]
        self.ridge_strength = ridge_strength
        self.coefficient_scale = coefficient_scale
        self.prediction_strength = prediction_strength
        self.prediction_scale = prediction_scale

    def components(self, parameters, corrected=None):
        theta = parameters[: self.p]
        logsigma = parameters[self.p : self.p + N_BINS]
        prediction = self.x.dot(theta)
        u_corrected = self.u - prediction if corrected is None else corrected
        residual = u_corrected[:, None] - self.means
        sigma2_by_bin = np.exp(2.0 * logsigma)
        sigma2 = sigma2_by_bin[self.b]
        variance = sigma2[:, None] + self.tracker
        logpdf = (
            T_NORMALIZATION
            - 0.5 * np.log(variance)
            - ((NU + 1.0) / 2.0)
            * np.log1p(residual * residual / (NU * variance))
        )
        return logpdf, residual, variance, sigma2, prediction

    def value_gradient(self, parameters):
        logpdf, residual, variance, sigma2, prediction = self.components(parameters)
        with np.errstate(divide="ignore"):
            log_mixture = np.logaddexp(
                np.log1p(-self.prior) + logpdf[:, 0],
                np.log(self.prior) + logpdf[:, 1],
            )
        responsibilities = np.exp(logpdf - log_mixture[:, None])
        weights = responsibilities * np.column_stack(
            [1.0 - self.prior, self.prior]
        )
        denominator = NU * variance + residual * residual
        correction_gradient = -np.sum(
            weights * (NU + 1.0) * residual / denominator, axis=1
        )
        sigma_gradient = np.sum(
            weights
            * (sigma2[:, None] / variance)
            * (1.0 - (NU + 1.0) * residual * residual / denominator),
            axis=1,
        )
        n = float(len(self.u))
        gradient_theta = self.x.T.dot(correction_gradient) / n
        gradient_sigma = (
            np.bincount(self.b, weights=sigma_gradient, minlength=N_BINS) / n
        )
        value = -float(np.mean(log_mixture))

        if self.p > 1:
            scaled = parameters[1 : self.p] / self.coefficient_scale
            value += 0.5 * self.ridge_strength * float(np.sum(scaled * scaled))
            gradient_theta[1:] += (
                self.ridge_strength
                * parameters[1 : self.p]
                / (self.coefficient_scale ** 2)
            )
        scaled_prediction = prediction / self.prediction_scale
        value += 0.5 * self.prediction_strength * float(
            np.mean(scaled_prediction * scaled_prediction)
        )
        gradient_theta += (
            self.prediction_strength
            / (self.prediction_scale ** 2)
            * self.x.T.dot(prediction)
            / n
        )
        return value, np.concatenate([gradient_theta, gradient_sigma])

    def event_nll(self, parameters, corrected=None):
        logpdf = self.components(parameters, corrected=corrected)[0]
        with np.errstate(divide="ignore"):
            return -np.logaddexp(
                np.log1p(-self.prior) + logpdf[:, 0],
                np.log(self.prior) + logpdf[:, 1],
            )


def fit_density(
    name,
    data,
    x,
    prior,
    fit_theta,
    initial_logsigma=None,
    starts=1,
):
    density = MixtureDensity(data, x, prior)
    p_features = x.shape[1]
    if initial_logsigma is None:
        initial_logsigma = np.full(N_BINS, math.log(0.010), dtype=float)
    base = np.r_[np.zeros(p_features), np.asarray(initial_logsigma, dtype=float)]
    bounds = [(-0.05, 0.05)] * p_features + [
        (math.log(0.0015), math.log(0.08))
    ] * N_BINS
    if not fit_theta:
        for j in range(p_features):
            bounds[j] = (0.0, 0.0)

    solutions, records = [], []
    for start in range(starts):
        initial = base.copy()
        if start == 1 and fit_theta:
            initial[0] = 0.003
            initial[p_features:] += math.log(1.15)
        result = minimize(
            density.value_gradient,
            initial,
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": 1200,
                "maxls": 60,
                "ftol": 1e-12,
                "gtol": 2e-6,
            },
        )
        record = {
            "start": start,
            "success": bool(result.success),
            "status": int(result.status),
            "objective": float(result.fun),
            "iterations": int(result.nit),
            "message": str(result.message),
        }
        records.append(record)
        log("fit %s %s" % (name, record))
        if (
            np.isfinite(result.fun)
            and np.all(np.isfinite(result.x))
            and (result.success or result.status == 1)
        ):
            solutions.append(result)
    if not solutions:
        raise RuntimeError("All optimizer starts failed for " + name)

    best = min(solutions, key=lambda item: item.fun)
    correction = x.dot(best.x[:p_features])
    return {
        "name": name,
        "parameters": best.x,
        "theta": best.x[:p_features],
        "logsigma": best.x[p_features:],
        "objective": float(best.fun),
        "optimizer": records,
        "correction_rms_training": rms(correction),
        "correction_p99_abs_training": float(
            np.percentile(np.abs(correction), 99.0)
        ),
    }


def gradient_check(data, x, prior):
    rows = np.arange(min(500, len(data["u"])))
    small = subset(data, rows)
    density = MixtureDensity(small, x[rows], prior[rows])
    parameters = np.r_[
        np.zeros(x.shape[1]), np.full(N_BINS, math.log(0.012))
    ]
    parameters[0] = 0.001
    value, analytic = density.value_gradient(parameters)
    indices = list(range(min(6, x.shape[1]))) + list(
        range(x.shape[1], x.shape[1] + N_BINS)
    )
    errors = []
    for index in indices:
        plus, minus = parameters.copy(), parameters.copy()
        step = 1e-6
        plus[index] += step
        minus[index] -= step
        numeric = (
            density.value_gradient(plus)[0]
            - density.value_gradient(minus)[0]
        ) / (2.0 * step)
        error = abs(analytic[index] - numeric) / (1.0 + abs(numeric))
        errors.append(error)
    maximum = float(max(errors))
    log(
        "gradient check value=%.8g max_relative_error=%.3g"
        % (value, maximum)
    )
    if maximum > 3e-5:
        raise RuntimeError("Gradient check failed")
    return maximum


def fit_oracle(name, data, x, regularization=20.0):
    truth = data["truth_valid"]
    xt = x[truth]
    target = data["u"][truth] - data["u_true"][truth]
    if len(target) < x.shape[1] * 5:
        log(
            "warning: oracle %s has only %d truth rows for %d features"
            % (name, len(target), x.shape[1])
        )
    penalty = np.full(x.shape[1], regularization, dtype=float)
    penalty[0] = 1e-6
    theta = np.linalg.solve(
        xt.T.dot(xt) + np.diag(penalty), xt.T.dot(target)
    )
    for _ in range(20):
        residual = target - xt.dot(theta)
        delta = max(0.004, 2.0 * w68(residual))
        weight = np.minimum(
            1.0, delta / np.maximum(np.abs(residual), 1e-12)
        )
        weighted_x = xt * weight[:, None]
        updated = np.linalg.solve(
            xt.T.dot(weighted_x) + np.diag(penalty),
            xt.T.dot(weight * target),
        )
        if np.max(np.abs(updated - theta)) < 1e-9:
            theta = updated
            break
        theta = updated
    correction = x.dot(theta)
    log(
        "oracle %s truth_rows=%d correction_rms=%.6g max_coeff=%.6g"
        % (name, len(target), rms(correction), np.max(np.abs(theta)))
    )
    return {
        "name": name,
        "theta": theta,
        "truth_rows": int(len(target)),
        "regularization": regularization,
        "correction_rms_training": rms(correction),
    }


def fit_logistic(x, label):
    label = label.astype(np.float64)

    def loss(parameters):
        score = x.dot(parameters)
        value = float(
            np.mean(np.logaddexp(0.0, score) - label * score)
        )
        probability = 1.0 / (
            1.0 + np.exp(-np.clip(score, -40.0, 40.0))
        )
        gradient = x.T.dot(probability - label) / len(label)
        value += 0.5e-3 * float(np.sum(parameters[1:] ** 2))
        gradient[1:] += 1e-3 * parameters[1:]
        return value, gradient

    result = minimize(
        loss,
        np.zeros(x.shape[1]),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 600, "ftol": 1e-11, "gtol": 1e-6},
    )
    if not result.success:
        log(
            "warning: logistic diagnostic did not fully converge: %s"
            % result.message
        )
    return result.x


def method_metrics(method, corrected_u, valid, heldout_nll):
    rows = []
    truth = (
        valid["truth_valid"]
        & np.isfinite(corrected_u)
        & (corrected_u > 0.0)
    )
    mass2 = (CHARGE * valid["R"]) ** 2 * (
        corrected_u * corrected_u - 1.0
    )
    correction = valid["u"] - corrected_u

    groups = [
        ("all", np.ones(len(corrected_u), dtype=bool)),
        ("C12", valid["label"] == 0),
        ("C13", valid["label"] == 1),
    ]
    for group_name, group_mask in groups:
        select = truth & group_mask
        du = corrected_u[select] - valid["u_true"][select]
        dbeta = (
            1.0 / corrected_u[select]
            - 1.0 / valid["u_true"][select]
        )
        rows.append(
            {
                "method": method,
                "group": group_name,
                "n": int(np.sum(select)),
                "median_du": float(np.median(du)),
                "w68_du": w68(du),
                "rms_du": rms(du),
                "median_dbeta": float(np.median(dbeta)),
                "w68_dbeta": w68(dbeta),
                "median_mass2": float(np.median(mass2[select])),
                "w68_mass2": w68(mass2[select]),
                "correction_rms": rms(correction[select]),
                "correction_median": float(np.median(correction[select])),
                "heldout_nll": (
                    float(heldout_nll)
                    if group_name == "all"
                    else float("nan")
                ),
                "mass2_auc": float("nan"),
                "mass2_separation": float("nan"),
                "improvement_vs_raw_percent": float("nan"),
            }
        )

    select = truth
    rows[0]["mass2_auc"] = auc(
        mass2[select], valid["label"][select]
    )
    light = select & (valid["label"] == 0)
    heavy = select & (valid["label"] == 1)
    denominator = math.sqrt(
        w68(mass2[light]) ** 2 + w68(mass2[heavy]) ** 2
    )
    rows[0]["mass2_separation"] = (
        (
            float(np.median(mass2[heavy]))
            - float(np.median(mass2[light]))
        )
        / denominator
        if denominator > 0
        else float("nan")
    )
    return rows, mass2


def bin_metrics(method, corrected_u, mass2, valid):
    rows = []
    truth = (
        valid["truth_valid"]
        & np.isfinite(corrected_u)
        & (corrected_u > 0.0)
    )
    for b in range(N_BINS):
        in_bin = truth & (valid["bin"] == b)
        groups = [
            ("all", np.ones(len(corrected_u), dtype=bool)),
            ("C12", valid["label"] == 0),
            ("C13", valid["label"] == 1),
        ]
        for group_name, group_mask in groups:
            select = in_bin & group_mask
            du = corrected_u[select] - valid["u_true"][select]
            rows.append(
                {
                    "method": method,
                    "bin": b,
                    "r_low": R_EDGES[b],
                    "r_high": R_EDGES[b + 1],
                    "group": group_name,
                    "n": int(np.sum(select)),
                    "median_du": (
                        float(np.median(du))
                        if len(du)
                        else float("nan")
                    ),
                    "w68_du": w68(du),
                    "median_mass2": (
                        float(np.median(mass2[select]))
                        if np.sum(select)
                        else float("nan")
                    ),
                    "w68_mass2": w68(mass2[select]),
                    "mass2_auc": (
                        auc(mass2[select], valid["label"][select])
                        if group_name == "all"
                        else float("nan")
                    ),
                }
            )
    return rows


def write_tsv(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(str(path), "w") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=keys,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)



def main():
    destination = Path(__file__).resolve().parent
    started = time.time()

    c12_files = choose_files(C12_PATTERN, FILES_PER_ROLE)
    c13_files = choose_files(C13_PATTERN, C13_FILES_PER_ROLE)
    manifest_rows = []
    for role in (0, 1):
        for isotope, files in [
            ("C12", c12_files[role]),
            ("C13", c13_files[role]),
        ]:
            for path in files:
                manifest_rows.append(
                    {"role": role, "isotope": isotope, "path": path}
                )
    write_tsv(destination / "manifest.tsv", manifest_rows)

    extraction_audit = []
    pools = {0: {}, 1: {}}
    for role in (0, 1):
        pools[role][0] = extract_pool(
            c12_files[role],
            0,
            67,
            MASSES[0],
            role,
            extraction_audit,
        )
        pools[role][1] = extract_pool(
            c13_files[role],
            1,
            117,
            MASSES[1],
            role,
            extraction_audit,
        )
    write_tsv(destination / "extraction_audit.tsv", extraction_audit)

    train, train_mix_audit = controlled_mix(
        pools[0],
        TRAIN_FRACTIONS,
        TRAIN_TARGET_PER_BAG_BIN,
        SEED,
    )
    valid, valid_mix_audit = controlled_mix(
        pools[1],
        VALID_FRACTIONS,
        VALID_TARGET_PER_BIN,
        SEED + 1,
    )
    mix_rows = []
    for role_name, rows in [
        ("train", train_mix_audit),
        ("validation", valid_mix_audit),
    ]:
        for row in rows:
            item = dict(row)
            item["role"] = role_name
            mix_rows.append(item)
    write_tsv(destination / "mixture_audit.tsv", mix_rows)
    log(
        "controlled samples train=%d validation=%d"
        % (len(train["u"]), len(valid["u"]))
    )

    train_ids = set(
        zip(
            train["run"].tolist(),
            train["event"].tolist(),
            train["label"].tolist(),
        )
    )
    valid_ids = set(
        zip(
            valid["run"].tolist(),
            valid["event"].tolist(),
            valid["label"].tolist(),
        )
    )
    if len(train_ids) != len(train["u"]) or len(valid_ids) != len(valid["u"]):
        raise RuntimeError("Duplicate selected event")
    if train_ids.intersection(valid_ids):
        raise RuntimeError("Train/validation event overlap")

    states = {}
    matrices_train, matrices_valid = {}, {}
    for mode in ("geometry", "no_edep", "full"):
        state = make_feature_state(train, mode)
        states[mode] = state
        matrices_train[mode] = apply_feature_state(train, state)
        matrices_valid[mode] = apply_feature_state(valid, state)
        log(
            "feature map %s dimensions=%d"
            % (mode, matrices_train[mode].shape[1])
        )

    # Explicit allowlist: weak objectives cannot access labels or MC truth.
    weak_train = {
        key: train[key] for key in ["u", "R", "relR", "bin"]
    }
    weak_valid = {
        key: valid[key] for key in ["u", "R", "relR", "bin"]
    }
    gradient_error = gradient_check(
        weak_train,
        matrices_train["full"],
        train["prior_c13"],
    )

    raw_model = fit_density(
        "raw",
        weak_train,
        matrices_train["geometry"],
        train["prior_c13"],
        False,
    )
    initial_sigma = raw_model["logsigma"]
    geometry_model = fit_density(
        "weak_geometry",
        weak_train,
        matrices_train["geometry"],
        train["prior_c13"],
        True,
        initial_sigma,
        starts=2,
    )
    no_edep_model = fit_density(
        "weak_no_edep",
        weak_train,
        matrices_train["no_edep"],
        train["prior_c13"],
        True,
        initial_sigma,
        starts=2,
    )
    full_model = fit_density(
        "weak_full",
        weak_train,
        matrices_train["full"],
        train["prior_c13"],
        True,
        initial_sigma,
        starts=2,
    )
    wrong_prior = np.clip(
        train["prior_c13"] + 0.10, 0.01, 0.99
    )
    wrong_model = fit_density(
        "weak_wrong_fraction",
        weak_train,
        matrices_train["full"],
        wrong_prior,
        True,
        full_model["logsigma"],
        starts=1,
    )

    weak_models = {
        "raw": raw_model,
        "weak_geometry": geometry_model,
        "weak_no_edep": no_edep_model,
        "weak_full": full_model,
        "weak_wrong_fraction": wrong_model,
    }
    frozen = {
        "truth_policy": (
            "Per-event isotope labels and mevmom1[:,10] were excluded "
            "from every weak objective. Labels only constructed bags "
            "with known fractions. Truth arrays were first consulted "
            "after this weak-model file was written."
        ),
        "masses_GeV": MASSES,
        "rigidity_edges_GV": R_EDGES,
        "training_fractions": TRAIN_FRACTIONS,
        "validation_fractions": VALID_FRACTIONS,
        "gradient_check_max_relative_error": gradient_error,
        "feature_states": states,
        "models": weak_models,
        "script_sha256": sha256_file(Path(__file__)),
    }
    (destination / "weak_frozen_models.json").write_text(
        json.dumps(jsonable(frozen), indent=2, sort_keys=True)
    )
    log("weak models frozen before truth evaluation")

    # Truth is consulted from this point onward only.
    oracle_models = {}
    for mode, regularization in [
        ("geometry", 5.0),
        ("no_edep", 15.0),
        ("full", 20.0),
    ]:
        oracle_models["oracle_" + mode] = fit_oracle(
            "oracle_" + mode,
            train,
            matrices_train[mode],
            regularization,
        )
    (destination / "oracle_models.json").write_text(
        json.dumps(jsonable(oracle_models), indent=2, sort_keys=True)
    )

    # Post-freeze source/isotope-label diagnostics.
    diagnostic = {}
    for mode in ("geometry", "no_edep", "full"):
        classifier = fit_logistic(
            matrices_train[mode], train["label"]
        )
        diagnostic[mode + "_label_auc"] = auc(
            matrices_valid[mode].dot(classifier),
            valid["label"],
        )
    log("post-freeze label AUC diagnostics %s" % diagnostic)

    corrected = {"raw": valid["u"].copy()}
    corrected["weak_geometry"] = (
        valid["u"]
        - matrices_valid["geometry"].dot(geometry_model["theta"])
    )
    corrected["weak_no_edep"] = (
        valid["u"]
        - matrices_valid["no_edep"].dot(no_edep_model["theta"])
    )
    corrected["weak_full"] = (
        valid["u"] - matrices_valid["full"].dot(full_model["theta"])
    )
    corrected["weak_wrong_fraction"] = (
        valid["u"] - matrices_valid["full"].dot(wrong_model["theta"])
    )
    corrected["oracle_geometry"] = (
        valid["u"]
        - matrices_valid["geometry"].dot(
            oracle_models["oracle_geometry"]["theta"]
        )
    )
    corrected["oracle_no_edep"] = (
        valid["u"]
        - matrices_valid["no_edep"].dot(
            oracle_models["oracle_no_edep"]["theta"]
        )
    )
    corrected["oracle_full"] = (
        valid["u"]
        - matrices_valid["full"].dot(
            oracle_models["oracle_full"]["theta"]
        )
    )

    nominal_densities = {
        "geometry": MixtureDensity(
            weak_valid,
            matrices_valid["geometry"],
            valid["prior_c13"],
        ),
        "no_edep": MixtureDensity(
            weak_valid,
            matrices_valid["no_edep"],
            valid["prior_c13"],
        ),
        "full": MixtureDensity(
            weak_valid,
            matrices_valid["full"],
            valid["prior_c13"],
        ),
    }
    heldout_nll = {
        "raw": float(
            np.mean(
                nominal_densities["geometry"].event_nll(
                    raw_model["parameters"]
                )
            )
        ),
        "weak_geometry": float(
            np.mean(
                nominal_densities["geometry"].event_nll(
                    geometry_model["parameters"]
                )
            )
        ),
        "weak_no_edep": float(
            np.mean(
                nominal_densities["no_edep"].event_nll(
                    no_edep_model["parameters"]
                )
            )
        ),
        "weak_full": float(
            np.mean(
                nominal_densities["full"].event_nll(
                    full_model["parameters"]
                )
            )
        ),
        "weak_wrong_fraction": float(
            np.mean(
                nominal_densities["full"].event_nll(
                    wrong_model["parameters"]
                )
            )
        ),
    }
    # Oracle NLL uses the frozen weak-full response nuisance only.
    for name in (
        "oracle_geometry",
        "oracle_no_edep",
        "oracle_full",
    ):
        heldout_nll[name] = float(
            np.mean(
                nominal_densities["full"].event_nll(
                    full_model["parameters"],
                    corrected=corrected[name],
                )
            )
        )

    metric_rows, per_bin_rows, mass2_by_method = [], [], {}
    order = [
        "raw",
        "weak_geometry",
        "weak_no_edep",
        "weak_full",
        "weak_wrong_fraction",
        "oracle_geometry",
        "oracle_no_edep",
        "oracle_full",
    ]
    for name in order:
        rows, mass2 = method_metrics(
            name,
            corrected[name],
            valid,
            heldout_nll[name],
        )
        metric_rows.extend(rows)
        per_bin_rows.extend(
            bin_metrics(name, corrected[name], mass2, valid)
        )
        mass2_by_method[name] = mass2

    raw_width = {
        row["group"]: row["w68_du"]
        for row in metric_rows
        if row["method"] == "raw"
    }
    for row in metric_rows:
        baseline = raw_width[row["group"]]
        row["improvement_vs_raw_percent"] = 100.0 * (
            1.0 - row["w68_du"] / baseline
        )
    write_tsv(destination / "metrics.tsv", metric_rows)
    write_tsv(destination / "bin_metrics.tsv", per_bin_rows)

    predictions = {
        "run": valid["run"],
        "event": valid["event"],
        "label": valid["label"],
        "version": valid["version"],
        "bin": valid["bin"],
        "R": valid["R"],
        "relR": valid["relR"],
        "truth_valid": valid["truth_valid"],
        "u_true": valid["u_true"],
        "R_true": valid["R_true"],
        "prior_c13": valid["prior_c13"],
    }
    for name, values in corrected.items():
        predictions["u_" + name] = values
        predictions["mass2_" + name] = mass2_by_method[name]
    np.savez_compressed(
        str(destination / "predictions.npz"), **predictions
    )

    all_rows = {
        row["method"]: row
        for row in metric_rows
        if row["group"] == "all"
    }
    denominator = (
        raw_width["all"] - all_rows["oracle_full"]["w68_du"]
    )
    summary = {
        "train_events": int(len(train["u"])),
        "validation_events": int(len(valid["u"])),
        "train_truth_events": int(np.sum(train["truth_valid"])),
        "validation_truth_events": int(
            np.sum(valid["truth_valid"])
        ),
        "versions": {
            "C12": sorted(
                set(
                    int(v)
                    for v in np.unique(
                        np.r_[
                            pools[0][0]["version"],
                            pools[1][0]["version"],
                        ]
                    )
                )
            ),
            "C13": sorted(
                set(
                    int(v)
                    for v in np.unique(
                        np.r_[
                            pools[0][1]["version"],
                            pools[1][1]["version"],
                        ]
                    )
                )
            ),
        },
        "label_diagnostics": diagnostic,
        "overall_metrics": all_rows,
        "weak_full_recovery_fraction_of_oracle": (
            (
                raw_width["all"]
                - all_rows["weak_full"]["w68_du"]
            )
            / denominator
            if denominator != 0.0
            else float("nan")
        ),
        "elapsed_seconds": time.time() - started,
        "script_sha256": sha256_file(Path(__file__)),
    }
    (destination / "results.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )

    readme = """# C12/C13 isotope-proportion closure

This directory is an MC closure for the mass-peak/known-proportion route. It does not modify the production estimator.

- Inputs: pure C12 B1308 and pure C13 B1315 amsd69n ntuples, split by whole source file.
- Weak supervision: two synthetic training bags per reconstructed-rigidity bin with fixed known C13 fractions. Per-event isotope labels and `mevmom1[:,10]` are absent from every weak objective.
- Rigidity contract: `tk_rigidityk[1,1]`, with `tk_erigidity1[2,1]` and `tk_chis1[2,2,1]`, matching the earlier mass-study definition.
- Target for closure only: inverse beta from MC momentum at checkpoint 10 and the GEANT masses C12=11.17793 GeV, C13=12.11221 GeV.
- Full detector feature map: TOF geometry/path/orthogonal timing residuals, four TOF energy deposits, seven Tracker `tk_ql` dE/dx proxies, and per-layer TOF bar identities. No rigidity value, baseline beta, absolute time, isotope label, MC momentum, or production version enters the correction.
- Controls: raw, legacy 8-dimensional geometry weak fit, high-capacity no-energy-loss weak fit, full weak fit, +0.10 wrong-fraction weak fit, and truth-supervised oracle fits with the same feature maps.

Important limitation: C12 is reconstruction version B1308 while C13 is B1315. The post-freeze label diagnostics quantify how readily detector features distinguish the two source classes, but cannot fully separate physical isotope response from production-version response. A same-production C12/C13 pair is required before transfer claims.

Primary outputs: `metrics.tsv`, `bin_metrics.tsv`, `results.json`, `weak_frozen_models.json`, `oracle_models.json`, `predictions.npz`, and the source/mixture audit tables.

Reproduce on sdiat:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 results/c_isotope_proportion_closure/run_closure.py
```
"""
    (destination / "README.md").write_text(readme)
    log(
        "DONE elapsed=%.1fs results=%s"
        % (time.time() - started, destination / "results.json")
    )


if __name__ == "__main__":
    main()
