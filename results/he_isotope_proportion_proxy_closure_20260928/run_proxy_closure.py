#!/usr/bin/env python3
"""Phase-0 He isotope proportion-constrained beta closure.

This deliberately uses a He4-derived kinematic He3 proxy. It tests only whether
an externally supplied selected isotope fraction can identify a shared, low-
dimensional TOF inverse-beta calibration. It is not a replacement for a full
He3 detector simulation.

Truth branches are used only to construct the proxy and to evaluate frozen
models. The fit objective sees only reconstructed R, TOF-derived u, the eight
pre-existing geometry/timing-residual features, tracker uncertainty, known
isotope masses, and a bag-level isotope proportion.
"""
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

import numpy as np
import uproot
from scipy.optimize import minimize

C_CM_PER_NS = 29.9792458
M3 = 2.80839160743
M4 = 3.7273794066
MASSES = np.array([M3, M4], dtype=float)  # components: He3, He4
EDGES = np.array([2.5, 3.0, 4.0, 5.0, 6.0], dtype=float)
B = len(EDGES) - 1
P = 8
NU = 6.0
LOG_NORM = math.lgamma((NU + 1.0) / 2.0) - math.lgamma(NU / 2.0) - 0.5 * math.log(NU * math.pi)
BASE = Path('/lustre02/user/yachen/public/hxwu/amsd69n_mc_z1_z2_z6_z8_20260903/Z2_He4')
MASS_STUDY = Path('/nas01/home/hxwu/AMSSoft/betaEstimator/results/mass_study')
AMS_TABLE = Path('/tmp/ams_he_isotope_tables_smi_smxxi.csv')
OUT = Path('/tmp/beta_estimator_he_proxy_closure_20260928')
N_FILES = 40
MAX_PER_SPLIT = 12000
SEED = 20260928

BRANCHES = [
    'run', 'event', 'mmom', 'mmass', 'mch',
    'tof_tl', 'tof_leng', 'tof_pos', 'tof_edep', 'tof_barid', 'tof_bpatt',
    'tk_qin', 'tk_q', 'tk_edep', 'tk_pos', 'tk_rigidity',
    'tk_rigidityk', 'tk_erigidity1', 'tk_chis1',
]


def w68(values: np.ndarray) -> float:
    if len(values) == 0:
        return float('nan')
    lo, hi = np.percentile(values, [16.0, 84.0])
    return float((hi - lo) / 2.0)


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    labels = labels.astype(np.int8)
    n1 = int(np.sum(labels == 1))
    n0 = int(np.sum(labels == 0))
    if n0 == 0 or n1 == 0:
        return float('nan')
    order = np.argsort(scores, kind='mergesort')
    s = scores[order]
    y = labels[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    ends = np.r_[starts[1:], len(s)]
    ranks = np.repeat((starts + 1 + ends) / 2.0, ends - starts)
    return float((np.sum(ranks[y == 1]) - n1 * (n1 + 1) / 2.0) / (n0 * n1))


def flight_bin_fractions() -> np.ndarray:
    counts = np.zeros(B, dtype=float)
    with (MASS_STUDY / 'bins.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if int(row['Z']) != 2 or row['method'] != 'raw':
                continue
            lo = float(row['R_low'])
            idx = int(np.where(np.isclose(EDGES[:-1], lo))[0][0])
            counts[idx] = float(row['N_holdout'])
    if np.any(counts <= 0):
        raise RuntimeError(f'Invalid flight bin counts: {counts}')
    return counts / counts.sum()


def official_f3_by_bin() -> Tuple[np.ndarray, List[dict]]:
    if not AMS_TABLE.is_file():
        raise FileNotFoundError(AMS_TABLE)
    rows: List[dict] = []
    with AMS_TABLE.open(newline='') as handle:
        for row in csv.DictReader(handle):
            if int(row['bartels_rotation_start']) != 2426 or int(row['bartels_rotation_end']) != 2429:
                continue
            try:
                rows.append({
                    'lo': float(row['rigidity_min GV']),
                    'hi': float(row['rigidity_max GV']),
                    'phi3': float(row['helium3_flux m^-2sr^-1s^-1GV^-1']),
                    'phi4': float(row['helium4_flux m^-2sr^-1s^-1GV^-1']),
                })
            except (TypeError, ValueError):
                continue
    fractions = np.zeros(B, dtype=float)
    provenance: List[dict] = []
    for b, (blo, bhi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
        i3 = 0.0
        i4 = 0.0
        pieces = []
        for row in rows:
            width = max(0.0, min(bhi, row['hi']) - max(blo, row['lo']))
            if width <= 0:
                continue
            i3 += row['phi3'] * width
            i4 += row['phi4'] * width
            pieces.append(dict(row, overlap=width))
        if i3 <= 0 or i4 <= 0:
            raise RuntimeError(f'Missing official flux for bin {blo}-{bhi}')
        fractions[b] = i3 / (i3 + i4)
        provenance.append({'R_low': blo, 'R_high': bhi, 'integral_phi3': i3,
                           'integral_phi4': i4, 'f3': fractions[b], 'pieces': pieces})
    return fractions, provenance


def extract(files: List[Path], label: str) -> Dict[str, np.ndarray]:
    chunks: List[Dict[str, np.ndarray]] = []
    specs = [str(path) + ':amstreea' for path in files]
    raw_total = 0
    kept_total = 0
    for arrays in uproot.iterate(specs, BRANCHES, step_size='120 MB', library='np'):
        n = len(arrays['mmom'])
        raw_total += n
        qin = arrays['tk_qin'][:, 0, 2]
        q = np.where(qin < 2.5, arrays['tk_q'][:, 1], qin)
        z = np.floor(q + 0.5).astype(np.int16)
        good = np.isfinite(q) & (z == 2)
        r0 = arrays['tk_rigidity']
        good &= np.isfinite(r0) & (r0 >= 1.0) & (r0 <= 12.0) & (arrays['tof_bpatt'] == 4444)
        good &= np.all(np.isfinite(arrays['tk_edep'][:, 1:8]) & (arrays['tk_edep'][:, 1:8] > 0), axis=1)
        tk_pos = arrays['tk_pos']
        good &= np.all(np.isfinite(tk_pos[:, 1:8, :]), axis=(1, 2))
        good &= (tk_pos[:, 4, 2] > 0) & (tk_pos[:, 5, 2] < 0)

        tl = arrays['tof_tl'].astype(float)
        leng = arrays['tof_leng'].astype(float)
        tof_pos = arrays['tof_pos'].astype(float)
        edep = arrays['tof_edep']
        barid = arrays['tof_barid']
        tof_good = np.all(np.isfinite(tl), axis=1)
        tof_good &= np.all(np.isfinite(leng), axis=1)
        tof_good &= np.all(np.isfinite(tof_pos[:, :, :2]), axis=(1, 2))
        tof_good &= np.all(np.isfinite(edep) & (edep > 0), axis=1)
        tof_good &= np.all((barid >= 0) & (barid < 10), axis=1)
        good &= tof_good

        path_ns = leng / C_CM_PER_NS
        xc = path_ns - np.mean(path_ns, axis=1, keepdims=True)
        tc = tl - np.mean(tl, axis=1, keepdims=True)
        den = np.sum(xc * xc, axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            u = np.sum(xc * tc, axis=1) / den
        span = np.max(xc, axis=1) - np.min(xc, axis=1)
        good &= np.isfinite(u) & (u > 0) & np.isfinite(den) & (den > 1e-12) & (span > 0)

        r_rec = arrays['tk_rigidityk'][:, 1, 1].astype(float)
        rel_r = np.abs(arrays['tk_erigidity1'][:, 2, 1].astype(float) * r_rec)
        cx = arrays['tk_chis1'][:, 2, 2, 1, 0].astype(float)
        cy = arrays['tk_chis1'][:, 2, 2, 1, 1].astype(float)
        good &= np.isfinite(r_rec) & (r_rec >= EDGES[0]) & (r_rec < EDGES[-1])
        good &= np.isfinite(rel_r) & (rel_r > 0) & (rel_r < 0.20)
        good &= np.isfinite(cx) & np.isfinite(cy) & (cx > 0) & (cy > 0) & (cx <= 10) & (cy <= 10)
        r_true = arrays['mmom'].astype(float) / np.abs(arrays['mch'].astype(float))
        good &= np.isfinite(r_true) & (r_true > 0)

        residual = tc - u[:, None] * xc
        features = np.column_stack([
            np.ones(n),
            np.sum(tof_pos[:, :, 0] / 240.0, axis=1),
            np.sum(tof_pos[:, :, 1] / 240.0, axis=1),
            np.square(np.sum(tof_pos[:, :, 0] / 240.0, axis=1)),
            np.square(np.sum(tof_pos[:, :, 1] / 240.0, axis=1)),
            np.log(np.maximum(span, 1e-300) / 4.0),
            (residual[:, 0] - residual[:, 1]) / span,
            (residual[:, 2] - residual[:, 3]) / span,
        ])
        good &= np.all(np.isfinite(features), axis=1)
        b = np.searchsorted(EDGES, r_rec, side='right') - 1
        keep = np.flatnonzero(good)
        kept_total += len(keep)
        chunks.append({
            'run': arrays['run'][keep].astype(np.int64),
            'event': arrays['event'][keep].astype(np.int64),
            'source': np.zeros(len(keep), dtype=np.int32),
            'u4_measured': u[keep],
            'features': features[keep],
            'R': r_rec[keep],
            'Rtrue': r_true[keep],
            'relR': rel_r[keep],
            'bin': b[keep].astype(np.int8),
            # Retained only for extraction audit; not passed to the fit objective.
            'truth_mass_input': arrays['mmass'][keep].astype(float),
            'truth_charge_input': arrays['mch'][keep].astype(float),
        })
    print(f'extract {label}: files={len(files)} raw={raw_total} selected={kept_total}', flush=True)
    if not chunks:
        raise RuntimeError(f'No chunks for {label}')
    return {key: np.concatenate([chunk[key] for chunk in chunks]) for key in chunks[0]}


def stratify_and_assign(data: Dict[str, np.ndarray], target: np.ndarray, f3: np.ndarray,
                        seed: int, maximum: int) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    available = np.bincount(data['bin'], minlength=B).astype(int)
    total = int(min(maximum, math.floor(np.min(available / target) * 0.97)))
    requested = np.floor(total * target).astype(int)
    while requested.sum() < total:
        residual = total * target - requested
        requested[int(np.argmax(residual))] += 1
    indices = []
    isotope = []  # 0 He3, 1 He4
    for b in range(B):
        pool = np.flatnonzero(data['bin'] == b)
        chosen = rng.choice(pool, requested[b], replace=False)
        order = rng.permutation(len(chosen))
        n3 = int(round(f3[b] * len(chosen)))
        labels = np.ones(len(chosen), dtype=np.int8)
        labels[order[:n3]] = 0
        indices.append(chosen)
        isotope.append(labels)
    rows = np.concatenate(indices)
    labels = np.concatenate(isotope)
    perm = rng.permutation(len(rows))
    rows = rows[perm]
    labels = labels[perm]
    out = {key: value[rows] for key, value in data.items()}
    out['isotope'] = labels
    u3_true = np.sqrt(1.0 + (M3 / (2.0 * out['Rtrue'])) ** 2)
    u4_true = np.sqrt(1.0 + (M4 / (2.0 * out['Rtrue'])) ** 2)
    out['u_true'] = np.where(out['isotope'] == 0, u3_true, u4_true)
    out['u'] = out['u4_measured'] + np.where(out['isotope'] == 0, u3_true - u4_true, 0.0)
    return out


def feature_state(train: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    center = np.median(train['features'], axis=0)
    q25, q75 = np.percentile(train['features'], [25, 75], axis=0)
    scale = np.maximum(1e-6, 0.741301109 * (q75 - q25))
    center[0] = 0.0
    scale[0] = 1.0
    return center, scale


def standardized(data: Dict[str, np.ndarray], center: np.ndarray, scale: np.ndarray) -> np.ndarray:
    x = np.clip((data['features'] - center) / scale, -4.0, 4.0)
    x[:, 0] = 1.0
    return x


class MixtureDensity:
    def __init__(self, data: Dict[str, np.ndarray], x: np.ndarray, fixed_f4: Optional[np.ndarray]):
        self.u = data['u'].astype(float)
        self.x = x.astype(float)
        self.b = data['bin'].astype(np.int64)
        self.fixed_f4 = None if fixed_f4 is None else np.asarray(fixed_f4, dtype=float)
        self.means = np.sqrt(1.0 + (MASSES[None, :] / (2.0 * data['R'][:, None])) ** 2)
        sensitivity = (self.means ** 2 - 1.0) / self.means
        self.tracker = (4.0 / 6.0) * sensitivity ** 2 * data['relR'][:, None] ** 2

    def components(self, parameters: np.ndarray, corrected: Optional[np.ndarray] = None):
        u = self.u - self.x.dot(parameters[:P]) if corrected is None else corrected
        r = u[:, None] - self.means
        s2 = np.exp(2.0 * parameters[P:P + B])[self.b]
        v = s2[:, None] + self.tracker
        lp = LOG_NORM - 0.5 * np.log(v) - ((NU + 1.0) / 2.0) * np.log1p(r * r / (NU * v))
        return lp, r, v, s2

    def fraction(self, parameters: np.ndarray) -> np.ndarray:
        f4_bin = self.fixed_f4 if self.fixed_f4 is not None else parameters[P + B:P + 2 * B]
        return f4_bin[self.b]

    def __call__(self, parameters: np.ndarray):
        lp, r, v, s2 = self.components(parameters)
        f4 = self.fraction(parameters)
        with np.errstate(divide='ignore'):
            lm = np.logaddexp(np.log1p(-f4) + lp[:, 0], np.log(f4) + lp[:, 1])
        density_ratio = np.exp(lp - lm[:, None])
        posterior = density_ratio * np.column_stack((1.0 - f4, f4))
        denominator = NU * v + r * r
        gc = -np.sum(posterior * (NU + 1.0) * r / denominator, axis=1)
        gs = np.sum(posterior * (s2[:, None] / v) *
                    (1.0 - (NU + 1.0) * r * r / denominator), axis=1)
        gradient = [self.x.T.dot(gc), np.bincount(self.b, weights=gs, minlength=B)]
        if self.fixed_f4 is None:
            gf = density_ratio[:, 0] - density_ratio[:, 1]
            gradient.append(np.bincount(self.b, weights=gf, minlength=B))
        g = np.concatenate(gradient) / len(self.u)
        value = -float(np.mean(lm)) + 0.0005 * np.sum((parameters[1:P] / 0.03) ** 2)
        g[1:P] += 0.001 * parameters[1:P] / 0.03 ** 2
        return value, g

    def score(self, parameters: np.ndarray, corrected: np.ndarray) -> np.ndarray:
        lp = self.components(parameters, corrected)[0]
        return lp[:, 1] - lp[:, 0]


def fit_model(data: Dict[str, np.ndarray], x: np.ndarray, fixed_f4: Optional[np.ndarray],
              mode: str, initial: Optional[np.ndarray] = None, starts: int = 2) -> dict:
    density = MixtureDensity(data, x, fixed_f4)
    npar = P + B + (0 if fixed_f4 is not None else B)
    if initial is None:
        p0 = np.zeros(npar, dtype=float)
        p0[P:P + B] = math.log(0.022)
        if fixed_f4 is None:
            p0[P + B:] = 0.87
    else:
        p0 = np.zeros(npar, dtype=float)
        ncopy = min(npar, len(initial))
        p0[:ncopy] = initial[:ncopy]
        if fixed_f4 is None and len(initial) < npar:
            p0[P + B:] = 0.87
    bounds = [(-0.08, 0.08)] * P + [(math.log(0.002), math.log(0.12))] * B
    if fixed_f4 is None:
        bounds += [(0.001, 0.999)] * B
    if mode == 'raw':
        for j in range(P):
            bounds[j] = (0.0, 0.0)
            p0[j] = 0.0
    solutions = []
    records = []
    for start in range(starts):
        q0 = p0.copy()
        if start:
            q0[P:P + B] += math.log(1.35)
            if fixed_f4 is None:
                q0[P + B:] = 0.5
            if mode != 'raw':
                q0[0] = -0.008
        result = minimize(density, q0, jac=True, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter': 2500, 'maxls': 80, 'ftol': 1e-13, 'gtol': 2e-7})
        record = {'start': start, 'success': bool(result.success), 'status': int(result.status),
                  'objective': float(result.fun), 'iterations': int(result.nit), 'message': str(result.message)}
        records.append(record)
        print('fit', mode, 'fixed' if fixed_f4 is not None else 'free', record, flush=True)
        if result.success and np.isfinite(result.fun) and np.all(np.isfinite(result.x)):
            solutions.append(result)
    if not solutions:
        raise RuntimeError(f'No successful fit for {mode}, fixed={fixed_f4 is not None}')
    best = min(solutions, key=lambda result: result.fun)
    f4 = np.asarray(fixed_f4 if fixed_f4 is not None else best.x[P + B:P + 2 * B])
    return {'parameters': best.x, 'objective': float(best.fun), 'optimizer': records,
            'f4': f4, 'fixed_fraction': fixed_f4 is not None}


def evaluate(name: str, data: Dict[str, np.ndarray], x: np.ndarray, model: dict) -> dict:
    parameters = model['parameters']
    corrected = data['u'] - x.dot(parameters[:P])
    density = MixtureDensity(data, x, model['f4'])
    score = density.score(parameters[:P + B], corrected)
    residual_u = corrected - data['u_true']
    residual_beta = 1.0 / corrected - 1.0 / data['u_true']
    mass2 = (2.0 * data['R']) ** 2 * (corrected ** 2 - 1.0)
    row = {
        'method': name,
        'N': int(len(corrected)),
        'u_bias_median': float(np.median(residual_u)),
        'u_W68': w68(residual_u),
        'u_RMS': rms(residual_u),
        'beta_bias_median': float(np.median(residual_beta)),
        'beta_W68': w68(residual_beta),
        'AUC_He4': auc(score, data['isotope']),
        'mass2_He3_median': float(np.median(mass2[data['isotope'] == 0])),
        'mass2_He3_W68': w68(mass2[data['isotope'] == 0]),
        'mass2_He4_median': float(np.median(mass2[data['isotope'] == 1])),
        'mass2_He4_W68': w68(mass2[data['isotope'] == 1]),
    }
    row['mass2_separation_W68'] = ((row['mass2_He4_median'] - row['mass2_He3_median']) /
                                    math.sqrt(row['mass2_He3_W68'] ** 2 + row['mass2_He4_W68'] ** 2))
    per_bin = []
    for b in range(B):
        q = data['bin'] == b
        per_bin.append({
            'R_low': float(EDGES[b]), 'R_high': float(EDGES[b + 1]), 'N': int(q.sum()),
            'f4_truth': float(np.mean(data['isotope'][q] == 1)),
            'u_bias_median': float(np.median(residual_u[q])), 'u_W68': w68(residual_u[q]),
            'AUC_He4': auc(score[q], data['isotope'][q]),
        })
    row['per_bin'] = per_bin
    return row


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [serializable(item) for item in value]
    return value


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(BASE.glob('*.root'))[:N_FILES]
    if len(paths) != N_FILES:
        raise RuntimeError(f'Need {N_FILES} files, found {len(paths)}')
    train_raw = extract(paths[:N_FILES // 2], 'train-files')
    valid_raw = extract(paths[N_FILES // 2:], 'validation-files')
    target = flight_bin_fractions()
    f3_nominal, flux_provenance = official_f3_by_bin()
    f4_nominal = 1.0 - f3_nominal
    train = stratify_and_assign(train_raw, target, f3_nominal, SEED + 1, MAX_PER_SPLIT)
    valid = stratify_and_assign(valid_raw, target, f3_nominal, SEED + 2, MAX_PER_SPLIT)
    center, scale = feature_state(train)
    x_train = standardized(train, center, scale)
    x_valid = standardized(valid, center, scale)

    # Explicit fit allowlist: truth arrays below are not passed to MixtureDensity.
    fit_train = {key: train[key] for key in ['u', 'features', 'R', 'relR', 'bin']}
    fit_valid = {key: valid[key] for key in ['u', 'features', 'R', 'relR', 'bin']}
    del fit_valid  # Evaluation uses truth, but fitting does not.

    raw = fit_model(fit_train, x_train, f4_nominal, 'raw', starts=2)
    constrained = fit_model(fit_train, x_train, f4_nominal, 'full', raw['parameters'], starts=2)
    free = fit_model(fit_train, x_train, None, 'full', constrained['parameters'], starts=3)
    flight_raw_f4 = np.array([0.7613551049585812, 0.76140866223888,
                              0.7882802305806085, 0.7997684354299699])
    wrong = fit_model(fit_train, x_train, flight_raw_f4, 'full', constrained['parameters'], starts=2)

    injected_train = dict(fit_train)
    injected_valid = dict(valid)
    delta_train = 0.010 * x_train[:, 1] + 0.006 * x_train[:, 6]
    delta_valid = 0.010 * x_valid[:, 1] + 0.006 * x_valid[:, 6]
    injected_train['u'] = fit_train['u'] + delta_train
    injected_valid['u'] = valid['u'] + delta_valid
    injected = fit_model(injected_train, x_train, f4_nominal, 'full', constrained['parameters'], starts=2)

    metrics = [
        evaluate('raw_fixed_fraction', valid, x_valid, raw),
        evaluate('constrained_official_fraction', valid, x_valid, constrained),
        evaluate('free_fraction', valid, x_valid, free),
        evaluate('wrong_flight_raw_fraction', valid, x_valid, wrong),
        evaluate('injected_then_constrained', injected_valid, x_valid, injected),
    ]
    nominal_corrected = valid['u'] - x_valid.dot(constrained['parameters'][:P])
    injected_corrected = injected_valid['u'] - x_valid.dot(injected['parameters'][:P])
    closure_difference = injected_corrected - nominal_corrected
    theta_delta = injected['parameters'][:P] - constrained['parameters'][:P]
    expected_delta = np.array([0.0, 0.010, 0.0, 0.0, 0.0, 0.0, 0.006, 0.0])

    actual_train_f4 = [float(np.mean(train['isotope'][train['bin'] == b] == 1)) for b in range(B)]
    actual_valid_f4 = [float(np.mean(valid['isotope'][valid['bin'] == b] == 1)) for b in range(B)]
    result = {
        'scope': 'Phase-0 kinematic proxy only; not full He3 detector simulation',
        'truth_policy': 'mmom/mmass/mch used only to construct proxy and evaluate frozen fits',
        'fit_allowlist': ['u', 'features[8]', 'R_reconstructed', 'relative_R_uncertainty',
                          'R_bin', 'known_masses', 'bag_level_f4'],
        'input': {
            'files_total': N_FILES,
            'train_files': [str(path) for path in paths[:N_FILES // 2]],
            'validation_files': [str(path) for path in paths[N_FILES // 2:]],
            'train_selected': len(train['u']), 'validation_selected': len(valid['u']),
            'flight_bin_fraction': target,
            'official_f3_flux_fraction': f3_nominal,
            'official_f4_flux_fraction': f4_nominal,
            'actual_train_f4': actual_train_f4,
            'actual_validation_f4': actual_valid_f4,
            'flux_provenance': flux_provenance,
        },
        'feature_center': center,
        'feature_scale': scale,
        'models': {
            'raw': raw, 'constrained': constrained, 'free': free,
            'wrong_fraction': wrong, 'injected': injected,
        },
        'metrics': metrics,
        'injection_closure': {
            'expected_theta_delta': expected_delta,
            'fitted_theta_delta': theta_delta,
            'theta_delta_error': theta_delta - expected_delta,
            'corrected_difference_median': float(np.median(closure_difference)),
            'corrected_difference_W68': w68(closure_difference),
            'corrected_difference_RMS': rms(closure_difference),
        },
    }
    result = serializable(result)
    (OUT / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    with (OUT / 'metrics.tsv').open('w') as handle:
        columns = [key for key in metrics[0] if key != 'per_bin']
        handle.write('\t'.join(columns) + '\n')
        for row in metrics:
            handle.write('\t'.join(str(row[key]) for key in columns) + '\n')
    with (OUT / 'bin_metrics.tsv').open('w') as handle:
        columns = ['method', 'R_low', 'R_high', 'N', 'f4_truth', 'u_bias_median', 'u_W68', 'AUC_He4']
        handle.write('\t'.join(columns) + '\n')
        for row in metrics:
            for item in row['per_bin']:
                record = dict(item, method=row['method'])
                handle.write('\t'.join(str(record[key]) for key in columns) + '\n')
    print('RESULT_DIR', OUT, flush=True)
    print('OFFICIAL_F4', f4_nominal, flush=True)
    print('FREE_F4', free['f4'], flush=True)
    for row in metrics:
        print('METRIC', row['method'], 'u_bias', row['u_bias_median'], 'u_W68', row['u_W68'],
              'AUC', row['AUC_He4'], 'mass_sep', row['mass2_separation_W68'], flush=True)
    print('INJECTION_THETA_DELTA', theta_delta, flush=True)
    print('INJECTION_CLOSURE_RMS', rms(closure_difference), 'W68', w68(closure_difference), flush=True)


if __name__ == '__main__':
    main()
