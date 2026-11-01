"""TOF-only correction trained against an independent measured reference.

Python 3.6 / NumPy 1.19 compatible; no scikit-learn dependency.

The caller supplies ``target = u0 - u_reference`` to ``fit``.  Neither this
module nor its saved model receives the reference measurement as an input.
``predict`` returns ``(u0 - correction, correction)``.  Every correction
feature is invariant under ``t_i -> t_i + a + delta * (length_i / c)``.

An external run-separated CV must choose mode and alpha.  Alpha uses the
conventional *mean* squared loss: 1/(2 N) sum w e**2 + alpha/2 ||coef||**2;
the intercept is unpenalized.  There is no cut on target values, clipping of
the reconstructed speed, or reading of simulated truth fields.
"""
from __future__ import division

import json
import numpy as np


C_CM_PER_NS = 29.9792458
VERSION = 1
ALLOWED_FIELDS = (
    "tof_tl", "tof_leng", "position", "tof_barid", "tof_edep",
    "tk_edep", "direction", "tof_bpatt", "charge",
)
_H = np.array([
    [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(6.0), 1.0 / np.sqrt(12.0)],
    [-1.0 / np.sqrt(2.0), 1.0 / np.sqrt(6.0), 1.0 / np.sqrt(12.0)],
    [0.0, -2.0 / np.sqrt(6.0), 1.0 / np.sqrt(12.0)],
    [0.0, 0.0, -3.0 / np.sqrt(12.0)],
], dtype=np.float64)


def _n_events(d):
    shape = np.shape(d["tof_tl"])
    if len(shape) != 2 or shape[1] != 4:
        raise ValueError("tof_tl must have shape (N, 4)")
    return shape[0]


def _field(d, name, n, trailing, default=np.nan):
    if name not in d:
        return np.full((n,) + trailing, default, dtype=np.float64)
    value = np.asarray(d[name], dtype=np.float64)
    expected = (n,) + trailing
    if value.shape == () and not trailing:
        return np.full(n, float(value), dtype=np.float64)
    if value.shape != expected:
        raise ValueError("{} must have shape {}, received {}".format(
            name, expected, value.shape))
    return value


def _slice_data(d, rows):
    # Slice before casting to float64, so a batch never casts the full dataset.
    out = {}
    for name in ALLOWED_FIELDS:
        if name in d:
            value = np.asarray(d[name])
            out[name] = value[rows] if value.ndim else value
    return out


def select_data(d, rows):
    """Slice the explicit input allowlist; unrelated fields are never copied."""
    return _slice_data(d, rows)


def kinematics(d):
    """Compute the unweighted TOF slope and an orthonormal residual basis.

    A two- or three-layer event can have a raw slope but is never given a
    learned four-layer correction.  Fewer than two valid distinct path
    coordinates produce NaN.  If supplied, tof_bpatt must equal decimal 4444
    for model use: it is a per-layer status string, NOT a bit mask.  Other
    status codes are deliberately not decoded.  Coordinates are
    x = tof_leng / c, in nanoseconds.
    """
    n = _n_events(d)
    t = _field(d, "tof_tl", n, (4,))
    x = _field(d, "tof_leng", n, (4,)) / C_CM_PER_NS
    layer_valid = np.isfinite(t) & np.isfinite(x)
    count = np.sum(layer_valid, axis=1)
    safe_count = np.maximum(count, 1)
    tm = np.sum(np.where(layer_valid, t, 0.0), axis=1) / safe_count
    xm = np.sum(np.where(layer_valid, x, 0.0), axis=1) / safe_count
    tc = np.where(layer_valid, t - tm[:, None], 0.0)
    xc = np.where(layer_valid, x - xm[:, None], 0.0)
    denom = np.sum(xc * xc, axis=1)
    valid_baseline = (count >= 2) & np.isfinite(denom) & (denom > 1e-16)
    u0 = np.full(n, np.nan)
    u0[valid_baseline] = (np.sum(xc * tc, axis=1)[valid_baseline]
                          / denom[valid_baseline])
    xmax = np.max(np.where(layer_valid, x, -np.inf), axis=1)
    xmin = np.min(np.where(layer_valid, x, np.inf), axis=1)
    span = np.where(count > 0, xmax - xmin, 0.0)
    valid_model = valid_baseline & (count == 4) & (span > 1e-8)
    if "tof_bpatt" in d:
        patt = _field(d, "tof_bpatt", n, ())
        valid_model &= np.isfinite(patt) & (patt == 4444)
    residual = tc - np.where(valid_baseline, u0, 0.0)[:, None] * xc

    # In the 3D subspace perpendicular to (1,1,1,1), a Householder reflection
    # rotates the path direction onto the first axis.  Its other two columns
    # give a deterministic orthonormal basis perpendicular to both 1 and x.
    norm = np.sqrt(np.maximum(denom, 1e-16))
    w = np.dot(xc / norm[:, None], _H)
    h = w.copy()
    h[:, 0] += np.where(w[:, 0] >= 0.0, 1.0, -1.0)
    h /= np.sqrt(np.maximum(np.sum(h * h, axis=1), 1e-16))[:, None]
    hh = np.dot(h, _H.T)
    q0 = _H[:, 1][None, :] - 2.0 * hh * h[:, 1:2]
    q1 = _H[:, 2][None, :] - 2.0 * hh * h[:, 2:3]
    orth = np.column_stack((np.sum(residual * q0, axis=1),
                            np.sum(residual * q1, axis=1)))
    orth[~valid_model] = np.nan
    return {
        "u0": u0, "orthogonal_residuals": orth,
        "path_span": span, "centered_x": xc,
        "valid_baseline": valid_baseline, "valid_model": valid_model,
        "layer_valid": layer_valid,
    }


def baseline(d, batch_size=8192):
    """The raw-time OLS estimate of inverse beta, with partial-layer fallback."""
    n = _n_events(d)
    out = np.empty(n, dtype=np.float64)
    for first in range(0, n, batch_size):
        sl = slice(first, min(n, first + batch_size))
        out[sl] = kinematics(_slice_data(d, sl))["u0"]
    return out


def _charge(d, n, state):
    if state["charge_mode"] == "event":
        z = _field(d, "charge", n, ())
        z = np.where(np.isfinite(z) & (z > 0.0), z,
                     state["charge_fallback"])
    else:
        z = np.full(n, state["charge_fallback"], dtype=np.float64)
    return z


def _raw_base(d, state):
    n = _n_events(d)
    kin = kinematics(d)
    z = _charge(d, n, state)
    pos = _field(d, "position", n, (4, 2)) / 100.0
    direction = _field(d, "direction", n, (3,))
    dn = np.sqrt(np.sum(direction * direction, axis=1))
    direction = direction / np.where(dn > 0.0, dn, np.nan)[:, None]
    span = np.maximum(kin["path_span"], 1e-8)
    shape = kin["centered_x"] / span[:, None]
    log_span = np.log(span)
    tof = _field(d, "tof_edep", n, (4,)) / (z * z)[:, None]
    tk = _field(d, "tk_edep", n, (9,)) / (z * z)[:, None]
    # log1p admits measured zero deposits without an artificial log floor;
    # negative/missing deposits are imputed from training data below.
    tof = np.where(tof >= 0.0, np.log1p(np.maximum(tof, 0.0)), np.nan)
    tk = np.where(tk >= 0.0, np.log1p(np.maximum(tk, 0.0)), np.nan)
    residual = kin["orthogonal_residuals"] / span[:, None]
    pieces = [pos.reshape(n, 8), direction, shape, log_span[:, None],
              np.log(z)[:, None], tof, tk, residual]
    names = (["pos{}_{}".format(i + 1, a) for i in range(4) for a in ("x", "y")]
             + ["dir_x", "dir_y", "dir_z"]
             + ["path_shape{}".format(i + 1) for i in range(4)]
             + ["log_path_span", "log_charge"]
             + ["log_tof_edep_z2_{}".format(i + 1) for i in range(4)]
             + ["log_tk_edep_z2_{}".format(i + 1) for i in range(9)]
             + ["orth_residual1_over_span", "orth_residual2_over_span"])
    bars = _field(d, "tof_barid", n, (4,))
    return np.column_stack(pieces), names, bars, kin


def _robust_column_state(a, quantiles):
    p = a.shape[1]
    med = np.empty(p)
    lo = np.empty(p)
    hi = np.empty(p)
    for j in range(p):
        finite = a[np.isfinite(a[:, j]), j]
        if finite.size:
            med[j] = np.median(finite)
            lo[j], hi[j] = np.percentile(finite, quantiles)
        else:
            med[j], lo[j], hi[j] = 0.0, 0.0, 0.0
    return med, lo, hi


def _base_transform(raw, state):
    med = np.asarray(state["base_median"])
    lo = np.asarray(state["base_low"])
    hi = np.asarray(state["base_high"])
    filled = np.where(np.isfinite(raw), raw, med[None, :])
    filled = np.minimum(np.maximum(filled, lo[None, :]), hi[None, :])
    center = np.asarray(state["base_center"])
    scale = np.asarray(state["base_scale"])
    return (filled - center[None, :]) / scale[None, :]


def _expand(base, bars, state, mode, return_names=False):
    """Small explicit bases; no hidden networks or data-dependent target bases."""
    names0 = state["base_names"]
    columns = []
    names = []

    def add(v, name):
        columns.append(v)
        if return_names:
            names.append(name)

    # 30 primary invariant geometry/deposit features (residuals are separate).
    for j in range(30):
        add(base[:, j], names0[j])
    # Modest curvature in impact position, direction, span, and deposit.
    square_ids = list(range(11)) + [15] + list(range(17, 30))
    for j in square_ids:
        add(base[:, j] ** 2, names0[j] + "^2")
    for layer in range(4):
        ix, iy, ie = 2 * layer, 2 * layer + 1, 17 + layer
        add(base[:, ix] * base[:, iy], names0[ix] + "*" + names0[iy])
        for ip in (ix, iy):
            add(base[:, ip] * base[:, ie], names0[ip] + "*" + names0[ie])
            add(base[:, ip] ** 2 * base[:, ie], names0[ip] + "^2*" + names0[ie])
    for i in range(17, 21):
        for j in range(i + 1, 21):
            add(base[:, i] * base[:, j], names0[i] + "*" + names0[j])
    tof_mean = np.mean(base[:, 17:21], axis=1)
    for j in range(21, 30):
        add(base[:, j] * tof_mean, names0[j] + "*mean_log_tof")
    for j in (8, 9, 10, 15):
        add(base[:, j] * tof_mean, names0[j] + "*mean_log_tof")
    bar_features = []
    for layer, vocab in enumerate(state["bar_values"]):
        seen = np.zeros(len(base), dtype=bool)
        for value in vocab:
            indicator = bars[:, layer] == value
            seen |= indicator
            term = indicator.astype(np.float64)
            label = "bar{}={}".format(layer + 1, value)
            add(term, label)
            bar_features.append((term, label))
        term = (~seen).astype(np.float64)
        label = "bar{}=unknown".format(layer + 1)
        add(term, label)
        bar_features.append((term, label))

    if mode == "full":
        for ir in (30, 31):
            r = base[:, ir]
            add(r, names0[ir])
            for j in range(30):
                add(r * base[:, j], names0[ir] + "*" + names0[j])
            for bar_term, label in bar_features:
                add(r * bar_term, names0[ir] + "*" + label)
            add(r * r, names0[ir] + "^2")
            add(r * np.abs(r), names0[ir] + "*abs")
        add(base[:, 30] * base[:, 31], "orth_residual1*orth_residual2")
    elif mode != "geometry":
        raise ValueError("mode must be 'geometry' or 'full'")
    matrix = np.column_stack(columns)
    return (matrix, names) if return_names else matrix


def _feature_matrix(d, state, mode):
    raw, unused_names, bars, kin = _raw_base(d, state)
    expanded = _expand(_base_transform(raw, state), bars, state, mode)
    expanded -= np.asarray(state["feature_center"])[None, :]
    expanded /= np.asarray(state["feature_scale"])[None, :]
    return expanded, kin


def features(d, model):
    """Return standardized feature matrix and names, for small-batch audits."""
    matrix, unused_kin = _feature_matrix(d, model["feature_state"], model["mode"])
    return matrix, list(model["feature_state"]["feature_names"])


def _mad_scale(residual):
    residual = residual[np.isfinite(residual)]
    if not residual.size:
        return 1e-10
    median = np.median(residual)
    return max(float(1.482602218505602 * np.median(np.abs(residual - median))),
               1e-10)


def fit(d, z, target, alpha, mode="full", batch_size=8192,
        max_iter=4, huber_k=1.5, feature_sample=50000,
        feature_quantiles=(0.25, 99.75)):
    """Fit a serializable correction using an independent *measured* target.

    Parameters
    ----------
    d : mapping of measured input arrays
    z : positive scalar charge, or N measured charges
        With an N-vector, inference requires ``d['charge']``; scalar-charge
        models retain the scalar in their state and use it consistently.
    target : N-vector
        Caller-computed ``baseline(d) - 1 / independent_measured_beta``.
        The caller alone is responsible for independent reference selection,
        run-separated folds, and keeping MC truth out of model selection.
    alpha : nonnegative float
        Conventional mean-loss ridge; fixed candidate set .01, .1, 1, 10.
    mode : 'geometry' or 'full'
        'geometry' excludes *all* time residual features and interactions.
    max_iter : integer from 1 to 5
        1 is ordinary ridge; 4 (default) uses three Huber reweightings.

    Invalid kinematics and nonfinite labels cannot be fitted.  There is no
    magnitude-based label rejection.  Label values never choose the feature
    basis, clipping thresholds, imputation, or standardization.
    """
    n = _n_events(d)
    y = np.asarray(target, dtype=np.float64)
    if y.shape != (n,):
        raise ValueError("target must have shape (N,)")
    if not np.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    if mode not in ("geometry", "full"):
        raise ValueError("mode must be 'geometry' or 'full'")
    if not 1 <= max_iter <= 5:
        raise ValueError("max_iter must be between 1 and 5")
    if batch_size < 1 or feature_sample < 1 or huber_k <= 0:
        raise ValueError("batch_size, feature_sample and huber_k must be positive")
    zarr = np.asarray(z, dtype=np.float64)
    train_d = dict((key, d[key]) for key in ALLOWED_FIELDS if key in d)
    if zarr.ndim == 0:
        if not np.isfinite(zarr) or zarr <= 0.0:
            raise ValueError("scalar z must be positive and finite")
        charge_mode, charge_fallback = "scalar", float(zarr)
    elif zarr.shape == (n,):
        if np.any(~np.isfinite(zarr) | (zarr <= 0.0)):
            raise ValueError("all event charges must be positive and finite")
        charge_mode, charge_fallback = "event", float(np.median(zarr))
        train_d["charge"] = zarr
    else:
        raise ValueError("z must be a scalar or N-vector")
    state = {"charge_mode": charge_mode, "charge_fallback": charge_fallback}
    valid = np.zeros(n, dtype=bool)
    for first in range(0, n, batch_size):
        sl = slice(first, min(n, first + batch_size))
        valid[sl] = kinematics(_slice_data(train_d, sl))["valid_model"]
    rows = np.flatnonzero(valid & np.isfinite(y))
    if len(rows) < 10:
        raise ValueError("at least 10 finite four-layer labelled events are required")

    # Sample only for robust feature quantiles.  Feature scales below use every
    # training row.  Evenly spread deterministic sampling uses no label values.
    if len(rows) > feature_sample:
        sample_rows = rows[np.linspace(0, len(rows) - 1, feature_sample).astype(np.int64)]
    else:
        sample_rows = rows
    raw, names, bars, unused_kin = _raw_base(_slice_data(train_d, sample_rows), state)
    med, low, high = _robust_column_state(raw, feature_quantiles)
    filled = np.minimum(np.maximum(np.where(np.isfinite(raw), raw, med[None, :]),
                                   low[None, :]), high[None, :])
    center = np.mean(filled, axis=0)
    scale = np.std(filled, axis=0)
    scale = np.where(scale > 1e-10, scale, 1.0)
    state.update({"base_names": names, "base_median": med.tolist(),
                  "base_low": low.tolist(), "base_high": high.tolist(),
                  "base_center": center.tolist(), "base_scale": scale.tolist(),
                  "feature_quantiles": list(feature_quantiles)})
    vocab = []
    for layer in range(4):
        v = bars[:, layer]
        v = v[np.isfinite(v) & (v >= 0) & (v == np.floor(v))]
        values, counts = np.unique(v.astype(np.int64), return_counts=True)
        # At most 16 bar indicators per layer bounds even malformed inputs.
        if len(values) > 16:
            values = values[np.argsort(-counts, kind="mergesort")[:16]]
        vocab.append(sorted(int(value) for value in values))
    state["bar_values"] = vocab
    small_matrix, fnames = _expand(_base_transform(raw, state), bars, state, mode, True)
    p = small_matrix.shape[1]
    del raw, bars, filled, small_matrix
    sums = np.zeros(p)
    squares = np.zeros(p)
    for first in range(0, len(rows), batch_size):
        ids = rows[first:first + batch_size]
        raw, unused_names, bars, unused_kin = _raw_base(_slice_data(train_d, ids), state)
        f = _expand(_base_transform(raw, state), bars, state, mode)
        sums += np.sum(f, axis=0)
        squares += np.sum(f * f, axis=0)
    fc = sums / len(rows)
    fs = np.sqrt(np.maximum(squares / len(rows) - fc * fc, 0.0))
    fs = np.where(fs > 1e-10, fs, 1.0)
    state.update({"feature_center": fc.tolist(), "feature_scale": fs.tolist(),
                  "feature_names": fnames})

    coef = np.zeros(p)
    intercept = float(np.median(y[rows]))
    scale_h = _mad_scale(y[rows] - intercept)
    history = []
    resid = np.empty(len(rows))
    for iteration in range(max_iter):
        gram = np.zeros((p + 1, p + 1))
        rhs = np.zeros(p + 1)
        sw = 0.0
        for first in range(0, len(rows), batch_size):
            ids = rows[first:first + batch_size]
            f, unused_kin = _feature_matrix(_slice_data(train_d, ids), state, mode)
            error = y[ids] - intercept - np.dot(f, coef)
            if iteration == 0:
                weight = np.ones(len(ids))
            else:
                weight = np.minimum(1.0, huber_k * scale_h / np.maximum(np.abs(error), 1e-300))
            root_w = np.sqrt(weight)
            a = np.empty((len(ids), p + 1))
            a[:, 0] = root_w
            a[:, 1:] = f * root_w[:, None]
            gram += np.dot(a.T, a)
            rhs += np.dot(a.T, y[ids] * root_w)
            sw += float(np.sum(weight))
        gram /= len(rows)
        rhs /= len(rows)
        gram.flat[(p + 1) + 1::(p + 2)] += float(alpha)
        # alpha=0 may leave redundant one-hot columns; use a least-squares
        # solve only when the direct normal-equation solve is singular.
        try:
            solution = np.linalg.solve(gram, rhs)
        except np.linalg.LinAlgError:
            solution = np.linalg.lstsq(gram, rhs, rcond=1e-12)[0]
        intercept = float(solution[0])
        coef = solution[1:]
        for first in range(0, len(rows), batch_size):
            ids = rows[first:first + batch_size]
            f, unused_kin = _feature_matrix(_slice_data(train_d, ids), state, mode)
            resid[first:first + len(ids)] = y[ids] - intercept - np.dot(f, coef)
        scale_h = _mad_scale(resid)
        history.append({"iteration": iteration + 1, "residual_mad_sigma": scale_h,
                        "residual_rms": float(np.sqrt(np.mean(resid * resid))),
                        "mean_weight": sw / len(rows)})
    model = {
        "version": VERSION, "model_type": "independent_measured_reference_tof",
        "baseline": "OLS with intercept, x=tof_leng/29.9792458 ns",
        "target_definition": "OLS_inverse_beta - independent_measured_inverse_beta",
        "mode": mode, "alpha": float(alpha), "huber_k": float(huber_k),
        "intercept": intercept, "coef": coef.tolist(), "feature_state": state,
        "training": {"n_input": n, "n_fit": len(rows), "n_features": p,
                     "n_feature_sample": len(sample_rows), "iterations": history},
    }
    # Fail at training time if a numerical failure produced an invalid artifact.
    json.dumps(model, allow_nan=False)
    return model


def predict(d, model, batch_size=8192):
    """Return (inverse_beta, correction), with zero correction for partial TOF.

    Four valid layers permit a learned correction.  Missing geometry or
    deposits use the stored training medians; unseen bars use their explicit
    unknown indicator.  A baseline that cannot be reconstructed remains NaN.
    The output inverse beta is never clipped or fed back into the features.
    """
    if model.get("version") != VERSION:
        raise ValueError("unsupported reference model version")
    n = _n_events(d)
    if model["feature_state"]["charge_mode"] == "event" and "charge" not in d:
        raise ValueError("event-charge model requires measured d['charge']")
    coef = np.asarray(model["coef"], dtype=np.float64)
    u = np.empty(n, dtype=np.float64)
    corr = np.zeros(n, dtype=np.float64)
    for first in range(0, n, batch_size):
        sl = slice(first, min(n, first + batch_size))
        f, kin = _feature_matrix(_slice_data(d, sl), model["feature_state"], model["mode"])
        if f.shape[1] != len(coef):
            raise ValueError("model coefficient and feature dimensions disagree")
        value = float(model["intercept"]) + np.dot(f, coef)
        apply = kin["valid_model"] & np.isfinite(value)
        corr[sl] = np.where(apply, value, 0.0)
        u[sl] = kin["u0"] - corr[sl]
    return u, corr


def reconstruct(d, model=None, batch_size=8192, return_details=False):
    """Safe reconstruction with the producer's fallback for ineligible events.

    By default returns only inverse beta.  With return_details=True returns
    (inverse_beta, details).  Missing/invalid models fall back to raw OLS.
    Events ineligible for four-layer correction (including bpatt != 4444)
    then use 1/tof_betah, when that producer measurement is available.  This
    fallback branch is the module's ONLY permitted reader of tof_betah.
    If the producer beta is zero/nonfinite, its fallback is NaN.  When the
    producer field is absent, the available raw OLS estimate is retained.
    No beta value can affect any learned correction or feature.
    """
    reason = None
    if model is None:
        reason = "no_model"
    else:
        try:
            u, corr = predict(d, model, batch_size=batch_size)
        except (ValueError, KeyError, TypeError) as error:
            reason = "{}: {}".format(type(error).__name__, error)
    if reason is not None:
        u = baseline(d, batch_size=batch_size)
        corr = np.zeros(len(u), dtype=np.float64)
    producer_fallback = np.zeros(len(u), dtype=bool)
    if "tof_betah" in d:
        for first in range(0, len(u), batch_size):
            sl = slice(first, min(len(u), first + batch_size))
            kin = kinematics(_slice_data(d, sl))
            bad = ~kin["valid_model"]
            if np.any(bad):
                # Slice only these fallback rows, never put beta in train_d,
                # ALLOWED_FIELDS, a feature matrix, or a fitted model state.
                beta = np.asarray(d["tof_betah"])[sl][bad].astype(np.float64)
                fallback_u = np.full(len(beta), np.nan)
                good_beta = np.isfinite(beta) & (beta != 0.0)
                fallback_u[good_beta] = 1.0 / beta[good_beta]
                ub = u[sl]
                cb = corr[sl]
                ub[bad] = fallback_u
                cb[bad] = 0.0
                producer_fallback[sl] = bad
    if return_details:
        return u, {"correction": corr, "fallback_reason": reason,
                   "producer_fallback": producer_fallback}
    return u
