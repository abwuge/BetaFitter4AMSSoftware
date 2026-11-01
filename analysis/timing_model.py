"""Self-supervised TOF timing model with no beta or mass training target.

Only two orthogonal timing contrasts enter the likelihood. The learned mean
is applied in that subspace, preserving the original affine calibration.
A conditional diagonal layer-noise model determines GLS weights.
"""
import numpy as np
from scipy.optimize import minimize

C = 29.9792458  # cm/ns
NU = 6.0
RIDGE = 0.01


def kinematics(d):
    x = d['tof_leng'].astype(float) / C
    x -= x.mean(axis=1)[:, None]
    t = d['tof_tl'].astype(float)
    t -= t.mean(axis=1)[:, None]
    with np.errstate(invalid='ignore', divide='ignore'):
        u = (x * t).sum(axis=1) / (x * x).sum(axis=1)
    return x, t, u


def clean(d):
    x, t, u = kinematics(d)
    return ((d['tof_bpatt'] == 4444) & np.all(np.isfinite(x) & np.isfinite(t), axis=1)
            & (np.ptp(x, axis=1) > 0) & np.isfinite(u) & (u > 0)
            & np.all(np.isfinite(d['position']), axis=(1, 2))
            & np.all(np.isfinite(d['tof_edep']) & (d['tof_edep'] > 0), axis=1)
            & np.all((d['tof_barid'] >= 0) & (d['tof_barid'] < 10), axis=1))


def contrasts(x):
    # Batched QR is unavailable in the host NumPy. Orthogonalize two rows
    # of the rank-two projector, choosing the largest norm at each step.
    xc = x - x.mean(axis=1)[:, None]
    p = np.eye(4)[None, :, :] - .25
    p = p - xc[:, :, None]*xc[:, None, :]/np.sum(xc*xc, axis=1)[:, None, None]
    rows = np.arange(len(x))
    first = p[rows, np.argmax(np.diagonal(p, axis1=1, axis2=2), axis=1)]
    first /= np.linalg.norm(first, axis=1)[:, None]
    p = p - first[:, :, None]*first[:, None, :]
    second = p[rows, np.argmax(np.diagonal(p, axis1=1, axis2=2), axis=1)]
    second /= np.linalg.norm(second, axis=1)[:, None]
    return np.stack([first, second], axis=1)


def raw_features(d, z):
    # No measured time, fitted beta, rigidity, or simulation quantity is a feature.
    e = np.log(d['tof_edep'].astype(float) / float(z * z))
    p = d['position'].astype(float) / 60.
    px, py = p[:, :, 0], p[:, :, 1]
    f = [np.ones_like(e), e, e*e, px, py, px*px, py*py, px*py, e*px, e*py]
    f += [(d['tof_barid'] == b).astype(float) for b in range(10)]
    return np.stack(f, axis=2)


def prepare_features(d, z, state=None):
    f = raw_features(d, z)
    if state is None:
        lo, hi = np.percentile(f, [1, 99], axis=0)
        lo[:, 10:], hi[:, 10:] = 0, 1
        center = np.mean(np.clip(f, lo, hi), axis=0)
        scale = np.std(np.clip(f, lo, hi), axis=0)
        scale[scale < .05] = 1
        center[:, 0], scale[:, 0] = 0, 1
        state = dict(low=lo.tolist(), high=hi.tolist(), center=center.tolist(), scale=scale.tolist())
    f = (np.clip(f, state['low'], state['high']) - state['center']) / state['scale']
    return f, state


def likelihood(parameters, f, q, y, scale, regularize=True):
    """Student-t contrast likelihood and analytic gradient (2 residual dof)."""
    k = f.shape[2]
    a = parameters[:4*k].reshape(4, k)
    b = parameters[4*k:].reshape(4, k)
    mean = scale * np.einsum('nik,ik->ni', f, a)
    logv = np.einsum('nik,ik->ni', f, b)
    v = scale**2 * np.exp(4 * np.tanh(logv / 4))
    residual = y - np.einsum('nri,ni->nr', q, mean)
    cov = np.einsum('nri,ni,nsi->nrs', q, v, q)
    det = cov[:, 0, 0]*cov[:, 1, 1] - cov[:, 0, 1]**2
    inv = np.empty_like(cov)
    inv[:, 0, 0], inv[:, 1, 1] = cov[:, 1, 1]/det, cov[:, 0, 0]/det
    inv[:, 0, 1] = inv[:, 1, 0] = -cov[:, 0, 1]/det
    w = np.einsum('nrs,ns->nr', inv, residual)
    chi = np.sum(residual*w, axis=1)
    robust = (NU + 2)/(NU + chi)
    n = len(f)
    value = np.mean(.5*np.log(det) + .5*(NU+2)*np.log1p(chi/NU))
    dmean = -scale * np.einsum('nri,nr->ni', q, w) * robust[:, None]
    dcov = .5 * (inv - robust[:, None, None] * w[:, :, None]*w[:, None, :])
    dv = np.einsum('nri,nrs,nsi->ni', q, dcov, q)
    dlogv = dv*v*(1-np.tanh(logv/4)**2)
    grad = np.concatenate([np.einsum('ni,nik->ik', dmean, f).ravel(),
                           np.einsum('ni,nik->ik', dlogv, f).ravel()])/n
    if regularize:
        # Fix mean gauge by minimum norm and shrink poorly identified variances.
        value += .5*RIDGE*np.sum(parameters**2)
        grad += RIDGE*parameters
    return value, grad


def fit(d, z, max_iter=200):
    x, t, _ = kinematics(d)
    q = contrasts(x)
    y = np.einsum('nri,ni->nr', q, t)
    scale = float(np.median(np.abs(y - np.median(y, axis=0)))/.67448975)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Nonpositive timing scale')
    f, state = prepare_features(d, z)
    initial = np.zeros(8*f.shape[2])
    result = minimize(likelihood, initial, args=(f, q, y, scale), jac=True,
                      method='L-BFGS-B', options=dict(maxiter=max_iter, ftol=1e-9, gtol=1e-5))
    if not result.success or not np.all(np.isfinite(result.x)):
        raise RuntimeError('Timing model did not converge: ' + str(result.message))
    return dict(kind='tof_contrast_student_t', charge=z, scale_ns=scale, features=state,
                parameters=result.x.tolist(), nu=NU, ridge=RIDGE,
                training=dict(entries=len(x), iterations=int(result.nit),
                              objective=float(result.fun), converged=bool(result.success)))


def predict(d, model):
    x, t, _ = kinematics(d)
    q = contrasts(x)
    f, _ = prepare_features(d, model['charge'], model['features'])
    k = f.shape[2]
    p = np.asarray(model['parameters'])
    raw_mean = model['scale_ns'] * np.einsum('nik,ik->ni', f, p[:4*k].reshape(4,k))
    # Never extrapolate an unidentifiable mean along beta/t0 directions.
    mean = np.einsum('nri,nr->ni', q, np.einsum('nri,ni->nr', q, raw_mean))
    logv = np.einsum('nik,ik->ni', f, p[4*k:].reshape(4,k))
    v = model['scale_ns']**2 * np.exp(4*np.tanh(logv/4))
    weights = 1/v
    xc = x - (weights*x).sum(axis=1)[:,None]/weights.sum(axis=1)[:,None]
    u = np.sum(weights*xc*(t-mean), axis=1)/np.sum(weights*xc*xc, axis=1)
    return u, mean, v


def score(d, model):
    # Keep the full flight holdout, but bound feature/covariance memory.
    n = len(d['tof_tl'])
    if not n:
        raise ValueError('Empty validation sample')
    totals = np.zeros(2)
    residuals = {name: [[] for _ in range(4)] for name in ('equal_weight', 'conditional')}
    p = np.asarray(model['parameters'])
    for start in range(0, n, 100000):
        block = {k: v[start:start+100000] for k,v in d.items()}
        x, t, _ = kinematics(block)
        q = contrasts(x)
        f, _ = prepare_features(block, model['charge'], model['features'])
        y = np.einsum('nri,ni->nr', q, t)
        totals[0] += len(x)*likelihood(p, f, q, y, model['scale_ns'], regularize=False)[0]
        totals[1] += len(x)*likelihood(np.zeros_like(p), f, q, y, model['scale_ns'], regularize=False)[0]
        _, mean, variance = predict(block, model)
        for name, mu, var in [('equal_weight', np.zeros_like(t), np.ones_like(t)),
                               ('conditional', mean, variance)]:
            for layer in range(4):
                weights = 1/var
                weights[:, layer] = 0
                xb = np.sum(weights*x, axis=1)/weights.sum(axis=1)
                tb = np.sum(weights*(t-mu), axis=1)/weights.sum(axis=1)
                xc = x-xb[:, None]
                slope = np.sum(weights*xc*(t-mu), axis=1)/np.sum(weights*xc*xc, axis=1)
                r = t[:, layer]-mu[:, layer]-tb-slope*(x[:, layer]-xb)
                residuals[name][layer].append(r)
    loo = {}
    for name, layers in residuals.items():
        quantiles = [np.percentile(np.concatenate(values), [16,50,84]) for values in layers]
        loo[name] = dict(central68_half_width_ns=[float((v[2]-v[0])/2) for v in quantiles],
                         median_ns=[float(v[1]) for v in quantiles])
    new, baseline = totals/n
    return dict(conditional_nll=float(new), pooled_equal_noise_nll=float(baseline),
                nll_improvement=float(baseline-new), entries=n, leave_one_layer_out=loo)


def reconstruct(d, model):
    valid = clean(d)
    with np.errstate(invalid='ignore', divide='ignore'):
        u = 1/d['tof_betah'].astype(float)
    indices = np.flatnonzero(valid)
    for start in range(0, len(indices), 100000):
        block = indices[start:start+100000]
        part = {k: v[block] for k,v in d.items()}
        candidate, _, _ = predict(part, model)
        good = np.isfinite(candidate) & (candidate > 0)
        u[block[good]] = candidate[good]
        valid[block[~good]] = False
    return u, valid
