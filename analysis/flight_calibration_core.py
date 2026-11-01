"""Linearized global calibration with free event time and slope.

Templates and noise scales must be independently justified. In this linear
step they may not depend on the same event's measured timing slope. This is
an inference primitive, not an implementation of AMS material transport.
"""
import numpy as np


def project(x, values, sigma):
    """Whiten and remove each event's intercept and timing slope."""
    x = np.asarray(x, float)
    sigma = np.broadcast_to(np.asarray(sigma, float), x.shape)
    if x.ndim != 2 or not np.isfinite(x).all() or (np.ptp(x, axis=1) <= 0).any():
        raise ValueError('Finite paths with nonzero span required')
    if not np.isfinite(sigma).all() or (sigma <= 0).any():
        raise ValueError('Positive finite noise scales required')
    a = np.stack((np.ones_like(x), x), axis=-1)/sigma[:, :, None]
    inv = np.linalg.inv(np.einsum('nmi,nmj->nij', a, a))
    whitened = values/(sigma if values.ndim == 2 else sigma[:, :, None])
    if values.ndim == 2:
        rhs = np.einsum('nmi,nm->ni', a, whitened)
        local = np.einsum('nij,nj->ni', inv, rhs)
        return whitened-np.einsum('nmi,ni->nm', a, local)
    rhs = np.einsum('nmi,nmp->nip', a, whitened)
    local = np.einsum('nij,njp->nip', inv, rhs)
    return whitened-np.einsum('nmi,nip->nmp', a, local)


def solve_global(x, times, design, sigma, relative_cut=1e-10,
                 penalty=None, allow_rank_deficient=False):
    """Fit global columns after profiling two local event parameters.

    Rank uses data alone. A smoothing penalty is not measured information.
    Sampling covariance uses the sandwich form when a penalty is present;
    the inverse penalized precision is returned separately. Neither object
    certifies unidentified parameters or accounts for model bias.
    """
    x, times, design = [np.asarray(v, float) for v in (x, times, design)]
    if times.shape != x.shape or design.shape[:2] != x.shape or design.ndim != 3:
        raise ValueError('Incompatible input dimensions')
    if not np.isfinite(times).all() or not np.isfinite(design).all():
        raise ValueError('Finite measurements and design required')
    bp, yp = project(x, design, sigma), project(x, times, sigma)
    fisher = np.einsum('nmi,nmj->ij', bp, bp)
    rhs = np.einsum('nmi,nm->i', bp, yp)
    raw_norm2 = np.sum((design/np.broadcast_to(np.asarray(sigma, float), x.shape)[:, :, None])**2, axis=(0, 1))
    visible = np.diag(fisher) > 1e-24*np.maximum(raw_norm2, 1e-300)
    fisher[~visible, :] = 0.
    fisher[:, ~visible] = 0.
    rhs[~visible] = 0.
    scales = np.where(visible, np.sqrt(np.maximum(np.diag(fisher), 1e-300)), 1.)
    scaled = fisher/scales[:, None]/scales[None, :]
    eigen, vectors = np.linalg.eigh(scaled)
    keep = eigen > max(float(eigen.max()), 1.)*relative_cut
    if not keep.all() and not allow_rank_deficient:
        raise ValueError('Data identify %d of %d global combinations' % (keep.sum(), len(eigen)))
    precision = scaled.copy()
    if penalty is not None:
        penalty = np.asarray(penalty, float)
        if penalty.shape != fisher.shape or not np.allclose(penalty, penalty.T):
            raise ValueError('Invalid penalty matrix')
        if np.linalg.eigvalsh(penalty).min() < -1e-9:
            raise ValueError('Penalty must be positive semidefinite')
        precision += penalty/scales[:, None]/scales[None, :]
    inverse = np.linalg.pinv(precision, rcond=relative_cut)
    theta = inverse.dot(rhs/scales)/scales
    precision_inverse = inverse/scales[:, None]/scales[None, :]
    covariance = inverse.dot(scaled).dot(inverse)/scales[:, None]/scales[None, :]
    residual = yp-np.einsum('nmp,p->nm', bp, theta)
    return dict(theta=theta, covariance_on_estimable_subspace=covariance, fisher=fisher,
                penalized_precision_inverse=precision_inverse,
                covariance_valid_all_parameters=bool(keep.all()),
                scaled_eigenvalues=eigen, data_rank=int(keep.sum()),
                weak_vectors_scaled=vectors[:, ~keep], column_scales=scales,
                residual_sum_squares=float(np.sum(residual*residual)),
                residual_dof=int(len(x)*(x.shape[1]-2)-keep.sum()))


def reconstruct(x, times, design, sigma, theta):
    """GLS with an unrestricted measured inverse speed."""
    x, times, design = [np.asarray(v, float) for v in (x, times, design)]
    sigma = np.broadcast_to(np.asarray(sigma, float), x.shape)
    corrected = times-np.einsum('nmp,p->nm', design, np.asarray(theta, float))
    weights = 1/sigma**2
    xc = x-np.sum(weights*x, axis=1)[:, None]/weights.sum(axis=1)[:, None]
    denominator = np.sum(weights*xc*xc, axis=1)
    coefficients = weights*xc/denominator[:, None]
    u = np.sum(coefficients*corrected, axis=1)
    t0 = np.sum(weights*(corrected-u[:, None]*x), axis=1)/weights.sum(axis=1)
    return dict(u=u, t0=t0, timing_variance=1/denominator, slope_weights=coefficients)


def impulse_kernel(x, position, center=0.):
    """Time shift in ns per unit inverse-speed jump; path/c is in ns."""
    x = np.asarray(x, float)
    position = np.broadcast_to(np.asarray(position, float), (len(x),))
    center = np.broadcast_to(np.asarray(center, float), (len(x),))
    return (np.maximum(x-position[:, None], 0)
            - np.maximum(center-position, 0)[:, None]
            - (x-center[:, None])*(center >= position)[:, None])


def time_penalty(periods, features, strength):
    """First differences of period-major coefficient vectors."""
    if periods < 2 or features < 1 or strength < 0:
        raise ValueError('Invalid smoothing dimensions or strength')
    d = np.diff(np.eye(periods), axis=0)
    return strength*np.kron(d.T.dot(d), np.eye(features))
