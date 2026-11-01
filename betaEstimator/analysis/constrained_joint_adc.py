"""ADC-assisted time slopes with explicit signal and gain response constraints."""
import numpy as np
import joint_adc_likelihood as joint
import flight_conditional_likelihood as base

VARIANTS = ('speed_null', 'speed_gain_null')


def constrained_weights(precision, direction_time, direction_adc, common_gain=False):
    """Minimize w.T Sigma w subject to unit time and null ADC-signal response.

    A third constraint can remove a common additive shift in log ADC, which
    corresponds to a common multiplicative gain of all associated end ADCs.
    """
    q = np.asarray(precision, dtype=float)
    t = np.asarray(direction_time, dtype=float)
    adc = np.asarray(direction_adc, dtype=float)
    if q.shape != (7, 7) or t.ndim != 2 or t.shape[1] != 3 or adc.shape != (len(t), 4):
        raise ValueError('Expected a 7-D precision, n x 3 time directions and n x 4 ADC directions')
    if not all(np.isfinite(x).all() for x in (q, t, adc)):
        raise ValueError('Nonfinite constrained-weight input')
    if not np.allclose(q, q.T, rtol=1e-12, atol=1e-12) or np.linalg.eigvalsh(q)[0] <= 0:
        raise ValueError('Precision must be symmetric positive definite')
    nconstraints = 3 if common_gain else 2
    c = np.zeros((len(t), 7, nconstraints))
    c[:, :3, 0] = t
    c[:, 3:, 1] = adc
    if common_gain:
        c[:, 3:, 2] = 1.
    qc = np.einsum('ij,njk->nik', q, c)
    gram = np.einsum('nik,nil->nkl', c, qc)
    target = np.zeros((len(t), nconstraints))
    target[:, 0] = 1.
    if np.any(np.linalg.eigvalsh(gram)[:, 0] <= 0):
        raise ValueError('Dependent response constraints')
    multiplier = np.linalg.solve(gram, target)
    w = np.einsum('nik,nk->ni', qc, multiplier)
    response = np.einsum('ni,nik->nk', w, c)
    if np.max(abs(response - target)) > 2e-10:
        raise FloatingPointError('Constrained response failed')
    scale_variance = np.einsum('ni,ij,nj->n', w, np.linalg.inv(q), w)
    return dict(weights=w, constraint_matrix=c, response=response, scale_variance=scale_variance,
                gram_condition=np.linalg.cond(gram))


def reconstruct(theta, x, y, features, variant='speed_null'):
    if variant not in VARIANTS:
        raise ValueError('Unknown ADC time-constraint variant')
    unconstrained = joint.reconstruct(theta, x, y, features)
    _, q = joint.precision(theta)
    result = constrained_weights(q, unconstrained['direction'][:, :3], unconstrained['direction'][:, 3:],
                                 common_gain=variant == 'speed_gain_null')
    w = result['weights']
    u = np.sum(w * unconstrained['centered'], axis=1)
    result.update(u=u, time_weights=w[:, :3].dot(base.D), centered=unconstrained['centered'],
                  time_u=unconstrained['time_u'], joint_u=unconstrained['joint_u'],
                  model_width_gain_vs_time=1 - np.sqrt(result['scale_variance'] * unconstrained['time_information']))
    if not np.isfinite(u).all():
        raise FloatingPointError('Nonfinite constrained inverse speed')
    return result
