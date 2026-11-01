"""One-way material transport with an explicit unrestricted measurement extension.

Physical subluminal states use the pre-existing continuous-mass mean-loss
propagator. Nonphysical measured inverse speeds u<=1 are retained with zero
loss correction, not clipped to a mass peak or removed from evaluation.
"""
import numpy as np
import rich_transport as material


def subset_geometry(geometry, select):
    n = len(geometry['valid'])
    output = dict(geometry)
    for key in ['column', 'valid', 'target_z', 'NaF_predicted']:
        if key in geometry:
            if len(geometry[key]) != n:
                raise ValueError('Inconsistent event geometry')
            output[key] = np.asarray(geometry[key])[select]
    return output


def subdivide(geometry, parts=2):
    if not isinstance(parts, int) or parts < 1:
        raise ValueError('Positive integer subdivision required')
    result = dict(geometry)
    result['column'] = np.repeat(np.asarray(geometry['column'], float) / parts, parts, axis=1)
    result['segment'] = np.repeat(geometry['segment'], parts)
    return result


def propagate(u, rigidity, geometry, material_scale=1.):
    u = np.asarray(u, float); rigidity = np.asarray(rigidity, float)
    if u.ndim != 1 or rigidity.shape != u.shape or len(geometry['valid']) != len(u):
        raise ValueError('Incompatible inverse-speed, rigidity or geometry input')
    if not np.isfinite(u).all() or np.any(u == 0) or not np.isfinite(rigidity).all() or np.any(rigidity <= 0):
        raise ValueError('Finite nonzero inverse speed and positive rigidity required')
    if not np.isfinite(material_scale) or material_scale < 0:
        raise ValueError('Finite nonnegative material scale required')
    if not np.asarray(geometry['valid'], bool).all():
        raise ValueError('Invalid measured geometry must be resolved before estimator comparison')
    column = np.asarray(geometry['column'], float)
    if column.shape[0] != len(u) or column.shape[2] != 9 or not np.isfinite(column).all() or np.any(column < 0):
        raise ValueError('Finite nonnegative nine-element material columns required')
    if len(geometry['segment']) != column.shape[1]:
        raise ValueError('Material column and segment lengths differ')
    out = u.copy(); beta = 1 / u; valid = np.ones(len(u), bool)
    losses = np.zeros((len(u), 7)); positive = u > 1
    # The inherited Bethe routine is only used where it does not floor beta.
    physical = positive & (u <= 20)
    outside = positive & ~physical
    valid[outside] = False; out[outside] = np.nan; beta[outside] = np.nan; losses[outside] = np.nan
    if physical.any() and material_scale > 0:
        selected = subset_geometry(geometry, physical)
        result = material.propagate(u[physical], rigidity[physical], selected, material_scale=material_scale)
        indices = np.flatnonzero(physical); accepted = result['valid']
        valid[indices] = accepted
        beta[indices] = result['beta_planes'][:, -1]
        out[indices] = 1 / beta[indices]
        losses[indices] = result['loss_GeV']
    if material_scale == 0:
        valid[:] = True; out[:] = u; beta[:] = 1 / u; losses[:] = 0
    if np.any(valid & (~np.isfinite(out) | ~np.isfinite(beta))):
        raise FloatingPointError('Accepted propagation produced nonfinite output')
    if np.any(valid & physical & (out < u - 1e-12)):
        raise FloatingPointError('Positive mean energy loss increased physical beta')
    return dict(u=out, beta=beta, valid=valid, physical_input=positive, zero_correction_extension=~positive,
                outside_Bethe_domain=outside, loss_GeV=losses)
