"""Period-dependent raw-ADC response, with no new additive timing calibration.

Flight uses declared UTC eras. MC uses an explicitly labelled run-hash negative
control, not a fabricated physical detector date. Both are calibrated separately.
"""
import calendar
import hashlib
from datetime import datetime
from pathlib import Path
import numpy as np
import joint_adc_likelihood as original

ROOT = original.ROOT
OUT = ROOT / 'results/flight_calibration_research/continuation_20261013/periodic_adc'
VARIANTS = ['static', 'step', 'spline']
YEARS = [2011, 2014, 2017, 2020, 2023, 2027]
EDGES = np.asarray([calendar.timegm(datetime(year, 1, 1).timetuple()) for year in YEARS], float)
KNOTS = .5 * (EDGES[:-1] + EDGES[1:])


def period_coordinates(run, dataset):
    values = np.asarray(run)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError('Finite run identifiers required')
    if dataset == 'flight':
        coordinate = values.astype(float)
        if np.any(coordinate < EDGES[0]) or np.any(coordinate >= EDGES[-1]):
            raise ValueError('Flight run outside declared 2011-2027 response eras')
        return coordinate
    if dataset == 'mc':
        unique, inverse = np.unique(values, return_inverse=True)
        ranks = np.asarray([int.from_bytes(hashlib.sha256(('periodic_adc_negative_control:%d' % int(run)).encode()).digest()[:8], 'big') / float(2 ** 64)
                            for run in unique])
        return (KNOTS[0] + ranks * (KNOTS[-1] - KNOTS[0]))[inverse]
    raise ValueError('Unknown measurement dataset')


def period_basis(run, dataset, variant):
    if variant not in VARIANTS:
        raise ValueError('Unknown period response variant')
    n = len(run)
    if variant == 'static':
        return np.empty((n, 0))
    coordinate = period_coordinates(run, dataset)
    if variant == 'step':
        group = np.searchsorted(EDGES, coordinate, side='right') - 1
        full = (group[:, None] == np.arange(5)[None, :]).astype(float)
    else:
        full = np.column_stack([np.interp(coordinate, KNOTS, np.eye(5)[:, j]) for j in range(5)])
    if not np.allclose(full.sum(1), 1, atol=1e-14) or np.any(full < 0):
        raise ValueError('Period basis is not a partition of unity')
    return full[:, 1:]


def features(data, train, dataset, variant, normalization=None):
    train = np.asarray(train, bool)
    if train.shape != (len(data['run']),):
        raise ValueError('Training mask shape differs')
    base_norm = None if normalization is None else normalization['geometry']
    baseline, geometry_norm = original.response_features(data, train, base_norm)
    raw = period_basis(data['run'], dataset, variant)
    if normalization is None:
        center = raw[train].mean(0) if raw.shape[1] else np.empty(0)
        scale = raw[train].std(0) if raw.shape[1] else np.empty(0)
        if np.any(scale < 1e-6):
            raise ValueError('Insufficient training-era variation')
        normalization = dict(dataset=dataset, variant=variant, geometry=geometry_norm,
            period_center=center.tolist(), period_scale=scale.tolist(), clip=3.,
            utc_year_boundaries=YEARS, knots_unix=KNOTS.tolist(),
            mapping='Observed UTC run time' if dataset == 'flight' else 'Deterministic run-hash negative control; not real instrument time')
    if normalization['dataset'] != dataset or normalization['variant'] != variant:
        raise ValueError('Period normalization belongs to another calibration')
    if raw.shape[1]:
        values = np.clip((raw - np.asarray(normalization['period_center'])) / np.asarray(normalization['period_scale']), -3., 3.)
        baseline = np.concatenate([baseline, np.broadcast_to(values[:, None, :], (len(values), 4, values.shape[1]))], axis=2)
    if not np.isfinite(baseline).all():
        raise ValueError('Invalid period-dependent response features')
    return baseline, normalization


def lift(theta, nfeatures):
    theta = np.asarray(theta, float)
    old_features = (len(theta) - 39) // 4
    if len(theta) != 39 + 4 * old_features or old_features > nfeatures:
        raise ValueError('Invalid period-model initialization')
    output = np.zeros(39 + 4 * nfeatures)
    output[:39] = theta[:39]
    output[39:].reshape(4, nfeatures)[:, :old_features] = theta[39:].reshape(4, old_features)
    return output
