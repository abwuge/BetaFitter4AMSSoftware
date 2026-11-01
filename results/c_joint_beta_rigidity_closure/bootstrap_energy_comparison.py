#!/usr/bin/env python3
"""Paired, isotope-stratified source-file bootstrap with frozen predictions."""
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr
from plot_beta_residual_trends import EDGES

HERE = Path(__file__).resolve().parent / 'expanded_validation'
OUT = HERE / 'energy_comparison'
METHODS = ['raw', 'joint_anchored', 'beta_only_frozen', 'anchored_energy']


def fit_core(residual):
    """Same conditional likelihood using sufficient statistics and a simplex fit."""
    low, median, high = np.percentile(residual, [16, 50, 84])
    width = (high - low) / 2
    left, right = median - 2.5 * width, median + 2.5 * width
    x = residual[(residual >= left) & (residual <= right)]
    average, variance = np.mean(x), np.var(x)

    def nll(p):
        mu, log_sigma = p
        if not (left <= mu <= right and np.log(width / 5) <= log_sigma <= np.log(width * 5)):
            return np.inf
        sigma = np.exp(log_sigma)
        norm = ndtr((right - mu) / sigma) - ndtr((left - mu) / sigma)
        if norm <= 0:
            return np.inf
        return log_sigma + np.log(norm) + .5 * (variance + (average - mu)**2) / sigma**2

    initial = np.array([median, np.log(width)])
    simplex = np.array([initial, initial + [.1 * width, 0], initial + [0, .05]])
    fit = minimize(nll, initial, method='Nelder-Mead',
                   options={'xatol': 1e-6, 'fatol': 1e-10, 'maxiter': 3000,
                            'initial_simplex': simplex})
    if not fit.success:
        raise RuntimeError('%s; n=%d median=%g width=%g fit=%s fun=%g' %
                           (fit.message, len(x), median, width, fit.x, fit.fun))
    return {'mean': fit.x[0], 'sigma': np.exp(fit.x[1])}


def replicate(seed):
    rng = np.random.RandomState(seed)
    idx = np.concatenate([groups[l][j] for l in (0, 1)
                          for j in rng.randint(len(groups[l]), size=len(groups[l]))])
    widths = np.diff(np.percentile(du[:, idx], [16, 84], axis=1), axis=0)[0] / 2
    fits = np.empty((4, 15, 2))
    for b in range(15):
        selected = idx[bins[idx] == b]
        for m in range(4):
            try:
                f = fit_core(dbeta[m, selected])
            except RuntimeError as error:
                raise RuntimeError('seed=%d method=%d bin=%d: %s' % (seed, m, b, error))
            fits[m, b] = [f['mean'], f['sigma']]
    return widths, fits


def main():
    global groups, du, dbeta, bins
    with np.load(HERE / 'predictions.npz') as f:
        d = dict(f)
    with np.load(OUT / 'energy_predictions.npz') as f:
        for key in ['run', 'event', 'label', 'source']:
            assert np.array_equal(d[key], f[key]), key
        d['u_anchored_energy'] = f['u_anchored_energy']
    good = d['truth_valid'] & np.isfinite(d['u_true']) & (d['u_true'] > 0)
    truth = d['u_true'][good]
    u = np.array([d['u_' + m][good] for m in METHODS])
    du = u - truth
    dbeta = 1000 * (1 / u - 1 / truth)
    bins = np.searchsorted(EDGES, 1 / truth, side='right') - 1
    saved = json.loads((OUT / 'probe_results.json').read_text())
    for row in saved['expanded']['beta_trend']:
        m = METHODS.index(row['method'])
        b = int(np.argmin(abs((EDGES[:-1] + EDGES[1:]) / 2 - row['beta_center'])))
        fit = fit_core(dbeta[m, bins == b])
        for metric in ['mean', 'sigma']:
            assert abs(fit[metric] - row[metric]) < 1e-4, (m, b, metric)
    labels, sources = d['label'][good], d['source'][good]
    groups = {l: [np.flatnonzero((labels == l) & (sources == s))
                  for s in np.unique(sources[labels == l])] for l in (0, 1)}
    values = []
    with mp.Pool(6) as pool:
        for i, value in enumerate(pool.imap(replicate, range(292900, 293200)), 1):
            values.append(value)
            if i % 25 == 0:
                print('bootstrap %d/300' % i, flush=True)
    w = np.array([v[0] for v in values])
    fits = np.array([v[1] for v in values])
    gain = 100 * (1 - fits[:, :, :, 1] / fits[:, 0:1, :, 1])
    overall_gain = 100 * (1 - w / w[:, 0:1])
    np.savez_compressed(str(OUT / 'bootstrap_replicates.npz'),
                        methods=np.array(METHODS), w68=w, fits=fits,
                        sigma_gain=gain, overall_gain=overall_gain)
    result = json.loads((OUT / 'probe_results.json').read_text())
    for r in result['expanded']['beta_trend']:
        m = METHODS.index(r['method'])
        b = int(np.argmin(abs((EDGES[:-1] + EDGES[1:]) / 2 - r['beta_center'])))
        r['mean_bootstrap_error'] = float(np.std(fits[:, m, b, 0], ddof=1))
        r['sigma_bootstrap_error'] = float(np.std(fits[:, m, b, 1], ddof=1))
        r['sigma_gain_error'] = float(np.std(gain[:, m, b], ddof=1))
    summary = {}
    for m, method in enumerate(METHODS):
        summary[method] = {
            'gain_lo95': float(np.percentile(overall_gain[:, m], 2.5)),
            'gain_hi95': float(np.percentile(overall_gain[:, m], 97.5)),
            'gain_error': float(np.std(overall_gain[:, m], ddof=1)),
            'mean_bin_sigma_gain': float(np.mean([
                100 * (1 - r['sigma'] / raw['sigma'])
                for r, raw in zip([r for r in result['expanded']['beta_trend'] if r['method'] == method],
                                  [r for r in result['expanded']['beta_trend'] if r['method'] == 'raw'])])),
        }
    paired = {}
    for m in [1, 2]:
        v = 100 * (1 - w[:, 3] / w[:, m])
        paired[METHODS[m]] = dict(zip(['lo95', 'median', 'hi95'],
                                      np.percentile(v, [2.5, 50, 97.5]).tolist()))
    result['bootstrap'] = {'replicates': 300, 'seed_start': 292900,
                           'unit': 'source file, stratified by isotope; paired across methods',
                           'trend_errors': 'one bootstrap standard deviation; core windows refitted',
                           'model': 'fixed; training uncertainty excluded',
                           'summary': summary, 'energy_gain_vs_other': paired}
    (OUT / 'probe_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['bootstrap'], indent=2), flush=True)


if __name__ == '__main__':
    main()
