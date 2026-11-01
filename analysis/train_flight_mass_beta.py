#!/usr/bin/env python3
"""Fit the carbon mass-constrained beta correction using flight measurements.

Training reads a measured-branch allowlist. RICH is opened only by the
evaluation stage, after all nominal and fraction-sensitivity fits are frozen.
"""
import argparse
import json
import math
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import uproot

import rich_transport as rt
import evaluate_rich_transport as ev

rc, er = rt.rc, rt.er
DEST = rt.ROOT / 'results/flight_mass_beta'
CACHE = Path('/tmp/betaEstimator-flight-mass-beta')
TRAIN_BRANCHES = ['run', 'event', 'tof_tl', 'tof_leng', 'tof_pos',
                  'tof_edep', 'tof_barid', 'tk_ql', 'tk_rigidityk', 'tk_pos',
                  'tk_qin', 'tk_q', 'tof_bpatt', 'tk_erigidity1', 'tk_chis1']
TRAIN_KEYS = {'run', 'event', 'source', 'u', 'R', 'relR', 'bin', 'continuous', 'bars'}


def sources():
    records = []
    for line in (rt.ROOT / 'config/data_sources.conf').read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            role, _, path = line.split()
            records.append((len(records), role, path))
    return records


def extract_part(task):
    source, role, path = task
    if role != 'train':
        raise ValueError('Training extraction received a validation source')
    part = CACHE / 'parts' / ('%04d.npz' % source)
    note = part.with_suffix('.json')
    signature = {'path': path, 'source': source, 'branches': TRAIN_BRANCHES,
                 'code_sha256': ev.digest(Path(__file__)),
                 'selection_sha256': ev.digest(Path(rt.__file__))}
    if part.exists() and note.exists():
        audit = json.loads(note.read_text())
        if audit['signature'] == signature and audit['complete']:
            with np.load(str(part), allow_pickle=False) as f:
                if set(f.files) != TRAIN_KEYS:
                    raise ValueError('Unexpected training fields')
                return {k: f[k] for k in f.files}, audit
    with uproot.open(path) as f:
        t = f['amstreea']
        a = t.arrays(TRAIN_BRANCHES, library='np')
        rows = np.flatnonzero(rt.core_selection(a))
        d = rt.center_features(ev.sub(a, rows))
        d.update(run=a['run'][rows], event=a['event'][rows],
                 source=np.full(len(rows), source, np.int32),
                 relR=np.abs(a['tk_erigidity1'][rows, 2, 1] * d['R']))
        audit = {'signature': signature, 'raw': int(t.num_entries),
                 'selected': len(rows), 'complete': True}
    temporary = part.with_suffix('.partial.npz')
    np.savez_compressed(str(temporary), **d)
    temporary.replace(part)
    ev.dump(note, audit)
    return d, audit


def extract(workers):
    DEST.mkdir(parents=True, exist_ok=True)
    (CACHE / 'parts').mkdir(parents=True, exist_ok=True)
    tasks = [s for s in sources() if s[1] == 'train']
    pieces, audit = [], []
    with mp.Pool(workers) as pool:
        for j, (d, record) in enumerate(pool.imap(extract_part, tasks), 1):
            pieces.append(d)
            audit.append(record)
            if j % 20 == 0 or j == len(tasks):
                print('training sources %d/%d; carbon events %d' %
                      (j, len(tasks), sum(a['selected'] for a in audit)), flush=True)
    data = ev.concatenate(pieces)
    ids = np.rec.fromarrays([data['run'], data['event']], names='run,event')
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicate training events')
    validation_runs = {int(Path(s[2]).name.split('_')[0]) for s in sources() if s[1] == 'validation'}
    if set(data['run'].tolist()) & validation_runs:
        raise ValueError('Training/validation run overlap')
    np.savez_compressed(str(CACHE / 'training.npz'), **data)
    ev.dump(DEST / 'extraction.json', {'complete': True, 'sources': audit,
            'source_list_sha256': ev.digest(rt.ROOT / 'config/data_sources.conf'),
            'n': len(ids), 'bins': np.bincount(data['bin'], minlength=4),
            'branches': TRAIN_BRANCHES, 'unique_events': True, 'disjoint_runs': True})


class RawMixture:
    """Fit four C13 fractions and four TOF widths from the raw training spectrum."""
    def __init__(self, data):
        self.data = data
        self.zero = np.zeros((len(data['u']), 1))
        self.base = rc.MixtureDensity(data, self.zero, np.full(len(data['u']), .1))

    def value_gradient(self, params):
        sigma, fraction = params[:4], params[4:]
        p = fraction[self.data['bin']]
        self.base.prior = p
        v, g = self.base.value_gradient(np.r_[0., sigma])
        logpdf = self.base.components(np.r_[0., sigma])[0]
        lm = np.logaddexp(np.log1p(-p) + logpdf[:, 0], np.log(p) + logpdf[:, 1])
        df = np.exp(logpdf[:, 0] - lm) - np.exp(logpdf[:, 1] - lm)
        gf = np.bincount(self.data['bin'], weights=df, minlength=4) / len(p)
        return v, np.r_[g[1:], gf]


def fit_raw(data):
    density = RawMixture(data)
    bounds = [(math.log(.0015), math.log(.08))] * 4 + [(.001, .999)] * 4
    solutions = []
    for fraction in (.05, .3, .7):
        fit = minimize(density.value_gradient, np.r_[np.full(4, math.log(.012)), np.full(4, fraction)],
                       jac=True, method='L-BFGS-B', bounds=bounds,
                       options={'maxiter': 2000, 'maxls': 100, 'ftol': 1e-13, 'gtol': 1e-7})
        if not fit.success:
            raise RuntimeError('Raw-spectrum fit failed: ' + str(fit.message))
        solutions.append(fit)
    best = min(solutions, key=lambda x: x.fun)
    return {'logsigma': best.x[:4], 'fraction': best.x[4:], 'objective': best.fun,
            'starts': [{'objective': f.fun, 'fraction': f.x[4:], 'iterations': f.nit} for f in solutions]}


def feature_matrix(data, state):
    base = rc.apply_feature_state(data, state['base_state'])
    energy = er.apply_standardizer(er.energy_raw(data, 'fitted', state['response']), state['standardizer'])
    derived = np.column_stack([energy, energy * energy])
    derived = np.clip((derived - state['derived_center']) / state['derived_scale'], -5, 5)
    raw = np.column_stack([base, derived])[:, 1:]
    anchor = state['anchor']
    projected = (raw - anchor['bin_centers'][data['bin']]).dot(anchor['null_projection']) / anchor['scale']
    return np.column_stack([np.zeros(len(raw)), projected])


def features(data):
    state = {'base_state': rc.make_feature_state(data, 'no_edep'),
             'response': er.response_state(data)}
    energy = er.energy_raw(data, 'fitted', state['response'])
    state['standardizer'] = er.standardizer(energy)
    z = er.apply_standardizer(energy, state['standardizer'])
    derived = np.column_stack([z, z * z])
    state['derived_center'] = np.mean(derived, axis=0)
    state['derived_scale'] = np.std(derived, axis=0)
    state['derived_scale'][state['derived_scale'] < 1e-6] = 1.
    derived = np.clip((derived - state['derived_center']) / state['derived_scale'], -5, 5)
    raw = np.column_stack([rc.apply_feature_state(data, state['base_state']), derived])[:, 1:]
    centers = np.array([np.mean(raw[data['bin'] == b], axis=0) for b in range(4)])
    centered = raw - centers[data['bin']]
    scale = np.std(centered, axis=0)
    keep = scale > 1e-8
    state['anchor'] = {'bin_centers': centers, 'null_projection': np.eye(raw.shape[1])[:, keep],
                       'scale': scale[keep], 'name': 'flight_per_bin_zero_mean',
                       'note': 'No synthetic isotope bags or known flight bag proportions are available.'}
    return state, feature_matrix(data, state)


def fit_correction(data, matrix, raw, fraction):
    prior = fraction[data['bin']]
    density = rc.MixtureDensity(data, matrix, prior)
    p = matrix.shape[1]
    fit = minimize(density.value_gradient, np.r_[np.zeros(p), raw['logsigma']],
                   jac=True, method='L-BFGS-B',
                   bounds=[(0., 0.)] + [(-.05, .05)] * (p - 1) +
                          [(math.log(.0015), math.log(.08))] * 4,
                   options={'maxiter': 5000, 'maxls': 100, 'maxcor': 30, 'ftol': 1e-13, 'gtol': 1e-6})
    if not fit.success:
        raise RuntimeError('Beta fit failed: ' + str(fit.message))
    correction = matrix.dot(fit.x[:p])
    return fit.x, {'success': bool(fit.success), 'iterations': fit.nit, 'objective': fit.fun,
                   'gradient_inf': np.max(np.abs(density.value_gradient(fit.x)[1][1:])),
                   'fraction': fraction, 'rms_correction': np.sqrt(np.mean(correction**2)),
                   'p99_abs_correction': np.percentile(np.abs(correction), 99),
                   'mean_correction_per_bin': [np.mean(correction[data['bin'] == b]) for b in range(4)]}


def dependency_hashes():
    paths = [Path(__file__), Path(rt.__file__), Path(ev.__file__),
             rc.__file__, er.__file__, rt.ROOT / 'config/data_sources.conf']
    return {str(Path(p).resolve().relative_to(rt.ROOT)): ev.digest(p) for p in paths}


def train():
    DEST.mkdir(parents=True, exist_ok=True)
    if (DEST / 'freeze.json').exists():
        raise RuntimeError('Models already frozen; evaluate without retraining')
    with np.load(str(CACHE / 'training.npz'), allow_pickle=False) as f:
        if set(f.files) != TRAIN_KEYS:
            raise ValueError('Unexpected training fields')
        data = {k: f[k] for k in f.files}
    extraction = json.loads((DEST / 'extraction.json').read_text())
    if not extraction['complete'] or extraction['n'] != len(data['u']):
        raise ValueError('Incomplete training extraction')
    raw = fit_raw(data)
    print('raw-spectrum fractions', raw['fraction'], flush=True)
    state, matrix = features(data)
    models = {'nominal': raw['fraction'],
              'fraction_minus_0.05': np.clip(raw['fraction'] - .05, .001, .999),
              'fraction_plus_0.05': np.clip(raw['fraction'] + .05, .001, .999)}
    summary = {}
    for name, fraction in models.items():
        params, record = fit_correction(data, matrix, raw, fraction)
        model = dict(state, parameters=params)
        ev.dump(DEST / (name + '_model.json'), model)
        summary[name] = record
        print('fitted', name, json.dumps(record, default=lambda x: x.tolist()), flush=True)
    ev.dump(DEST / 'training.json', {'n': len(data['u']), 'raw_spectrum': raw, 'models': summary,
            'features': matrix.shape[1], 'training_fields': sorted(TRAIN_KEYS),
            'supervision': 'known C12/C13 masses; training-only raw-spectrum fraction estimates',
            'RICH_used': False, 'MC_used': False,
            'fraction_status': 'model-dependent nuisance estimates, not independent abundance measurements'})
    freeze_models()


def freeze_models():
    summary = json.loads((DEST / 'training.json').read_text())
    models = summary['models']
    for name, record in models.items():
        if not record['success']:
            raise RuntimeError('Cannot freeze a failed fit')
        state = json.loads((DEST / (name + '_model.json')).read_text())
        if not np.isfinite(state['parameters']).all():
            raise RuntimeError('Non-finite fitted parameters')
    executed = DEST / 'executed_training_source.py'
    ev.dump(DEST / 'freeze.json', {'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'code_sha256': dependency_hashes(),
            'executed_training_source_sha256': ev.digest(executed) if executed.exists() else ev.digest(Path(__file__)),
            'model_sha256': {name: ev.digest(DEST / (name + '_model.json')) for name in models},
            'training_cache_sha256': ev.digest(CACHE / 'training.npz'),
            'material_sha256': ev.digest(rt.MAP),
            'reference': 'existing v2 RICH quality; 4<Q<8; 5<Q<7 sensitivity',
            'evaluation': 'existing 98 validation sources, never used to fit or select candidates',
            'primary': 'nominal', 'sensitivity': ['fraction_minus_0.05', 'fraction_plus_0.05']})
    print('FROZEN', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['extract', 'train', 'freeze'])
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.stage == 'extract':
        extract(args.workers)
    elif args.stage == 'train':
        train()
    else:
        freeze_models()


if __name__ == '__main__':
    main()
