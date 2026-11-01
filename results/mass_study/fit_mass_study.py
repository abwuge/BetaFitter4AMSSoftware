#!/usr/bin/env python3
"""Flight-only constrained mass calibration. ROOT extraction/plots are retained.

SciPy L-BFGS-B is used because boundary isotope weights made the ROOT Minuit
exploration fail its convergence checks. No MC file or truth column is read.
"""
import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import uproot

P, B = 8, 4
EDGES = np.array([2.5, 3., 4., 5., 6.])
MASSES = {1: np.array([.938272, 1.875613]), 2: np.array([2.808392, 3.727379])}
NORMALIZATION = math.lgamma(3.5) - math.lgamma(3.) - .5 * math.log(6 * math.pi)
FIELDS = ['run', 'event', 'role', 'source', 'Z', 'u', 'richU', 'richOK', 'NaF',
          'features', 'rigidityK', 'erigidity', 'chis']
NAMES = ['raw', 'offset', 'mass_calibrated', 'rich_reference']


def subset(d, rows):
    return {k: v[rows] for k, v in d.items()}


def read(path):
    chunks = []
    with uproot.open(str(path)) as f:
        for a in f['flight'].iterate(FIELDS, step_size=100000, library='np'):
            r = a['rigidityK'][:, 8]
            relative = np.abs(a['erigidity'][:, 15] * r)
            cx, cy = a['chis'][:, 171], a['chis'][:, 172]
            keep = np.isin(a['Z'], [1, 2]) & np.isfinite(r) & (r >= 2.5) & (r < 6)
            keep &= np.isfinite(relative) & (relative > 0) & (relative < .20)
            keep &= np.isfinite(cx) & np.isfinite(cy) & (cx > 0) & (cy > 0) & (cx <= 10) & (cy <= 10)
            rr, ur = a['rigidityK'][:, 15], a['richU']
            rich = (a['richOK'] != 0) & (a['NaF'] != 0) & np.isfinite(rr) & (rr > 0)
            rich &= np.isfinite(ur) & (ur > 1/.97) & (ur < 1/.78)
            mr = np.full(len(r), np.nan)
            mr[rich] = a['Z'][rich] * rr[rich] * np.sqrt(ur[rich]**2 - 1)
            tag = np.full(len(r), -1, dtype=np.int32)
            tag[rich & (a['Z'] == 1) & (mr > .78) & (mr < 1.10)] = 0
            tag[rich & (a['Z'] == 1) & (mr > 1.65) & (mr < 2.05)] = 1
            tag[rich & (a['Z'] == 2) & (mr > 2.60) & (mr < 3.02)] = 0
            tag[rich & (a['Z'] == 2) & (mr > 3.52) & (mr < 3.94)] = 1
            d = {k: a[k][keep] for k in ['run', 'event', 'role', 'source', 'Z', 'u', 'features']}
            d.update(R=r[keep], Rref=rr[keep], relR=relative[keep], richU=ur[keep],
                     richOK=rich[keep], massRich=mr[keep], tag=tag[keep],
                     bin=(np.searchsorted(EDGES, r[keep], side='right') - 1).astype(np.int32))
            chunks.append(d)
    d = {k: np.concatenate([v[k] for v in chunks]) for k in chunks[0]}
    r0, r1 = np.unique(d['run'][d['role'] == 0]), np.unique(d['run'][d['role'] == 1])
    if np.intersect1d(r0, r1).size:
        raise ValueError('Train/validation run overlap')
    ids = np.rec.fromarrays([d['run'], d['event']])
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicate selected event')
    print('Selected={} train_runs={} validation_runs={}'.format(len(ids), len(r0), len(r1)), flush=True)
    return d


def w68(v):
    return float(np.diff(np.percentile(v, [16, 84]))[0]/2) if len(v) else float('nan')


def rms(v):
    return float(np.sqrt(np.mean(np.square(v))))


def state(d):
    center = np.median(d['features'], axis=0)
    q = np.percentile(d['features'], [25, 75], axis=0)
    scale = np.maximum(1e-6, .741301109 * (q[1]-q[0]))
    center[0], scale[0] = 0, 1
    return {'center': center.tolist(), 'scale': scale.tolist()}


def features(d, model):
    x = np.clip((d['features']-model['center'])/model['scale'], -4, 4)
    x[:, 0] = 1
    return x


class Density:
    def __init__(self, d, feature_state, charge):
        # Explicit allowlist. RICH/labels/run/source never enter this objective.
        self.x = features({'features': d['features']}, feature_state)
        self.u = d['u'].astype(float)
        self.b = d['bin']
        self.means = np.sqrt(1 + (MASSES[charge][None, :]/(charge*d['R'][:, None]))**2)
        sensitivity = (self.means**2 - 1)/self.means
        self.tracker = (4./6) * sensitivity**2 * d['relR'][:, None]**2

    def log_components(self, parameters, corrected=None):
        u = self.u - self.x.dot(parameters[:P]) if corrected is None else corrected
        r = u[:, None] - self.means
        s2 = np.exp(2*parameters[P:P+B])[self.b]
        v = s2[:, None] + self.tracker
        lp = NORMALIZATION - .5*np.log(v) - 3.5*np.log1p(r*r/(6*v))
        return lp, r, v, s2

    def __call__(self, parameters):
        lp, r, v, s2 = self.log_components(parameters)
        f = parameters[P+B:][self.b]
        with np.errstate(divide='ignore'):
            lm = np.logaddexp(np.log1p(-f)+lp[:, 0], np.log(f)+lp[:, 1])
        dr = np.exp(lp-lm[:, None])
        weights = dr * np.column_stack((1-f, f))
        denominator = 6*v + r*r
        gc = -np.sum(weights*7*r/denominator, axis=1)
        gs = np.sum(weights*(s2[:, None]/v)*(1-7*r*r/denominator), axis=1)
        gf = dr[:, 0]-dr[:, 1]
        g = np.concatenate((self.x.T.dot(gc), np.bincount(self.b, weights=gs, minlength=B),
                            np.bincount(self.b, weights=gf, minlength=B))) / len(self.u)
        value = -float(np.mean(lm)) + .0005*np.sum((parameters[1:P]/.03)**2)
        g[1:P] += .001*parameters[1:P]/.03**2
        return value, g

    def event_nll(self, parameters, corrected):
        lp = self.log_components(parameters, corrected)[0]
        f = parameters[P+B:][self.b]
        with np.errstate(divide='ignore'):
            return -np.logaddexp(np.log1p(-f)+lp[:, 0], np.log(f)+lp[:, 1])

    def score(self, parameters, corrected):
        lp = self.log_components(parameters, corrected)[0]
        return lp[:, 1]-lp[:, 0]


def fit(d, feature_state, z, mode, initial=None, starts=1):
    density = Density(d, feature_state, z)
    x0 = np.zeros(P+2*B) if initial is None else np.asarray(initial, dtype=float).copy()
    if initial is None:
        x0[P:P+B] = math.log(.04 if z == 1 else .025)
        x0[P+B:] = .02 if z == 1 else .85
    bounds = [(-.08, .08)]*P + [(math.log(.003), math.log(.15))]*B + [(0., 1.)]*B
    for j in range(P):
        if mode == 'raw' or (mode == 'offset' and j > 0):
            bounds[j] = (0., 0.)
            x0[j] = 0
        elif mode == 'fixed':
            bounds[j] = (x0[j], x0[j])
    records, solutions = [], []
    for start in range(starts):
        start_x = x0.copy()
        if start:
            start_x[P:P+B] += math.log(1.2)
            start_x[P+B:] = .15 if z == 1 else .5
            if mode == 'full':
                start_x[:P] = 0
                start_x[0] = -.015
        result = minimize(density, start_x, jac=True, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter': 2000, 'maxls': 60, 'ftol': 1e-12, 'gtol': 1e-6})
        row = {'start': start, 'success': bool(result.success), 'status': int(result.status),
               'objective': float(result.fun), 'iterations': int(result.nit), 'message': str(result.message)}
        records.append(row)
        print('fit Z={} mode={} {}'.format(z, mode, row), flush=True)
        if result.success and np.isfinite(result.fun) and np.all(np.isfinite(result.x)):
            solutions.append(result)
    if not solutions:
        raise RuntimeError('All optimizer starts failed for Z{} {}'.format(z, mode))
    best = min(solutions, key=lambda r: r.fun)
    p = best.x
    boundary = {'calibration': np.flatnonzero(np.abs(p[:P]) > .079).tolist(),
                'fraction': np.flatnonzero((p[P+B:] < 1e-5) | (p[P+B:] > 1-1e-5)).tolist()}
    print('chosen Z={} {} boundaries={}'.format(z, mode, boundary), flush=True)
    return dict(feature_state, parameters=p.tolist(), objective=float(best.fun),
                optimizer=records, boundary=boundary, n_training=len(d['u']))


def ref_u0(d):
    return np.sqrt(1+(d['Rref']/d['R'])**2*(d['richU']**2-1))


def fit_reference(d, feature_state, z):
    q = subset(d, d['richOK'])
    if len(q['u']) < 500:
        raise ValueError('Too few RICH references')
    x = features(q, feature_state)
    target = q['u'] - ref_u0(q)
    delta = 3*w68(target)
    def loss(p):
        r = target-x.dot(p)
        a = np.abs(r)
        value = np.mean(np.where(a < delta, .5*r*r, delta*(a-.5*delta)))
        grad = -x.T.dot(np.clip(r, -delta, delta))/len(r)
        value += .0005*np.sum(p[1:]**2)
        grad[1:] += .001*p[1:]
        return value, grad
    result = minimize(loss, np.zeros(P), jac=True, method='L-BFGS-B', bounds=[(-.08, .08)]*P,
                      options={'maxiter': 2000, 'maxls': 60, 'ftol': 1e-14, 'gtol': 1e-9})
    if not result.success:
        raise RuntimeError('RICH reference optimization failed: '+str(result.message))
    print('RICH fit Z={} N={} objective={:.9g}'.format(z, len(target), result.fun), flush=True)
    return result.x, len(target)


def auc(score, tag):
    select = tag >= 0
    s, y = score[select], tag[select]
    n1, n0 = np.sum(y == 1), np.sum(y == 0)
    if not n1 or not n0:
        return float('nan')
    order = np.argsort(s, kind='mergesort')
    s, y = s[order], y[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1])+1]
    ends = np.r_[starts[1:], len(s)]
    ranks = np.repeat((starts+1+ends)/2., ends-starts)
    return float((np.sum(ranks[y == 1]) - n1*(n1+1)/2.)/(float(n0)*n1))


def gradient_check(d, feature_state, z):
    small = subset(d, np.arange(min(4000, len(d['u']))))
    density = Density(small, feature_state, z)
    p = np.r_[np.zeros(P), np.full(B, math.log(.04 if z == 1 else .025)), np.full(B, .3)]
    g = density(p)[1]
    errors = []
    for j in range(len(p)):
        plus, minus = p.copy(), p.copy()
        h = 1e-6
        plus[j] += h
        minus[j] -= h
        fd = (density(plus)[0]-density(minus)[0])/(2*h)
        errors.append(abs(g[j]-fd)/(1+abs(fd)))
    maximum = float(max(errors))
    print('Gradient check Z={} max_relative_error={:.4g}'.format(z, maximum), flush=True)
    if maximum > 2e-4:
        raise ValueError('Gradient check failed')
    return maximum


def write_tsv(path, rows):
    with path.open('w') as f:
        keys = list(rows[0])
        f.write('\t'.join(keys)+'\n')
        for row in rows:
            f.write('\t'.join(str(row[k]) for k in keys)+'\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    dest = args.directory
    d = read(dest/'flight.root')
    train = {}; valid = {}; all_models = {}; checks = {}; reference_counts = {}
    # Fit every model before any holdout metric is evaluated.
    for z in [1, 2]:
        original = subset(d, (d['Z'] == z) & (d['role'] == 0))
        va = subset(d, (d['Z'] == z) & (d['role'] == 1))
        idx = np.random.RandomState(2981+z).permutation(len(original['u']))[:60000]
        tr = subset(original, idx)
        # Explicitly strip reference measurements, tags and identities from fit inputs.
        tr = {k: tr[k] for k in ['features', 'u', 'R', 'relR', 'bin']}
        st = state(tr)
        checks[str(z)] = gradient_check(tr, st, z)
        m0 = fit(tr, st, z, 'raw', starts=2)
        m1 = fit(tr, st, z, 'offset', m0['parameters'])
        m2 = fit(tr, st, z, 'full', m1['parameters'], starts=2)
        refcoef, nr = fit_reference(original, st, z)
        p = np.array(m0['parameters']); p[:P] = refcoef
        mr = fit(tr, st, z, 'fixed', p)
        ti = dict(tr); x = features(tr, st)
        ti['u'] = tr['u']+.010*x[:, 1]+.006*x[:, 6]
        mi = fit(ti, st, z, 'full', m2['parameters'])
        ri = dict(tr); factor = np.exp(.03*x[:, 1])
        ri['R'], ri['relR'] = tr['R']*factor, tr['relR']*factor
        mri = fit(ri, st, z, 'full', m2['parameters'])
        all_models[str(z)] = dict(raw=m0, offset=m1, mass_calibrated=m2,
                                  rich_reference=mr, injection=mi, rigidity_injection=mri)
        reference_counts[str(z)] = nr
        train[z], valid[z] = tr, va
    freeze = {'truth_access': False, 'frozen_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'optimizer': 'SciPy L-BFGS-B, analytic gradient, physical fraction bounds [0,1]',
              'gradient_checks': checks, 'rich_training_counts': reference_counts,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_audit_sha256': hashlib.sha256((dest/'flight.root.sources.tsv').read_bytes()).hexdigest(),
              'models': all_models}
    (dest/'frozen_models.json').write_text(json.dumps(freeze, indent=2, allow_nan=False))
    print('FROZEN all models, before holdout metrics', flush=True)
    model_rows = []
    for zs, models in all_models.items():
        for name, m in models.items():
            for j, p in enumerate(m['parameters']):
                model_rows.append(dict(Z=int(zs), model=name, parameter=j, value=p,
                                       feature_center=m['center'][j] if j < P else 0,
                                       feature_scale=m['scale'][j] if j < P else 1,
                                       fit_status=0, objective=m['objective']))
    write_tsv(dest/'models.tsv', model_rows)
    stats, bins, closures, rclosures, boot, preds = [], [], [], [], [], []
    for z in [1, 2]:
        va = valid[z]; models = all_models[str(z)]; x = features(va, models['raw'])
        density = Density(va, models['raw'], z)
        uv, mass2, scores, nll = {}, {}, {}, {}
        ru0 = np.full(len(va['u']), np.nan)
        ru0[va['richOK']] = ref_u0(subset(va, va['richOK']))
        for name in NAMES:
            p = np.array(models[name]['parameters'])
            uv[name] = va['u']-x.dot(p[:P])
            mass2[name] = (z*va['R'])**2*(uv[name]**2-1)
            scores[name] = density.score(p, uv[name])
            nll[name] = density.event_nll(p, uv[name])
            res = (uv[name]-ru0)[va['richOK']]
            light, heavy = mass2[name][va['tag'] == 0], mass2[name][va['tag'] == 1]
            sep = (np.median(heavy)-np.median(light))/math.sqrt(w68(light)**2+w68(heavy)**2)
            stats.append(dict(Z=z, method=name, N_train=len(train[z]['u']), N_holdout=len(va['u']),
                              N_rich=len(res), N_tag_light=len(light), N_tag_heavy=len(heavy),
                              mean_nll=float(np.mean(nll[name])), rich_residual_median=float(np.median(res)),
                              rich_residual_W68=w68(res), tag_AUC=auc(scores[name], va['tag']),
                              light_mass2_median=float(np.median(light)), light_mass2_W68=w68(light),
                              heavy_mass2_median=float(np.median(heavy)), heavy_mass2_W68=w68(heavy),
                              separation_W68=float(sep)))
            for b in range(B):
                q = va['bin'] == b
                bins.append(dict(Z=z, method=name, R_low=EDGES[b], R_high=EDGES[b+1],
                                 N_holdout=int(q.sum()), mean_nll=float(np.mean(nll[name][q])),
                                 train_width_scale=math.exp(p[P+b]), train_heavy_fraction=float(p[P+B+b]),
                                 N_tag_light=int(np.sum(va['tag'][q] == 0)), N_tag_heavy=int(np.sum(va['tag'][q] == 1)),
                                 tag_score_AUC=auc(scores[name][q], va['tag'][q]),
                                 tag_mass2_AUC=auc(mass2[name][q], va['tag'][q])))
        injection = .010*x[:, 1]+.006*x[:, 6]
        restored = va['u']+injection-x.dot(np.array(models['injection']['parameters'][:P]))
        error = restored-uv['mass_calibrated']
        closures.append(dict(Z=z, N_holdout=len(error), injected_RMS=rms(injection),
                             recovered_difference_RMS=rms(error), residual_W68=w68(error)))
        factor = np.exp(.03*x[:, 1]); ir = va['R']*factor
        iu = va['u']-x.dot(np.array(models['rigidity_injection']['parameters'][:P]))
        before = mass2['mass_calibrated']*factor**2
        after = (z*ir)**2*(iu**2-1)
        e0, e1 = before-mass2['mass_calibrated'], after-mass2['mass_calibrated']
        rclosures.append(dict(Z=z, N_holdout=len(e0), R_relative_injection_RMS=rms(factor-1),
                              mass2_change_RMS_before=rms(e0), mass2_change_RMS_after=rms(e1),
                              mass2_change_W68_before=w68(e0), mass2_change_W68_after=w68(e1)))
        forced = MASSES[z][np.argmin(np.abs(mass2['raw'][:, None]-MASSES[z][None, :]**2), axis=1)]
        pred = dict(Z=va['Z'].astype(np.int32), source=va['source'], run=va['run'], event=va['event'],
                    R=va['R'], Rref=va['Rref'], rawU=uv['raw'], offsetU=uv['offset'],
                    calibratedU=uv['mass_calibrated'], referenceU=uv['rich_reference'], richU=va['richU'],
                    richU0=ru0, massRich=va['massRich'], richOK=va['richOK'].astype(np.int32), tag=va['tag'],
                    bin=va['bin'], mass2Raw=mass2['raw'], mass2Offset=mass2['offset'],
                    mass2Calibrated=mass2['mass_calibrated'], mass2Reference=mass2['rich_reference'],
                    scoreRaw=scores['raw'], scoreCalibrated=scores['mass_calibrated'],
                    scoreReference=scores['rich_reference'], nearestMass=forced,
                    nearestU=np.sqrt(1+(forced/(z*va['R']))**2), injectedU=va['u']+injection,
                    restoredU=restored, Rinjected=ir, mass2RigidityInjectionBefore=before,
                    mass2RigidityInjectionAfter=after)
        preds.append(pred)
        runs = np.unique(va['run'])
        reference_rows = {int(r): np.flatnonzero((va['run'] == r) & va['richOK']) for r in runs}
        tag_rows = {int(r): np.flatnonzero((va['run'] == r) & (va['tag'] >= 0)) for r in runs}
        rng = np.random.RandomState(8831+z)
        for rep in range(200):
            draw = rng.choice(runs, len(runs), replace=True)
            qr = np.concatenate([reference_rows[int(r)] for r in draw])
            qt = np.concatenate([tag_rows[int(r)] for r in draw])
            narrow = 1-w68((uv['mass_calibrated']-ru0)[qr])/w68((uv['raw']-ru0)[qr])
            da = auc(scores['mass_calibrated'][qt], va['tag'][qt])-auc(scores['raw'][qt], va['tag'][qt])
            boot.append(dict(Z=z, replicate=rep, rich_width_narrowing=narrow, auc_difference=da))
        print('Z={} holdout={} RICH={} tags={}/{} W68 raw={:.6g} mass={:.6g} RICH={:.6g}'.format(
            z, len(va['u']), int(va['richOK'].sum()), int((va['tag'] == 0).sum()), int((va['tag'] == 1).sum()),
            stats[-4]['rich_residual_W68'], stats[-2]['rich_residual_W68'], stats[-1]['rich_residual_W68']), flush=True)
    with uproot.recreate(str(dest/'predictions.root')) as f:
        f['predictions'] = {k: np.concatenate([p[k] for p in preds]) for k in preds[0]}
    for name, rows in [('metrics', stats), ('bins', bins), ('injection', closures),
                       ('rigidity_injection', rclosures), ('bootstrap', boot)]:
        write_tsv(dest/(name+'.tsv'), rows)
    print('DONE mass study; no MC truth accessed', flush=True)


if __name__ == '__main__':
    main()
