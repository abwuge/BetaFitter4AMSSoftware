#!/usr/bin/env python3
"""Post-freeze response and feature-group audit on the existing RICH sample.

This does not retrain, select a replacement model, or open a new blind sample.
Group removal is an inference-time diagnostic, not a retrained ablation.
"""
import json
from pathlib import Path
import numpy as np
import time_local_beta as t


def quantiles(x):
    return np.percentile(x, [5, 25, 50, 75, 95]).tolist()


def stats(x):
    q = np.percentile(x, [16, 84])
    return dict(n=len(x), mean=float(np.mean(x)), median=float(np.median(x)),
                w68=float((q[1]-q[0])/2), std=float(np.std(x)),
                rmse=float(np.sqrt(np.mean(x*x))))


def centered(x, code):
    count = np.bincount(code)
    mean = np.bincount(code, weights=x)/np.maximum(count, 1)
    return x-mean[code]


def main():
    t.w.check_freeze()
    out = t.OUT
    p = t.load(out/'predictions.npz')
    c = t.load(t.OLD/'diagnostic_covariates.npz')
    for key in ('run', 'event', 'source'):
        assert np.array_equal(p[key], c[key]), key
    plan = json.loads((out/'plan.json').read_text())
    n = len(p['run'])
    group = {k: np.full(n, np.nan) for k in ('geometry', 'timing', 'energy')}
    response = {eps: np.full(n, np.nan) for eps in (1e-4, .001, .01)}
    sigma_raw = np.full(n, np.nan)
    fraction = np.full(n, np.nan)
    fold_fractions = [[], [], []]
    replay_max = 0.
    for j, f in enumerate(plan['folds']):
        rows = np.flatnonzero(p['fold'] == f['fold'])
        dest = out/'folds'/('%04d' % f['fold'])
        freeze = json.loads((dest/'freeze.json').read_text())
        training = json.loads((dest/'training.json').read_text())
        for bi, edges in enumerate(t.w.BLOCKS):
            model_path = dest/('block_%d_model.json' % bi)
            assert t.w.ev.digest(model_path) == freeze['model_sha256'][str(bi)]
            s = t.w.rt.arrays(json.loads(model_path.read_text()))
            record = training['blocks'][bi]
            fold_fractions[bi].extend(record['fraction'])
            ix = rows[(p['R'][rows] >= edges[0]) & (p['R'][rows] < edges[-1])]
            if not len(ix):
                continue
            t.w.set_edges(edges)
            d = dict(u=p['raw_u'][ix], R=p['R'][ix],
                     continuous=c['continuous'][ix].astype(float), bars=c['bars'][ix])
            d['bin'] = np.searchsorted(edges, d['R'], side='right')-1
            base = t.w.rt.rc.apply_feature_state(d, s['base_state'])
            z = t.w.rt.er.apply_standardizer(t.w.rt.er.energy_raw(d, 'fitted', s['response']), s['standardizer'])
            derived = np.clip((np.column_stack([z, z*z])-s['derived_center'])/s['derived_scale'], -5, 5)
            raw = np.column_stack([base, derived])[:, 1:]
            a = s['anchor']
            nf = len(a['scale'])+1
            coef = s['parameters'][:nf]
            raw_weights = a['null_projection'].dot(coef[1:]/a['scale'])
            names = list(s['base_state']['feature_names'][1:])+['energy_tof', 'energy_tracker', 'energy_tof_square', 'energy_tracker_square']
            assert len(names) == raw.shape[1]
            x = raw-a['bin_centers'][d['bin']]
            labels = np.array(['energy' if name.startswith('energy_') else 'timing' if 'time_residual' in name else 'geometry' for name in names])
            for key in group:
                select = labels == key
                group[key][ix] = x[:, select].dot(raw_weights[select])
            correction = sum(group[key][ix] for key in group)
            difference = np.max(np.abs(d['u']-correction-p['local_corrected_u'][ix]))
            replay_max = max(replay_max, difference)
            assert difference < 1e-11, difference
            for eps in response:
                shifted = dict(d, u=d['u']+eps)
                xm = t.tr.feature_matrix(shifted, s)
                xp = t.tr.feature_matrix(dict(d, u=d['u']-eps), s)
                response[eps][ix] = 1-(xm.dot(coef)-xp.dot(coef))/(2*eps)
            sigma_raw[ix] = np.exp(np.asarray(record['raw_spectrum']['logsigma']))[d['bin']]
            fraction[ix] = np.asarray(record['fraction'])[d['bin']]
        if j % 50 == 0:
            print('replayed fold', j, flush=True)
    assert all(np.isfinite(v).all() for v in group.values())
    assert all(np.isfinite(v).all() for v in response.values())
    methods = dict(raw=p['raw_beta_rich'], local=p['local_beta_rich'], constant_one=np.ones(n))
    for name, included in [('geometry_only', ['geometry']), ('timing_only', ['timing']),
                           ('energy_only', ['energy']), ('no_energy', ['geometry','timing']),
                           ('no_timing', ['geometry','energy']), ('no_geometry', ['timing','energy'])]:
        correction = sum(group[k] for k in included)
        methods[name] = 1/(p['raw_u']-correction+p['delay_u_C12'])
    r12, r13 = [np.sqrt(1+(m/(6*p['R']))**2) for m in t.w.rt.rc.MASSES]
    separation = r13-r12
    sensitivity = (r12*r12-1)/r12
    tracker_u_sd = sensitivity*c['relR']
    combined_sd = np.sqrt(1.5*sigma_raw*sigma_raw+tracker_u_sd*tracker_u_sd)
    corr = p['raw_u']-p['local_corrected_u']
    error0 = 1/p['raw_beta_rich']-1/p['rich_beta']
    error1 = 1/p['local_beta_rich']-1/p['rich_beta']
    rows = []
    bins = list(zip(t.w.R_EDGES[:-1], t.w.R_EDGES[1:]))
    selections = [('block', lo, hi) for lo,hi in [(2.5,6), (6,20), (20,200)]]+ [('bin',float(lo),float(hi)) for lo,hi in bins]
    for kind, lo, hi in selections:
        m = p['reference_mask'] & (p['R'] >= lo) & (p['R'] < hi)
        if not np.any(m):
            continue
        ref = p['rich_beta'][m]
        record = dict(kind=kind, low=lo, high=hi, n=int(m.sum()),
                      metrics={k: stats(v[m]-ref) for k,v in methods.items()},
                      time_slope_injection={str(eps): quantiles(v[m]) for eps,v in response.items()},
                      relR_quantiles=quantiles(c['relR'][m]),
                      isotope_u_separation_quantiles=quantiles(separation[m]),
                      tracker_predicted_u_sd_quantiles=quantiles(tracker_u_sd[m]),
                      combined_raw_model_u_sd_quantiles=quantiles(combined_sd[m]),
                      isotope_separation_over_combined_sd_quantiles=quantiles(separation[m]/combined_sd[m]),
                      mean_u_correction=float(corr[m].mean()),
                      group_mean_u_correction={k: float(v[m].mean()) for k,v in group.items()},
                      fraction_reference_quantiles=quantiles(fraction[m]),
                      C13_vs_C12_delay_difference=stats((p['delay_u_C13']-p['delay_u_C12'])[m]))
        for val in record['metrics'].values():
            val['gain_beta_pct'] = 100*(1-val['w68']/record['metrics']['raw']['w68'])
        # Approximate conditional response. RICH noise remains in the regressor;
        # compare to raw on the identical sample, not to an assumed slope of 1.
        # No cut on raw residual or on model agreement is used here.
        rb = np.searchsorted(np.geomspace(lo,hi,33), p['R'][m], side='right')-1
        code = ((p['fold'][m]*33+rb)*2+p['NaF'][m].astype(int)).astype(int)
        uref = centered(1/ref, code)
        denom = float(uref.dot(uref))
        slopes = {k: float(centered(1/v[m], code).dot(uref)/denom) for k,v in methods.items()}
        record['conditional_u_response'] = dict(slopes=slopes, variance_ref=float(np.mean(uref*uref)),
            local_over_raw=slopes['local']/slopes['raw'] if abs(slopes['raw'])>1e-12 else None,
            conditioning='target 15-day period, 32 log-spaced measured-R bins within each range, radiator; full selected sample')
        e, d = error0[m], corr[m]
        cov = float(np.mean((e-e.mean())*(d-d.mean())))
        record['inverse_beta_variance'] = dict(raw=float(np.var(e)), local=float(np.var(error1[m])),
            correction=float(np.var(d)), twice_covariance=2*cov,
            identity_error=float(np.var(error1[m])-(np.var(e)+np.var(d)-2*cov)))
        rows.append(record)
        if kind == 'block':
            print(json.dumps({k: record[k] for k in ['low','high','n','metrics','conditional_u_response','isotope_separation_over_combined_sd_quantiles','group_mean_u_correction']}), flush=True)
    quality = []
    for lo,hi in [(2.5,6),(6,20),(20,200)]:
        for radiator in [None, False, True]:
            for rlo,rhi in [(0,.20),(0,.08),(.08,.12),(.12,.20)]:
                m = p['reference_mask']&(p['R']>=lo)&(p['R']<hi)&(c['relR']>=rlo)&(c['relR']<rhi)
                if radiator is not None:
                    m &= p['NaF'] == radiator
                if m.sum()<100:
                    continue
                quality.append(dict(low=lo,high=hi,NaF=radiator,relR_low=rlo,relR_high=rhi,n=int(m.sum()),
                    metrics={k:stats(methods[k][m]-p['rich_beta'][m]) for k in ['raw','local','no_energy','energy_only']}))
    fraction_summary = []
    for bi, values in enumerate(fold_fractions):
        values = np.asarray(values)
        fraction_summary.append(dict(block=bi,number=len(values),lower_bound_count=int((values<=.001+1e-8).sum()),quantiles=quantiles(values)))
    result = dict(complete=True,posthoc=True,models_retrained=False,new_blind_test=False,
        n=n,folds=len(plan['folds']),replay_max_abs_u=replay_max,
        predictions_sha256=t.w.ev.digest(out/'predictions.npz'),
        covariates_sha256=t.w.ev.digest(t.OLD/'diagnostic_covariates.npz'),
        script_sha256=t.w.ev.digest(Path(__file__)),fraction_summary=fraction_summary,
        time_injection='t_i -> t_i +/- epsilon*(path_i-mean_path)/c; fixed geometry and measured deposits; the least-squares u shifts by +/- epsilon and orthogonal residuals are invariant. This tests a timing-mode response, not a full physical change of speed including dE/dx.',
        group_removal='Fixed frozen coefficients and training anchors; contributions add exactly to full correction. Not separately trained models and not independent causal attribution.',
        conditional_response_limit='RICH has finite resolution and selected tails; slopes are errors-in-variables diagnostics. High-R absolute values cannot alone establish true-speed fidelity.',
        rows=rows,quality=quality)
    t.w.ev.dump(out/'response_diagnostic.json',result)
    np.savez_compressed(str(out/'response_components.npz'),run=p['run'],event=p['event'],
        geometry_correction=group['geometry'],timing_correction=group['timing'],energy_correction=group['energy'],
        time_mode_response=response[1e-4])
    print('COMPLETE',out/'response_diagnostic.json',flush=True)


if __name__ == '__main__':
    main()
