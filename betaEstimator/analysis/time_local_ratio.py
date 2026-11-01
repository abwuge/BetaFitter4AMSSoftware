#!/usr/bin/env python3
"""Adjacent-period fits with the supplied carbon energy-spectrum ratio."""
import argparse
import copy
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
import time_local_beta as local

ROOT = local.w.rt.ROOT
OLD = local.OUT
OUT = ROOT/'results/time_local_ratio'
MASSES = (11.17793, 12.11221)
CACHE = None


def dump(path, value):
    local.w.ev.dump(path, value)


def energy(r, i):
    return (np.hypot(6*np.asarray(r), MASSES[i])-MASSES[i])/(12+i)


def rigidity(e, i):
    t = (12+i)*np.asarray(e)
    return np.sqrt(t*(t+2*MASSES[i]))/6


def mapping(r):
    e = energy(r, 1)
    s = rigidity(e, 0)
    jac = (12/13.) * np.hypot(6*s, MASSES[0])/np.hypot(6*r, MASSES[1]) * r/s
    return s, jac


def ratio_fraction(data, reference, edges):
    """Infer an R-only C12 density, then integrate both isotope densities.

    The external ratio is interpreted as a selected-event dN/dE ratio.
    Equal relative selection efficiency is assumed. Measured R approximates
    latent R; this conversion does not unfold Tracker resolution.
    """
    x = np.asarray(reference['x']); y = np.asarray(reference['ratio'])
    maximum = rigidity(x[-1], 0)
    hist_edges = np.geomspace(2.5, 200., 97)
    counts, _ = np.histogram(data['R'], hist_edges)
    centers = np.sqrt(hist_edges[:-1]*hist_edges[1:])
    log_density = gaussian_filter1d(np.log((counts+.5)/np.diff(hist_edges)), 1.5, mode='nearest')
    slope = float(np.clip(np.polyfit(np.log(centers[:6]), log_density[:6], 1)[0], -8, 4))
    grid = np.geomspace(2.5, maximum, 2400)
    total = np.exp(np.interp(np.log(grid), np.log(centers), log_density))
    total[grid < centers[0]] = np.exp(log_density[0])*(grid[grid < centers[0]]/centers[0])**slope
    mapped, jac = mapping(grid)
    q = np.interp(energy(grid, 1), x, y)
    density12 = np.empty_like(grid)
    # The unobserved interval immediately below 2.5 GV uses the fitted local
    # power-law slope. Its normalization solves the constant-slope boundary.
    boundary = total[0]/(1+q[0]*jac[0]*(mapped[0]/grid[0])**slope)
    for j in range(len(grid)):
        lower = boundary*(mapped[j]/grid[0])**slope if mapped[j] < grid[0] else np.interp(mapped[j], grid[:j], density12[:j])
        density12[j] = total[j]-q[j]*jac[j]*lower
    if not np.all(density12 > 0):
        raise RuntimeError('External ratio incompatible with positive R-only density')
    density13 = total-density12
    fraction = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        z = np.r_[lo, grid[(grid>lo)&(grid<hi)], hi]
        a = np.trapz(np.interp(z, grid, density13), z)
        b = np.trapz(np.interp(z, grid, total), z)
        fraction.append(a/b)
    return np.asarray(fraction), dict(R_grid=grid, total_density=total, C12_density=density12,
        C13_density=density13, boundary_slope=slope, boundary_C12_density=boundary,
        fraction_per_bin=fraction,
        assumptions=['External dN13/dE divided by dN12/dE for the selected population',
            'Equal relative isotope selection efficiency at the same Ekn',
            'Measured R used as latent-R approximation; no unfolding',
            'R-only density: 96 log bins, log-density Gaussian smoothing sigma=1.5 bins',
            'Below-2.5-GV boundary: first-six-bin power-law continuation'])


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'protocol.json').exists():
        raise FileExistsError('Existing external-ratio protocol')
    reference = json.loads((ROOT/'results/material_comparison/ratio_reference.json').read_text())
    upper = float(rigidity(reference['x'][-1], 0))
    blocks = [list(x) for x in local.w.BLOCKS[:2]] + [[20.,30.,40.,50.,upper]]
    plan = json.loads((OLD/'plan.json').read_text())
    dump(OUT/'protocol.json',dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        reference=reference, blocks=blocks, original_plan_sha256=local.w.ev.digest(OLD/'plan.json'),
        training_sha256=plan['training_sha256'], script_sha256=local.w.ev.digest(__file__),
        RICH_used_for_training=False, MC_used_for_training=False,
        fixed_ratio='Central external ratio, linearly interpolated in Ekn; no extrapolation',
        uncertainty='Graph errors retained; central-value refit only, no ratio-error propagation',
        coverage_GV=[2.5,upper], time_definition='Original 333 target periods, original guards and neighbour windows',
        low_mid_features='Reuse each original frozen feature state and fitted parameters as initialization',
        high_features='Same feature construction with four bins over the covered 20--62.55 GV domain',
        original_artifacts_overwritten=False))
    (OUT/'folds').mkdir(exist_ok=True)
    (OUT/'condor').mkdir(exist_ok=True)
    for worker in range(50):
        (OUT/'condor'/('fit_%03d.txt'%worker)).write_text('\n'.join(str(x['fold']) for x in plan['folds'][worker::50])+'\n')


def check():
    p = json.loads((OUT/'protocol.json').read_text())
    if local.w.ev.digest(__file__) != p['script_sha256']:
        raise RuntimeError('External-ratio training implementation changed')
    if local.w.ev.digest(OLD/'plan.json') != p['original_plan_sha256']:
        raise RuntimeError('Original time plan changed')
    return p


def fit_one(fold_id):
    global CACHE
    p = check(); dest = OUT/'folds'/('%04d'%fold_id);dest.mkdir(exist_ok=True)
    if (dest/'freeze.json').exists():
        return True
    if CACHE is None:
        if local.w.ev.digest(OLD/'training.npz') != p['training_sha256']:
            raise RuntimeError('Training input changed')
        CACHE = local.load(OLD/'training.npz')
    plan = json.loads((OLD/'plan.json').read_text())
    record = next(x for x in plan['folds'] if x['fold']==fold_id)
    maxdays = max(b['training_days'] for b in record['blocks'])
    train = local.training_mask(CACHE['run'],record['start'],record['stop'],maxdays)
    dump(dest/'training_runs.json', np.unique(CACHE['run'][train]))
    summaries = []
    for i, edge_list in enumerate(p['blocks']):
        edges = np.asarray(edge_list);local.w.set_edges(edges)
        train = local.training_mask(CACHE['run'],record['start'],record['stop'],record['blocks'][i]['training_days'])
        population = local.w.ev.sub(CACHE, train)
        fraction, conversion = ratio_fraction(population,p['reference'],edges)
        d = local.w.ev.sub(population,(population['R']>=edges[0])&(population['R']<edges[-1]))
        d['bin'] = np.searchsorted(edges,d['R'],side='right')-1
        counts = np.bincount(d['bin'],minlength=4)
        if counts.min() < 100:
            raise RuntimeError('Insufficient covered training bin '+str(counts))
        oldmodel = local.w.rt.arrays(json.loads((OLD/'folds'/('%04d'%fold_id)/('block_%d_model.json'%i)).read_text()))
        if i < 2:
            state = copy.deepcopy(oldmodel);matrix = local.tr.feature_matrix(d,state)
            initial = oldmodel['parameters'].copy()
        else:
            state,matrix = local.tr.features(d)
            initial = np.r_[np.zeros(matrix.shape[1]),oldmodel['parameters'][-4:]]
        density = local.w.rt.rc.MixtureDensity(d,matrix,fraction[d['bin']])
        n = matrix.shape[1]
        bounds = [(0.,0.)]+[(-.05,.05)]*(n-1)+[(math.log(.0015),math.log(.08))]*4
        fit = minimize(density.value_gradient,initial,jac=True,method='L-BFGS-B',bounds=bounds,
            options=dict(maxiter=20000,maxls=100,maxcor=30,ftol=1e-13,gtol=1e-6))
        if not fit.success:
            raise RuntimeError('External-ratio fit failed '+str(fit.message))
        correction = matrix.dot(fit.x[:n])
        means = [float(correction[d['bin']==b].mean()) for b in range(4)]
        if max(abs(x) for x in means)>1e-10:
            raise RuntimeError('Training zero-mean anchor failed')
        state['parameters']=fit.x
        state['anchor']['note']='External Ekn isotope ratio converted through training-only R density'
        dump(dest/('block_%d_model.json'%i),dict(state,rigidity_edges_GV=edges))
        summaries.append(dict(block=i,edges=edges,n=len(d['u']),counts=counts,fraction=fraction,
            objective=fit.fun,iterations=fit.nit,success=bool(fit.success),mean_correction_per_bin=means,
            conversion=conversion))
        print('FIT',fold_id,i,len(d['u']),fit.nit,fraction,flush=True)
    dump(dest/'training.json',dict(complete=True,fold=record,blocks=summaries))
    dump(dest/'freeze.json',dict(complete=True,protocol_sha256=local.w.ev.digest(OUT/'protocol.json'),
        model_sha256={str(i):local.w.ev.digest(dest/('block_%d_model.json'%i)) for i in range(3)},
        training_sha256=local.w.ev.digest(dest/'training.json')))
    if (dest/'failure.json').exists(): (dest/'failure.json').unlink()
    return True


def safe_fit(fold):
    try:
        return fit_one(fold)
    except Exception:
        dest=OUT/'folds'/('%04d'%fold);dest.mkdir(exist_ok=True)
        dump(dest/'failure.json',dict(complete=False,traceback=traceback.format_exc()))
        print('FAILED',fold,traceback.format_exc(),flush=True)
        return False


def evaluate():
    p=check();plan=json.loads((OLD/'plan.json').read_text());freezes={}
    for fold in plan['folds']:
        dest=OUT/'folds'/('%04d'%fold['fold']);freeze=json.loads((dest/'freeze.json').read_text())
        if freeze['protocol_sha256']!=local.w.ev.digest(OUT/'protocol.json'):raise RuntimeError('Protocol changed')
        for i,h in freeze['model_sha256'].items():
            if local.w.ev.digest(dest/('block_%s_model.json'%i))!=h:raise RuntimeError('Model changed')
        freezes[str(fold['fold'])]=local.w.ev.digest(dest/'freeze.json')
    dump(OUT/'freeze.json',dict(complete=True,folds=freezes,RICH_used_for_training=False,MC_used_for_training=False))
    a=local.load(OLD/'predictions.npz');features=local.load(local.OLD/'diagnostic_covariates.npz',['run','event','source','continuous','bars','relR'])
    for k in ['run','event','source']:
        if not np.array_equal(a[k],features[k]):raise RuntimeError('Evaluation join mismatch')
    covered=(a['R']>=p['blocks'][0][0])&(a['R']<p['blocks'][-1][-1])
    correction=np.full(len(a['run']),np.nan);seen=np.zeros(len(a['run']),np.int8)
    for fold in plan['folds']:
        rows=np.flatnonzero(covered&(a['fold']==fold['fold']));dest=OUT/'folds'/('%04d'%fold['fold'])
        train_runs=np.asarray(json.loads((dest/'training_runs.json').read_text()))
        if np.intersect1d(a['run'][rows],train_runs).size:raise RuntimeError('Target run in training')
        for i,edge_list in enumerate(p['blocks']):
            edges=np.asarray(edge_list);local.w.set_edges(edges);idx=rows[(a['R'][rows]>=edges[0])&(a['R'][rows]<edges[-1])]
            d=dict(u=a['raw_u'][idx],R=a['R'][idx],continuous=features['continuous'][idx].astype(float),bars=features['bars'][idx])
            d['bin']=np.searchsorted(edges,d['R'],side='right')-1
            state=local.w.rt.arrays(json.loads((dest/('block_%d_model.json'%i)).read_text()))
            x=local.tr.feature_matrix(d,state);correction[idx]=d['u']-x.dot(state['parameters'][:x.shape[1]])
        seen[rows]+=1
    if not np.all(seen[covered]==1) or not np.isfinite(correction[covered]).all():raise RuntimeError('Prediction coverage failed')
    a['ratio_corrected_u']=correction;a['reference_mask'] &= covered
    np.savez_compressed(str(OUT/'predictions.npz'),**a)
    from finish_material_comparison import stats
    metrics={}
    for naf,name in [(True,'NaF'),(False,'AGL')]:
        s=a['reference_mask']&(a['NaF'].astype(bool)==naf);ref=1/a['rich_beta'][s]-a['delay_u_C12'][s]
        metrics[name]={key:stats(a[value][s]-ref) for key,value in [('linear','raw_u'),('old_neighbour','local_corrected_u'),('external_ratio','ratio_corrected_u')]}
    dump(OUT/'verification.json',dict(passed=True,folds=len(freezes),covered_reference_n=int(covered.sum()),
        excluded_outside_external_coverage=int((~covered).sum()),target_run_overlap=0,each_event_predicted_once=True,
        prediction_sha256=local.w.ev.digest(OUT/'predictions.npz'),metrics=metrics))
    print('EXTERNAL_RATIO_VERIFIED',metrics,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','fit','evaluate'])
    parser.add_argument('--batch',type=Path);parser.add_argument('--fold',type=int);parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='evaluate':evaluate()
    else:
        ids=[args.fold] if args.fold is not None else [int(x) for x in args.batch.read_text().split()]
        with mp.Pool(min(args.workers,len(ids))) as pool:
            result=list(pool.imap_unordered(safe_fit,ids))
        if not all(result):raise SystemExit(1)
