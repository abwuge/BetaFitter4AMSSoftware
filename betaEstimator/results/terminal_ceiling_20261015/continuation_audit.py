#!/usr/bin/env python3
"""Resume audit from saved event predictions; never transfer event files."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json, hashlib, time
from pathlib import Path
import numpy as np
from audit_common import HERE, metric, write_json

def load_json(name):
    return json.loads((HERE/name).read_text())

def wm(x,w):
    ok=np.isfinite(x)&np.isfinite(w)&(w>0)
    x=x[ok];w=w[ok];order=np.argsort(x)
    q=np.interp([.16,.84],(np.cumsum(w[order])-.5*w[order])/w.sum(),x[order])
    return dict(n=len(x),w68=float((q[1]-q[0])/2),mean=float(np.average(x,weights=w)),rmse=float(np.sqrt(np.average(x*x,weights=w))))

start=time.monotonic()
out={'completed_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'scope':'Independent re-reading of saved predictions and source manifests; no retraining; no global upper-bound claim.'}
native=load_json('native_counter_analysis_summary.json')
out['native_covariance']={s:native['covariance'][s] for s in ['all_beta_gt06','R2p5_6','beta06_09']}
out['native_geometry']={s:native['geometry'][s] for s in ['all_beta_gt06','R2p5_6']}
out['native_extra_results']={s:{k:native['results'][s][k] for k in ['oracle_geometry_only','oracle_loss_only','transport_prediction_only']} for s in ['all_beta_gt06','R2p5_6']}
man=load_json('mixed_manifest.json')
out['mixed_manifest_counts']={'events':man['events'],'primary_valid':man['primary_valid'],'sources':len(man['sources']),'C13_example':next(s for s in man['sources'] if s.get('isotope')=='C13'),'C12_example':next(s for s in man['sources'] if s.get('isotope')=='C12')}
out['broad_by_true_beta']={}
for z in [1,2,6,8]:
    fit=load_json('fit_Z%d_measured_all_t640_l15_m160_f100_core.json'%z)
    out['broad_by_true_beta'][str(z)]={s:{k:r[k] for k in ['n','w68','std','rmse','mean','gain_vs_hybrid_pct','relative_beta_w68','previous_hybridBeta']} for s,r in fit['test'].items() if s.startswith('beta_')}
a=np.load(str(HERE/'extended_errors.npz'))
mask=a['role']==1;lab=a['label'][mask];weights=a['weights'][mask]
methods=['ols','measured_tof','simple_unrestricted','all_noRich','all_noRich_capacity','all_with_RICH','all_noRich_known_isotope']
out['test_composition_sensitivity']={}
for ratio in [.05,.1,.2,.5]:
    w=weights*np.where(lab==1,2*ratio,2*(1-ratio));res={}
    base=wm(a['ols'][mask],w)
    for name in methods:
        m=wm(a[name][mask],w);m['gain_w68_pct']=100*(1-m['w68']/base['w68']);m['gain_rmse_pct']=100*(1-m['rmse']/base['rmse']);res[name]=m
    out['test_composition_sensitivity'][str(ratio)]=res
out['extended_by_species']={}
for isotope in [0,1]:
    sel=mask&(a['label']==isotope);w=a['weights'][sel]
    out['extended_by_species']['C%d'%(12+isotope)]={name:wm(a[name][sel],w) for name in methods}
out['extended_true_beta_bins']={}
beta=1/a['u_true']
for lo,hi in zip([.6,.7,.8,.9,.95],[.7,.8,.9,.95,1.000001]):
    sel=mask&(beta>=lo)&(beta<hi)
    if sel.sum()<40: continue
    out['extended_true_beta_bins']['%.2f_%.2f'%(lo,hi)]={name:wm(a[name][sel],a['weights'][sel]) for name in methods}
# Re-evaluate direct beta errors, not only inverse-beta, on identical mixed events.
out['relative_beta']={name:wm(a['u_true'][mask]/(a['u_true'][mask]+a[name][mask])-1,weights) for name in methods}
# Verify disjoint source and run identities and comparison-replay counts.
B=np.load(str(HERE/'broad_Z6.npz'));N=np.load(str(HERE/'native_counter_collected.npz'));E=np.load(str(HERE/'native_counter_analysis_errors.npz'))
from audit_common import partitions
tr,va,te=partitions({k:B[k] for k in ['role','run','source']})
id0=B['run'].astype(np.uint64)*np.uint64(2**32)+B['event'].astype(np.uint64)
out['replay_checks']={'C12_unique_events':bool(len(np.unique(id0))==len(id0)),'C12_train_test_runs_disjoint':bool(not np.intersect1d(B['run'][tr|va],B['run'][te]).size),'native_event_order_equal':bool(np.array_equal(B['run'],N['run']) and np.array_equal(B['event'],N['event'])),'native_ols_error_replay_max':float(np.max(np.abs((B['u'][E['index']]-N['u_true'][E['index']])-E['ols']))),'native_decomposition_max':float(np.max(np.abs(E['ols']-E['oracle_transport']-E['oracle_noise_removed'])))}
out['scalar_baselines']={}
for z in [1,2,6,8]:
    x=load_json('bootstrap_broad_Z%d.json'%z)
    out['scalar_baselines'][str(z)]={s:dict(n=r['n'],sources=r['sources'],baseline=r['methods']['hybrid']) for s,r in x.items()}
out['beta_relative_previous']=load_json('beta_relative_summary.json')
out['seconds']=time.monotonic()-start
write_json(HERE/'continuation_audit.json',out)
print('CONTINUATION_AUDIT_COMPLETE',round(out['seconds'],2),flush=True)
print('CHECKS',json.dumps(out['replay_checks']),flush=True)
print('LOWR_COVARIANCE',json.dumps(out['native_covariance']['R2p5_6']),flush=True)
print('NATIVE_EXTRA_LOW',json.dumps(out['native_extra_results']['R2p5_6']),flush=True)
print('RELATIVE_BETA',json.dumps(out['relative_beta']),flush=True)
print('C13_FIRST',json.dumps(out['mixed_manifest_counts']['C13_example']),flush=True)
print('SHIFTS',json.dumps({r:{k:{f:m[f] for f in ['gain_w68_pct','mean']} for k,m in methods.items()} for r,methods in out['test_composition_sensitivity'].items()}),flush=True)
