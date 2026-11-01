#!/usr/bin/env python3
"""Version-2 common-event timing study after explicit sentinel lineage audit.

Old diagnostic outputs are intentionally retained. Only these v2 directories
are the validity-corrected run-separated comparisons.
"""
import argparse
import json
import numpy as np
import research_flight_timing_redundancy as r
ROOT=r.ROOT
PARENT=r.OUT
LINEAGE=ROOT/'results/flight_calibration_research/timing_lineage'


def prepare(baseline):
    out=ROOT/'results/flight_calibration_research'/('valid_timing_'+baseline)
    out.mkdir(exist_ok=True)
    with np.load(str(PARENT/'measured.npz'),allow_pickle=False) as f:data={k:f[k] for k in f.files}
    with np.load(str(LINEAGE/'measured.npz'),allow_pickle=False) as f:
        runs=f['run'];events=f['event'];new=f['tof_tl_NewCalib'][:,0,0,0,0,0,:].astype(float)
    lookup={key:i for i,key in enumerate(zip(runs,events))}
    if len(lookup)!=len(runs):raise ValueError('Duplicate lineage identities')
    ix=np.array([lookup[key] for key in zip(data['run'],data['event'])],int);new=new[ix]
    native=data['t'].copy()
    valid_native=np.isfinite(native).all(axis=1)&(np.abs(native)<1000).all(axis=1)&(native!=-1).all(axis=1)&(native!=0).all(axis=1)
    valid_new=np.isfinite(new).all(axis=1)&(np.abs(new)<1000).all(axis=1)&(new!=0).all(axis=1)&(new!=-1).all(axis=1)
    common=valid_native&valid_new
    data['t']=native if baseline=='native' else new
    data=r.subset(data,common);np.savez_compressed(str(out/'measured.npz'),**data)
    r.write(out/'protocol.json',dict(version=2,baseline=baseline,
        branch='tof_tl' if baseline=='native' else 'tof_tl_NewCalib[:,0,0,0,0,0,:]',
        event_set='Exact same valid-event intersection for native and NewCalib LT comparisons; preserve original three whole-run cohorts',
        validity='Four times finite, abs<1000 ns and unequal to missing sentinels -1/0; no inverse-speed or residual quality cut added',
        sentinel_evidence='amsd69n/analysis_amsd69n.C line3359 initializes tof_tl=-1; only IsGoodTime() overwrites it at3399-3400. NewCalib has zero initial/missing entries.',
        baseline_semantics='native means already standard-calibrated MIT input, not raw uncalibrated time; newcalib means separate custom LT calibration',
        parent_events=len(common),native_invalid=int((~valid_native).sum()),newcalib_invalid=int((~valid_new).sum()),common_events=int(common.sum()),
        unsupported_claims=['No beta-resolution measurement in flight','No full upstream-beta-dependence injection test','No absolute energy or passive-loss calibration'],
        parent_cache_sha256=r.digest(PARENT/'measured.npz'),new_time_cache_sha256=r.digest(LINEAGE/'measured.npz'),
        script_sha256=r.digest(__file__),RICH_read=False,MC_coefficients_used=False))
    p,w=r.m.projector(data['x']);old=native[common];nn=new[common];du=np.sum(w*(nn-old),axis=1)
    comparison=dict(common_events=int(common.sum()),u_change=r.stats(du),
        abs_u_change_above_point1=int((np.abs(du)>.1).sum()),abs_u_change_above_one=int((np.abs(du)>1).sum()),
        layer_contrasts=[dict(layer=l,native=r.stats(np.einsum('nij,nj->ni',p,old)[:,l]),newcalib=r.stats(np.einsum('nij,nj->ni',p,nn)[:,l])) for l in range(4)])
    r.write(out/'baseline_comparison.json',comparison)
    print('VALID_TIMING_PREPARED',baseline,int(common.sum()),'events, invalid native/new',int((~valid_native).sum()),int((~valid_new).sum()),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',choices=['native','newcalib']);p.add_argument('stage',choices=['prepare','fit','evaluate','bootstrap','all']);p.add_argument('--replicates',type=int,default=32);a=p.parse_args()
    r.OUT=ROOT/'results/flight_calibration_research'/('valid_timing_'+a.baseline)
    if a.stage in ['prepare','all']:prepare(a.baseline)
    if a.stage in ['fit','all']:r.freeze()
    if a.stage in ['evaluate','all']:r.evaluate()
    if a.stage in ['bootstrap','all']:r.bootstrap(a.replicates)
