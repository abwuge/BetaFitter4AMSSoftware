#!/usr/bin/env python3
"""Independent AMS MC closure of the readout-timing module, not energy transport.

Truth is opened only by the evaluation stage after a serialized measured-only
fit exists. Coefficients remain in this MC directory and never enter flight.
"""
import argparse
import hashlib
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import flight_timing_redundancy as m
import research_flight_timing_redundancy as r
ROOT=r.ROOT;OUT=ROOT/'results/flight_calibration_research/timing_mc_closure'
MCROOT=Path('/lustre02/user/yachen/public/hxwu/newtofcalib_mc_fullsize_mitdst/amsd69nMCC12')
FIELDS=r.FIELDS+['tk_qin','tk_q']


def measured(record):
    with uproot.open(record['path']) as f:a=f['amstreea'].arrays(FIELDS,entry_stop=3000,library='np')
    n=len(a['run']);q=np.where(a['tk_qin'][:,0,2]<2.5,a['tk_q'][:,1],a['tk_qin'][:,0,2])
    pos=a['tk_pos'].astype(float);direction=a['tk_dir'].astype(float)
    valid=np.isfinite(pos).all(axis=2)&np.isfinite(direction).all(axis=2)&(np.abs(direction[:,:,2])>.05)
    valid&=(np.max(np.abs(pos),axis=2)<2000)&(np.linalg.norm(direction,axis=2)>.1)
    xy=np.full((n,4,2),np.nan);good=np.zeros((n,4),bool)
    for l in range(4):
        z=a['tof_pos'][:,l,2].astype(float);distance=np.where(valid,np.abs(pos[:,:,2]-z[:,None]),np.inf)
        nearest=np.argmin(distance,axis=1);pp=pos[np.arange(n),nearest];d=direction[np.arange(n),nearest]
        with np.errstate(invalid='ignore',divide='ignore'):
            xy[:,l]=pp[:,:2]+(z-pp[:,2])[:,None]*d[:,:2]/d[:,2,None]
        good[:,l]=np.isfinite(xy[:,l]).all(axis=1)&(np.abs(xy[:,l])<150).all(axis=1)&(np.min(distance,axis=1)<80)&(np.abs(z)<100)
    adc=a['tof_adca'].astype(float);good&=np.isfinite(adc).all(axis=2)&(adc>0).all(axis=2)&(adc<3700).all(axis=2)
    bar=a['tof_barid'].astype(int);good&=(bar>=0)&(bar<10)
    t=a['tof_tl'].astype(float);x=a['tof_leng'].astype(float)/29.9792458
    event_good=np.isfinite(t).all(axis=1)&(np.abs(t)<1000).all(axis=1)&np.isfinite(x).all(axis=1)&(np.ptp(x,axis=1)>1)
    event_good&=((bar>=0)&(bar<10)).all(axis=1)&np.isfinite(q)&(q>=5.5)&(q<6.5)
    event_good&=np.isfinite(a['tof_edep']).all(axis=1)&(a['tof_edep']>0).all(axis=1)
    safe=x-x.mean(axis=1)[:,None]
    with np.errstate(invalid='ignore',divide='ignore'):u=(safe*t).sum(axis=1)/(safe*safe).sum(axis=1)
    event_good&=np.isfinite(u)&(u>0)
    data=dict(run=a['run'],event=a['event'],entry=np.arange(n),source=np.full(n,record['source'],int),epoch=np.zeros(n,int),
        raw=np.stack([a[k] for k in ['tof_rtime','tof_rtime_ht','tof_rtime_sht']],axis=-1).astype(float),
        x=x,t=t,xy=xy,adc=adc,bar=bar,good=good,edep=a['tof_edep'].astype(float))
    return r.subset(data,event_good)


def prepare(workers):
    OUT.mkdir(exist_ok=True);paths=sorted(MCROOT.glob('*/*.root'))
    chosen=np.unique(np.linspace(0,len(paths)-1,min(48,len(paths))).astype(int))
    records=[dict(source=int(i),path=str(paths[j])) for i,j in enumerate(chosen)]
    r.write(OUT/'protocol.json',dict(records=records,measured_fields=FIELDS,entry_stop=3000,
        sample='C12 production files uniformly spaced by sorted path; observed carbon charge selection; no truth-based event cut before freeze',
        train_split='Whole MC run SHA256 modulo5: response1,2; correction3,4; evaluation0',
        MC_conditions='One stationary response epoch; no flight response, noise or correction coefficients reused',
        truth_target='GLS timing slope of stored tof_tl_truth using the same measured path; not L5 beta or RICH beta',
        reference_caveat='TofCalib producer selects earliest matching-bar TOF MC cluster and adds a common offset; no primary-track label in this stored truth field',
        RICH_read=False,truth_read_at_training=False))
    pieces=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i,a in enumerate(pool.map(measured,records),1):
            pieces.append(a)
            if i%12==0:print('MC_MEASURED_READ',i,'/',len(records),flush=True)
    data={k:np.concatenate([a[k] for a in pieces]) for k in pieces[0]}
    if len(set(zip(data['run'],data['event'])))!=len(data['run']):raise ValueError('Duplicate MC identity')
    lookup={int(run):int(hashlib.sha256(str(int(run)).encode()).hexdigest(),16)%5 for run in np.unique(data['run'])}
    data['fold']=np.array([lookup[int(run)] for run in data['run']],np.int8)
    np.savez_compressed(str(OUT/'measured.npz'),**data)
    print('MC_PREPARED',len(data['run']),'events',len(lookup),'runs',flush=True)


def load():
    with np.load(str(OUT/'measured.npz'),allow_pickle=False) as f:return {k:f[k] for k in f.files}


def fit():
    data=load();models,fitted=r.fit(data)
    r.write(OUT/'freeze.json',dict(models=models,candidates=fitted,measured_sha256=r.digest(OUT/'measured.npz'),
        train_runs=np.unique(data['run'][data['fold']!=0]),test_runs=np.unique(data['run'][data['fold']==0]),
        module_sha256=r.digest(ROOT/'analysis/flight_timing_redundancy.py'),data_coefficients_imported=False,truth_opened=False))
    print('MC_FROZEN',len(models),'response models',flush=True)


def evaluate():
    frozen=json.loads((OUT/'freeze.json').read_text());data=load();data=r.subset(data,data['fold']==0)
    if set(frozen['train_runs'])&set(frozen['test_runs']):raise ValueError('MC run leakage')
    records=json.loads((OUT/'protocol.json').read_text())['records'];truth=np.full_like(data['t'],np.nan)
    for record in records:
        take=data['source']==record['source'];ix=np.flatnonzero(take)
        if not len(ix):continue
        with uproot.open(record['path']) as f:a=f['amstreea'].arrays(['run','event','tof_tl_truth'],entry_stop=3000,library='np')
        entries=data['entry'][ix]
        if not np.array_equal(a['event'][entries],data['event'][ix]) or not np.array_equal(a['run'][entries],data['run'][ix]):raise ValueError('MC truth identity mismatch')
        truth[ix]=a['tof_tl_truth'][entries]
    good=np.isfinite(truth).all(axis=1)&(np.abs(truth)<1000).all(axis=1)&(truth!=0).all(axis=1)
    data=r.subset(data,good);truth=truth[good];p,w=m.projector(data['x']);u_truth=np.sum(w*truth,axis=1)
    h,supported,_,_=m.apply_response(data,frozen['models']);channel=data['bar']+10*np.arange(4)[None,:]
    u0=np.sum(w*data['t'],axis=1);base=r.stats(u0-u_truth);out=dict(complete=True,events=len(u0),runs=len(np.unique(data['run'])),
        missing_truth=int((~good).sum()),truth_target='Effective TOF timing slope at measured path, not L5 beta',baseline=base,candidates={},
        freeze_sha256=r.digest(OUT/'freeze.json'),MC_parameters_transferred_to_flight=False)
    predictions=dict(run=data['run'],event=data['event'],baseline_u=u0,truth_effective_u=u_truth)
    for name,nf in [('end_difference',1),('all_thresholds',3)]:
        rec=m.reconstruct(data['x'],data['t'],h[:,:,:nf],channel,np.asarray(frozen['candidates'][name]['theta']))
        predictions[name+'_u']=rec['u']
        metric=r.stats(rec['u']-u_truth);row=dict(metric=metric,w68_fractional_gain=1-metric['w68']/base['w68'],
            rmse_fractional_gain=1-metric['rmse']/base['rmse'],bias_change=metric['mean']-base['mean'],bins=[])
        for lo,hi in [(0,1.02),(1.02,1.1),(1.1,1.3),(1.3,10)]:
            mask=(u_truth>=lo)&(u_truth<hi)
            row['bins'].append(dict(u_truth_range=[lo,hi],baseline=r.stats((u0-u_truth)[mask]),corrected=r.stats((rec['u']-u_truth)[mask])))
        out['candidates'][name]=row
    r.write(OUT/'evaluation.json',out);np.savez_compressed(str(OUT/'predictions.npz'),**predictions)
    print(json.dumps(r.native(out),indent=2),flush=True)


def bootstrap():
    with np.load(str(OUT/'predictions.npz'),allow_pickle=False) as f:a={k:f[k] for k in f.files}
    rng=np.random.RandomState(72919);rows={name:[] for name in ['end_difference','all_thresholds']}
    for rep in range(512):
        weight=r.run_weights(a['run'],rng).astype(int);ix=np.repeat(np.arange(len(weight)),weight)
        old=r.stats((a['baseline_u']-a['truth_effective_u'])[ix])
        for name in rows:
            new=r.stats((a[name+'_u']-a['truth_effective_u'])[ix])
            rows[name].append([1-new['w68']/old['w68'],1-new['rmse']/old['rmse'],new['mean']-old['mean']])
    result=dict(replicates=512,method='Paired resampling of held-out MC runs, frozen estimator; excludes retraining and simulation-model uncertainty',
        candidates={name:dict(columns=['fractional_W68_gain','fractional_RMSE_gain','bias_change'],percentiles_2p5_50_97p5=np.percentile(values,[2.5,50,97.5],axis=0)) for name,values in rows.items()})
    r.write(OUT/'evaluation_bootstrap.json',result);print(json.dumps(r.native(result),indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','fit','evaluate','bootstrap','all']);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    if a.stage in ['prepare','all']:prepare(a.workers)
    if a.stage in ['fit','all']:fit()
    if a.stage in ['evaluate','all']:evaluate()
    if a.stage in ['bootstrap','all']:bootstrap()
