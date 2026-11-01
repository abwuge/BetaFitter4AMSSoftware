#!/usr/bin/env python3
"""Measured-only run-separated timing redundancy experiment.

The reported timing-contrast widths are not beta-resolution measurements.
No MC parameters, isotope labels, RICH measurements or beta targets enter.
"""
import argparse
import hashlib
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import flight_timing_redundancy as m

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research/timing_redundancy'
FIELDS=['run','event','tof_tl','tof_leng','tof_barid','tof_adca','tof_rtime',
        'tof_rtime_ht','tof_rtime_sht','tof_pos','tk_pos','tk_dir','tof_edep']


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def native(v):
    if isinstance(v,dict):return {k:native(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [native(x) for x in v]
    if isinstance(v,np.ndarray):return native(v.tolist())
    if isinstance(v,(np.integer,np.bool_)):return v.item()
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    return v

def write(path,value):path.write_text(json.dumps(native(value),indent=2,allow_nan=False)+'\n')

def stats(x):
    x=np.asarray(x,float);x=x[np.isfinite(x)]
    if not len(x):return dict(n=0)
    q=np.percentile(x,[16,50,84]);return dict(n=len(x),mean=x.mean(),median=q[1],w68=(q[2]-q[0])/2,rmse=np.sqrt(np.mean(x*x)))

def subset(a,mask):return {k:v[mask] for k,v in a.items()}

def read_source(record):
    source=int(record['source']);part=ROOT/'results/time_local_beta/parts'/('%05d'%source)
    with np.load(str(part/'measured.npz'),allow_pickle=False) as f:
        wanted=set(zip(f['run'].astype(int),f['event'].astype(int)))
    with uproot.open(record['path']) as f:a=f['amstreea'].arrays(FIELDS,library='np')
    ix=np.array([i for i,k in enumerate(zip(a['run'].astype(int),a['event'].astype(int))) if k in wanted],int)
    a=subset(a,ix);n=len(ix)
    pos=a['tk_pos'].astype(float);direction=a['tk_dir'].astype(float)
    valid=np.isfinite(pos).all(axis=2)&np.isfinite(direction).all(axis=2)&(np.abs(direction[:,:,2])>.05)
    valid&=(np.max(np.abs(pos),axis=2)<2000)&(np.linalg.norm(direction,axis=2)>.1)
    xy=np.full((n,4,2),np.nan);good=np.zeros((n,4),bool)
    for layer in range(4):
        z=a['tof_pos'][:,layer,2].astype(float)
        distance=np.where(valid,np.abs(pos[:,:,2]-z[:,None]),np.inf)
        nearest=np.argmin(distance,axis=1);p=pos[np.arange(n),nearest];d=direction[np.arange(n),nearest]
        with np.errstate(invalid='ignore',divide='ignore'):
            xy[:,layer]=p[:,:2]+(z-p[:,2])[:,None]*d[:,:2]/d[:,2,None]
        good[:,layer]=np.isfinite(xy[:,layer]).all(axis=1)&(np.abs(xy[:,layer])<150).all(axis=1)&(np.min(distance,axis=1)<80)&(np.abs(z)<100)
    adc=a['tof_adca'].astype(float)
    good&=np.isfinite(adc).all(axis=2)&(adc>0).all(axis=2)&(adc<3700).all(axis=2)
    bar=a['tof_barid'].astype(int);good&=(bar>=0)&(bar<10)
    t=a['tof_tl'].astype(float);x=a['tof_leng'].astype(float)/29.9792458
    event_good=np.isfinite(t).all(axis=1)&(np.abs(t)<1000).all(axis=1)&np.isfinite(x).all(axis=1)&(np.ptp(x,axis=1)>1)
    event_good&=((bar>=0)&(bar<10)).all(axis=1)
    out=dict(run=a['run'],event=a['event'],source=np.full(n,source,int),epoch=(a['run'].astype(np.int64)//(15*86400)),
        raw=np.stack([a[k] for k in ['tof_rtime','tof_rtime_ht','tof_rtime_sht']],axis=-1).astype(float),
        x=x,t=t,xy=xy,adc=adc,bar=bar,good=good,edep=a['tof_edep'].astype(float))
    return subset(out,event_good)


def prepare(workers):
    OUT.mkdir(exist_ok=True)
    config=json.loads((ROOT/'results/time_local_beta/configuration.json').read_text());indices=np.asarray(config['file_indices'],int)
    inv=json.loads((ROOT/'results/flight_calibration_research/input_inventory.json').read_text())['records']
    anchors=[int(np.argmin(abs(indices-r['source']))) for r in inv]
    chosen=sorted(set(int(indices[j]) for c in anchors for j in range(max(0,c-12),min(len(indices),c+13))))
    records=json.loads((ROOT/'results/tofcalib_validation/preparation.json').read_text())['records']
    records=sorted([r for r in records if int(r['source']) in chosen],key=lambda r:int(r['source']))
    write(OUT/'protocol.json',dict(fields=FIELDS,sources=chosen,source_count=len(chosen),
        source_selection='Each of nine inventory anchors plus/minus 12 positions in the pre-existing 4096-source grid; endpoint windows clipped',
        event_selection='Reuse detector-only carbon training identities, including inherited charge and positive TOF slope cuts; require finite four-layer times and geometry. No RICH cut or new residual cut.',
        response_split='SHA256(decimal run) modulo 5 in {1,2}',correction_split='modulo 5 in {3,4}',evaluation_split='modulo 5 equals 0',
        epochs='floor(Unix run/(15 days)); diagnostic grid, not established hardware-calibration release dates',
        features=['LT end difference','mean HT minus LT over two ends','mean SHT minus LT over two ends'],
        response_model='Per bar: cubic independent Tracker coordinate, inverse square-root ADC and square terms, position interaction, epoch intercepts',
        correction_model='Shared per-bar contrast coefficients after profiling free event t0 and u; no new time intercepts; equal layer metric with robust event weighting',
        candidates=['end_difference','all_thresholds'],baseline='Original stored NEWTOF layer times; unconstrained equal-weight timing slope',
        RICH_read=False,MC_read=False,physical_loss_model=False,ADC_upper_bound=3700,ADC_bound_source='Native TofRecPar::PUXMXCH',
        bootstrap_seed=204817,script_sha256=digest(__file__),module_sha256=digest(ROOT/'analysis/flight_timing_redundancy.py')))
    pieces=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for j,piece in enumerate(pool.map(read_source,records),1):
            pieces.append(piece)
            if j%20==0 or j==len(records):print('READ',j,'/',len(records),'events',sum(len(p['run']) for p in pieces),flush=True)
    data={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    identity=list(zip(data['run'],data['event']))
    if len(set(identity))!=len(identity):raise ValueError('Duplicate event identity')
    runs=np.unique(data['run']);lookup={int(r):int(hashlib.sha256(str(int(r)).encode()).hexdigest(),16)%5 for r in runs}
    data['fold']=np.array([lookup[int(r)] for r in data['run']],np.int8)
    np.savez_compressed(str(OUT/'measured.npz'),**data)
    write(OUT/'input_summary.json',dict(events=len(identity),runs=len(runs),epochs=len(np.unique(data['epoch'])),
        cohorts={str(i):dict(events=int((data['fold']==i).sum()),runs=len(np.unique(data['run'][data['fold']==i]))) for i in range(5)},
        usable_ADC_hit_fraction=data['good'].mean(),cache_sha256=digest(OUT/'measured.npz')))
    print('PREPARATION_COMPLETE',len(identity),'events',len(runs),'runs',flush=True)


def load():
    with np.load(str(OUT/'measured.npz'),allow_pickle=False) as f:return {k:f[k] for k in f.files}


def fit(data,weights=None):
    if weights is not None:
        # Explicit replication keeps robust medians and PCA identical to a
        # genuine whole-run bootstrap, not only its least-squares weights.
        counts=np.asarray(weights,int)
        if np.any(counts<0) or not np.allclose(counts,weights):raise ValueError('Integer bootstrap multiplicities required')
        data=subset(data,np.repeat(np.arange(len(counts)),counts));weights=None
    response=np.isin(data['fold'],[1,2]);training=np.isin(data['fold'],[3,4])
    models=m.fit_response(data,response,weights)
    h,supported,raw_res,pooled=m.apply_response(data,models)
    channel=data['bar']+10*np.arange(4)[None,:]
    fitted={}
    for name,nf in [('end_difference',1),('all_thresholds',3)]:
        active=training&np.any(supported[:,:,:nf],axis=(1,2))
        fit=m.fit_correction(data['x'][active],data['t'][active],h[active,:,:nf],channel[active],
                             None if weights is None else weights[active])
        fitted[name]=dict(theta=fit['theta'],active=fit['active'],eigenvalues=fit['eigenvalues'],scale=fit['scale'],iterations=fit['iterations'],training_events=int(active.sum()))
    return models,fitted


def freeze():
    data=load();models,fitted=fit(data)
    result=dict(models=models,candidates=fitted,training_runs=np.unique(data['run'][data['fold']!=0]),
        evaluation_runs=np.unique(data['run'][data['fold']==0]),response_runs=np.unique(data['run'][np.isin(data['fold'],[1,2])]),
        correction_runs=np.unique(data['run'][np.isin(data['fold'],[3,4])]),
        module_sha256=digest(ROOT/'analysis/flight_timing_redundancy.py'),script_sha256=digest(__file__),input_sha256=digest(OUT/'measured.npz'))
    write(OUT/'freeze.json',result)
    print('FROZEN',len(models),'response models', {k:int(v['active'].sum()) for k,v in fitted.items()},flush=True)


def evaluate():
    frozen=json.loads((OUT/'freeze.json').read_text());data=load();test=data['fold']==0
    if set(frozen['training_runs'])&set(frozen['evaluation_runs']):raise ValueError('Run leakage')
    if set(frozen['response_runs'])&set(frozen['correction_runs']):raise ValueError('Response/correction run leakage')
    data=subset(data,test);h,support,raw_res,pooled=m.apply_response(data,frozen['models']);channel=data['bar']+10*np.arange(4)[None,:]
    p,w=m.projector(data['x']);base=np.einsum('nij,nj->ni',p,data['t']);u0=np.sum(w*data['t'],axis=1)
    response=[]
    for layer in range(4):
        for feature in range(3):
            good=support[:,layer,feature]
            response.append(dict(layer=layer,feature=feature,pooled=stats(pooled[good,layer,feature]),epoch_conditioned=stats(raw_res[good,layer,feature])))
    result=dict(complete=True,events=len(u0),runs=len(np.unique(data['run'])),response=response,candidates={},freeze_sha256=digest(OUT/'freeze.json'),
                beta_resolution_measured=False,absolute_energy_calibration=False,physical_loss_model=False)
    predictions=dict(run=data['run'],event=data['event'],epoch=data['epoch'],baseline_u=u0,baseline_residual=base,support=support,bar=data['bar'])
    rng=np.random.RandomState(338);a=rng.normal(size=len(u0));b=rng.uniform(-.15,.15,len(u0));delta=a[:,None]+b[:,None]*data['x']
    # Missing thresholds remain missing under a physical timing injection.
    # Shifting a stored zero would otherwise fabricate a newly detected pulse.
    valid_raw=np.isfinite(data['raw'])&(np.abs(data['raw'])<1000)&(data['raw']!=0)
    shifted=dict(data);shifted['raw']=np.where(valid_raw,data['raw']-delta[:,:,None,None],data['raw'])
    hh,_,_,_=m.apply_response(shifted,frozen['models'])
    for name,nf in [('end_difference',1),('all_thresholds',3)]:
        theta=np.asarray(frozen['candidates'][name]['theta'])
        rec=m.reconstruct(data['x'],data['t'],h[:,:,:nf],channel,theta)
        check=m.reconstruct(data['x'],data['t']+delta,hh[:,:,:nf],channel,theta)
        response_error=np.max(np.abs(check['u']-rec['u']-b))
        if response_error>1.e-10:raise ValueError('Coherent timing response failed')
        row=dict(timing_contrasts=[dict(layer=l,baseline=stats(base[:,l]),corrected=stats(rec['residual'][:,l])) for l in range(4)],
            inverse_speed_change=stats(rec['u']-u0),maximum_signal_response_error=response_error,
            baseline_beta_above_one=int(((u0>0)&(u0<1)).sum()),corrected_beta_above_one=int(((rec['u']>0)&(rec['u']<1)).sum()),
            affected_events=int(np.any(support[:,:,:nf],axis=(1,2)).sum()),epochs=[])
        for epoch in np.unique(data['epoch']):
            ix=data['epoch']==epoch
            row['epochs'].append(dict(epoch=int(epoch),n=int(ix.sum()),baseline=stats(base[ix].ravel()),corrected=stats(rec['residual'][ix].ravel()),u_shift=stats(rec['u'][ix]-u0[ix])))
        result['candidates'][name]=row
        predictions[name+'_u']=rec['u'];predictions[name+'_residual']=rec['residual']
        print(name,'layer contrast W68:',[(r['baseline']['w68'],r['corrected']['w68']) for r in row['timing_contrasts']],'u shift',row['inverse_speed_change'],flush=True)
    write(OUT/'evaluation.json',result);np.savez_compressed(str(OUT/'predictions.npz'),**predictions)


def run_weights(runs,rng):
    unique,inverse=np.unique(runs,return_inverse=True);draw=rng.randint(0,len(unique),len(unique));count=np.bincount(draw,minlength=len(unique))
    return count[inverse].astype(float)


def bootstrap(replicates):
    data=load();test=data['fold']==0;evaluation=subset(data,test);channel=evaluation['bar']+10*np.arange(4)[None,:]
    p,w=m.projector(evaluation['x']);base=np.einsum('nij,nj->ni',p,evaluation['t']);u0=np.sum(w*evaluation['t'],axis=1)
    rng=np.random.RandomState(204817);results=[];all_u={name:[] for name in ['end_difference','all_thresholds']}
    for rep in range(replicates):
        weight=np.ones(len(test))
        for cohort in [[1,2],[3,4]]:
            ix=np.isin(data['fold'],cohort);weight[ix]=run_weights(data['run'][ix],rng)
        try:
            models,fitted=fit(data,weight)
            h,support,_,_=m.apply_response(evaluation,models)
            eval_weight=run_weights(evaluation['run'],rng).astype(int);ix=np.repeat(np.arange(len(u0)),eval_weight)
            row=dict(replicate=rep,complete=True,candidates={})
            for name,nf in [('end_difference',1),('all_thresholds',3)]:
                rec=m.reconstruct(evaluation['x'],evaluation['t'],h[:,:,:nf],channel,fitted[name]['theta'])
                changes=[]
                for l in range(4):
                    old=stats(base[ix,l]);new=stats(rec['residual'][ix,l]);changes.append(1-new['w68']/old['w68'])
                row['candidates'][name]=dict(layer_contrast_fractional_gain=changes,mean_u_shift=float(np.mean(rec['u'][ix]-u0[ix])))
                all_u[name].append(rec['u'])
        except Exception as exc:
            row=dict(replicate=rep,complete=False,error=repr(exc))
        results.append(row)
        if (rep+1)%4==0:print('BOOTSTRAP',rep+1,'/',replicates,'success',sum(r['complete'] for r in results),flush=True)
        write(OUT/'bootstrap_progress.json',dict(replicates=results))
    uncertainty={}
    for name,values in all_u.items():
        if len(values)>1:
            sd=np.std(np.array(values),axis=0,ddof=1)
            good=[r['candidates'][name]['layer_contrast_fractional_gain'] for r in results if r['complete']]
            uncertainty[name]=dict(successful_refits=len(values),calibration_u_sd=stats(sd),layer_contrast_gain_percentiles=np.percentile(good,[2.5,50,97.5],axis=0))
    write(OUT/'bootstrap.json',dict(complete=True,requested_replicates=replicates,successful_replicates=sum(r['complete'] for r in results),
        retrained='Both response and correction using whole-run resampling; independently resample evaluation runs for reported contrast intervals',
        caveat='Only defined response family and sampled runs; not a detector-model or energy-loss systematic uncertainty',uncertainty=uncertainty,replicates=results))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','fit','evaluate','bootstrap','all']);p.add_argument('--workers',type=int,default=4);p.add_argument('--replicates',type=int,default=32);a=p.parse_args()
    if a.stage in ['prepare','all']:prepare(a.workers)
    if a.stage in ['fit','all']:freeze()
    if a.stage in ['evaluate','all']:evaluate()
    if a.stage in ['bootstrap','all']:bootstrap(a.replicates)
