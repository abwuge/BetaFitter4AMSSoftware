#!/usr/bin/env python3
"""Whole-training-run recalibration for geometry/readout noise, no truth imports."""
import argparse, json, time, socket
from pathlib import Path
import numpy as np
import feature_conditional as feature
import flight_conditional_likelihood as base
from research_amplitude_noise_20261013 import OUT,ROOT,inputs,feature_arrays
from research_feature_noise import solve
from research_flight_conditional import write,digest
from research_conditional_mc import stamp

BOOT=OUT/'bootstrap'
CONFIGS=[(name,nu,variant) for name,nu in [('native_gaussianR',0.),('native_studentR5',5.)] for variant in ['geometry','readout']]


def rebase(theta,old,new):
    """Affine standardization change is exact wherever neither feature is clipped."""
    out=np.asarray(theta,float).copy();n=len(new['center']);B=out[13:].reshape(6,n).copy()
    om=np.asarray(old['center']);os=np.asarray(old['scale']);nm=np.asarray(new['center']);ns=np.asarray(new['scale'])
    out[2:8]+=B.dot((nm-om)/os);out[13:]=(B*(ns/os)[None,:]).ravel()
    return out


def prepare(replicates):
    BOOT.mkdir(exist_ok=True)
    if (BOOT/'protocol.json').exists():raise FileExistsError('Bootstrap protocol already frozen')
    folder,a,g,adc,domain,train,test=inputs('mc');p=json.loads((OUT/'mc/protocol.json').read_text());freezes={}
    for name,nu,variant in CONFIGS:
        label=name+'_'+variant;d=OUT/'mc'/label;meta=json.loads((d/'measurement_evaluation.json').read_text())
        if not meta['complete'] or digest(d/'freeze.json')!=meta['freeze_sha256'] or digest(d/'predictions_before_truth.npz')!=meta['output_sha256']:raise ValueError('Changed base fit')
        freezes[label]=dict(model=meta['freeze_sha256'],output=meta['output_sha256'])
    write(BOOT/'protocol.json',dict(created_utc=stamp(),replicates=replicates,seed=4420261013 % (2**32-1),configs=CONFIGS,
        sampling='Resample whole measured-training run IDs with replacement; all events and duplicates retained; evaluation events fixed',
        normalization='Recalculate feature means/scales on every resampled training set including duplicate runs; clip +/-3 after standardization',
        starts='Rebased original measured-only optimum and independent physical kappa1/mixture0.5 zero-feature start, each chosen only by penalized resampled training likelihood',
        feature_rebase='Exactly transforms the unclipped linear predictor; clipped tails may differ and both starts are fully optimized',
        train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),base_freezes=freezes,
        input_protocol_sha256=digest(OUT/'mc/protocol.json'),script_sha256=digest(Path(__file__)),
        measurement_sha256=digest(folder/'measured.npz'),geometry_sha256=digest(folder/'moments.npz'),supplement_sha256=digest(OUT/'mc_supplement.npz'),
        nll_convergence_tolerance=1e-3,MC_truth_read=False,RICH_read=False,evaluation_reference_used=False))
    print('AMPLITUDE_BOOTSTRAP_PROTOCOL',replicates,flush=True)


def replicate(rep,a,g,adc,domain,train,test,protocol,config_protocol):
    note=BOOT/('%04d.json'%rep);output=note.with_suffix('.npz')
    if note.exists():
        old=json.loads(note.read_text())
        if old['complete'] and digest(output)==old['output_sha256']:return old
        raise ValueError('Incomplete bootstrap output exists')
    start=time.monotonic();rng=np.random.RandomState(protocol['seed']+rep);runs=np.asarray(protocol['train_runs']);draw=rng.choice(runs,len(runs),replace=True)
    index=np.concatenate([np.flatnonzero(train&(a['run']==r)) for r in draw]);result=[];pred={k:a[k][test] for k in ['run','event']}
    x=g['x'][test];t=a['t_native'][test]
    for name,nu,variant in CONFIGS:
        label=name+'_'+variant;spec=config_protocol['variants'][variant];directory=OUT/'mc'/label
        original=json.loads((directory/'freeze.json').read_text())['model'];theta=np.array(original['theta'])
        ff,norm=feature_arrays(a,adc,domain,index,variant)
        initial=rebase(theta,spec['normalization'],norm)
        # The rebased start can lie just outside a declared coefficient bound.
        # Project only the starting vector; final output is never beta-clipped.
        bounds=np.asarray(spec['bounds'],float);initial=np.maximum(bounds[:,0],np.minimum(bounds[:,1],initial))
        obj=feature.FeatureConditional(a,g,index,ff,penalty=config_protocol['penalty'],nodes=config_protocol['nodes'],rigidity_nu=nu,threads=4)
        starts=[initial,np.r_[base.initial(1,.5),np.zeros(6*spec['nfeatures'])]]
        fits=[solve(obj,point,spec['bounds']) for point in starts];best=min(fits,key=lambda v:v['penalized_nll'])
        if not best['success'] or best['projected_gradient_max']>protocol['nll_convergence_tolerance']:
            fits.append(solve(obj,np.asarray(best['theta']),spec['bounds']));best=min(fits,key=lambda v:v['penalized_nll'])
        w,Q=feature.precision_weights(np.asarray(best['theta']),x,ff[test]);u=np.sum(w*t,axis=1)
        if not np.isfinite(u).all():raise ValueError('Nonfinite bootstrap slope')
        pred[label+'_u']=u
        result.append(dict(label=label,normalization=norm,model=best,starts=fits,
            coherent_constraint_max_error=float(max(np.max(abs(w.sum(1))),np.max(abs(np.sum(w*x,axis=1)-1))))))
        print('AMPLITUDE_BOOT_MODEL',rep,label,best['penalized_nll'],best['projected_gradient_max'],flush=True)
    np.savez_compressed(output,**pred)
    record=dict(complete=True,created_utc=stamp(),replicate=rep,drawn_runs=draw,training_events=len(index),models=result,
        elapsed_seconds=time.monotonic()-start,host=socket.gethostname(),output_sha256=digest(output),protocol_sha256=digest(BOOT/'protocol.json'),MC_truth_read=False,RICH_read=False)
    write(note,record);print('AMPLITUDE_BOOT_COMPLETE',rep,record['elapsed_seconds'],flush=True)
    return record


def worker(index,workers):
    p=json.loads((BOOT/'protocol.json').read_text());cp=json.loads((OUT/'mc/protocol.json').read_text());folder,a,g,adc,domain,train,test=inputs('mc')
    if digest(Path(__file__))!=p['script_sha256'] or digest(OUT/'mc/protocol.json')!=p['input_protocol_sha256']:raise ValueError('Bootstrap protocol/source changed')
    for name,expected in [('measured.npz','measurement_sha256'),('moments.npz','geometry_sha256')]:
        if digest(folder/name)!=p[expected]:raise ValueError('Bootstrap measurements changed')
    if digest(OUT/'mc_supplement.npz')!=p['supplement_sha256']:raise ValueError('ADC supplement changed')
    for label,h in p['base_freezes'].items():
        if digest(OUT/'mc'/label/'freeze.json')!=h['model'] or digest(OUT/'mc'/label/'predictions_before_truth.npz')!=h['output']:raise ValueError('Bootstrap original model changed')
    for rep in range(index,p['replicates'],workers):replicate(rep,a,g,adc,domain,train,test,p,cp)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','worker']);parser.add_argument('--replicates',type=int,default=128);parser.add_argument('--worker',type=int,default=0);parser.add_argument('--workers',type=int,default=32);args=parser.parse_args()
    if args.stage=='prepare':prepare(args.replicates)
    else:worker(args.worker,args.workers)
