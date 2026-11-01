#!/usr/bin/env python3
"""Fit geometry-dependent time noise using separately normalized measurement sets."""
import argparse
import json
import time
import hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import feature_conditional as feature
import flight_conditional_likelihood as base
from research_flight_conditional import load,write,digest
from research_conditional_mc import OUT as MCOUT,stamp

OUT=feature.OUT
FLIGHT=base.ROOT/'results/flight_calibration_research/mean_loss_likelihood'
CONFIGS={'mc':[('native_gaussianR','t_native',0.),('native_studentR5','t_native',5.)],
         'flight':[('native_gaussianR','t_native',0.),('newcalib_gaussianR','t',0.),('native_studentR5','t_native',5.)]}


def inputs(dataset):
    directory=MCOUT if dataset=='mc' else FLIGHT
    a=load(directory/'measured.npz');g=load(directory/'moments.npz')
    if not np.array_equal(a['run'],g['run']) or not np.array_equal(a['event'],g['event']):raise ValueError('Measurement/geometry identity mismatch')
    train=a['fold']!=0;test=~train
    if set(a['run'][train])&set(a['run'][test]):raise ValueError('Run leakage')
    return directory,a,g,train,test


def prepare(dataset):
    folder=OUT/dataset;folder.mkdir(exist_ok=True)
    if (folder/'protocol.json').exists():raise FileExistsError('Feature noise protocol exists')
    directory,a,g,train,test=inputs(dataset);ff,norm=feature.prepare_features(a,train)
    frozen={}
    for name,field,nu in CONFIGS[dataset]:
        source=(directory if dataset=='mc' else directory/'conditional')/name/'freeze.json';frozen[name]=dict(path=str(source),sha256=digest(source))
    bounds=base.BOUNDS+[(-.75,.75)]*(6*ff.shape[1])
    write(folder/'protocol.json',dict(created_utc=stamp(),dataset=dataset,feature_normalization=norm,bounds=bounds,
        features=feature.FEATURE_NAMES,parameter_count=len(bounds),penalty=.03,penalty_formula='0.5 * 0.03 * sum(feature_precision_coefficients ** 2), added to mean conditional NLL',
        train_events=int(train.sum()),test_events=int(test.sum()),train_runs=len(np.unique(a['run'][train])),test_runs=len(np.unique(a['run'][test])),
        train_nodes=128,evaluation_nodes=256,starts=2,configs=CONFIGS[dataset],qrange=[.02,.7],
        start_definition='Matching dataset baseline optimum with zero feature coefficients, plus coefficient-free physical initialization with zero feature coefficients',
        selection='Smallest penalized training NLL, no reference beta or evaluation-NLL choice',
        base_freezes=frozen,measured_sha256=digest(directory/'measured.npz'),geometry_sha256=digest(directory/'moments.npz'),
        implementation_sha256={n:digest(base.ROOT/'analysis'/n) for n in ['feature_conditional.py','feature_conditional_likelihood.cpp','research_feature_noise.py']},
        output='Free linear time slope with event precision depending only on measured Tracker position and direction; no event isotope posterior, no velocity clipping',
        additive_time_calibration_opened=False,RICH_read=False,MC_truth_read=False,flight_to_MC_parameters=False,
        beta_evaluation_domain='Already declared conditional-MC run holdout; new candidate fixed before its truth scoring, source pool is a development pool'))
    print('FEATURE_PROTOCOL',dataset,int(train.sum()),int(test.sum()),flush=True)


def solve(obj,start,bounds):
    begin=time.monotonic();result=minimize(obj.evaluate,np.asarray(start),jac=True,method='L-BFGS-B',bounds=bounds,
        options={'maxiter':900,'maxls':40,'ftol':2e-12,'gtol':1e-7,'maxcor':20})
    value,grad=obj.evaluate(result.x);pg=grad.copy()
    for i,(low,high) in enumerate(bounds):
        if result.x[i]<=low+1e-8:pg[i]=min(pg[i],0)
        if result.x[i]>=high-1e-8:pg[i]=max(pg[i],0)
    return dict(theta=result.x,penalized_nll=value,success=bool(result.success),message=str(result.message),iterations=int(result.nit),
        projected_gradient_max=float(np.max(abs(pg))),elapsed_seconds=time.monotonic()-begin,
        unregularized_nll=obj.evaluate(result.x,regularized=False)[0])


def fit(dataset,name):
    directory,a,g,train,test=inputs(dataset);parent=OUT/dataset;p=json.loads((parent/'protocol.json').read_text())
    if digest(directory/'measured.npz')!=p['measured_sha256'] or digest(directory/'moments.npz')!=p['geometry_sha256']:raise ValueError('Feature inputs changed')
    folder=parent/name;folder.mkdir(exist_ok=True)
    if (folder/'freeze.json').exists():raise FileExistsError('Feature model already frozen')
    config=[c for c in p['configs'] if c[0]==name][0];_,field,nu=config;ff,norm=feature.prepare_features(a,train,p['feature_normalization'])
    old=p['base_freezes'][name]
    if digest(old['path'])!=old['sha256']:raise ValueError('Dataset baseline changed')
    baseline=json.loads(Path(old['path']).read_text());theta=np.asarray(baseline['model']['theta'])
    model=feature.FeatureConditional(a,g,train,ff,penalty=p['penalty'],nodes=p['train_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    starts=[np.r_[theta,np.zeros(36)],np.r_[base.initial(1.,.5),np.zeros(36)]]
    results=[]
    for j,start in enumerate(starts):
        path=folder/('start_%d.json'%j)
        if path.exists():result=json.loads(path.read_text())
        else:result=solve(model,start,p['bounds']);write(path,result)
        results.append(result);print('FEATURE_FIT',dataset,name,j,result['penalized_nll'],result['projected_gradient_max'],result['elapsed_seconds'],flush=True)
    best=min(results,key=lambda r:r['penalized_nll']);chosen=np.asarray(best['theta'])
    frozen=dict(complete=True,created_utc=stamp(),model=best,starts=results,protocol_sha256=digest(parent/'protocol.json'),normalization=norm,
        measurement_only=True,reference_read=False,baseline_sha256=old['sha256'])
    write(folder/'freeze.json',frozen)
    output(dataset,name)


def output(dataset,name):
    directory,a,g,train,test=inputs(dataset);parent=OUT/dataset;folder=parent/name;p=json.loads((parent/'protocol.json').read_text())
    frozen=json.loads((folder/'freeze.json').read_text())
    if frozen['protocol_sha256']!=digest(parent/'protocol.json') or digest(directory/'measured.npz')!=p['measured_sha256'] or digest(directory/'moments.npz')!=p['geometry_sha256']:
        raise ValueError('Frozen feature inputs changed')
    _,field,nu=[c for c in p['configs'] if c[0]==name][0];ff,_=feature.prepare_features(a,train,p['feature_normalization'])
    old=p['base_freezes'][name]
    if digest(old['path'])!=old['sha256']:raise ValueError('Baseline changed after feature freeze')
    theta=np.asarray(json.loads(Path(old['path']).read_text())['model']['theta']);chosen=np.asarray(frozen['model']['theta'])
    if (folder/'measurement_evaluation.json').exists():raise FileExistsError('Feature output already complete')
    evaluation=feature.FeatureConditional(a,g,test,ff,penalty=0,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    newloss,_,per=evaluation.evaluate(chosen,events=True)
    oldmodel=base.ConditionalModel(a,g,test,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4);oldloss,_,oldper=oldmodel.evaluate(theta,events=True)
    xx=g['x'][test];tt=a[field][test];w,precision=feature.precision_weights(chosen,xx,ff[test]);u=np.sum(w*tt,axis=1)
    injection=[]
    for shift in [-.02,.02]:injection.append(np.sum(w*(tt+shift*xx),axis=1)-u-shift)
    if np.max(np.abs(injection))>2e-12 or np.max(abs(w.sum(1)))>2e-12 or np.max(abs(np.sum(w*xx,axis=1)-1))>2e-12:raise ValueError('Time equivariance failed')
    np.savez_compressed(folder/'predictions_before_truth.npz',run=a['run'][test],event=a['event'][test],u=u,weights=w,
        new_nll=per,baseline_nll=oldper,R=a['R'][test],features=ff[test])
    runids=np.unique(a['run'][test]);groups=[np.flatnonzero(a['run'][test]==r) for r in runids];rng=np.random.RandomState(1020261012)
    diffs=[np.concatenate([per[groups[i]]-oldper[groups[i]] for i in rng.randint(0,len(groups),len(groups))]).mean() for _ in range(512)]
    write(folder/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),mean_nll=newloss,baseline_nll=oldloss,delta_nll=newloss-oldloss,
        paired_run_delta_nll_interval=np.percentile(diffs,[2.5,50,97.5]),
        injection_max_error=float(np.max(np.abs(injection))),weight_median=np.median(w,axis=0),beta_above_one_retained=int(np.sum(u<1)),
        freeze_sha256=digest(folder/'freeze.json'),predictions_sha256=digest(folder/'predictions_before_truth.npz'),
        output_script_sha256=digest(Path(__file__)),coefficients_refitted_during_output=False,
        scope='Measurement likelihood and frozen output only; beta resolution is evaluated separately',MC_truth_read=False,RICH_read=False))
    print('FEATURE_FREEZE',dataset,name,'delta NLL',newloss-oldloss,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','fit','output']);parser.add_argument('--dataset',choices=['mc','flight'],required=True);parser.add_argument('--name',default='native_gaussianR');args=parser.parse_args()
    if args.stage=='prepare':prepare(args.dataset)
    elif args.stage=='fit':fit(args.dataset,args.name)
    else:output(args.dataset,args.name)
