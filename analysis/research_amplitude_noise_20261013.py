#!/usr/bin/env python3
"""Measured-only amplitude-conditioned noise; independent MC/flight parameters.

Scoring is a separate explicitly invoked stage. No legacy files are modified.
"""
import argparse, json, time, hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import flight_conditional_likelihood as base
import feature_conditional as feature
from research_feature_noise import solve
from research_flight_conditional import load,write,digest
from research_conditional_mc import OUT as MCOUT, stamp

ROOT=base.ROOT
OUT=ROOT/'results/flight_calibration_research/continuation_20261013/amplitude_noise'
FLIGHT=ROOT/'results/flight_calibration_research/mean_loss_likelihood'
CONFIGS={'mc':[('native_gaussianR','t_native',0.),('native_studentR5','t_native',5.)],
         'flight':[('native_gaussianR','t_native',0.),('newcalib_gaussianR','t',0.),('native_studentR5','t_native',5.)]}


def supplement():
    OUT.mkdir(parents=True,exist_ok=True)
    target=OUT/'mc_supplement.npz'
    if target.exists():raise FileExistsError(target)
    a=load(MCOUT/'measured.npz');p=json.loads((MCOUT/'protocol.json').read_text());n=len(a['run'])
    result={k:a[k] for k in ['run','event','source','entry']};result['adc']=np.full((n,4,2),np.nan);result['available']=np.zeros(n,bool)
    records=[]
    def one(r):
        ix=np.flatnonzero(a['source']==r['source'])
        if not len(ix):return None
        with uproot.open(r['path']) as s:
            tree=s['amstreea'];keys=set(tree.keys());field='tof_Aadc' if 'tof_Aadc' in keys else ('tof_adca' if 'tof_adca' in keys else None)
            if field is None:return (ix,None,dict(source=r['source'],events=len(ix),field=None,reason='No stored end ADC branch'))
            fields=['run','event','itrtrack','ibetah','tof_barid',field]
            raw=tree.arrays(fields,entry_stop=r['entry_stop'],library='np');j=a['entry'][ix]
            for k in ['run','event','itrtrack','ibetah']:
                if not np.array_equal(raw[k][j],a[k][ix]):raise ValueError('Supplement identity mismatch '+k)
            if not np.array_equal(raw['tof_barid'][j],a['bar'][ix]):raise ValueError('Bar association mismatch')
            adc=raw[field][j].astype(float)
        return (ix,adc,dict(source=r['source'],path=r['path'],events=len(ix),field=field,fields_read=fields,identity_match=True))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i,item in enumerate(pool.map(one,p['records'])):
            if item is None:continue
            ix,adc,note=item;records.append(note)
            if adc is not None:result['adc'][ix]=adc;result['available'][ix]=True
            if i%40==0:print('ADC_SUPPLEMENT',i,'/',len(p['records']),flush=True)
    positive=np.isfinite(result['adc']).all((1,2))&(result['adc']>0).all((1,2));result['valid']=positive
    np.savez_compressed(target,**result)
    write(OUT/'supplement_audit.json',dict(complete=True,created_utc=stamp(),events=n,available=int(result['available'].sum()),valid=int(positive.sum()),
        unavailable=int((~result['available']).sum()),nonpositive_available=int((result['available']&~positive).sum()),records=records,
        cache_original_tof_adca_finite=int(np.isfinite(a['tof_adca']).sum()),source_sha256=digest(ROOT/'../amsd69n/analysis_amsd69n.C'),
        lineage='analysis_amsd69n.C: tof_Aadc[layer][side] = associated TofClusterHR::Aadc[side]; no beta argument',
        output_sha256=digest(target),measured_sha256=digest(MCOUT/'measured.npz'),MC_truth_read=False,RICH_read=False))
    print('ADC_COMPLETE',n,int(positive.sum()),flush=True)


def inputs(dataset):
    folder=MCOUT if dataset=='mc' else FLIGHT;a=load(folder/'measured.npz');g=load(folder/'moments.npz')
    for k in ['run','event']:
        if not np.array_equal(a[k],g[k]):raise ValueError('Material identity differs')
    if dataset=='mc':
        b=load(OUT/'mc_supplement.npz')
        for k in ['run','event']:
            if not np.array_equal(a[k],b[k]):raise ValueError('ADC identity differs')
        adc=b['adc'];domain=b['valid']
    else:
        adc=a['adc'];domain=np.isfinite(adc).all((1,2))&(adc>0).all((1,2))
    domain&=np.isfinite(a['edep']).all(1)&(a['edep']>0).all(1)
    train=domain&(a['fold']!=0);test=domain&(a['fold']==0)
    if np.intersect1d(a['run'][train],a['run'][test]).size:raise ValueError('Run leakage')
    return folder,a,g,adc,domain,train,test


def feature_arrays(a,adc,domain,train,variant,normalization=None):
    geo=feature.raw_features(a);names=list(feature.FEATURE_NAMES)
    # Invalid/missing rows are outside the predeclared common measured domain.
    safe=np.where(np.isfinite(adc)&(adc>0),adc,1.);dep=np.where(a['edep']>0,a['edep'],1.)
    blocks=[geo]
    if variant=='readout':
        blocks.extend([np.log(safe.sum(2)),np.log(safe[:,:,0]/safe[:,:,1]),np.log(dep)])
        names+=['log_adc_sum_%d'%l for l in range(4)]+['log_adc_end_ratio_%d'%l for l in range(4)]+['log_edep_%d'%l for l in range(4)]
    raw=np.column_stack(blocks)
    if normalization is None:
        center=raw[train].mean(0);scale=raw[train].std(0)
        if np.any(scale<1e-9):raise ValueError('Uninformative feature')
        normalization=dict(names=names,center=center.tolist(),scale=scale.tolist(),clip=3.,source='Own dataset, own training runs only')
    if names!=normalization['names']:raise ValueError('Feature order changed')
    ff=np.clip((raw-np.array(normalization['center']))/np.array(normalization['scale']),-3.,3.)
    if not np.isfinite(ff[domain]).all():raise ValueError('Nonfinite features')
    ff[~domain]=0.
    return ff,normalization


def prepare(dataset):
    parent=OUT/dataset;parent.mkdir(exist_ok=True)
    if (parent/'protocol.json').exists():raise FileExistsError('Protocol already frozen')
    folder,a,g,adc,domain,train,test=inputs(dataset);variants=['geometry','readout'] if dataset=='mc' else ['readout']
    specs={}
    for v in variants:
        ff,norm=feature_arrays(a,adc,domain,train,v);specs[v]=dict(normalization=norm,nfeatures=ff.shape[1],bounds=base.BOUNDS+[(-.75,.75)]*(6*ff.shape[1]))
    seeds={}
    for name,field,nu in CONFIGS[dataset]:
        p=(folder if dataset=='mc' else folder/'conditional')/name/'freeze.json';seeds[name]=dict(path=str(p),sha256=digest(p))
    protocol=dict(created_utc=stamp(),dataset=dataset,configs=CONFIGS[dataset],variants=specs,penalty=.03,nodes=128,evaluation_nodes=256,
        domain='Existing measured selection intersect finite positive ADC at both ends of all four layers and positive stored GetEdep; same common events for geometry/readout controls',
        domain_events=int(domain.sum()),train_events=int(train.sum()),test_events=int(test.sum()),train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),
        missing_ADC_policy='Excluded by detector-schema availability, never imputed from isotope; no C13 claim where ADC absent',
        initialization='Matching own-dataset constant-covariance model with zero slopes; physical kappa1 and mixture0.5 with zero feature slopes',
        selection='Minimum penalized training likelihood of two starts; no truth or holdout likelihood selection',
        output='Free GLS slope only, event precision uses geometry, direct end ADC and stored GetEdep; no isotope/source posterior in output',
        invariance='Exact at stored-time interface when covariates fixed; not a claim about rerunning upstream beta-fed reconstruction',
        energy_caveat='GetEdep is calibrated signal, not newly inferred visible MeV; no additive timing or deposited-energy response fit',
        MC_truth_read=False,RICH_read=False,flight_MC_parameter_transfer=False,baseline_seeds=seeds,
        measurement_sha256=digest(folder/'measured.npz'),geometry_sha256=digest(folder/'moments.npz'),supplement_sha256=digest(OUT/'mc_supplement.npz') if dataset=='mc' else None,
        script_sha256=digest(Path(__file__)),feature_library_sha256=digest(feature.OUT/'libfeature.so'))
    write(parent/'protocol.json',protocol);print('AMPLITUDE_PROTOCOL',dataset,protocol['domain_events'],protocol['train_events'],protocol['test_events'],flush=True)


def fit(dataset,name,variant):
    folder,a,g,adc,domain,train,test=inputs(dataset);parent=OUT/dataset;p=json.loads((parent/'protocol.json').read_text());spec=p['variants'][variant]
    if digest(folder/'measured.npz')!=p['measurement_sha256'] or digest(folder/'moments.npz')!=p['geometry_sha256']:raise ValueError('Inputs changed')
    if digest(Path(__file__))!=p['script_sha256']:raise ValueError('Source changed after protocol')
    ff,norm=feature_arrays(a,adc,domain,train,variant,spec['normalization']);_,field,nu=[x for x in p['configs'] if x[0]==name][0]
    dest=parent/(name+'_'+variant);dest.mkdir(exist_ok=True)
    if (dest/'freeze.json').exists():raise FileExistsError('Already frozen')
    seed=p['baseline_seeds'][name]
    if digest(seed['path'])!=seed['sha256']:raise ValueError('Seed changed')
    initial=np.array(json.loads(Path(seed['path']).read_text())['model']['theta'])
    obj=feature.FeatureConditional(a,g,train,ff,penalty=p['penalty'],nodes=p['nodes'],time_field=field,rigidity_nu=nu,threads=4)
    starts=[np.r_[initial,np.zeros(6*spec['nfeatures'])],np.r_[base.initial(1,.5),np.zeros(6*spec['nfeatures'])]];results=[]
    for i,x in enumerate(starts):
        path=dest/('start_%d.json'%i)
        result=json.loads(path.read_text()) if path.exists() else solve(obj,x,spec['bounds'])
        write(path,result);results.append(result);print('AMPLITUDE_FIT',dataset,name,variant,i,result['penalized_nll'],result['projected_gradient_max'],result['elapsed_seconds'],flush=True)
    best=min(results,key=lambda x:x['penalized_nll'])
    if not best['success'] or best['projected_gradient_max']>1e-3:raise RuntimeError('Fit did not converge')
    frozen=dict(complete=True,created_utc=stamp(),model=best,starts=results,protocol_sha256=digest(parent/'protocol.json'),MC_truth_read=False,RICH_read=False)
    write(dest/'freeze.json',frozen)
    theta=np.array(best['theta']);xx=g['x'][test];tt=a[field][test];w,Q=feature.precision_weights(theta,xx,ff[test]);u=np.sum(w*tt,axis=1)
    c=np.linspace(-2,2,len(u));b=np.linspace(-.1,.1,len(u));err=np.sum(w*(tt+c[:,None]+b[:,None]*xx),axis=1)-u-b
    if np.max(abs(err))>2e-12:raise ValueError('Affine response failed')
    out={k:a[k][test] for k in ['run','event','R']};out.update(u=u,weights=w,domain_indices=np.flatnonzero(test))
    evalobj=feature.FeatureConditional(a,g,test,ff,penalty=0,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    loss,grad,per=evalobj.evaluate(theta,events=True);out['nll']=per
    old=base.ConditionalModel(a,g,test,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    oloss,_,oper=old.evaluate(initial,events=True);out['original_baseline_nll']=oper
    np.savez_compressed(dest/'predictions_before_truth.npz',**out)
    write(dest/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),events=len(u),mean_nll=loss,original_baseline_nll=oloss,delta_nll=loss-oloss,
        injection_max_error=float(np.max(abs(err))),freeze_sha256=digest(dest/'freeze.json'),output_sha256=digest(dest/'predictions_before_truth.npz'),
        truth_used=False,RICH_read=False,beta_above_one_retained=int((u<1).sum()),median_weights=np.median(w,axis=0)))
    print('AMPLITUDE_FROZEN',dataset,name,variant,'delta NLL',loss-oloss,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['supplement','prepare','fit']);parser.add_argument('--dataset',choices=['mc','flight']);parser.add_argument('--name',default='native_gaussianR');parser.add_argument('--variant',choices=['geometry','readout'],default='readout');args=parser.parse_args()
    if args.stage=='supplement':supplement()
    elif args.stage=='prepare':prepare(args.dataset)
    else:fit(args.dataset,args.name,args.variant)
