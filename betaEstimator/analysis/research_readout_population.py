#!/usr/bin/env python3
"""Ablate conditional latent population versus conditional time-noise precision."""
import argparse,json,time
from pathlib import Path
import numpy as np
import conditional_readout_population as population
import flight_conditional_likelihood as base
from research_amplitude_noise_20261013 import OUT as AMPOUT,ROOT,inputs,feature_arrays,CONFIGS
from research_feature_noise import solve
from research_flight_conditional import write,load,digest
from research_conditional_mc import stamp

OUT=population.OUT
VARIANTS={'mc':['population_only','joint'],'flight':['joint']}


def lift(theta,nfeatures,geometry_only=False):
    theta=np.asarray(theta,float);result=np.zeros(13+10*nfeatures);result[:13]=theta[:13]
    slopes=result[13:].reshape(10,nfeatures)
    if geometry_only:slopes[:6,:6]=theta[13:].reshape(6,6)
    else:slopes[:6,:]=theta[13:].reshape(6,nfeatures)
    return result


def prepare(dataset):
    directory=OUT/dataset;directory.mkdir(exist_ok=True)
    if (directory/'protocol.json').exists():raise FileExistsError('Population protocol already frozen')
    folder,a,g,adc,domain,train,test=inputs(dataset);old=json.loads((AMPOUT/dataset/'protocol.json').read_text())
    norm=old['variants']['readout']['normalization'];ff,norm=feature_arrays(a,adc,domain,train,'readout',norm);nf=ff.shape[1]
    variants={}
    for variant in VARIANTS[dataset]:
        bounds=base.BOUNDS+[(-.75,.75)]*(10*nf)
        if variant=='population_only':
            for h in range(6):
                for j in range(6,nf):bounds[13+h*nf+j]=(0.,0.)
        variants[variant]=dict(bounds=bounds,nparameters=len(bounds),free_parameters=sum(x!=y for x,y in bounds))
    seeds={}
    for name,field,nu in CONFIGS[dataset]:
        seeds[name]={}
        for variant in VARIANTS[dataset]:
            old_variant='geometry' if variant=='population_only' else 'readout';path=AMPOUT/dataset/(name+'_'+old_variant)/'freeze.json'
            seeds[name][variant]=dict(path=str(path),sha256=digest(path),geometry_only=variant=='population_only')
    write(directory/'protocol.json',dict(created_utc=stamp(),dataset=dataset,variants=variants,configs=CONFIGS[dataset],feature_normalization=norm,nfeatures=nf,
        penalty=.03,penalty_formula='0.5 * 0.03 * sum(all feature coefficient squares), both precision and latent-population coefficients',
        nodes=128,evaluation_nodes=256,qrange=[.02,.7],seeds=seeds,
        hypothesis='Amplitude encodes the physical speed population as well as detector noise; allow q-spectrum and isotope log-odds coefficients to depend on the same measured covariates',
        conditional_density='p(Dt | r, ADC, GetEdep, geometry), normalized in all three Dt dimensions; latent population is conditional, not a measured cosmic isotope abundance',
        parameterization='Six precision Cholesky parameters and four latent-spectrum/isotope parameters have linear predictors in 18 standardized measured features',
        ablation='population_only retains geometric covariance and conditions only latent population on readout; joint also conditions covariance on readout',
        fixed_baseline='No additive timing function is fitted; original timing baseline, geometry and two material groups retained',
        output='Free GLS time slope; discard all event population coefficients and all isotope posteriors at reconstruction',
        optimization='Two starts selected only by penalized training NLL; continue unsuccessful starts up to three 900-iteration rounds before accepting failure',
        evaluation='All fixed MC candidates are reported in the already-used development pool; not a fresh blind discovery test',
        train_events=int(train.sum()),test_events=int(test.sum()),train_runs=np.unique(a['run'][train]),test_runs=np.unique(a['run'][test]),
        MC_truth_read=False,RICH_read=False,MC_parameters_transferred_to_flight=False,
        amplitude_protocol_sha256=digest(AMPOUT/dataset/'protocol.json'),measurement_sha256=digest(folder/'measured.npz'),geometry_sha256=digest(folder/'moments.npz'),
        implementation_sha256={name:digest(ROOT/'analysis'/name) for name in ['conditional_readout_population.cpp','conditional_readout_population.py','research_readout_population.py']},
        library_sha256=digest(OUT/'libpopulation.so')))
    print('POPULATION_PROTOCOL',dataset,'events',int(train.sum()),int(test.sum()),flush=True)


def fit(dataset,name,variant):
    folder,a,g,adc,domain,train,test=inputs(dataset);parent=OUT/dataset;protocol=json.loads((parent/'protocol.json').read_text())
    for f,h in protocol['implementation_sha256'].items():
        if digest(ROOT/'analysis'/f)!=h:raise ValueError('Population source changed '+f)
    if digest(OUT/'libpopulation.so')!=protocol['library_sha256']:raise ValueError('Library changed')
    if digest(folder/'measured.npz')!=protocol['measurement_sha256'] or digest(folder/'moments.npz')!=protocol['geometry_sha256']:raise ValueError('Measured input changed')
    spec=protocol['variants'][variant];seed=protocol['seeds'][name][variant]
    if digest(seed['path'])!=seed['sha256']:raise ValueError('Population seed changed')
    dest=parent/(name+'_'+variant);dest.mkdir(exist_ok=True)
    if (dest/'freeze.json').exists():raise FileExistsError('Population already frozen')
    ff,norm=feature_arrays(a,adc,domain,train,'readout',protocol['feature_normalization']);_,field,nu=[c for c in protocol['configs'] if c[0]==name][0]
    original=np.array(json.loads(Path(seed['path']).read_text())['model']['theta']);initial=lift(original,protocol['nfeatures'],seed['geometry_only'])
    obj=population.PopulationConditional(a,g,train,ff,penalty=protocol['penalty'],nodes=protocol['nodes'],time_field=field,rigidity_nu=nu,threads=4)
    starts=[initial,np.r_[base.initial(1,.5),np.zeros(10*protocol['nfeatures'])]];allfits=[]
    for i,point in enumerate(starts):
        path=dest/('start_%d.json'%i)
        if path.exists():result=json.loads(path.read_text())
        else:
            history=[]
            for repeat in range(3):
                result=solve(obj,point,spec['bounds']);history.append({k:v for k,v in result.items() if k!='theta'})
                print('POPULATION_FIT',dataset,name,variant,i,repeat,result['penalized_nll'],result['projected_gradient_max'],result['elapsed_seconds'],flush=True)
                if result['success'] and result['projected_gradient_max']<1e-3:break
                point=np.asarray(result['theta'])
            result['optimization_rounds']=history;write(path,result)
        allfits.append(result)
    best=min(allfits,key=lambda v:v['penalized_nll'])
    if not best['success'] or best['projected_gradient_max']>=1e-3:raise RuntimeError('Population fit failed to converge')
    write(dest/'freeze.json',dict(complete=True,created_utc=stamp(),model=best,starts=allfits,protocol_sha256=digest(parent/'protocol.json'),MC_truth_read=False,RICH_read=False))
    theta=np.array(best['theta']);w,Q=population.precision_weights(theta,g['x'][test],ff[test]);u=np.sum(w*a[field][test],axis=1)
    injection=np.linspace(-.1,.1,len(u));offset=np.linspace(-3,3,len(u))
    changed=a[field][test]+offset[:,None]+injection[:,None]*g['x'][test]
    response=np.sum(w*changed,axis=1)-u-injection
    if np.max(abs(response))>2e-12:raise ValueError('Population output response failed')
    ev=population.PopulationConditional(a,g,test,ff,penalty=0,nodes=protocol['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    value,grad,per=ev.evaluate(theta,events=True);ev128=population.PopulationConditional(a,g,test,ff,penalty=0,nodes=128,time_field=field,rigidity_nu=nu,threads=4)
    v128,_,p128=ev128.evaluate(theta,events=True)
    data={k:a[k][test] for k in ['run','event','R']};data.update(u=u,weights=w,nll=per,nll128=p128,features=ff[test])
    np.savez_compressed(dest/'predictions_before_truth.npz',**data)
    write(dest/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),events=len(u),mean_nll=value,quadrature_mean128_minus256=v128-value,
        quadrature_rms128_minus256=float(np.sqrt(np.mean((p128-per)**2))),quadrature_max_abs=float(np.max(abs(p128-per))),
        injection_max_error=float(np.max(abs(response))),freeze_sha256=digest(dest/'freeze.json'),output_sha256=digest(dest/'predictions_before_truth.npz'),
        MC_truth_read=False,RICH_read=False,median_weights=np.median(w,axis=0),beta_above_one_retained=int((u<1).sum())))
    print('POPULATION_FROZEN',dataset,name,variant,'holdout_nll',value,'quadrature_rms',float(np.sqrt(np.mean((p128-per)**2))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','fit']);p.add_argument('--dataset',choices=['mc','flight'],required=True);p.add_argument('--name',default='native_gaussianR');p.add_argument('--variant',choices=['population_only','joint'],default='joint');args=p.parse_args()
    if args.stage=='prepare':prepare(args.dataset)
    else:fit(args.dataset,args.name,args.variant)
