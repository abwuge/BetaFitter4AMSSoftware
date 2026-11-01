#!/usr/bin/env python3
"""Fit exogenous geometry/era time-noise variations after ADC-mean selection."""
import argparse
import json
from pathlib import Path
import numpy as np
import conditional_time_noise as noise
import periodic_joint_adc as period
import joint_adc_likelihood as joint
import flight_conditional_likelihood as timing
import free_timing_transport as center
from research_amplitude_noise_20261013 import inputs
from research_flight_conditional import ROOT,write,load,digest
from research_conditional_mc import stamp
from joint_adc_scaled_optimizer import solve_scaled

OUT=noise.OUT
SOURCES=['research_conditional_noise.py','conditional_time_noise.py','conditional_time_noise.cpp','periodic_joint_adc.py',
    'joint_adc_likelihood.py','joint_adc_scaled_optimizer.py','free_timing_transport.py']
FAMILIES={'flight':('newcalib_gaussianR','t'),'mc':('native_gaussianR','t_native')}


def plan():
    if (OUT/'protocol.json').exists():raise FileExistsError('Conditional-noise plan already exists')
    selected=json.loads((period.OUT/'measurement_comparison/evaluation.json').read_text())
    if not selected['complete']:raise ValueError('Finish the measurement-only ADC period selection')
    datasets={}
    for dataset,(family,field) in FAMILIES.items():
        folder,a,g,adc,domain,train,test=inputs(dataset)
        variant=selected['selected'][dataset][family];source=period.OUT/dataset/(family+'_'+variant)/'freeze.json'
        fr=json.loads(source.read_text());features,_=period.features(a,train,dataset,variant,fr['normalization'])
        variants={}
        for label in noise.VARIANTS:
            hf,normalization=noise.features(a,train,dataset,label)
            variants[label]=dict(normalization=normalization,nnoise=hf.shape[1],
                bounds=joint.bounds(features.shape[2])+[[-.35,.35]]*(3*hf.shape[1]))
        datasets[dataset]=dict(family=family,time_field=field,ADC_period_variant=variant,mean_freeze=str(source.relative_to(ROOT)),
            mean_freeze_sha256=digest(source),mean_normalization=fr['normalization'],initial_theta=fr['model']['theta'],
            data_sha256=digest(folder/'measured.npz'),moment_sha256=digest(folder/'moments.npz'),
            train_events=int(train.sum()),test_events=int(test.sum()),variants=variants)
    write(OUT/'protocol.json',dict(created_utc=stamp(),datasets=datasets,source_sha256={name:digest(ROOT/'analysis'/name) for name in SOURCES},
        library_sha256=digest(OUT/'libconditionalnoise.so'),measurement_selection_sha256=digest(period.OUT/'measurement_comparison/evaluation.json'),
        model='Scale the first three columns of the seven-dimensional precision Cholesky factor by exp(Gamma * exogenous features). ADC marginal noise and conditional time-on-ADC regression are preserved for fixed base parameters.',
        features='Only measured L5 geometry and declared UTC era (MC run-hash negative-control era). No event ADC, TOF speed, isotope posterior, RICH or truth enters the noise features.',
        conditioning='Given geometry and era the joint Student density remains normalized; the inherited hit/charge selection and source-response assumptions are not promoted to a complete raw-readout likelihood.',
        response='Continuous GLS has exact unit full time-plus-ADC mean-signal response because noise features are unchanged by that injection.',
        time_offsets_added=False,MC_coefficients_to_flight=False,reference_values_for_fitting=False,
        variants=noise.VARIANTS,constant_model='Reuse the selected ADC-mean model after numerical identity validation; do not refit an identical zero-feature density.',
        new_starts_per_variant=2,initialization='Own-dataset selected ADC-mean parameters, with Gamma=0 or deterministic N(0,.04) Gamma. All original nuisance parameters remain free.',
        mean_penalty=.01,noise_penalty=.1,train_nodes=256,evaluate_nodes=512,maximum_rounds=3,projected_gradient_max=1e-3,
        optimizer='Existing curvature-scaled L-BFGS-B, 6000 iterations per round, unchanged unscaled objective and bounds.',
        selection='Within each declared family, minimum held-out joint NLL among constant, geometry and geometry_period. Never select using beta-reference performance.',
        prospective_beta_status='The already exposed MC and RICH pools remain development data; no blind acceptance claim from these pools.'))
    print('CONDITIONAL_NOISE_PLAN',flush=True)


def context(dataset,variant):
    p=json.loads((OUT/'protocol.json').read_text())
    for name,fingerprint in p['source_sha256'].items():
        if digest(ROOT/'analysis'/name)!=fingerprint:raise ValueError('Noise implementation changed '+name)
    if digest(OUT/'libconditionalnoise.so')!=p['library_sha256']:raise ValueError('Noise binary changed')
    item=p['datasets'][dataset];folder,a,g,adc,domain,train,test=inputs(dataset)
    if digest(folder/'measured.npz')!=item['data_sha256'] or digest(folder/'moments.npz')!=item['moment_sha256']:
        raise ValueError('Noise measurements changed')
    if digest(ROOT/item['mean_freeze'])!=item['mean_freeze_sha256']:raise ValueError('ADC-mean seed changed')
    features,_=period.features(a,train,dataset,item['ADC_period_variant'],item['mean_normalization'])
    hf,_=noise.features(a,train,dataset,variant,item['variants'][variant]['normalization'])
    return p,item,a,g,adc,train,test,features,hf


def obj(dataset,variant,select=None,nodes=None):
    p,item,a,g,adc,train,test,f,h=context(dataset,variant)
    model=noise.ConditionalNoise(a,g,train if select is None else select,f,adc,h,nodes=nodes or p['train_nodes'],
        time_field=item['time_field'],rigidity_nu=0,threads=4,penalty=p['mean_penalty'],noise_penalty=p['noise_penalty'])
    return model,(p,item,a,g,adc,train,test,f,h)


def verify():
    if (OUT/'numerical_verification.json').exists():raise FileExistsError('Noise numerical verification already exists')
    records=[];rng=np.random.RandomState(1580)
    for dataset in ['flight','mc']:
        p,item,a,g,adc,train,test,f,h=context(dataset,'geometry_period');choose=np.flatnonzero(train)[:8]
        legacy=joint.JointADC(a,g,choose,f,adc,nodes=64,time_field=item['time_field'],rigidity_nu=0,threads=1)
        candidate=noise.ConditionalNoise(a,g,choose,f,adc,h,nodes=64,time_field=item['time_field'],rigidity_nu=0,threads=1)
        initial=np.asarray(item['initial_theta']);zero=np.r_[initial,np.zeros(3*h.shape[1])]
        old=legacy.evaluate(initial,events=True);new=candidate.evaluate(zero,events=True)
        constant_error=max(abs(old[0]-new[0]),np.max(abs(old[1]-new[1][:len(initial)])),np.max(abs(old[2]-new[2])))
        theta=zero.copy();theta[len(initial):]=rng.normal(0,.035,len(theta)-len(initial))
        value,gradient,per=candidate.evaluate(theta,events=True)
        numpy_error=max(abs(per[i]+candidate.numpy_event(theta,i)) for i in range(len(choose)))
        finite=np.empty(len(theta))
        for j in range(len(theta)):
            step=2e-6*max(1,abs(theta[j]));hi=theta.copy();lo=theta.copy();hi[j]+=step;lo[j]-=step
            finite[j]=(candidate.evaluate(hi)[0]-candidate.evaluate(lo)[0])/(2*step)
        gradient_error=float(np.max(abs(gradient-finite)))
        if constant_error>2e-10 or numpy_error>2e-10 or gradient_error>5e-6:raise ValueError('Conditional-noise density or gradient check failed')
        rec=noise.reconstruct(theta,g['x'][choose],candidate.y,f[choose],h[choose]);delta=np.linspace(-.03,.03,len(choose))
        shifted=noise.reconstruct(theta,g['x'][choose],candidate.y+delta[:,None]*rec['direction'],f[choose],h[choose])
        response=float(np.max(abs(shifted['joint_u']-rec['joint_u']-delta)))
        if response>2e-11:raise ValueError('Noise signal-response identity failed')
        records.append(dict(dataset=dataset,parameters=len(theta),old_model_identity_max_error=float(constant_error),
            numpy_density_max_error=float(numpy_error),all_parameter_gradient_max_error=gradient_error,signal_injection_max_error=response))
    write(OUT/'numerical_verification.json',dict(complete=True,created_utc=stamp(),records=records,protocol_sha256=digest(OUT/'protocol.json')))
    print('CONDITIONAL_NOISE_NUMERICAL',json.dumps(records),flush=True)


def start(dataset,variant,number):
    if variant=='constant':raise ValueError('The constant noise model is the already fitted selected mean response')
    check=json.loads((OUT/'numerical_verification.json').read_text())
    if not check['complete'] or check['protocol_sha256']!=digest(OUT/'protocol.json'):raise ValueError('Verify the current noise protocol')
    model,c=obj(dataset,variant);p,item,a,g,adc,train,test,f,h=c
    folder=OUT/dataset/variant;folder.mkdir(parents=True,exist_ok=True);path=folder/('start_%d.json'%number)
    if path.exists():raise FileExistsError('Noise start already saved')
    if number not in [0,1]:raise ValueError('Only two starts declared')
    theta=np.r_[item['initial_theta'],np.zeros(3*h.shape[1])]
    if number:theta[-3*h.shape[1]:]=np.random.RandomState(802+len(theta)+(0 if dataset=='flight' else 10)).normal(0,.04,3*h.shape[1])
    rounds=[]
    for round_index in range(p['maximum_rounds']):
        result=solve_scaled(model,theta,item['variants'][variant]['bounds']);rounds.append(result);theta=np.asarray(result['theta'])
        print('CONDITIONAL_NOISE_ROUND',dataset,variant,number,round_index,result['penalized_nll'],result['projected_gradient_max'],flush=True)
        if result['projected_gradient_max']<p['projected_gradient_max']:break
    best=min(rounds,key=lambda r:r['penalized_nll'])
    record=dict(best);record.update(dataset=dataset,variant=variant,start=number,rounds=rounds,created_utc=stamp(),
        numerically_converged=best['projected_gradient_max']<p['projected_gradient_max'],protocol_sha256=digest(OUT/'protocol.json'))
    write(path,record)


def freeze(dataset,variant):
    p,item,a,g,adc,train,test,f,h=context(dataset,variant);folder=OUT/dataset/variant;folder.mkdir(parents=True,exist_ok=True)
    if (folder/'freeze.json').exists():raise FileExistsError('Noise output already frozen')
    if variant=='constant':
        source=json.loads((ROOT/item['mean_freeze']).read_text());best=dict(source['model']);source_paths=[ROOT/item['mean_freeze']]
    else:
        source_paths=[folder/('start_%d.json'%i) for i in range(2)];records=[json.loads(path.read_text()) for path in source_paths]
        for record in records:
            if record['protocol_sha256']!=digest(OUT/'protocol.json'):raise ValueError('Noise start protocol changed')
        best=min(records,key=lambda record:record['penalized_nll'])
        if not best['numerically_converged']:raise ValueError('Best noise likelihood start is not converged')
    theta=np.asarray(best['theta'])
    write(folder/'freeze.json',dict(created_utc=stamp(),dataset=dataset,variant=variant,model=best,
        mean_normalization=item['mean_normalization'],noise_normalization=item['variants'][variant]['normalization'],
        protocol_sha256=digest(OUT/'protocol.json'),source_sha256={str(path.relative_to(ROOT)):digest(path) for path in source_paths},
        RICH_read=False,MC_truth_read=False))
    full,_=obj(dataset,variant,select=test,nodes=p['evaluate_nodes']);value,gradient,per=full.evaluate(theta,events=True,regularized=False)
    coarse,_=obj(dataset,variant,select=test,nodes=p['train_nodes']);_,_,low=coarse.evaluate(theta,events=True,regularized=False)
    selected=np.flatnonzero(test);rec=noise.reconstruct(theta,g['x'][selected],full.y,f[selected],h[selected])
    raw_x=g['x'][selected];xx=raw_x-raw_x.mean(1)[:,None];weight=xx/np.sum(xx*xx,axis=1)[:,None]
    branches={'ols':(np.sum(weight*a[item['time_field']][selected],axis=1),weight),
        'time':(rec['time_u'],rec['time_weights']),'joint':(rec['joint_u'],rec['weights'][:,:3].dot(timing.D))}
    for label in ['speed_null','speed_gain_null']:branches[label]=(rec[label]['u'],rec[label]['time_weights'])
    output={key:a[key][selected] for key in ['run','event','R']};output['nll']=per
    for label,(u,w) in branches.items():
        transported=center.solve_center(u,a['R'][selected],w,g['moment'][selected],theta[:2])
        if not transported['converged'].all():raise ValueError('Conditional-noise material inverse failed')
        output[label+'_u']=u;output[label+'_transport_u']=transported['u']
    output['time_only_gain']=rec['time_injection_gain'];output['information_width_gain']=1-np.sqrt(rec['time_information']/rec['information'])
    np.savez_compressed(folder/'predictions_before_reference.npz',**output)
    write(folder/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),dataset=dataset,variant=variant,events=len(selected),
        heldout_nll=float(value),quadrature_mean256_minus512=float((low-per).mean()),quadrature_rms256_minus512=float(np.sqrt(np.mean((low-per)**2))),
        material_kappa=theta[:2],ADC_lambda=np.exp(theta[35:39]),
        time_only_gain_percentiles=np.percentile(rec['time_injection_gain'],[0,5,50,95,100]),
        information_width_gain_percentiles=np.percentile(output['information_width_gain'],[0,5,50,95,100]),
        freeze_sha256=digest(folder/'freeze.json'),predictions_sha256=digest(folder/'predictions_before_reference.npz')))
    print('CONDITIONAL_NOISE_FROZEN',dataset,variant,value,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['plan','verify','start','freeze']);parser.add_argument('--dataset',choices=['flight','mc'])
    parser.add_argument('--variant',choices=noise.VARIANTS);parser.add_argument('--start',type=int,default=0);args=parser.parse_args()
    if args.stage in ['plan','verify']:globals()[args.stage]()
    elif args.dataset is None or args.variant is None:parser.error('dataset and variant required')
    elif args.stage=='start':start(args.dataset,args.variant,args.start)
    else:freeze(args.dataset,args.variant)
