#!/usr/bin/env python3
"""Optimization-only amendment for the frozen joint-ADC statistical model."""
import argparse,json
from pathlib import Path
import numpy as np
import joint_adc_likelihood as model
import flight_conditional_likelihood as base
import free_timing_transport as transport
from joint_adc_scaled_optimizer import solve_scaled
from research_joint_adc import OUT,ROOT,CONFIGS,inputs
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp

REFINE=OUT/'optimizer_refinement'


def prepare():
    REFINE.mkdir(exist_ok=True)
    if (REFINE/'protocol.json').exists():raise FileExistsError('Joint optimizer amendment already frozen')
    records=[]
    for dataset in ['mc','flight']:
        p=json.loads((OUT/dataset/'protocol.json').read_text())
        for name,field,nu in CONFIGS[dataset]:
            d=OUT/dataset/name
            if (d/'freeze.json').exists():raise FileExistsError('Never replace an already frozen model')
            prior=[]
            for f in sorted(d.glob('start_*.json')):
                a=json.loads(f.read_text());prior.append(dict(path=str(f),sha256=digest(f),penalized_nll=a['penalized_nll'],success=a['success'],projected_gradient_max=a['projected_gradient_max']))
            seed=min(prior,key=lambda a:a['penalized_nll']) if prior else None
            records.append(dict(dataset=dataset,family=name,time_field=field,nu=nu,original_protocol_sha256=digest(OUT/dataset/'protocol.json'),prior_attempts=prior,warm_start=seed))
    write(REFINE/'protocol.json',dict(created_utc=stamp(),records=records,
        reason='Original unscaled memory20 L-BFGS-B exhausted repeated iteration limits; no outcome scores were opened for these candidates',
        statistical_model_changed=False,bounds_changed=False,penalty_changed=False,data_changed=False,
        optimizer='Diagonal observed-curvature coordinate scaling, memory80 L-BFGS-B, maxiter6000, ftol2e-13, gtol1e-8, maxls50',
        preconditioning='Central derivative of analytic gradient with step1e-4*max(1,abs(theta)); 1/sqrt(abs(diagonal)) clipped to [0.001,10]. Only a coordinate transformation; not a covariance/uncertainty estimate.',
        starts='Lowest training-NLL complete retained original attempt, otherwise original own-dataset conditional-time initialization with ADC lambda1; plus independent physical initialization with ADC lambda3',
        acceptance='Optimizer success AND original-coordinate projected gradient max <1e-4; select lowest penalized training NLL among the two prescribed starts',
        output_location='Previously empty original joint model freeze/output paths, with explicit optimizer_amendment_sha256; retain all original failure files',
        interruption_sha256=digest(OUT/'initial_optimizer_interruption.json'),
        source_sha256={name:digest(ROOT/'analysis'/name) for name in ['joint_adc_scaled_optimizer.py','refine_joint_adc_optimizer.py']},MC_truth_read=False,RICH_read=False))
    print('JOINT_OPTIMIZER_AMENDMENT',len(records),flush=True)


def check(dataset,name):
    amendment=json.loads((REFINE/'protocol.json').read_text());p=json.loads((OUT/dataset/'protocol.json').read_text())
    row=[r for r in amendment['records'] if r['dataset']==dataset and r['family']==name][0]
    if digest(OUT/dataset/'protocol.json')!=row['original_protocol_sha256']:raise ValueError('Original joint protocol changed')
    for source,expected in list(amendment['source_sha256'].items())+list(p['implementation_sha256'].items()):
        if digest(ROOT/'analysis'/source)!=expected:raise ValueError('Joint optimizer/model source changed '+source)
    if digest(OUT/'libjointadc.so')!=p['library_sha256']:raise ValueError('Joint likelihood library changed')
    for attempt in row['prior_attempts']:
        if digest(attempt['path'])!=attempt['sha256']:raise ValueError('Original retained attempt changed')
    return amendment,p,row


def fit(dataset,name):
    amendment,p,row=check(dataset,name);folder,a,g,adc,domain,train,test=inputs(dataset)
    if digest(folder/'measured.npz')!=p['measurement_sha256'] or digest(folder/'moments.npz')!=p['geometry_sha256']:raise ValueError('Joint input changed')
    destination=OUT/dataset/name;destination.mkdir(exist_ok=True)
    if (destination/'freeze.json').exists():raise FileExistsError('Refinement never overwrites a frozen model')
    store=REFINE/dataset/name;store.mkdir(parents=True,exist_ok=True)
    ff,norm=model.response_features(a,train,p['feature_normalization']);field=row['time_field'];nu=row['nu']
    obj=model.JointADC(a,g,train,ff,adc,penalty=p['penalty'],nodes=p['nodes'],time_field=field,rigidity_nu=nu,threads=4)
    if row['warm_start'] is not None:warm=np.asarray(json.loads(Path(row['warm_start']['path']).read_text())['theta'])
    else:
        seed=p['seeds'][name]
        if digest(seed['path'])!=seed['sha256']:raise ValueError('Conditional initialization changed')
        own=np.asarray(json.loads(Path(seed['path']).read_text())['model']['theta']);warm=model.parameters_initial(ff,adc,train,1.,own)
    points=[warm,model.parameters_initial(ff,adc,train,3.)];fits=[]
    for i,point in enumerate(points):
        path=store/('start_%d.json'%i)
        result=json.loads(path.read_text()) if path.exists() else solve_scaled(obj,point,p['bounds'])
        write(path,result);fits.append(result)
        print('JOINT_REFINED_START',dataset,name,i,result['success'],result['penalized_nll'],result['projected_gradient_max'],result['iterations'],result['elapsed_seconds'],flush=True)
    best=min(fits,key=lambda v:v['penalized_nll'])
    if not best['success'] or best['projected_gradient_max']>=1e-4:raise RuntimeError('Scaled joint optimizer has not met the stricter criterion')
    write(destination/'freeze.json',dict(complete=True,created_utc=stamp(),model=best,starts=fits,protocol_sha256=digest(OUT/dataset/'protocol.json'),
        optimizer_amendment_sha256=digest(REFINE/'protocol.json'),retained_prior_attempts=row['prior_attempts'],MC_truth_read=False,RICH_read=False))
    export(dataset,name,a,g,adc,test,ff,p,row,best)


def export(dataset,name,a,g,adc,test,ff,p,row,best):
    dest=OUT/dataset/name;theta=np.asarray(best['theta']);field=row['time_field'];nu=row['nu']
    ev=model.JointADC(a,g,test,ff,adc,penalty=0,nodes=p['evaluation_nodes'],time_field=field,rigidity_nu=nu,threads=4)
    value,grad,per=ev.evaluate(theta,events=True);result=model.reconstruct(theta,g['x'][test],ev.y,ff[test])
    data={k:a[k][test] for k in ['run','event','R']}
    for k in ['joint_u','time_u','adc_u','weights','time_weights','adc_weights','time_injection_gain','adc_injection_gain','information','time_information','additional_adc_information']:data[k]=result[k]
    inverse={}
    for label,weight in [('time',result['time_weights']),('joint',result['weights'][:,:3].dot(base.D))]:
        solved=transport.solve_center(result[label+'_u'],a['R'][test],weight,g['moment'][test],theta[:2])
        if not solved['converged'].all():raise ValueError('Joint loss inverse did not converge')
        data[label+'_transport_u']=solved['u'];inverse[label]=dict(max_residual=float(np.max(abs(solved['residual']))),median_correction=float(np.median(solved['correction'])),max_abs_correction=float(np.max(abs(solved['correction']))))
    delta=np.linspace(-.1,.1,len(per));changed=model.reconstruct(theta,g['x'][test],ev.y+delta[:,None]*result['direction'],ff[test])
    error=changed['joint_u']-result['joint_u']-delta
    if np.max(abs(error))>2e-12:raise ValueError('Joint coherent response failed')
    sample=model.JointADC(a,g,test,ff,adc,penalty=0,nodes=128,time_field=field,rigidity_nu=nu,threads=4);v128,_,p128=sample.evaluate(theta,events=True)
    data['nll']=per;data['nll128']=p128;np.savez_compressed(dest/'predictions_before_truth.npz',**data)
    write(dest/'measurement_evaluation.json',dict(complete=True,created_utc=stamp(),events=len(per),mean_joint_nll=value,
        quadrature_mean128_minus256=v128-value,quadrature_rms128_minus256=float(np.sqrt(np.mean((p128-per)**2))),quadrature_max_abs128_minus256=float(np.max(abs(p128-per))),
        coherent_joint_injection_max_error=float(np.max(abs(error))),time_only_injection_gain_percentiles=np.percentile(result['time_injection_gain'],[0,5,50,95,100]),
        model_information_width_gain_percentiles=np.percentile(1-np.sqrt(result['time_information']/result['information']),[0,5,50,95,100]),
        model_information_scope='Exact linear mean and correctly specified Student scale only, not measured truth performance. Student covariance is 5/3 times this scale; ratios are unchanged.',
        ADC_lambda=np.exp(theta[35:39]),material_kappa=theta[:2],loss_inverse=inverse,
        beta_above_one_retained={k:int((v<1).sum()) for k,v in data.items() if k.endswith('_u')},freeze_sha256=digest(dest/'freeze.json'),
        output_sha256=digest(dest/'predictions_before_truth.npz'),optimizer_amendment_sha256=digest(REFINE/'protocol.json'),MC_truth_read=False,RICH_read=False))
    print('JOINT_REFINED_FROZEN',dataset,name,'NLL',value,'lambda',np.exp(theta[35:39]),'timegain',np.median(result['time_injection_gain']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','fit']);parser.add_argument('--dataset',choices=['mc','flight']);parser.add_argument('--name',default='native_gaussianR');args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:fit(args.dataset,args.name)
