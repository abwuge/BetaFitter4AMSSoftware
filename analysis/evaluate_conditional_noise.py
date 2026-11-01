#!/usr/bin/env python3
"""Measured-only noise selection followed by explicitly developmental beta checks."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import conditional_time_noise as noise
import periodic_joint_adc as period
import joint_adc_likelihood as joint
import flight_conditional_likelihood as timing
import free_timing_transport as center
import continuous_rich_transport as forward
from research_conditional_noise import ROOT,OUT,context
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics,matrix_difference,METRICS,DELTA

DEST=OUT/'evaluation'
POOLS={dataset:period.OUT.parent/name for dataset,name in [('flight','prospective_flight'),('mc','additional_mc')]}
BRANCHES=['ols','time','joint','speed_null','speed_gain_null']
SOURCES=['evaluate_conditional_noise.py','conditional_time_noise.py','joint_adc_likelihood.py','periodic_joint_adc.py',
    'free_timing_transport.py','continuous_rich_transport.py','evaluate_additional_mc_20261013.py']


def plan():
    DEST.mkdir(exist_ok=True)
    if (DEST/'protocol.json').exists():raise FileExistsError('Noise evaluation protocol already exists')
    write(DEST/'protocol.json',dict(created_utc=stamp(),source_sha256={name:digest(ROOT/'analysis'/name) for name in SOURCES},
        fit_protocol_sha256=digest(OUT/'protocol.json'),noise_variants=noise.VARIANTS,branches=BRANCHES,material_modes=['fit','unit'],
        pools={dataset:{name:digest(folder/name) for name in ['measured.npz','moments.npz']} for dataset,folder in POOLS.items()},
        flight_forward_geometry_sha256=digest(POOLS['flight']/'rich_geometry.npz'),
        bootstrap_replicates=4096,measurement_seed=139132026,beta_seed=140132026,ordinary_quantiles=[2.5,50,97.5],two_comparison_quantiles=[1.25,50,98.75],
        measured_selection='Per dataset, select the smallest held-out 512-node joint NLL among constant, geometry and geometry_period. Beta values are not selection inputs.',
        beta_primary_pairs=['selected joint versus its matched OLS','selected joint versus constant joint'],
        beta_secondary_pairs=['selected speed_null versus its matched OLS','selected speed_gain_null versus its matched OLS'],
        reference_scope='Previously opened development pools only. No new blind reference, recalibration interval or production acceptance.',
        no_beta_reference_training=True,no_reference_residual_cuts=True,no_MC_coefficients_to_flight=True,
        coverage='Preserve every measured event. A failed or nonfinite prediction aborts the fixed-cohort output rather than removing that event.',
        beta_above_one='Retain original free-slope outputs and zero-loss continuation for inverse beta <=1.',
        endpoint={'flight':'Existing nominal RICH emission plane with unit forward material scale','mc':'Audited primary L5 reference'},
        uncertainty='Paired whole evaluation-run bootstrap only. Training refit uncertainty and reference/systematic model error are not included.',
        bins=dict(R=[2.5,3.,3.5,4.5,6.],reference_beta=[.78,.82,.86,.9,.94,.98,1.02],year=period.YEARS)))
    print('NOISE_EVALUATION_PROTOCOL',flush=True)


def check():
    p=json.loads((DEST/'protocol.json').read_text())
    for name,fingerprint in p['source_sha256'].items():
        if digest(ROOT/'analysis'/name)!=fingerprint:raise ValueError('Noise evaluation code changed '+name)
    if digest(OUT/'protocol.json')!=p['fit_protocol_sha256']:raise ValueError('Noise fitting protocol changed')
    return p


def finalize():
    p=check();status=[]
    for dataset in ['flight','mc']:
        for variant in noise.VARIANTS:
            folder=OUT/dataset/variant
            if not (folder/'measurement_evaluation.json').exists():
                paths=[folder/('start_%d.json'%i) for i in range(2)]
                if not all(path.exists() for path in paths):
                    status.append(dict(dataset=dataset,variant=variant,complete=False,starts=sum(path.exists() for path in paths)));continue
                if (folder/'freeze.json').exists():raise ValueError('Incomplete existing noise freeze requires inspection')
                log=folder/'freeze.log'
                if log.exists():raise FileExistsError('Earlier noise freeze attempt requires inspection')
                with log.open('w') as stream:
                    result=subprocess.run([sys.executable,'analysis/research_conditional_noise.py','freeze','--dataset',dataset,'--variant',variant],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
                if result.returncode:raise RuntimeError('Noise freeze failed '+str(log))
            status.append(dict(dataset=dataset,variant=variant,complete=True))
    print('NOISE_FREEZE_STATUS',json.dumps(status),flush=True)
    if all(row['complete'] for row in status) and not (DEST/'all_frozen.json').exists():
        write(DEST/'all_frozen.json',dict(complete=True,created_utc=stamp(),models=status,
            freeze_sha256={str(path.relative_to(ROOT)):digest(path) for path in OUT.glob('*/*/freeze.json')}))


def compare():
    p=check()
    if (DEST/'measurement_comparison.json').exists():raise FileExistsError('Noise comparison already saved')
    completed=json.loads((DEST/'all_frozen.json').read_text())
    for name,fingerprint in completed['freeze_sha256'].items():
        if digest(ROOT/name)!=fingerprint:raise ValueError('Noise freeze changed')
    rows=[];selected={};rng=np.random.RandomState(p['measurement_seed'])
    for dataset in ['flight','mc']:
        predictions=[];models=[]
        for variant in noise.VARIANTS:
            folder=OUT/dataset/variant;ev=json.loads((folder/'measurement_evaluation.json').read_text())
            if not ev['complete'] or digest(folder/'predictions_before_reference.npz')!=ev['predictions_sha256']:
                raise ValueError('Noise held-out prediction changed')
            v=load(folder/'predictions_before_reference.npz');predictions.append(v);models.append(ev)
        for v in predictions[1:]:
            for key in ['run','event']:
                if not np.array_equal(v[key],predictions[0][key]):raise ValueError('Noise comparison identity differs')
        runs=predictions[0]['run'];unique,inverse=np.unique(runs,return_inverse=True)
        scores=np.column_stack([v['nll'] for v in predictions]);delta=scores[:,1:]-scores[:,:1]
        sums=np.column_stack([np.bincount(inverse,weights=delta[:,k]) for k in range(2)]);counts=np.bincount(inverse)
        bootstrap=np.empty((p['bootstrap_replicates'],2))
        for i in range(len(bootstrap)):
            index=rng.randint(0,len(unique),len(unique));bootstrap[i]=sums[index].sum(0)/counts[index].sum()
        selected[dataset]=noise.VARIANTS[int(np.argmin(scores.mean(0)))];comparisons=[]
        for k,variant in enumerate(noise.VARIANTS[1:]):
            comparisons.append(dict(variant=variant,mean_delta_nll=float(delta[:,k].mean()),
                paired_run_interval=np.percentile(bootstrap[:,k],p['ordinary_quantiles']),
                two_comparison_interval=np.percentile(bootstrap[:,k],p['two_comparison_quantiles'])))
        rows.append(dict(dataset=dataset,events=len(runs),runs=len(unique),models=models,selected=selected[dataset],comparisons=comparisons))
        np.savez_compressed(DEST/(dataset+'_measurement_bootstrap.npz'),delta_nll=bootstrap)
        print('NOISE_MEASURED_COMPARISON',dataset,selected[dataset],json.dumps(comparisons,default=lambda a:a.tolist()),flush=True)
    write(DEST/'measurement_comparison.json',dict(complete=True,created_utc=stamp(),selected=selected,rows=rows,
        protocol_sha256=digest(DEST/'protocol.json'),all_frozen_sha256=digest(DEST/'all_frozen.json'),reference_values_read=False))


def predict(dataset):
    p=check();selection=json.loads((DEST/'measurement_comparison.json').read_text());folder=DEST/dataset;folder.mkdir(exist_ok=True)
    if (folder/'predictions.npz').exists():raise FileExistsError('Noise development predictions already saved')
    pool=POOLS[dataset]
    for name,fingerprint in p['pools'][dataset].items():
        if digest(pool/name)!=fingerprint:raise ValueError('Development measurements changed')
    a=load(pool/'measured.npz');g=load(pool/'moments.npz');n=len(a['run'])
    for key in ['run','event']:
        if not np.array_equal(a[key],g[key]):raise ValueError('Development geometry identity differs')
    item=json.loads((OUT/'protocol.json').read_text())['datasets'][dataset]
    f,_=period.features(a,np.ones(n,bool),dataset,item['ADC_period_variant'],item['mean_normalization'])
    y=joint.observation(a,a['adc'],item['time_field']);x=g['x'];xc=x-x.mean(1)[:,None];w=xc/np.sum(xc*xc,axis=1)[:,None]
    if dataset=='flight':
        if digest(pool/'rich_geometry.npz')!=p['flight_forward_geometry_sha256']:raise ValueError('Forward geometry changed')
        rich_geometry=load(pool/'rich_geometry.npz')
    values={key:a[key] for key in ['run','event','R']};models={};response_checks={}
    for variant in noise.VARIANTS:
        source=OUT/dataset/variant/'freeze.json';fr=json.loads(source.read_text());theta=np.asarray(fr['model']['theta'])
        h,_=noise.features(a,np.ones(n,bool),dataset,variant,fr['noise_normalization']);rec=noise.reconstruct(theta,x,y,f,h)
        branches={'ols':(np.sum(w*a[item['time_field']],axis=1),w),'time':(rec['time_u'],rec['time_weights']),
            'joint':(rec['joint_u'],rec['weights'][:,:3].dot(timing.D))}
        for name in ['speed_null','speed_gain_null']:branches[name]=(rec[name]['u'],rec[name]['time_weights'])
        delta=np.linspace(-.03,.03,n);shifted=noise.reconstruct(theta,x,y+delta[:,None]*rec['direction'],f,h)
        error=float(np.max(abs(shifted['joint_u']-rec['joint_u']-delta)))
        if error>2e-11:raise ValueError('Conditional-noise full-signal response failed')
        response_checks[variant]=dict(full_signal_max_error=error,time_only_gain_percentiles=np.percentile(rec['time_injection_gain'],[0,5,50,95,100]))
        models[variant]=digest(source)
        for branch,(slope,weight) in branches.items():
            for mode in p['material_modes']:
                loss=center.solve_center(slope,a['R'],weight,g['moment'],theta[:2] if mode=='fit' else np.ones(2))
                if not loss['converged'].all():raise ValueError('Noise material inverse failed; preserve fixed cohort')
                value=loss['u']
                if dataset=='flight':
                    forward_value=forward.propagate(value,a['R'],rich_geometry)
                    if not forward_value['valid'].all():raise ValueError('Noise forward transport failed; preserve cohort')
                    value=forward_value['u']
                if not np.isfinite(value).all() or np.any(value==0):raise ValueError('Invalid noise output; preserve cohort')
                values[variant+'_'+branch+'_'+mode]=value
        print('NOISE_DEVELOPMENT_PREDICTED',dataset,variant,n,flush=True)
    np.savez_compressed(folder/'predictions.npz',**values)
    write(folder/'output_freeze.json',dict(complete=True,created_utc=stamp(),events=n,selected=selection['selected'][dataset],models=models,
        prediction_sha256=digest(folder/'predictions.npz'),protocol_sha256=digest(DEST/'protocol.json'),
        measurement_selection_sha256=digest(DEST/'measurement_comparison.json'),response_checks=response_checks,
        beta_references_read_by_predictor=False,all_events_retained=True))


def score(dataset):
    p=check();folder=DEST/dataset;fr=json.loads((folder/'output_freeze.json').read_text())
    if (folder/'evaluation.json').exists():raise FileExistsError('Noise beta score already saved')
    if digest(folder/'predictions.npz')!=fr['prediction_sha256'] or digest(DEST/'measurement_comparison.json')!=fr['measurement_selection_sha256']:
        raise ValueError('Noise selection or predictions changed')
    for variant,fingerprint in fr['models'].items():
        if digest(OUT/dataset/variant/'freeze.json')!=fingerprint:raise ValueError('Noise frozen model changed')
    a=load(folder/'predictions.npz');ref=load(POOLS[dataset]/'reference/reference.npz')
    for key in ['run','event']:
        if not np.array_equal(a[key],ref[key]):raise ValueError('Noise beta reference identity differs')
    valid=ref['valid'].astype(bool);truth=ref['beta'][valid] if dataset=='flight' else 1/ref['u'][valid];runs=a['run'][valid]
    names=sorted(key for key in a if key not in ['run','event','R']);values=np.column_stack([a[name][valid] for name in names])
    errors={'beta':1/values-truth[:,None],'u':values-1/truth[:,None]};points={space:matrix_metrics(error) for space,error in errors.items()}
    selected=fr['selected'];pairs=[(selected+'_joint_fit',selected+'_ols_fit'),(selected+'_joint_fit','constant_joint_fit'),
        (selected+'_speed_null_fit',selected+'_ols_fit'),(selected+'_speed_gain_null_fit',selected+'_ols_fit')]
    needed=sorted(set(name for pair in pairs for name in pair));indices=[names.index(name) for name in needed]
    groups=[np.flatnonzero(runs==run) for run in np.unique(runs)];rng=np.random.RandomState(p['beta_seed']+(dataset=='flight'))
    reduced={space:error[:,indices] for space,error in errors.items()};bootstrap={space:np.empty((p['bootstrap_replicates'],len(needed),len(METRICS))) for space in errors}
    for i in range(p['bootstrap_replicates']):
        pick=np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))])
        for space,error in reduced.items():bootstrap[space][i]=matrix_metrics(error[pick])
    result=dict(complete=True,created_utc=stamp(),dataset=dataset,events=len(runs),runs=len(groups),selected=selected,
        metrics={},comparisons={},bins=[],development_only=True,retraining_uncertainty_included=False,production_replacement=False,
        output_freeze_sha256=digest(folder/'output_freeze.json'),reference_sha256=digest(POOLS[dataset]/'reference/reference.npz'))
    for space,error in errors.items():
        result['metrics'][space]={name:dict(zip(METRICS,points[space][k])) for k,name in enumerate(names)};comparison=[]
        for candidate,baseline in pairs:
            difference=matrix_difference(points[space][names.index(baseline)],points[space][names.index(candidate)])
            draws=matrix_difference(bootstrap[space][:,needed.index(baseline)],bootstrap[space][:,needed.index(candidate)])
            comparison.append(dict(candidate=candidate,baseline=baseline,point=dict(zip(DELTA,difference)),
                paired_run_interval={key:np.percentile(draws[:,k],p['ordinary_quantiles']) for k,key in enumerate(DELTA)},
                two_primary_interval={key:np.percentile(draws[:,k],p['two_comparison_quantiles']) for k,key in enumerate(DELTA)}))
        result['comparisons'][space]=comparison
    for coordinate,edges in p['bins'].items():
        if coordinate=='year' and dataset=='mc':continue
        xx=a['R'][valid] if coordinate=='R' else (truth if coordinate=='reference_beta' else runs);boundaries=period.EDGES if coordinate=='year' else edges
        for k,(low,high) in enumerate(zip(boundaries[:-1],boundaries[1:])):
            mask=(xx>=low)&(xx<high)
            if mask.sum()<30:continue
            mm={space:matrix_metrics(error[mask]) for space,error in errors.items()}
            result['bins'].append(dict(coordinate=coordinate,low=edges[k],high=edges[k+1],events=int(mask.sum()),
                metrics={space:{name:dict(zip(METRICS,mm[space][j])) for j,name in enumerate(names)} for space in errors}))
    np.savez_compressed(folder/'bootstrap_metrics.npz',names=np.asarray(needed),**bootstrap);result['bootstrap_sha256']=digest(folder/'bootstrap_metrics.npz')
    write(folder/'evaluation.json',result)
    print('NOISE_DEVELOPMENT_BETA',dataset,selected,json.dumps(result['comparisons']['beta'],default=lambda a:a.tolist()),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['plan','finalize','compare','predict','score']);parser.add_argument('--dataset',choices=['flight','mc']);args=parser.parse_args()
    if args.stage in ['plan','finalize','compare']:globals()[args.stage]()
    elif args.dataset is None:parser.error('dataset required')
    else:globals()[args.stage](args.dataset)
