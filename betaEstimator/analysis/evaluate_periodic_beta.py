#!/usr/bin/env python3
"""Development beta checks after measurement-only period-model selection."""
import argparse
import json
from pathlib import Path
import numpy as np
import periodic_joint_adc as period
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import flight_conditional_likelihood as timing
import free_timing_transport as center
import continuous_rich_transport as forward
from research_periodic_adc import ROOT, OUT, CONFIGS
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference, METRICS, DELTA

CONT = OUT.parent
DEST = OUT / 'development_beta'
POOLS = {'flight':CONT/'prospective_flight', 'mc':CONT/'additional_mc'}
BRANCHES = ['ols','time','joint','speed_null','speed_gain_null']
SOURCES = ['evaluate_periodic_beta.py','periodic_joint_adc.py','joint_adc_likelihood.py','constrained_joint_adc.py',
    'flight_conditional_likelihood.py','free_timing_transport.py','flight_mean_loss.py','continuous_rich_transport.py',
    'evaluate_additional_mc_20261013.py']


def label(family,variant,branch,mode='fit'):
    return '_'.join([family,variant,branch,mode])


def plan():
    DEST.mkdir(exist_ok=True)
    if (DEST/'protocol.json').exists():
        raise FileExistsError('Development beta protocol already exists')
    write(DEST/'protocol.json',dict(created_utc=stamp(),source_sha256={name:digest(ROOT/'analysis'/name) for name in SOURCES},
        selection_protocol_sha256=digest(OUT/'measurement_comparison/protocol.json'),
        pools={dataset:{name:digest(folder/name) for name in ['measured.npz','moments.npz']} for dataset,folder in POOLS.items()},
        flight_forward_geometry_sha256=digest(POOLS['flight']/'rich_geometry.npz'),
        variants=period.VARIANTS,branches=BRANCHES,material_modes=['fit','unit'],
        selection='The measured held-out NLL choice is fixed before evaluating these beta errors; no variant chosen from beta performance.',
        primary_family={'flight':'newcalib_gaussianR','mc':'native_gaussianR'},
        primary_pairs=['measurement-selected joint versus its matched OLS','measurement-selected joint versus retrained static joint'],
        secondary_pairs=['measurement-selected speed_null versus matched OLS','measurement-selected speed_gain_null versus matched OLS'],
        bootstrap=dict(replicates=4096,seed=101320269,unit='paired whole evaluation run',quantiles=[2.5,50,97.5]),
        reference_scope='Previously opened flight RICH and additional-MC primary-L5 pools. All resulting beta performance is development evidence, not a new blind test.',
        reference_in_training=False,MC_coefficients_to_flight=False,old_predictions_overwritten=False,
        endpoint={'flight':'Nominal NaF emission plane, same geometry and continuous transport as the original prospective experiment', 'mc':'Audited primary L5 inverse speed'},
        coordinate_limitation='The original nominal/native emission-coordinate equality failure remains; prior micrometer endpoint sensitivity does not validate the material map.',
        beta_above_one='Retain all finite nonzero outputs and the pre-existing zero-loss continuation for u<=1.',
        coverage='Any invalid output fails the complete-cohort check; do not remove events according to the reconstructed output.',
        metrics=METRICS,differences=DELTA,minimum_bin_events=30,
        bins=dict(R=[2.5,3.,3.5,4.5,6.],reference_beta=[.78,.82,.86,.9,.94,.98,1.02],year=period.YEARS),
        uncertainty='Evaluation-run variation only; no new-model training-run refit uncertainty or reference/model systematics in these intervals.',
        production_replacement=False))
    print('PERIOD_DEVELOPMENT_PROTOCOL',flush=True)


def check():
    p=json.loads((DEST/'protocol.json').read_text())
    for name,fingerprint in p['source_sha256'].items():
        if digest(ROOT/'analysis'/name)!=fingerprint:raise ValueError('Development implementation changed '+name)
    if digest(OUT/'measurement_comparison/protocol.json')!=p['selection_protocol_sha256']:raise ValueError('Selection protocol changed')
    return p


def predict(dataset):
    p=check();folder=DEST/dataset;folder.mkdir(exist_ok=True)
    if (folder/'output_freeze.json').exists() or (folder/'predictions_before_scoring.npz').exists():
        raise FileExistsError('Development predictions already saved')
    selection=json.loads((OUT/'measurement_comparison/evaluation.json').read_text())
    if not selection['complete']:raise ValueError('Complete measurement selection first')
    pool=POOLS[dataset]
    for name,fingerprint in p['pools'][dataset].items():
        if digest(pool/name)!=fingerprint:raise ValueError('Development input changed')
    a=load(pool/'measured.npz');g=load(pool/'moments.npz');n=len(a['run'])
    for key in ['run','event']:
        if not np.array_equal(a[key],g[key]):raise ValueError('Development geometry identity differs')
    if dataset=='flight':
        if digest(pool/'rich_geometry.npz')!=p['flight_forward_geometry_sha256']:raise ValueError('Forward geometry changed')
        rg=load(pool/'rich_geometry.npz')
    outputs={key:a[key] for key in ['run','event','R']};diagnostics={};models={}
    for family,field,_ in CONFIGS[dataset]:
        for variant in period.VARIANTS:
            source=OUT/dataset/(family+'_'+variant)/'freeze.json';fr=json.loads(source.read_text())
            theta=np.asarray(fr['model']['theta']);features,_=period.features(a,np.ones(n,bool),dataset,variant,fr['normalization'])
            y=joint.observation(a,a['adc'],field);rec=joint.reconstruct(theta,g['x'],y,features)
            x=g['x'];xc=x-x.mean(1)[:,None];w=xc/np.sum(xc*x,axis=1)[:,None]
            branches={'ols':(np.sum(w*a[field],axis=1),w),'time':(rec['time_u'],rec['time_weights']),
                'joint':(rec['joint_u'],rec['weights'][:,:3].dot(timing.D))}
            response={}
            for branch in ['speed_null','speed_gain_null']:
                r=constrained.reconstruct(theta,x,y,features,branch)
                branches[branch]=(r['u'],r['time_weights'])
                response[branch]=float(np.max(abs(r['response']-np.r_[1.,np.zeros(r['response'].shape[1]-1)])))
                if response[branch]>2e-10:raise ValueError('Constrained response failed')
            b=np.linspace(-.03,.03,n)
            shifted=joint.reconstruct(theta,x,y+b[:,None]*rec['direction'],features)
            response['joint_full_signal_max_error']=float(np.max(abs(shifted['joint_u']-rec['joint_u']-b)))
            if response['joint_full_signal_max_error']>2e-11:raise ValueError('Joint response failed')
            response['joint_time_only_gain_percentiles']=np.percentile(rec['time_injection_gain'],[0,5,50,95,100])
            models[family+'_'+variant]=digest(source);diagnostics[family+'_'+variant]=response
            for branch,(slope,weights) in branches.items():
                if not np.isfinite(slope).all() or np.any(slope==0):raise ValueError('Invalid slope; preserve cohort')
                for mode in ['fit','unit']:
                    loss=center.solve_center(slope,a['R'],weights,g['moment'],theta[:2] if mode=='fit' else np.ones(2))
                    if not loss['converged'].all():raise ValueError('Mean-loss inverse failed; preserve cohort')
                    value=loss['u'];name=label(family,variant,branch,mode)
                    if dataset=='flight':
                        transported=forward.propagate(value,a['R'],rg)
                        if not transported['valid'].all():raise ValueError('RICH propagation failed; preserve cohort')
                        value=transported['u']
                    if not np.isfinite(value).all() or np.any(value==0):raise ValueError('Invalid output; preserve cohort')
                    outputs[name]=value
            print('PERIOD_DEVELOPMENT_PREDICTED',dataset,family,variant,n,flush=True)
    np.savez_compressed(folder/'predictions_before_scoring.npz',**outputs)
    write(folder/'output_freeze.json',dict(complete=True,created_utc=stamp(),events=n,
        model_sha256=models,protocol_sha256=digest(DEST/'protocol.json'),
        measurement_selection_sha256=digest(OUT/'measurement_comparison/evaluation.json'),
        selected_variants=selection['selected'][dataset],prediction_sha256=digest(folder/'predictions_before_scoring.npz'),
        response_checks=diagnostics,reference_values_read_by_prediction=False,all_events_retained=True,
        reference_pool_was_previously_exposed=True))


def score(dataset):
    p=check();folder=DEST/dataset;fr=json.loads((folder/'output_freeze.json').read_text())
    if (folder/'evaluation.json').exists():raise FileExistsError('Development scoring already saved')
    if digest(folder/'predictions_before_scoring.npz')!=fr['prediction_sha256']:raise ValueError('Frozen predictions changed')
    if digest(OUT/'measurement_comparison/evaluation.json')!=fr['measurement_selection_sha256']:raise ValueError('Measurement selection changed')
    for name,fingerprint in fr['model_sha256'].items():
        if digest(OUT/dataset/name/'freeze.json')!=fingerprint:raise ValueError('Model changed')
    pool=POOLS[dataset];pred=load(folder/'predictions_before_scoring.npz');ref=load(pool/'reference/reference.npz')
    for key in ['run','event']:
        if not np.array_equal(pred[key],ref[key]):raise ValueError('Reference identity differs')
    cohort=ref['valid'].astype(bool);truth=ref['beta'][cohort] if dataset=='flight' else 1/ref['u'][cohort]
    names=sorted(name for name in pred if name not in ['run','event','R'])
    values=np.column_stack([pred[name][cohort] for name in names]);runs=pred['run'][cohort]
    errors={'beta':1/values-truth[:,None],'u':values-1/truth[:,None]}
    points={space:matrix_metrics(error) for space,error in errors.items()}
    family=p['primary_family'][dataset];chosen=fr['selected_variants'][family]
    pair_names=[(label(family,chosen,'joint'),label(family,chosen,'ols')),
        (label(family,chosen,'joint'),label(family,'static','joint')),
        (label(family,chosen,'speed_null'),label(family,chosen,'ols')),
        (label(family,chosen,'speed_gain_null'),label(family,chosen,'ols'))]
    required=sorted(set(name for pair in pair_names for name in pair));indices=[names.index(name) for name in required]
    groups=[np.flatnonzero(runs==run) for run in np.unique(runs)]
    rng=np.random.RandomState(p['bootstrap']['seed']+(1 if dataset=='flight' else 0))
    boot={space:np.empty((p['bootstrap']['replicates'],len(required),len(METRICS))) for space in errors}
    for i in range(p['bootstrap']['replicates']):
        pick=np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))])
        for space,error in errors.items():boot[space][i]=matrix_metrics(error[pick][:,indices])
        if (i+1)%1024==0:print('PERIOD_DEVELOPMENT_BOOTSTRAP',dataset,i+1,flush=True)
    result=dict(complete=True,created_utc=stamp(),dataset=dataset,events=len(runs),runs=len(groups),
        selected=chosen,primary_family=family,reference_sha256=digest(pool/'reference/reference.npz'),
        output_freeze_sha256=digest(folder/'output_freeze.json'),metrics={},comparisons={},bins=[],
        development_only=True,training_refit_uncertainty_included=False,production_replacement=False)
    for space,error in errors.items():
        result['metrics'][space]={name:dict(zip(METRICS,points[space][j])) for j,name in enumerate(names)}
        pairs=[]
        for candidate,baseline in pair_names:
            delta=matrix_difference(points[space][names.index(baseline)],points[space][names.index(candidate)])
            sample=matrix_difference(boot[space][:,required.index(baseline)],boot[space][:,required.index(candidate)])
            pairs.append(dict(candidate=candidate,baseline=baseline,point=dict(zip(DELTA,delta)),
                paired_run_interval={key:np.percentile(sample[:,k],p['bootstrap']['quantiles']) for k,key in enumerate(DELTA)}))
        result['comparisons'][space]=pairs
    for coordinate,edges in p['bins'].items():
        if coordinate=='year' and dataset=='mc':continue
        values_x=pred['R'][cohort] if coordinate=='R' else (truth if coordinate=='reference_beta' else runs)
        raw_edges=period.EDGES if coordinate=='year' else edges
        for index,(low,high) in enumerate(zip(raw_edges[:-1],raw_edges[1:])):
            mask=(values_x>=low)&(values_x<high)
            if mask.sum()<p['minimum_bin_events']:continue
            metrics={space:matrix_metrics(error[mask]) for space,error in errors.items()}
            result['bins'].append(dict(coordinate=coordinate,low=edges[index],high=edges[index+1],events=int(mask.sum()),
                metrics={space:{name:dict(zip(METRICS,metrics[space][j])) for j,name in enumerate(names)} for space in errors}))
    np.savez_compressed(folder/'bootstrap_metrics.npz',names=np.asarray(required),**boot)
    result['bootstrap_sha256']=digest(folder/'bootstrap_metrics.npz')
    write(folder/'evaluation.json',result)
    print('PERIOD_DEVELOPMENT_BETA',dataset,chosen,json.dumps(result['comparisons']['beta'],default=lambda a:a.tolist()),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['plan','predict','score']);parser.add_argument('--dataset',choices=['flight','mc'])
    args=parser.parse_args()
    if args.stage=='plan':plan()
    elif args.dataset is None:parser.error('dataset is required')
    else:globals()[args.stage](args.dataset)
