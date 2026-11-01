#!/usr/bin/env python3
"""Reference scoring kept separate from measured-only amplitude-noise fitting."""
import argparse, json
from pathlib import Path
import numpy as np
from research_amplitude_noise_20261013 import OUT,MCOUT,ROOT,inputs,CONFIGS
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
from evaluate_conditional_mc import metrics,pair_delta

SCORE=OUT/'scoring'
LABELS=[name+'_'+variant for name in ['native_gaussianR','native_studentR5'] for variant in ['geometry','readout']]


def identity(a):return (a['run'].astype(np.uint64)<<np.uint64(32))|a['event'].astype(np.uint64)


def join(left,right):
    a=identity(left);b=identity(right)
    if len(np.unique(a))!=len(a) or len(np.unique(b))!=len(b):raise ValueError('Duplicate identity')
    order=np.argsort(b);pos=np.searchsorted(b[order],a)
    if np.any(pos>=len(b)):raise ValueError('Missing identity')
    ix=order[pos]
    if not np.array_equal(a,b[ix]):raise ValueError('Missing identity')
    return ix


def prepare():
    SCORE.mkdir(exist_ok=True)
    if (SCORE/'protocol.json').exists():raise FileExistsError('Scoring protocol already declared')
    write(SCORE/'protocol.json',dict(created_utc=stamp(),candidates=LABELS,
        common_events='Intersection of the declared positive-ADC measured domain with existing primary-L5 audit; exact run/event association, no new reference-error cuts',
        baselines=['same-event OLS','stored native tof_betah','matched geometry-only model','constant beta=1 negative control'],
        metrics=['bias','median','W68=(q84-q16)/2','RMSE','absolute error tail 0.02 and 0.05'],spaces=['u','beta'],
        comparisons='Report all four candidates relative to OLS/native beta; amplitude versus same-family common-training-domain geometry',
        evaluation_bootstrap_replicates=2048,seed=1320261013,
        uncertainty='Paired whole evaluation-run bootstrap of frozen predictions; training-run refits reported separately',
        percentile_levels=[2.5,50,97.5],four_comparison_bonferroni_levels=[.625,50,99.375],
        bins=dict(measured_R=[2.5,3,3.5,4.5,6],reference_beta=[0,.8,.85,.9,.95,1.01]),
        selection='No estimator or hyperparameter is chosen using this reference evaluation; previous MC source pool is a development pool',
        reference_sha256=digest(MCOUT/'scoring/reference.npz'),input_protocol_sha256=digest(OUT/'mc/protocol.json'),
        script_sha256=digest(Path(__file__)),flight_RICH_read=False))
    print('AMPLITUDE_SCORING_PROTOCOL',flush=True)


def freeze_gate():
    protocol=json.loads((SCORE/'protocol.json').read_text())
    if digest(OUT/'mc/protocol.json')!=protocol['input_protocol_sha256']:raise ValueError('Measured protocol changed')
    if digest(Path(__file__))!=protocol['script_sha256']:raise ValueError('Scoring source changed')
    freeze={};pred={}
    for label in protocol['candidates']:
        directory=OUT/'mc'/label;meta=json.loads((directory/'measurement_evaluation.json').read_text())
        if not meta['complete'] or digest(directory/'freeze.json')!=meta['freeze_sha256'] or digest(directory/'predictions_before_truth.npz')!=meta['output_sha256']:
            raise ValueError('Candidate output changed')
        freeze[label]=dict(model=meta['freeze_sha256'],output=meta['output_sha256'])
        pred[label]=load(directory/'predictions_before_truth.npz')
    path=SCORE/'output_gate.json'
    record=dict(created_utc=stamp(),models=freeze,protocol_sha256=digest(SCORE/'protocol.json'),reference_opened_in_gate=False)
    if path.exists():
        old=json.loads(path.read_text())
        if old['models']!=freeze or old['protocol_sha256']!=record['protocol_sha256']:raise ValueError('Frozen scoring inputs changed')
    else:write(path,record)
    return protocol,pred


def score():
    protocol,pred=freeze_gate()
    if (SCORE/'evaluation.json').exists():raise FileExistsError('Scoring already complete')
    if digest(MCOUT/'scoring/reference.npz')!=protocol['reference_sha256']:raise ValueError('Reference changed')
    reference=load(MCOUT/'scoring/reference.npz');base=load(MCOUT/'outputs_before_truth.npz');a=load(MCOUT/'measured.npz')
    first=pred[protocol['candidates'][0]]
    for label,b in pred.items():
        if not np.array_equal(identity(b),identity(first)):raise ValueError('Candidate common events differ')
    ir=join(first,reference);ib=join(first,base);im=join(first,a)
    valid=reference['valid'][ir];truth=reference['u'][ir][valid];runs=first['run'][valid];n=len(truth)
    native=a['tof_betah'][im][valid]
    if not np.isfinite(native).all() or np.any(native==0):raise ValueError('Invalid native baseline; do not silently drop events')
    branches={'ols':base['ols_u'][ib][valid],'native':1/native,'constant':np.ones(n)}
    branches.update({label:b['u'][valid] for label,b in pred.items()})
    if not all(np.isfinite(v).all() and np.all(v!=0) for v in branches.values()):raise ValueError('Invalid prediction on common domain')
    unique=np.unique(runs);groups=[np.flatnonzero(runs==r) for r in unique];rng=np.random.RandomState(protocol['seed'])
    draws=[np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))]) for _ in range(protocol['evaluation_bootstrap_replicates'])]
    comparisons=[(label,baseline) for label in protocol['candidates'] for baseline in ['ols','native']]
    comparisons += [(label,label.replace('_readout','_geometry')) for label in protocol['candidates'] if label.endswith('_readout')]
    table={};bootstrap_values={};bins=[]
    for space in protocol['spaces']:
        errors={k:(v-truth if space=='u' else 1/v-1/truth) for k,v in branches.items()}
        table[space]=dict(metrics={k:metrics(e) for k,e in errors.items()},comparisons={})
        for label,baseline in comparisons:
            key=label+'_versus_'+baseline;rep=[pair_delta(errors[baseline][ix],errors[label][ix]) for ix in draws]
            record=dict(point=pair_delta(errors[baseline],errors[label]),interval={k:np.percentile([x[k] for x in rep],protocol['percentile_levels']) for k in rep[0]},
                four_comparison_bonferroni_interval={k:np.percentile([x[k] for x in rep],protocol['four_comparison_bonferroni_levels']) for k in rep[0]})
            table[space]['comparisons'][key]=record
            for k in rep[0]:bootstrap_values[space+'_'+key+'_'+k]=np.array([x[k] for x in rep])
        for field,edges in protocol['bins'].items():
            value=first['R'][valid] if field=='measured_R' else 1/truth
            for low,high in zip(edges[:-1],edges[1:]):
                select=(value>=low)&(value<high)
                if select.sum()<20:continue
                bins.append(dict(space=space,field=field,range=[low,high],events=int(select.sum()),
                    metrics={k:metrics(e[select]) for k,e in errors.items()},
                    comparisons={label+'_versus_'+old:pair_delta(errors[old][select],errors[label][select]) for label,old in comparisons}))
    np.savez_compressed(SCORE/'bootstrap_metrics.npz',**bootstrap_values)
    result=dict(complete=True,created_utc=stamp(),events=n,runs=len(unique),prediction_domain_events=len(first['run']),
        reference_selection_losses=len(first['run'])-n,table=table,bins=bins,
        native_vs_ols=pair_delta(branches['ols']-truth,branches['native']-truth),
        beta_above_one_retained={k:int(np.sum(v<1)) for k,v in branches.items()},
        output_gate_sha256=digest(SCORE/'output_gate.json'),reference_sha256=protocol['reference_sha256'],
        bootstrap_metrics_sha256=digest(SCORE/'bootstrap_metrics.npz'),scope='C12 primary L5 reference in an existing MC development pool; not isotope performance or new flight RICH validation',
        truth_used_for_fit=False,RICH_read=False)
    write(SCORE/'evaluation.json',result)
    for label in protocol['candidates']:
        for space in ['u','beta']:
            for old in ['ols','native']:
                r=table[space]['comparisons'][label+'_versus_'+old]
                print('AMPLITUDE_BETA',label,space,old,'gain',r['point']['w68_improvement'],'CI',r['interval']['w68_improvement'],'RMSE',r['point']['rmse_improvement'],flush=True)


def measured():
    path=SCORE/'measurement_comparison.json'
    if path.exists():raise FileExistsError('Measurement comparison already exists')
    rows=[]
    for dataset in ['mc','flight']:
        for name,field,nu in CONFIGS[dataset]:
            newdir=OUT/dataset/(name+'_readout');p=load(newdir/'predictions_before_truth.npz')
            if dataset=='mc':oldp=load(OUT/dataset/(name+'_geometry')/'predictions_before_truth.npz')
            else:oldp=load(ROOT/'results/flight_calibration_research/feature_noise/flight'/name/'predictions_before_truth.npz')
            ix=join(p,oldp);old=oldp['nll'][ix] if dataset=='mc' else oldp['new_nll'][ix]
            delta=p['nll']-old;run=p['run'];groups=[np.flatnonzero(run==r) for r in np.unique(run)];rng=np.random.RandomState(2620261013)
            means=[np.concatenate([delta[groups[j]] for j in rng.randint(0,len(groups),len(groups))]).mean() for _ in range(2048)]
            row=dict(dataset=dataset,family=name,events=len(delta),delta_nll=float(delta.mean()),paired_run_interval=np.percentile(means,[2.5,50,97.5]),
                comparison='Readout versus geometry covariance on exactly matched evaluation measurements',output_sha256=digest(newdir/'predictions_before_truth.npz'),
                RICH_read=False,MC_truth_read=False)
            rows.append(row);print('AMPLITUDE_MEASUREMENT',dataset,name,row['delta_nll'],row['paired_run_interval'],flush=True)
    write(path,dict(complete=True,created_utc=stamp(),comparisons=rows))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score','measured']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='score':score()
    else:measured()
