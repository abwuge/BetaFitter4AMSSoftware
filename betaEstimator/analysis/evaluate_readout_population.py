#!/usr/bin/env python3
"""All prespecified population ablations, scored after complete freezing."""
import argparse,json
from pathlib import Path
import numpy as np
from research_readout_population import OUT,AMPOUT,ROOT,CONFIGS
from research_conditional_mc import OUT as MCOUT,stamp
from research_flight_conditional import load,write,digest
from evaluate_amplitude_noise_20261013 import identity,join
from evaluate_conditional_mc import metrics,pair_delta

SCORE=OUT/'scoring'
CANDIDATES=[name+'_'+v for name in ['native_gaussianR','native_studentR5'] for v in ['population_only','joint']]


def prepare():
    SCORE.mkdir(exist_ok=True)
    if (SCORE/'protocol.json').exists():raise FileExistsError('Population scoring protocol exists')
    write(SCORE/'protocol.json',dict(created_utc=stamp(),candidates=CANDIDATES,reference_sha256=digest(MCOUT/'scoring/reference.npz'),
        measurement_protocol_sha256=digest(OUT/'mc/protocol.json'),script_sha256=digest(Path(__file__)),
        event_domain='Same ADC-valid events and original primary L5 audit as amplitude ablation; no new residual or reference cuts',
        metrics=['bias','median','W68','RMSE','tail02','tail05'],spaces=['u','beta'],replicates=2048,seed=170202610,
        comparisons='Each population-only/joint candidate against native beta and same-family geometry/noise-only controls; joint also against population-only',
        percentile_levels=[2.5,50,97.5],four_candidate_levels=[.625,50,99.375],
        uncertainty='Paired frozen-output evaluation-run sampling only; refit/structure error not included',
        bins=dict(R=[2.5,3,3.5,4.5,6],beta=[0,.8,.85,.9,.95,1.01]),
        model_selection='None using these metrics; report all prespecified outcomes, including worsening',RICH_read=False,
        scope='Sequential mechanistic diagnosis on a previously examined MC development pool, not independent discovery significance'))
    print('POPULATION_SCORING_PROTOCOL',flush=True)


def gather():
    p=json.loads((SCORE/'protocol.json').read_text())
    if p['measurement_protocol_sha256']!=digest(OUT/'mc/protocol.json') or p['script_sha256']!=digest(Path(__file__)):raise ValueError('Scoring source or protocol changed')
    output={};frozen={}
    for candidate in p['candidates']:
        d=OUT/'mc'/candidate;meta=json.loads((d/'measurement_evaluation.json').read_text())
        if not meta['complete'] or digest(d/'freeze.json')!=meta['freeze_sha256'] or digest(d/'predictions_before_truth.npz')!=meta['output_sha256']:raise ValueError('Unfrozen candidate '+candidate)
        output[candidate]=load(d/'predictions_before_truth.npz');frozen[candidate]=dict(model=meta['freeze_sha256'],output=meta['output_sha256'])
    gate=dict(created_utc=stamp(),models=frozen,scoring_protocol_sha256=digest(SCORE/'protocol.json'),reference_read=False)
    path=SCORE/'output_gate.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['models']!=frozen or old['scoring_protocol_sha256']!=gate['scoring_protocol_sha256']:raise ValueError('Candidate gate changed')
    else:write(path,gate)
    return p,output


def score():
    p,pred=gather()
    if (SCORE/'evaluation.json').exists():raise FileExistsError('Population scored already')
    if digest(MCOUT/'scoring/reference.npz')!=p['reference_sha256']:raise ValueError('Reference changed')
    first=pred[p['candidates'][0]];ref=load(MCOUT/'scoring/reference.npz');meas=load(MCOUT/'measured.npz');base=load(MCOUT/'outputs_before_truth.npz')
    ir=join(first,ref);im=join(first,meas);ib=join(first,base);valid=ref['valid'][ir];truth=ref['u'][ir][valid];runs=first['run'][valid]
    values={'ols':base['ols_u'][ib][valid],'native':1/meas['tof_betah'][im][valid],'constant':np.ones(len(truth))}
    for label,q in pred.items():
        if not np.array_equal(identity(first),identity(q)):raise ValueError('Candidate events differ')
        values[label]=q['u'][valid]
    for family in ['native_gaussianR','native_studentR5']:
        for variant in ['geometry','readout']:
            q=load(AMPOUT/'mc'/(family+'_'+variant)/'predictions_before_truth.npz');ix=join(first,q);values[family+'_'+variant]=q['u'][ix][valid]
    if not all(np.isfinite(x).all() and (x!=0).all() for x in values.values()):raise ValueError('Invalid prediction without allowed filtering')
    comparisons=[]
    for label in p['candidates']:
        family=label.split('_population_only')[0] if label.endswith('_population_only') else label.rsplit('_joint',1)[0]
        comparisons.extend([(label,'native'),(label,family+'_geometry'),(label,family+'_readout')])
        if label.endswith('_joint'):comparisons.append((label,family+'_population_only'))
    group=[np.flatnonzero(runs==r) for r in np.unique(runs)];rng=np.random.RandomState(p['seed'])
    draws=[np.concatenate([group[j] for j in rng.randint(0,len(group),len(group))]) for _ in range(p['replicates'])]
    result={};bins=[]
    for space in p['spaces']:
        errors={k:v-truth if space=='u' else 1/v-1/truth for k,v in values.items()};rows={}
        for label,baseline in comparisons:
            sample=[pair_delta(errors[baseline][ix],errors[label][ix]) for ix in draws]
            rows[label+'_versus_'+baseline]=dict(point=pair_delta(errors[baseline],errors[label]),
                interval={k:np.percentile([v[k] for v in sample],p['percentile_levels']) for k in sample[0]},
                four_candidate_interval={k:np.percentile([v[k] for v in sample],p['four_candidate_levels']) for k in sample[0]})
        result[space]=dict(metrics={k:metrics(v) for k,v in errors.items()},comparisons=rows)
        for key,edges in p['bins'].items():
            value=first['R'][valid] if key=='R' else 1/truth
            for low,high in zip(edges[:-1],edges[1:]):
                select=(value>=low)&(value<high)
                if select.sum()<20:continue
                bins.append(dict(space=space,key=key,range=[low,high],events=int(select.sum()),metrics={k:metrics(v[select]) for k,v in errors.items()},
                    comparisons={label+'_versus_'+old:pair_delta(errors[old][select],errors[label][select]) for label,old in comparisons}))
    write(SCORE/'evaluation.json',dict(complete=True,created_utc=stamp(),events=len(truth),runs=len(group),table=result,bins=bins,output_gate_sha256=digest(SCORE/'output_gate.json'),
        reference_sha256=p['reference_sha256'],MC_truth_fit=False,RICH_read=False))
    for space in ['u','beta']:
        for label in p['candidates']:
            row=result[space]['comparisons'][label+'_versus_native']
            print('POPULATION_BETA',space,label,row['point']['w68_improvement'],row['interval']['w68_improvement'],'RMSE',row['point']['rmse_improvement'],flush=True)


def measured():
    path=SCORE/'measurement_comparison.json'
    if path.exists():raise FileExistsError('Population measurement score exists')
    rows=[]
    for dataset in ['mc','flight']:
        for name,field,nu in CONFIGS[dataset]:
            old=load(AMPOUT/dataset/(name+'_readout')/'predictions_before_truth.npz')
            for variant in (['population_only','joint'] if dataset=='mc' else ['joint']):
                folder=OUT/dataset/(name+'_'+variant);new=load(folder/'predictions_before_truth.npz');ix=join(new,old);delta=new['nll']-old['nll'][ix]
                groups=[np.flatnonzero(new['run']==r) for r in np.unique(new['run'])];rng=np.random.RandomState(3320261013)
                means=[np.concatenate([delta[groups[j]] for j in rng.randint(0,len(groups),len(groups))]).mean() for _ in range(2048)]
                row=dict(dataset=dataset,family=name,variant=variant,events=len(delta),delta_nll=float(delta.mean()),interval=np.percentile(means,[2.5,50,97.5]),
                    quadrature=json.loads((folder/'measurement_evaluation.json').read_text()))
                rows.append(row);print('POPULATION_MEASUREMENT',dataset,name,variant,row['delta_nll'],row['interval'],flush=True)
    write(path,dict(complete=True,created_utc=stamp(),comparisons=rows,RICH_read=False,MC_truth_read=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score','measured']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='score':score()
    else:measured()
