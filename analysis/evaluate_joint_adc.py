#!/usr/bin/env python3
"""Evaluate frozen direct-ADC / timing reconstruction on the established L5 audit."""
import argparse,json
from pathlib import Path
import numpy as np
from research_joint_adc import OUT,ROOT,CONFIGS
from research_conditional_mc import OUT as MCOUT,stamp
from research_flight_conditional import load,write,digest
from evaluate_amplitude_noise_20261013 import identity,join
from evaluate_conditional_mc import metrics,pair_delta

SCORE=OUT/'scoring'
FAMILIES=['native_gaussianR','native_studentR5']
BRANCHES=['time_u','adc_u','joint_u','time_transport_u','joint_transport_u']


def prepare():
    SCORE.mkdir(exist_ok=True)
    if (SCORE/'protocol.json').exists():raise FileExistsError('Joint scoring protocol exists')
    write(SCORE/'protocol.json',dict(created_utc=stamp(),families=FAMILIES,branches=BRANCHES,
        primary='native_studentR5_joint_transport_u versus stored native tof_betah; fixed before opening these output errors, not selected by MC score',
        primary_caveat='This is a sequential development experiment on an already examined MC pool, not a prospective discovery confidence claim',
        reference_sha256=digest(MCOUT/'scoring/reference.npz'),measurement_protocol_sha256=digest(OUT/'mc/protocol.json'),script_sha256=digest(Path(__file__)),
        domain='Exactly the same positive-ADC domain and primary-L5 validity mask as previous amplitude and population experiments',
        metrics=['bias','median','W68','RMSE','tail02','tail05'],spaces=['u','beta'],replicates=2048,seed=229202610,
        comparisons='All ten outputs versus native baseline; each joint output also versus its own matching time-only marginal; report ADC-only as an intentionally noisy independent-measurement control',
        percentile_levels=[2.5,50,97.5],ten_candidate_levels=[.25,50,99.75],
        uncertainty='Paired frozen-estimator whole evaluation-run resampling; does not include retraining or response-model misspecification',
        bins=dict(R=[2.5,3,3.5,4.5,6],beta=[0,.8,.85,.9,.95,1.01]),
        invalid_prediction_policy='Do not silently reject failed, negative or beta>1 reconstructed outputs; finite nonzero u required for beta metrics, fail scoring otherwise',
        RICH_read=False,MC_truth_used_to_train=False))
    print('JOINT_ADC_SCORING_PROTOCOL',flush=True)


def gather():
    p=json.loads((SCORE/'protocol.json').read_text())
    if digest(Path(__file__))!=p['script_sha256'] or digest(OUT/'mc/protocol.json')!=p['measurement_protocol_sha256']:raise ValueError('Joint score source/protocol changed')
    pred={};frozen={}
    for family in p['families']:
        directory=OUT/'mc'/family;meta=json.loads((directory/'measurement_evaluation.json').read_text())
        if not meta['complete'] or digest(directory/'freeze.json')!=meta['freeze_sha256'] or digest(directory/'predictions_before_truth.npz')!=meta['output_sha256']:raise ValueError('Unfrozen joint model '+family)
        pred[family]=load(directory/'predictions_before_truth.npz');frozen[family]=dict(model=meta['freeze_sha256'],output=meta['output_sha256'])
    gate=dict(created_utc=stamp(),models=frozen,protocol_sha256=digest(SCORE/'protocol.json'),reference_read=False)
    path=SCORE/'output_gate.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['models']!=frozen or old['protocol_sha256']!=gate['protocol_sha256']:raise ValueError('Joint output gate changed')
    else:write(path,gate)
    return p,pred


def score():
    p,pred=gather()
    if (SCORE/'evaluation.json').exists():raise FileExistsError('Joint score already complete')
    if digest(MCOUT/'scoring/reference.npz')!=p['reference_sha256']:raise ValueError('Joint reference changed')
    first=pred[p['families'][0]];reference=load(MCOUT/'scoring/reference.npz');a=load(MCOUT/'measured.npz');base=load(MCOUT/'outputs_before_truth.npz')
    ir=join(first,reference);im=join(first,a);ib=join(first,base);valid=reference['valid'][ir];truth=reference['u'][ir][valid];runs=first['run'][valid]
    values={'native':1/a['tof_betah'][im][valid],'ols':base['ols_u'][ib][valid],'constant':np.ones(len(truth))}
    for family,q in pred.items():
        if not np.array_equal(identity(first),identity(q)):raise ValueError('Joint prediction identities differ')
        for branch in p['branches']:values[family+'_'+branch]=q[branch][valid]
    if not all(np.isfinite(v).all() and (v!=0).all() for v in values.values()):raise ValueError('Nonfinite or zero joint prediction; no outcome filtering allowed')
    comparisons=[(family+'_'+branch,'native') for family in p['families'] for branch in p['branches']]
    for family in p['families']:
        comparisons.extend([(family+'_joint_u',family+'_time_u'),(family+'_joint_transport_u',family+'_time_transport_u')])
    groups=[np.flatnonzero(runs==r) for r in np.unique(runs)];rng=np.random.RandomState(p['seed'])
    draws=[np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))]) for _ in range(p['replicates'])]
    result={};bins=[]
    for space in p['spaces']:
        errors={k:v-truth if space=='u' else 1/v-1/truth for k,v in values.items()};rows={}
        for label,old in comparisons:
            rep=[pair_delta(errors[old][ix],errors[label][ix]) for ix in draws]
            rows[label+'_versus_'+old]=dict(point=pair_delta(errors[old],errors[label]),
                interval={k:np.percentile([r[k] for r in rep],p['percentile_levels']) for k in rep[0]},
                ten_candidate_interval={k:np.percentile([r[k] for r in rep],p['ten_candidate_levels']) for k in rep[0]})
        result[space]=dict(metrics={k:metrics(v) for k,v in errors.items()},comparisons=rows)
        for key,edges in p['bins'].items():
            coordinate=first['R'][valid] if key=='R' else 1/truth
            for low,high in zip(edges[:-1],edges[1:]):
                take=(coordinate>=low)&(coordinate<high)
                if take.sum()<20:continue
                bins.append(dict(space=space,key=key,range=[low,high],events=int(take.sum()),metrics={k:metrics(v[take]) for k,v in errors.items()},
                    comparisons={label+'_versus_'+old:pair_delta(errors[old][take],errors[label][take]) for label,old in comparisons}))
    write(SCORE/'evaluation.json',dict(complete=True,created_utc=stamp(),events=len(truth),runs=len(groups),table=result,bins=bins,
        output_gate_sha256=digest(SCORE/'output_gate.json'),reference_sha256=p['reference_sha256'],
        predictions_beta_above_one={k:int((v<1).sum()) for k,v in values.items()},RICH_read=False,MC_truth_fit=False))
    for space in ['u','beta']:
        for family in p['families']:
            for branch in p['branches']:
                label=family+'_'+branch;row=result[space]['comparisons'][label+'_versus_native']
                print('JOINT_ADC_BETA',space,label,'widthgain',row['point']['w68_improvement'],'CI',row['interval']['w68_improvement'],'RMSEgain',row['point']['rmse_improvement'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:score()
