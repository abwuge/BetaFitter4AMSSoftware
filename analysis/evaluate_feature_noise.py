#!/usr/bin/env python3
"""Evaluate frozen geometry-noise outputs on the already audited MC identities."""
import argparse
import json
import numpy as np
import feature_conditional as feature
from research_conditional_mc import OUT as MCOUT,stamp
from research_flight_conditional import write,load,digest
from evaluate_conditional_mc import metrics,pair_delta

OUT=feature.OUT/'mc'
LABELS=['native_gaussianR','native_studentR5']


def prepare():
    if (OUT/'scoring_protocol.json').exists():raise FileExistsError('Feature scoring protocol exists')
    freezes={};outputs={}
    for label in LABELS:
        meta=json.loads((OUT/label/'measurement_evaluation.json').read_text())
        if digest(OUT/label/'freeze.json')!=meta['freeze_sha256'] or digest(OUT/label/'predictions_before_truth.npz')!=meta['predictions_sha256']:raise ValueError('Frozen feature output changed')
        freezes[label]=meta['freeze_sha256'];outputs[label]=meta['predictions_sha256']
    write(OUT/'scoring_protocol.json',dict(created_utc=stamp(),families=LABELS,freeze_sha256=freezes,output_sha256=outputs,
        reference_sha256=digest(MCOUT/'scoring/reference.npz'),selection='Exactly the established primary-verified common event set; no additional cuts',
        baseline='ols_u',secondary_baseline='same-family fitted constant-covariance free slope',
        bootstrap_replicates=1024,seed=920261011,bins=json.loads((MCOUT/'scoring/protocol.json').read_text())['bins'],
        model_selection='No selection by truth; report both predetermined rigidity-noise families',
        scope='MC development-pool closure, not blind flight RICH validation',frozen_parameter_retraining=False))
    print('FEATURE_SCORING_PROTOCOL',flush=True)


def score():
    protocol=json.loads((OUT/'scoring_protocol.json').read_text());reference=load(MCOUT/'scoring/reference.npz');old=load(MCOUT/'outputs_before_truth.npz')
    if digest(MCOUT/'scoring/reference.npz')!=protocol['reference_sha256']:raise ValueError('Reference changed')
    v=reference['valid'];truth=reference['u'][v];runs=old['run'][v];groups=[np.flatnonzero(runs==rr) for rr in np.unique(runs)]
    rng=np.random.RandomState(protocol['seed']);draws=[np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))]) for _ in range(protocol['bootstrap_replicates'])]
    results=[]
    for label in protocol['families']:
        if digest(OUT/label/'freeze.json')!=protocol['freeze_sha256'][label] or digest(OUT/label/'predictions_before_truth.npz')!=protocol['output_sha256'][label]:raise ValueError('Feature output changed')
        new=load(OUT/label/'predictions_before_truth.npz')
        if not np.array_equal(new['run'],old['run']) or not np.array_equal(new['event'],old['event']):raise ValueError('Feature identities differ')
        result=dict(family=label,events=int(v.sum()),runs=len(groups),spaces={},bins=[])
        for space in ['u','beta']:
            newerr=new['u'][v]-truth if space=='u' else 1/new['u'][v]-1/truth
            result['spaces'][space]=dict(metrics=metrics(newerr),comparisons={})
            for baseline,key in [('ols','ols_u'),('constant_covariance',label+'_fitted_u')]:
                olderr=old[key][v]-truth if space=='u' else 1/old[key][v]-1/truth
                bootstrap=[pair_delta(olderr[ix],newerr[ix]) for ix in draws]
                result['spaces'][space]['comparisons'][baseline]=dict(point=pair_delta(olderr,newerr),interval={k:np.percentile([r[k] for r in bootstrap],[2.5,50,97.5]) for k in bootstrap[0]})
            for field,edges in protocol['bins'].items():
                value=old['R'][v] if field=='measured_R' else 1/truth
                for low,high in zip(edges[:-1],edges[1:]):
                    sel=(value>=low)&(value<high)
                    if sel.sum()<20:continue
                    baseline=old['ols_u'][v]-truth if space=='u' else 1/old['ols_u'][v]-1/truth
                    result['bins'].append(dict(space=space,field=field,range=[low,high],events=int(sel.sum()),metrics=metrics(newerr[sel]),versus_ols=pair_delta(baseline[sel],newerr[sel])))
        results.append(result)
        print('FEATURE_BETA_SCORE',label,'u',result['spaces']['u'],'beta',result['spaces']['beta'],flush=True)
    write(OUT/'beta_evaluation.json',dict(complete=True,created_utc=stamp(),families=results,scoring_protocol_sha256=digest(OUT/'scoring_protocol.json'),
        uncertainty='Paired evaluation-run bootstrap of frozen feature models; whole-training-run refits are not included here',
        MC_truth_used_for_fitting=False,RICH_read=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:score()
