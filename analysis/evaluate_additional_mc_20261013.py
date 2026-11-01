#!/usr/bin/env python3
"""Vectorized common-event evaluation on the separately reserved MC runs."""
import argparse,json
from pathlib import Path
import numpy as np
from research_additional_mc_20261013 import OUT
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
from evaluate_conditional_mc import metrics,pair_delta
from reference_additional_mc_20261013 import gate

SCORE=OUT/'scoring'
METRICS=['bias','median','w68','rmse','tail02','tail05']
DELTA=['w68_improvement','rmse_improvement','bias_change','tail02_change','tail05_change']


def matrix_metrics(error):
    error=np.asarray(error,float)
    if error.ndim!=2 or not np.isfinite(error).all():raise ValueError('Finite matrix of errors required')
    q=np.percentile(error,[16,50,84],axis=0)
    return np.column_stack([error.mean(0),q[1],.5*(q[2]-q[0]),np.sqrt(np.mean(error*error,axis=0)),np.mean(abs(error)>.02,axis=0),np.mean(abs(error)>.05,axis=0)])


def matrix_difference(old,new):
    return np.stack([1-new[...,2]/old[...,2],1-new[...,3]/old[...,3],new[...,0]-old[...,0],new[...,4]-old[...,4],new[...,5]-old[...,5]],axis=-1)


def prepare():
    SCORE.mkdir(exist_ok=True)
    if (SCORE/'protocol.json').exists():raise FileExistsError('Additional score protocol already declared')
    original=json.loads((OUT/'protocol.json').read_text())
    write(SCORE/'protocol.json',dict(created_utc=stamp(),original_protocol_sha256=digest(OUT/'protocol.json'),script_sha256=digest(Path(__file__)),
        primary=original['evaluation_primary'],reference_rules=original['reference_conditions'],
        metrics=METRICS,spaces=original['evaluation_spaces'],replicates=original['evaluation_bootstrap'],seed=718202610,
        baseline='native_u',branches='Every output in the fixed prediction manifest; no ranking or pruning after the reference is opened',
        bootstrap='Paired whole evaluation-run resampling; vectorized simultaneous candidate metrics; training and response-model uncertainty excluded',
        quantiles=[2.5,50,97.5],multiplicity='Bonferroni per declared candidate, separately for each performance metric; does not correct all metrics and all historical experiments simultaneously',
        bins=original['bins'],no_reference_residual_cuts=True,no_model_output_cuts=True,RICH_read=False,
        exposure='Run-disjoint from prior conditional-MC measurements and no model refit here; entire project historical exposure not certified'))
    print('ADDITIONAL_SCORING_PROTOCOL',flush=True)


def score():
    frozen=gate();p=json.loads((SCORE/'protocol.json').read_text());audit=json.loads((OUT/'reference/audit.json').read_text())
    if digest(Path(__file__))!=p['script_sha256'] or digest(OUT/'protocol.json')!=p['original_protocol_sha256']:raise ValueError('Additional score protocol changed')
    if not audit['complete'] or digest(OUT/'reference/reference.npz')!=audit['reference_sha256']:raise ValueError('Additional reference is not complete/frozen')
    if (SCORE/'evaluation.json').exists():raise FileExistsError('Additional score already complete')
    data=load(OUT/'predictions_before_truth.npz');ref=load(OUT/'reference/reference.npz')
    for k in ['run','event']:
        if not np.array_equal(data[k],ref[k]):raise ValueError('Additional score identity mismatch')
    valid=ref['valid'];truth=ref['u'][valid];runs=data['run'][valid];branches=frozen['branches'];candidates=[k for k in branches if k not in ['native_u','ols_u','constant_u']]
    values=np.column_stack([data[k][valid] for k in branches]);baseindex=branches.index(p['baseline'])
    if not np.isfinite(values).all() or np.any(values==0):raise ValueError('Invalid prediction; do not silently change common events')
    errors=dict(u=values-truth[:,None],beta=1/values-1/truth[:,None]);point={s:matrix_metrics(e) for s,e in errors.items()}
    groups=[np.flatnonzero(runs==r) for r in np.unique(runs)];rng=np.random.RandomState(p['seed']);boot={s:np.zeros((p['replicates'],len(branches),len(METRICS))) for s in errors}
    for i in range(p['replicates']):
        index=np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))])
        for s,e in errors.items():boot[s][i]=matrix_metrics(e[index])
        if (i+1)%256==0:print('ADDITIONAL_BOOTSTRAP',i+1,'/',p['replicates'],flush=True)
    table={};multiplicity=[2.5/len(candidates),50,100-2.5/len(candidates)];bins=[]
    for s,e in errors.items():
        difference=matrix_difference(point[s][baseindex][None,:],point[s]);rep=matrix_difference(boot[s][:,baseindex:baseindex+1,:],boot[s]);rows={}
        for j,label in enumerate(branches):
            rows[label]=dict(metrics=dict(zip(METRICS,point[s][j])),versus_native=dict(zip(DELTA,difference[j])),
                paired_run_interval={k:np.percentile(rep[:,j,h],p['quantiles']) for h,k in enumerate(DELTA)},
                candidate_adjusted_interval={k:np.percentile(rep[:,j,h],multiplicity) for h,k in enumerate(DELTA)})
        extra={}
        for family in ['native_gaussianR','native_studentR5']:
            for ending in ['u','transport_u']:
                new=branches.index(family+'_joint_'+ending);old=branches.index(family+'_time_'+ending)
                diff=matrix_difference(point[s][old],point[s][new]);samples=matrix_difference(boot[s][:,old],boot[s][:,new]);key=family+'_joint_'+ending+'_versus_time'
                extra[key]=dict(point=dict(zip(DELTA,diff)),interval={k:np.percentile(samples[:,h],p['quantiles']) for h,k in enumerate(DELTA)})
        table[s]=dict(branches=rows,paired_joint_vs_time=extra)
        for key,edges in p['bins'].items():
            coordinate=data['R'][valid] if key=='R' else 1/truth
            for low,high in zip(edges[:-1],edges[1:]):
                select=(coordinate>=low)&(coordinate<high)
                if select.sum()<20:continue
                values=matrix_metrics(e[select]);differences=matrix_difference(values[baseindex][None,:],values)
                bins.append(dict(space=s,key=key,range=[low,high],events=int(select.sum()),branches={label:dict(metrics=dict(zip(METRICS,values[j])),versus_native=dict(zip(DELTA,differences[j]))) for j,label in enumerate(branches)}))
    np.savez_compressed(SCORE/'bootstrap_metrics.npz',**boot)
    result=dict(complete=True,created_utc=stamp(),events=len(truth),runs=len(groups),branches=branches,table=table,bins=bins,
        candidate_count=len(candidates),candidate_adjusted_quantiles=multiplicity,reference_sha256=audit['reference_sha256'],
        output_freeze_sha256=digest(OUT/'output_freeze.json'),scoring_protocol_sha256=digest(SCORE/'protocol.json'),bootstrap_metrics_sha256=digest(SCORE/'bootstrap_metrics.npz'),
        beta_above_one_retained={k:int((data[k][valid]<1).sum()) for k in branches},RICH_read=False,refit_on_additional_runs=False)
    write(SCORE/'evaluation.json',result)
    for s in ['u','beta']:
        for label in candidates:
            row=table[s]['branches'][label]
            print('ADDITIONAL_MC_RESULT',s,label,'width',row['versus_native']['w68_improvement'],'interval',row['paired_run_interval']['w68_improvement'],'RMSE',row['versus_native']['rmse_improvement'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','score']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:score()
