#!/usr/bin/env python3
"""Score frozen MC timing outputs against audited primary L5 acceptance records."""
import argparse
import json
from pathlib import Path
import numpy as np
import uproot
import awkward as ak
from research_conditional_mc import OUT,inputs,stamp
from research_conditional_mc_reference import REF,check_freeze
from research_flight_conditional import load,write,digest

SCORE=OUT/'scoring'
SEED=920261011


def metrics(error):
    e=np.asarray(error,float)
    if not len(e) or not np.isfinite(e).all():raise ValueError('Finite, nonempty errors required')
    q=np.percentile(e,[16,50,84])
    return dict(n=len(e),bias=float(e.mean()),median=float(q[1]),w68=float((q[2]-q[0])/2),
                rmse=float(np.sqrt(np.mean(e*e))),tail02=float(np.mean(abs(e)>.02)),tail05=float(np.mean(abs(e)>.05)))


def prepare():
    check_freeze();SCORE.mkdir(exist_ok=True)
    if (SCORE/'protocol.json').exists():raise FileExistsError('Scoring protocol already exists')
    freeze=json.loads((OUT/'output_freeze.json').read_text())
    p=load(OUT/'outputs_before_truth.npz')
    branches=[k for k in p if k.endswith('_u')]
    write(SCORE/'protocol.json',dict(created_utc=stamp(),branches=branches,output_freeze_sha256=digest(OUT/'output_freeze.json'),
        reference='Unique plane10 acceptance record with same primary particle, trackID, parentID, charge and mass as MCEventg[0]',
        plane_z_cm=1.698,plane_z_tolerance_cm=.05,association_tolerance_cm=1.,
        primary_requirements='Charge=6, positive finite momentum/mass, downward direction; native Tracker/BetaH association agrees with stored indices',
        primary_truth='sqrt(1+(native acceptance mass / native acceptance momentum)^2)',
        stored_reference='Audited equality to native last plane10 acceptance is checked, but primary reference is selected by identity, not storage order',
        scoring='Common accepted primary-reference events for every frozen output; no output or reference-residual cuts',
        missing_raw='Reported as missing native audit, excluded from primary-verified claim',
        bootstrap_replicates=1024,seed=SEED,bins=freeze['bins'],baseline='ols_u',
        uncertainty='Frozen-output paired evaluation-run bootstrap and joint independent training-run/evaluation-run bootstrap',
        truth_reused_for_fitting=False,RICH_read=False))
    print('SCORING_PROTOCOL_FROZEN',flush=True)


def collect(partial=False):
    check_freeze();protocol=json.loads((SCORE/'protocol.json').read_text())
    plan=json.loads((REF/'plan.json').read_text());p=load(OUT/'outputs_before_truth.npz');a,g=inputs();sel=a['fold']==0
    if not np.array_equal(p['run'],a['run'][sel]) or not np.array_equal(p['event'],a['event'][sel]):raise ValueError('Measured/prediction identities differ')
    stored=load(REF/'stored.npz')
    for k in ['run','event','source','entry']:
        if not np.array_equal(stored[k],p[k]):raise ValueError('Stored reference identity differs')
    n=len(p['run']);lookup={(int(rr),int(ee)):i for i,(rr,ee) in enumerate(zip(p['run'],p['event']))}
    native=np.zeros(n,bool);valid=np.zeros(n,bool);u=np.full(n,np.nan);mass=np.full(n,np.nan);mom=np.full(n,np.nan)
    reason=np.full(n,'raw_not_available',dtype='U32');distance=np.full(n,np.nan);z=np.full(n,np.nan)
    stored_match=np.zeros(n,bool);primary_matches_last=np.zeros(n,bool);mult=np.zeros(n,int);tk_du=np.full(n,np.nan)
    audits=[];missing=[];part_directory=REF/'parts';decoder='full_native';parent_summary=None
    if (REF/'parent_complete.json').exists():
        parent_summary=json.loads((REF/'parent_complete.json').read_text())
        if not parent_summary['complete'] or parent_summary['plan_sha256']!=digest(REF/'plan.json') or parent_summary['validation_sha256']!=digest(REF/'parent_validation.json'):
            raise ValueError('Selective decoder verification changed')
        validation=json.loads((REF/'parent_validation.json').read_text())
        if not validation['complete'] or validation['decoder_sha256']!=digest(REF/'parent_decoder.json'):
            raise ValueError('Selective decoder is not validated')
        part_directory=REF/'parent_parts';decoder='selective_native_parent'
    for r in plan['records']:
        if not r['raw']:continue
        note=part_directory/('%05d.json'%r['source']);path=note.with_suffix('.root')
        if not note.exists():missing.append(r['source']);continue
        meta=json.loads(note.read_text())
        if not meta.get('complete') or digest(path)!=meta.get('output_sha256'):raise ValueError('Invalid native audit '+str(path))
        if parent_summary is not None and (parent_summary['part_hashes'][str(r['source'])]!=meta['output_sha256'] or meta['decoder_sha256']!=parent_summary['decoder_sha256']):
            raise ValueError('Selective native part is not in the verified manifest')
        rows=ak.to_list(uproot.open(path)['reference'].arrays(library='ak'))
        if len(rows)!=r['events']:raise ValueError('Wrong native row count')
        for row in rows:
            i=lookup[(row['run'],row['event'])]
            if native[i]:raise ValueError('Duplicate native identity')
            native[i]=True;reason[i]='primary_record_missing'
            if not (row['found'] and row['track_exists'] and row['beta_exists'] and row['beta_track']==row['requested_track']):
                reason[i]='association_failure';continue
            if row['requested_track']!=a['itrtrack'][sel][i] or row['requested_beta']!=a['ibetah'][sel][i]:raise ValueError('Requested association changed')
            plane=np.asarray(row['plane']);candidates=np.flatnonzero(plane==10)
            if len(candidates):
                last=candidates[-1];pos=np.asarray(row['position']).reshape(-1,3)
                stored_match[i]=(stored['mevmom1'][i,10]==row['momentum'][last] and np.array_equal(stored['mevcoo1'][i,10],pos[last]))
                if not stored_match[i]:raise ValueError('Stored acceptance does not equal raw last plane10')
            selected=[j for j in candidates if row['particle'][j]==row['gen_particle'] and row['track'][j]==row['gen_track']
                      and row['parent'][j]==row['gen_parent'] and row['mass'][j]==row['gen_mass'] and row['charge'][j]==row['gen_charge']]
            mult[i]=len(selected)
            if len(selected)!=1:
                if len(selected)>1:reason[i]='multiple_primary_records'
                continue
            j=selected[0];pp=row['momentum'][j];mm=row['mass'][j];qq=row['charge'][j]
            xy=np.asarray(row['position'][3*j:3*j+3]);dd=np.asarray(row['direction'][3*j:3*j+3]);z[i]=xy[2]
            distance[i]=np.linalg.norm(xy-a['tk_pos'][sel][i,4]);primary_matches_last[i]=j==candidates[-1]
            if not np.isfinite([pp,mm,qq,*xy,*dd]).all() or pp<=0 or mm<=0 or qq!=6 or dd[2]>=0:
                reason[i]='invalid_primary_state';continue
            if abs(xy[2]-protocol['plane_z_cm'])>protocol['plane_z_tolerance_cm']:
                reason[i]='outside_L5_acceptance';continue
            if distance[i]>protocol['association_tolerance_cm']:
                reason[i]='primary_track_distance';continue
            valid[i]=True;reason[i]='accepted';u[i]=np.hypot(1,mm/pp);mass[i]=mm;mom[i]=pp
            kt=[j for j in range(len(row['tk_track'])) if row['tk_track'][j]==row['gen_track'] and row['tk_primary'][j]==1 and row['tk_particle'][j]==row['gen_particle'] and row['tk_momentum'][j]>0]
            if len(kt)==1:tk_du[i]=np.hypot(1,mm/row['tk_momentum'][kt[0]])-u[i]
        audits.append(dict(source=r['source'],events=len(rows),output_sha256=meta['output_sha256'],decoder=decoder))
    if missing and not partial:raise RuntimeError('Native sources pending: '+str(missing))
    counts={str(k):int(v) for k,v in zip(*np.unique(reason,return_counts=True))}
    summary=dict(created_utc=stamp(),complete=not missing,pending_sources=missing,events=n,native_events=int(native.sum()),accepted=int(valid.sum()),
        rejection_counts=counts,accepted_runs=len(np.unique(p['run'][valid])),native_source_audits=audits,
        stored_last_acceptance_equal=int(stored_match.sum()),primary_equal_last=int(primary_matches_last[valid].sum()),
        accepted_distance_cm_percentiles=np.percentile(distance[valid],[0,16,50,84,100]) if valid.any() else [],
        accepted_z_cm_percentiles=np.percentile(z[valid],[0,16,50,84,100]) if valid.any() else [],
        native_primary_masses_GeV=np.unique(mass[valid]),
        L5_tracker_cluster_u_difference=metrics(tk_du[np.isfinite(tk_du)]) if np.isfinite(tk_du).any() else None,
        protocol_sha256=digest(SCORE/'protocol.json'),decoder=decoder,
        parent_manifest_sha256=digest(REF/'parent_complete.json') if parent_summary is not None else None,
        selection_unchanged_after_decoder_validation=True)
    if partial:
        print(json.dumps({k:v for k,v in summary.items() if k!='native_source_audits'},default=lambda v:v.tolist(),indent=2));return
    np.savez_compressed(SCORE/'reference.npz',run=p['run'],event=p['event'],native=native,valid=valid,u=u,mass=mass,momentum=mom,reason=reason,distance=distance,z=z,tk_du=tk_du)
    summary['reference_sha256']=digest(SCORE/'reference.npz');write(SCORE/'reference_audit.json',summary)
    print('REFERENCE_AUDIT',counts,flush=True)


def pair_delta(old,new):
    aa=metrics(old);bb=metrics(new)
    return dict(w68_improvement=1-bb['w68']/aa['w68'],rmse_improvement=1-bb['rmse']/aa['rmse'],bias_change=bb['bias']-aa['bias'],tail02_change=bb['tail02']-aa['tail02'],tail05_change=bb['tail05']-aa['tail05'])


def score():
    check_freeze();protocol=json.loads((SCORE/'protocol.json').read_text());p=load(OUT/'outputs_before_truth.npz');r=load(SCORE/'reference.npz')
    if not np.array_equal(p['run'],r['run']) or not np.array_equal(p['event'],r['event']):raise ValueError('Reference/prediction identities differ')
    v=r['valid'];base=protocol['baseline'];branches=protocol['branches'];truth=r['u'][v]
    if not v.any():raise ValueError('No valid reference')
    runs=p['run'][v];ur=np.unique(runs);groups=[np.flatnonzero(runs==rr) for rr in ur]
    errors={space:{k:(p[k][v]-truth if space=='u' else 1/p[k][v]-1/truth) for k in branches} for space in ['u','beta']}
    if not all(np.isfinite(e).all() for vv in errors.values() for e in vv.values()):raise ValueError('Invalid frozen output; common-event scoring cannot silently drop it')
    table={space:{k:dict(metrics=metrics(e),versus_ols=pair_delta(errors[space][base],e)) for k,e in ee.items()} for space,ee in errors.items()}
    rng=np.random.RandomState(protocol['seed']);draws=[np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))]) for _ in range(protocol['bootstrap_replicates'])]
    cis={}
    for space,ee in errors.items():
        cis[space]={}
        for k,e in ee.items():
            values=[pair_delta(ee[base][ix],e[ix]) for ix in draws]
            cis[space][k]={m:np.percentile([x[m] for x in values],[2.5,50,97.5]) for m in values[0]}
    bins=[]
    for field,edges in protocol['bins'].items():
        values=p['R'][v] if field=='measured_R' else 1/truth
        for lo,hi in zip(edges[:-1],edges[1:]):
            select=(values>=lo)&(values<hi)
            if select.sum()<20:continue
            bins.append(dict(field=field,range=[lo,hi],events=int(select.sum()),
                metrics={s:{k:dict(metrics=metrics(e[select]),versus_ols=pair_delta(ee[base][select],e[select])) for k,e in ee.items()} for s,ee in errors.items()}))
    result=dict(complete=True,created_utc=stamp(),events=int(v.sum()),runs=len(ur),baseline=base,metrics=table,
        paired_evaluation_run_bootstrap=cis,bins=bins,output_freeze_sha256=digest(OUT/'output_freeze.json'),reference_sha256=digest(SCORE/'reference.npz'),
        scope='Existing MC development sources; primary L5 acceptance reference, not flight RICH validation',RICH_read=False,coefficient_retuning=False)
    write(SCORE/'evaluation.json',result)
    for k in branches:
        m=table['u'][k];print(k,'W68',m['metrics']['w68'],'gain',m['versus_ols']['w68_improvement'],'RMSE',m['metrics']['rmse'],'bias',m['metrics']['bias'],flush=True)


def joint_bootstrap():
    check_freeze();protocol=json.loads((SCORE/'protocol.json').read_text());bp=json.loads((OUT/'bootstrap/protocol.json').read_text())
    p=load(OUT/'outputs_before_truth.npz');r=load(SCORE/'reference.npz');v=r['valid'];truth=r['u'][v];runs=p['run'][v]
    groups=[np.flatnonzero(runs==rr) for rr in np.unique(runs)];results=[];parameters=[];failures=[]
    for rep in range(bp['replicates']):
        note=OUT/'bootstrap'/('%04d.json'%rep);path=note.with_suffix('.npz')
        if not note.exists():raise RuntimeError('Bootstrap pending '+str(rep))
        meta=json.loads(note.read_text())
        if not meta.get('complete') or digest(path)!=meta['output_sha256']:raise ValueError('Bootstrap failed '+str(rep))
        b=load(path)
        if not np.array_equal(b['run'],p['run']) or not np.array_equal(b['event'],p['event']):raise ValueError('Bootstrap identity mismatch')
        rng=np.random.RandomState(SEED+10000+rep);ix=np.concatenate([groups[j] for j in rng.randint(0,len(groups),len(groups))])
        for model in meta['models']:
            best=model['best'];parameters.append(dict(replicate=rep,label=model['label'],theta=best['theta'],projected_gradient_max=best['projected_gradient_max'],success=best['success']))
            if not best['success'] or best['projected_gradient_max']>1e-3:failures.append(dict(replicate=rep,label=model['label']))
        for k in [k for k in b if k.endswith('_u')]:
            for space in ['u','beta']:
                e=b[k][v]-truth if space=='u' else 1/b[k][v]-1/truth
                old=p['ols_u'][v]-truth if space=='u' else 1/p['ols_u'][v]-1/truth
                results.append(dict(replicate=rep,label=k,space=space,training_only=pair_delta(old,e),joint=pair_delta(old[ix],e[ix])))
    cis={}
    for space in ['u','beta']:
        cis[space]={}
        for label in sorted(set(x['label'] for x in results)):
            rr=[x for x in results if x['label']==label and x['space']==space]
            cis[space][label]={scope:{key:np.percentile([x[scope][key] for x in rr],[2.5,50,97.5]) for key in rr[0][scope]} for scope in ['training_only','joint']}
    write(SCORE/'joint_bootstrap.json',dict(complete=True,created_utc=stamp(),replicates=bp['replicates'],convergence_failures=failures,
        intervals=cis,parameters=parameters,replicate_metrics=results,reference_sha256=digest(SCORE/'reference.npz'),bootstrap_protocol_sha256=digest(OUT/'bootstrap/protocol.json')))
    print('JOINT_BOOTSTRAP_COMPLETE',bp['replicates'],'convergence failures',len(failures),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','collect','score','bootstrap']);parser.add_argument('--partial',action='store_true');args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='collect':collect(args.partial)
    elif args.stage=='score':score()
    else:joint_bootstrap()
