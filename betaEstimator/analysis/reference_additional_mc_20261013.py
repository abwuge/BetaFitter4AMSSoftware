#!/usr/bin/env python3
"""Native primary-L5 reference for additional runs, gated by output freezing."""
import argparse,json,subprocess,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import awkward as ak
from research_additional_mc_20261013 import OUT,ROOT,MCOUT,ENV
from research_flight_conditional import load,write,digest
from research_conditional_mc import stamp
from evaluate_conditional_mc import metrics

REF=OUT/'reference';DECODER=MCOUT/'reference'


def gate():
    frozen=json.loads((OUT/'output_freeze.json').read_text())
    if not frozen['complete'] or digest(OUT/'predictions_before_truth.npz')!=frozen['output_sha256']:raise ValueError('Additional predictions not frozen')
    if digest(OUT/'measured.npz')!=frozen['measured_sha256'] or digest(OUT/'moments.npz')!=frozen['moments_sha256']:raise ValueError('Additional frozen measurements changed')
    for label,row in frozen['models'].items():
        if digest(row['model_path'])!=row['model']:raise ValueError('Additional model changed '+label)
    decoder=json.loads((DECODER/'parent_decoder.json').read_text());validation=json.loads((DECODER/'parent_validation.json').read_text())
    if not validation['complete'] or validation['decoder_sha256']!=digest(DECODER/'parent_decoder.json') or digest(DECODER/'extract_parent')!=decoder['executable_sha256']:raise ValueError('Previously audited native decoder changed')
    return frozen


def plan():
    gate();REF.mkdir(exist_ok=True);(REF/'targets').mkdir(exist_ok=True);(REF/'parts').mkdir(exist_ok=True)
    if (REF/'plan.json').exists():raise FileExistsError('Additional native reference plan exists')
    p=json.loads((OUT/'protocol.json').read_text());a=load(OUT/'measured.npz');records=[]
    for record in p['records']:
        take=np.flatnonzero(a['source']==record['source'])
        if not len(take):continue
        path=REF/'targets'/('%05d.txt'%record['source']);values=np.column_stack([a['run'][take],a['event'][take],a['itrtrack'][take],a['ibetah'][take],a['bar'][take]])
        np.savetxt(path,values,fmt='%d');records.append(dict(record,events=len(take),targets=str(path),targets_sha256=digest(path)))
    write(REF/'plan.json',dict(created_utc=stamp(),records=records,events=len(a['run']),output_freeze_sha256=digest(OUT/'output_freeze.json'),
        native_decoder_sha256=digest(DECODER/'extract_parent'),native_validation_sha256=digest(DECODER/'parent_validation.json'),
        reference_selection=p['reference_conditions'],script_sha256=digest(Path(__file__)),MC_truth_read=False,RICH_read=False))
    print('ADDITIONAL_REFERENCE_PLAN',len(records),len(a['run']),flush=True)


def extract_one(record):
    tag='%05d'%record['source'];path=REF/'parts'/(tag+'.root');note=path.with_suffix('.json');log=path.with_suffix('.log')
    if note.exists():
        previous=json.loads(note.read_text())
        if previous['complete'] and digest(path)==previous['output_sha256']:return previous
        raise ValueError('Existing incomplete additional native part')
    if path.exists():raise FileExistsError('Unrecorded existing native part')
    if digest(record['targets'])!=record['targets_sha256']:raise ValueError('Native targets changed')
    start=time.monotonic();command=['bash','-c','source "$1" >/dev/null; exec "$2" "$3" "$4" "$5"','native_reference',str(ENV),str(DECODER/'extract_parent'),record['raw'],record['targets'],str(path)]
    with log.open('w') as stream:answer=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    if answer.returncode:raise RuntimeError('Additional native extraction failed '+tag)
    with uproot.open(path) as stream:a=stream['reference'].arrays(['run','event','found'],library='np')
    expected=np.loadtxt(record['targets'],dtype=np.int64,ndmin=2)
    if len(a['run'])!=record['events'] or not a['found'].all() or not np.array_equal(a['run'],expected[:,0]) or not np.array_equal(a['event'],expected[:,1]):raise ValueError('Additional native identity mismatch')
    result=dict(complete=True,created_utc=stamp(),source=record['source'],events=len(a['run']),elapsed_seconds=time.monotonic()-start,
        output_sha256=digest(path),log_sha256=digest(log),plan_sha256=digest(REF/'plan.json'))
    write(note,result);print('ADDITIONAL_NATIVE_SOURCE',tag,len(a['run']),result['elapsed_seconds'],flush=True)
    return result


def extract(workers):
    gate();plan=json.loads((REF/'plan.json').read_text())
    if digest(Path(__file__))!=plan['script_sha256'] or digest(DECODER/'extract_parent')!=plan['native_decoder_sha256']:raise ValueError('Additional native reference implementation changed')
    with ThreadPoolExecutor(max_workers=workers) as pool:parts=list(pool.map(extract_one,plan['records']))
    write(REF/'complete.json',dict(complete=True,created_utc=stamp(),events=sum(p['events'] for p in parts),sources=len(parts),
        part_hashes={str(p['source']):p['output_sha256'] for p in parts},plan_sha256=digest(REF/'plan.json'),output_freeze_sha256=digest(OUT/'output_freeze.json')))
    print('ADDITIONAL_NATIVE_COMPLETE',len(parts),sum(p['events'] for p in parts),flush=True)


def collect():
    gate();plan=json.loads((REF/'plan.json').read_text());complete=json.loads((REF/'complete.json').read_text())
    if not complete['complete'] or complete['plan_sha256']!=digest(REF/'plan.json'):raise ValueError('Additional native reference incomplete')
    if (REF/'reference.npz').exists():raise FileExistsError('Additional reference already collected')
    a=load(OUT/'measured.npz');p=load(OUT/'predictions_before_truth.npz');n=len(a['run']);lookup={(int(r),int(e)):i for i,(r,e) in enumerate(zip(a['run'],a['event']))}
    for k in ['run','event','source','entry']:
        if not np.array_equal(a[k],p[k]):raise ValueError('Additional measured prediction IDs differ')
    native=np.zeros(n,bool);valid=np.zeros(n,bool);truth=np.full(n,np.nan);mass=np.full(n,np.nan);mom=np.full(n,np.nan);distance=np.full(n,np.nan);z=np.full(n,np.nan);tkdu=np.full(n,np.nan)
    reason=np.full(n,'native_missing',dtype='U32');conditions=plan['reference_selection']
    for record in plan['records']:
        path=REF/'parts'/('%05d.root'%record['source'])
        if digest(path)!=complete['part_hashes'][str(record['source'])]:raise ValueError('Native additional part changed')
        with uproot.open(path) as stream:rows=ak.to_list(stream['reference'].arrays(library='ak'))
        for row in rows:
            i=lookup[(row['run'],row['event'])]
            if native[i]:raise ValueError('Duplicate native additional event')
            native[i]=True;reason[i]='primary_record_missing'
            if not (row['found'] and row['track_exists'] and row['beta_exists'] and row['beta_track']==row['requested_track']):reason[i]='association_failure';continue
            if row['requested_track']!=a['itrtrack'][i] or row['requested_beta']!=a['ibetah'][i]:raise ValueError('Requested additional association changed')
            candidates=np.flatnonzero(np.asarray(row['plane'])==conditions['plane'])
            choose=[j for j in candidates if row['particle'][j]==row['gen_particle'] and row['track'][j]==row['gen_track'] and row['parent'][j]==row['gen_parent'] and row['mass'][j]==row['gen_mass'] and row['charge'][j]==row['gen_charge']]
            if len(choose)!=1:
                if len(choose)>1:reason[i]='multiple_primary_records'
                continue
            j=choose[0];pp=row['momentum'][j];mm=row['mass'][j];qq=row['charge'][j];position=np.asarray(row['position'][3*j:3*j+3]);direction=np.asarray(row['direction'][3*j:3*j+3])
            z[i]=position[2];distance[i]=np.linalg.norm(position-a['tk_pos'][i,4])
            if not np.isfinite([pp,mm,qq,*position,*direction]).all() or pp<=0 or mm<=0 or qq!=conditions['charge'] or direction[2]>=0:reason[i]='invalid_primary_state';continue
            if abs(z[i]-conditions['z'])>conditions['z_tolerance']:reason[i]='outside_L5_acceptance';continue
            if distance[i]>conditions['track_distance']:reason[i]='primary_track_distance';continue
            valid[i]=True;reason[i]='accepted';truth[i]=np.hypot(1,mm/pp);mass[i]=mm;mom[i]=pp
            tracker=[j for j in range(len(row['tk_track'])) if row['tk_track'][j]==row['gen_track'] and row['tk_primary'][j]==1 and row['tk_particle'][j]==row['gen_particle'] and row['tk_momentum'][j]>0]
            if len(tracker)==1:tkdu[i]=np.hypot(1,mm/row['tk_momentum'][tracker[0]])-truth[i]
    if not native.all():raise ValueError('Missing additional native audit')
    np.savez_compressed(REF/'reference.npz',run=a['run'],event=a['event'],native=native,valid=valid,u=truth,mass=mass,momentum=mom,reason=reason,distance=distance,z=z,tk_du=tkdu)
    summary=dict(complete=True,created_utc=stamp(),events=n,accepted=int(valid.sum()),runs=len(np.unique(a['run'][valid])),
        rejection_counts={str(k):int(v) for k,v in zip(*np.unique(reason,return_counts=True))},reference_sha256=digest(REF/'reference.npz'),
        native_manifest_sha256=digest(REF/'complete.json'),output_freeze_sha256=digest(OUT/'output_freeze.json'),
        accepted_z_percentiles=np.percentile(z[valid],[0,16,50,84,100]),accepted_distance_percentiles=np.percentile(distance[valid],[0,16,50,84,100]),
        native_primary_mass_values=np.unique(mass[valid]),tracker_primary_u_difference=metrics(tkdu[np.isfinite(tkdu)]) if np.isfinite(tkdu).any() else None,
        selection='Same original primary-L5 audit conditions; no output-error or isotope-posterior cuts',reference_used_to_refit=False,RICH_read=False)
    write(REF/'audit.json',summary);print('ADDITIONAL_REFERENCE_ACCEPTED',int(valid.sum()),summary['rejection_counts'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['plan','extract','collect']);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    if args.stage=='extract':extract(args.workers)
    else:globals()[args.stage]()
