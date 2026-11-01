#!/usr/bin/env python3
"""Prepare and execute native-reference audits, gated by the output freeze."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import uproot
from research_conditional_mc import OUT,ROOT,stamp
from research_flight_conditional import write,digest,load

REF=OUT/'reference'
RAW_C12=Path('/lustre03/data/eos/ams/MC/AMS02/2022/C.B1308/c12.pl1.l1.612000.6_05')
ENV=Path('/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh')
FIELDS=['run','event','mpar','mmom','mch','mmass','mevmom1','mevcoo1','mevdir1',
        'mtrmom','mtrpri','mtrpar','mtrcoo','mtrcoo1','mtrz']


def check_freeze():
    frozen=json.loads((OUT/'output_freeze.json').read_text())
    if digest(OUT/'outputs_before_truth.npz')!=frozen['output_sha256']:raise ValueError('Changed prediction snapshot')
    for label,h in frozen['model_freeze_sha256'].items():
        if digest(OUT/label/'freeze.json')!=h:raise ValueError('Changed model freeze')
    return frozen


def prepare():
    check_freeze();REF.mkdir(exist_ok=True);(REF/'targets').mkdir(exist_ok=True);(REF/'parts').mkdir(exist_ok=True)
    if (REF/'plan.json').exists():raise FileExistsError('Reference audit plan already exists')
    data=load(OUT/'measured.npz');test=data['fold']==0
    protocol=json.loads((OUT/'protocol.json').read_text());records=[]
    for r in protocol['records']:
        rows=np.flatnonzero(test & (data['source']==r['source']))
        if not len(rows):continue
        p=Path(r['path']);raw=None
        if p.parent.name=='Z6_C12':
            candidate=RAW_C12/(p.stem.replace('_','.')+'.root')
            if candidate.exists():raw=str(candidate)
        target=REF/'targets'/('%05d.txt'%r['source'])
        values=np.column_stack([data['run'][rows],data['event'][rows],data['itrtrack'][rows],data['ibetah'][rows],data['bar'][rows]])
        np.savetxt(target,values,fmt='%d')
        records.append(dict(source=r['source'],path=r['path'],raw=raw,events=len(rows),entry_stop=r['entry_stop'],
                            targets=str(target),targets_sha256=digest(target)))
    value=dict(created_utc=stamp(),records=records,events=int(test.sum()),native_events=sum(r['events'] for r in records if r['raw']),
        native_sources=sum(bool(r['raw']) for r in records),missing_raw_sources=sum(not bool(r['raw']) for r in records),
        output_freeze_sha256=digest(OUT/'output_freeze.json'),extractor_sha256=digest(REF/'extract'),
        reference_fields=FIELDS,truth_read=False,source_status='C12 raw input path reconstructed from the original production task contract; C13 raw inputs not located')
    write(REF/'plan.json',value);print('REFERENCE_PLAN', {k:v for k,v in value.items() if k!='records'},flush=True)


def run_source(source):
    check_freeze();plan=json.loads((REF/'plan.json').read_text());r=[r for r in plan['records'] if r['source']==source][0]
    if r['raw'] is None:return
    if digest(REF/'extract')!=plan['extractor_sha256'] or digest(r['targets'])!=r['targets_sha256']:raise ValueError('Audit input changed')
    parent_path=REF/'parent_complete.json'
    if parent_path.exists():
        parent=json.loads(parent_path.read_text());validation=json.loads((REF/'parent_validation.json').read_text())
        candidate=REF/'parent_parts'/('%05d.root'%source);candidate_note=candidate.with_suffix('.json')
        if not parent['complete'] or parent['plan_sha256']!=digest(REF/'plan.json') or parent['validation_sha256']!=digest(REF/'parent_validation.json') or not validation['complete']:
            raise ValueError('Native parent reference manifest changed')
        if validation['decoder_sha256']!=digest(REF/'parent_decoder.json') or digest(candidate)!=parent['part_hashes'][str(source)]:
            raise ValueError('Native parent decoder or output changed')
        result=json.loads(candidate_note.read_text())
        if not result['complete'] or result['output_sha256']!=parent['part_hashes'][str(source)]:raise ValueError('Incomplete native parent output')
        print('NATIVE_REFERENCE_REUSED',source,'selective_parent',flush=True)
        return result
    output=REF/'parts'/('%05d.root'%source);note=output.with_suffix('.json');log=output.with_suffix('.log')
    if note.exists() and json.loads(note.read_text()).get('complete'):
        old=json.loads(note.read_text())
        if digest(output)!=old['output_sha256']:raise ValueError('Native audit output changed')
        return old
    if output.exists():raise FileExistsError('Incomplete native audit exists: '+str(output))
    command=['bash','-c','source "$1" >/dev/null; exec "$2" "$3" "$4" "$5"','_',str(ENV),str(REF/'extract'),r['raw'],r['targets'],str(output)]
    start=time.monotonic()
    with log.open('w') as stream:result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    value=dict(complete=result.returncode==0,created_utc=stamp(),source=source,events=r['events'],exit_code=result.returncode,
               elapsed_seconds=time.monotonic()-start,output_sha256=digest(output) if output.exists() else None,
               log_sha256=digest(log),plan_sha256=digest(REF/'plan.json'))
    if result.returncode==0:
        with uproot.open(output) as stream:
            a=stream['reference'].arrays(['run','event','found'],library='np')
        expected=np.loadtxt(r['targets'],dtype=np.int64,ndmin=2)
        if len(a['run'])!=r['events'] or not a['found'].all() or not np.array_equal(a['run'],expected[:,0]) or not np.array_equal(a['event'],expected[:,1]):
            value['complete']=False;value['identity_failure']=True
    write(note,value);print('NATIVE_REFERENCE',source,value,flush=True)
    if not value['complete']:raise RuntimeError('Native reference extraction failed')
    return value


def stored():
    check_freeze();plan=json.loads((REF/'plan.json').read_text());pred=load(OUT/'outputs_before_truth.npz')
    if (REF/'stored.npz').exists():raise FileExistsError('Stored reference snapshot already exists')
    data={k:pred[k] for k in ['run','event','source','entry']};seen=np.zeros(len(pred['run']),bool);arrays={};audit=[]
    for r in plan['records']:
        rows=np.flatnonzero(pred['source']==r['source']);entries=pred['entry'][rows]
        with uproot.open(r['path']) as stream:
            tree=stream['amstreea'];present=[k for k in FIELDS if k in tree.keys()]
            missing=[k for k in FIELDS if k not in present]
            if set(missing)-{'mmass'}:raise ValueError('Missing reference fields: '+str(missing))
            a=tree.arrays(present,entry_stop=r['entry_stop'],library='np')
            if 'mmass' not in present:a['mmass']=np.full(len(a['run']),np.nan,dtype=np.float32)
        if not np.array_equal(a['run'][entries],pred['run'][rows]) or not np.array_equal(a['event'][entries],pred['event'][rows]):raise ValueError('Stored reference identity mismatch')
        for k in FIELDS:
            if k in ['run','event']:continue
            if k not in arrays:arrays[k]=np.empty((len(pred['run']),)+a[k].shape[1:],a[k].dtype)
            arrays[k][rows]=a[k][entries]
        seen[rows]=True;audit.append(dict(source=r['source'],events=len(rows),fields_read=present,unavailable_fields=missing))
    if not seen.all():raise ValueError('Missing stored reference events')
    data.update(arrays);np.savez_compressed(REF/'stored.npz',**data)
    write(REF/'stored_extraction.json',dict(complete=True,created_utc=stamp(),events=len(seen),records=audit,
        fields_read=FIELDS,output_sha256=digest(REF/'stored.npz'),output_freeze_sha256=digest(OUT/'output_freeze.json'),RICH_read=False))
    print('STORED_REFERENCE_COMPLETE',len(seen),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','stored']);p.add_argument('--source',type=int)
    p.add_argument('--worker',type=int);p.add_argument('--workers',type=int,default=50);args=p.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='stored':stored()
    elif args.source is not None:run_source(args.source)
    else:
        tasks=[r['source'] for r in json.loads((REF/'plan.json').read_text())['records'] if r['raw']]
        if args.worker is None:raise ValueError('Source or worker is required')
        for source in tasks[args.worker::args.workers]:run_source(source)
