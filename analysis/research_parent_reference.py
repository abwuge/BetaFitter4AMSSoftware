#!/usr/bin/env python3
"""Validate a selective native decoder, then extract the frozen reference plan."""
import argparse
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import awkward as ak
from research_conditional_mc_reference import REF,check_freeze
from research_conditional_mc import stamp
from research_flight_conditional import write,digest

PARTS=REF/'parent_parts'
ENV='/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'


def compare(first,second):
    with uproot.open(first) as stream:a=stream['reference'].arrays(library='ak')
    with uproot.open(second) as stream:b=stream['reference'].arrays(library='ak')
    if set(a.fields)!=set(b.fields) or len(a)!=len(b):raise ValueError('Reference schema/length differs')
    checked=[]
    for field in a.fields:
        if field=='tk_position':continue
        for x,y in zip(ak.to_list(a[field]),ak.to_list(b[field])):
            xx=np.asarray(x);yy=np.asarray(y)
            if xx.shape!=yy.shape or not np.all((xx==yy)|(np.isnan(xx)&np.isnan(yy))):raise ValueError('Native decoder mismatch: '+field)
        checked.append(field)
    return dict(events=len(a),equal_fields=checked,deliberately_unread=['tk_position'],first_sha256=digest(first),second_sha256=digest(second))


def verify(path):
    check_freeze();decoder=json.loads((REF/'parent_decoder.json').read_text())
    if digest(REF/'extract_parent')!=decoder['executable_sha256']:raise ValueError('Parent decoder executable changed')
    record=compare(REF/'parts/00001.root',path)
    write(REF/'parent_validation.json',dict(complete=True,created_utc=stamp(),decoder_sha256=digest(REF/'parent_decoder.json'),
        original_source=1,comparison=record,scope='Exact native fields used by primary L5 reference and measured association; Tracker global MC-cluster positions are deliberately not recomputed'))
    print('PARENT_REFERENCE_VALIDATED',record['events'],'events',len(record['equal_fields']),'fields',flush=True)


def extract(record):
    tag='%05d'%record['source'];path=PARTS/(tag+'.root');note=PARTS/(tag+'.json');log=PARTS/(tag+'.log')
    if note.exists():
        previous=json.loads(note.read_text())
        if previous.get('complete') and digest(path)==previous['output_sha256']:return previous
        raise ValueError('Existing invalid parent part')
    if path.exists():raise FileExistsError('Incomplete parent output '+str(path))
    command=['bash','-c','source "$1" >/dev/null; exec "$2" "$3" "$4" "$5"','parent_extract',ENV,str(REF/'extract_parent'),record['raw'],str(REF/'targets'/(tag+'.txt')),str(path)]
    start=time.monotonic()
    with log.open('w') as stream:result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError('Parent extraction failed '+tag+' '+str(result.returncode))
    with uproot.open(path) as stream:
        a=stream['reference'].arrays(['run','event','found','track_exists','beta_exists','beta_track','requested_track'],library='np')
    if len(a['run'])!=record['events'] or not a['found'].all():raise ValueError('Incomplete parent reference '+tag)
    comparison=None;old=REF/'parts'/(tag+'.json')
    if old.exists() and json.loads(old.read_text()).get('complete'):comparison=compare(old.with_suffix('.root'),path)
    summary=dict(complete=True,created_utc=stamp(),source=record['source'],events=len(a['run']),elapsed_seconds=time.monotonic()-start,
        output_sha256=digest(path),log_sha256=digest(log),plan_sha256=digest(REF/'plan.json'),
        decoder_sha256=digest(REF/'parent_decoder.json'),validation_sha256=digest(REF/'parent_validation.json'),native_comparison=comparison)
    write(note,summary);print('PARENT_SOURCE',tag,len(a['run']),'seconds',summary['elapsed_seconds'],flush=True)
    return summary


def all_sources(workers):
    check_freeze();validation=json.loads((REF/'parent_validation.json').read_text());decoder=json.loads((REF/'parent_decoder.json').read_text())
    if not validation['complete'] or validation['decoder_sha256']!=digest(REF/'parent_decoder.json') or digest(REF/'extract_parent')!=decoder['executable_sha256']:raise ValueError('Parent decoder validation changed')
    PARTS.mkdir(exist_ok=True);plan=json.loads((REF/'plan.json').read_text());records=[r for r in plan['records'] if r['raw']]
    with ThreadPoolExecutor(max_workers=workers) as pool:parts=list(pool.map(extract,records))
    comparisons=[]
    for record in records:
        tag='%05d'%record['source'];note=REF/'parts'/(tag+'.json')
        if note.exists() and json.loads(note.read_text()).get('complete'):
            comparisons.append(dict(source=record['source'],**compare(note.with_suffix('.root'),PARTS/(tag+'.root'))))
    write(REF/'parent_complete.json',dict(complete=True,created_utc=stamp(),sources=len(parts),events=sum(r['events'] for r in parts),
        full_native_crosschecks=comparisons,part_hashes={str(r['source']):r['output_sha256'] for r in parts},
        validation_sha256=digest(REF/'parent_validation.json'),decoder_sha256=digest(REF/'parent_decoder.json'),plan_sha256=digest(REF/'plan.json')))
    print('PARENT_ALL_COMPLETE',len(parts),'sources',sum(r['events'] for r in parts),'events',len(comparisons),'full native comparisons',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['verify','extract']);parser.add_argument('--path',default='/tmp/conditional_mc_parent_82335416.root');parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if args.stage=='verify':verify(args.path)
    else:all_sources(args.workers)
