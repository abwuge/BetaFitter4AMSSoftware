#!/usr/bin/env python3
"""Prepare a charge-consistent RICH reference and the existing flight holdout.

Reference quality v2 was specified after the loose-reference MC audit, before
any flight beta score. It does not change either reconstruction component.
"""
import argparse
import csv
import json
import multiprocessing as mp
from pathlib import Path
import time
import numpy as np
import uproot
import rich_transport as rt
import evaluate_rich_transport as ev

DEST=rt.ROOT/'results/rich_transport'
CACHE=Path('/tmp/betaEstimator-rich-transport')
MC=None


def specification():
    p=DEST/'reference_quality_v2.json'
    value={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
           'main_RICH_charge_cut':'4 < rich_q[1] < 8; finite; in addition to frozen loose RICH cuts',
           'sensitivity_RICH_charge_cut':'5 < rich_q[1] < 7',
           'reason':'Reject low-charge rings attached to a reconstructed carbon track; beta agreement is never used',
           'timing':'specified following MC reference-failure audit, before flight beta evaluation',
           'data_scope':'all pre-existing validation-role files in config/data_sources.conf; no flight fitting',
           'first_all_flight_preparation':'cancelled after input-only reads; no flight predictions or scores were computed',
           'model_freeze_sha256':ev.digest(DEST/'freeze.json'),'preparation_code_sha256':ev.digest(Path(__file__).resolve())}
    if not p.exists():ev.dump(p,value)
    else:
        old=json.loads(p.read_text())
        for key in ['main_RICH_charge_cut','sensitivity_RICH_charge_cut','data_scope','model_freeze_sha256','preparation_code_sha256']:
            if value[key]!=old[key]:raise RuntimeError('Reference specification changed: '+key)


def flight_part(task):
    source,role,path=task
    part=CACHE/'flight_holdout_parts'/('source_%04d.npz'%source)
    note=part.with_suffix('.json')
    if part.exists() and note.exists():
        with np.load(str(part)) as f:a={k:f[k] for k in f.files}
        return a,json.loads(note.read_text())
    with uproot.open(path) as f:
        t=f['amstreea']
        a=t.arrays(ev.MEASURED+ev.SELECTION+['rich_q'],library='np')
        rows=np.flatnonzero(rt.core_selection(a))
        out={k:a[k][rows] for k in ev.MEASURED+['rich_q']}
        out['source']=np.full(len(rows),source,np.int32)
        out['outer_holdout']=np.ones(len(rows),bool)
        audit={'source':source,'role':role,'path':path,'raw':int(t.num_entries),'selected':len(rows)}
    np.savez_compressed(str(part),**out);ev.dump(note,audit)
    return out,audit


def mc_part(task):
    label,source,path=task
    rows=np.flatnonzero((MC['label']==label)&(MC['source']==source)&MC['rich_quality'])
    if not len(rows):return rows,np.empty((0,2),np.float32)
    with uproot.open(path) as f:
        a=f['amstreea'].arrays(['run','event','rich_q'],library='np')
    take=ev.match_rows(a['run'],a['event'],MC['run'][rows],MC['event'][rows])
    return rows,a['rich_q'][take]


def main():
    global MC
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['data','mc_charge']);ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args();ev.freeze(DEST);specification();CACHE.mkdir(exist_ok=True)
    if args.kind=='data':
        (CACHE/'flight_holdout_parts').mkdir(exist_ok=True)
        all_tasks=[]
        for line in (rt.ROOT/'config/data_sources.conf').read_text().splitlines():
            if line.strip() and not line.startswith('#'):
                role,charge,path=line.split();all_tasks.append((len(all_tasks),role,path))
        tasks=[t for t in all_tasks if t[1]=='validation']
        parts,audit=[],[]
        with mp.Pool(args.workers) as pool:
            for j,(a,note) in enumerate(pool.imap(flight_part,tasks),1):
                parts.append(a);audit.append(note)
                if j%12==0 or j==len(tasks):print('flight holdout',j,'/',len(tasks),'selected',sum(x['selected'] for x in audit),flush=True)
        a=ev.concatenate(parts)
        ids=np.rec.fromarrays([a['run'],a['event']],names='run,event')
        if len(np.unique(ids))!=len(ids):raise RuntimeError('Duplicated flight events')
        np.savez_compressed(str(CACHE/'data_evaluation_input.npz'),**a)
        ev.dump(DEST/'data_source_audit.json',audit)
        print('FLIGHT HOLDOUT READY',len(a['run']),len(tasks),flush=True)
    else:
        with np.load(str(DEST/'mc_predictions.npz')) as f:MC={k:f[k] for k in ['run','event','label','source','rich_quality']}
        with (ev.EXP/'manifest.tsv').open() as f:
            tasks=[({'C12':0,'C13':1}[r['isotope']],int(r['source']),r['path']) for r in csv.DictReader(f,delimiter='\t')]
        q=np.full((len(MC['run']),2),np.nan,np.float32);seen=np.zeros(len(q),bool)
        with mp.Pool(args.workers) as pool:
            for j,(rows,values) in enumerate(pool.imap(mc_part,tasks),1):
                q[rows]=values;seen[rows]=True
                if j%48==0 or j==len(tasks):print('MC RICH charge',j,'/',len(tasks),flush=True)
        if not np.array_equal(seen,MC['rich_quality']):raise RuntimeError('Incomplete charge-reference join')
        np.savez_compressed(str(DEST/'mc_rich_charge.npz'),run=MC['run'],event=MC['event'],rich_q=q)
        print('MC CHARGE READY',int(seen.sum()),flush=True)

if __name__=='__main__':main()
