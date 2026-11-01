#!/usr/bin/env python3
"""Prepare native-only jobs or validate their outputs on the analysis host."""
import argparse
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import uproot
import flight_beta_wide as w
from prepare_tofcalib_validation import OUT, EXTRA
from tofcalib_reference_preselection import reference_preselection


def original(source):
    dest=OUT/'parts'/('%05d'%source)
    with np.load(str(dest/'original.npz'),allow_pickle=False) as f:
        a={k:f[k] for k in f.files}
    candidates=len(a['run'])
    a=w.ev.sub(a,np.flatnonzero(reference_preselection(a)))
    return dest,a,candidates


def prepare_job(source):
    dest,a,candidates=original(source)
    n=len(a['run'])
    if not n:
        w.ev.dump(dest/'supplement.json',dict(complete=True,source=source,requested=0,
                  selected=0,candidates=candidates))
        return dict(source=source,requested=0,candidates=candidates)
    with (dest/'evaluation_targets.txt').open('w') as f:
        for i in range(n):
            q=a['tk_qin'][i,0,2] if a['tk_qin'][i,0,2]>8.5 else a['tk_q'][i,1]
            f.write('%d %d %d %d %d %.9g\n'%(a['run'][i],a['event'][i],a['itrtrack'][i],
                    a['ibetah'][i],a['irich'][i],q))
    runs=set(int(x) for x in a['run'])
    files=[x for x in (dest/'raw_files.txt').read_text().splitlines()
           if int(Path(x).name.split('.')[0]) in runs]
    (dest/'evaluation_raw_files.txt').write_text('\n'.join(files)+'\n')
    record=dict(source=source,requested=n,candidates=candidates,
                original_sha256=w.ev.digest(dest/'original.npz'))
    w.ev.dump(dest/'native_job.json',record)
    return record


def finalize_job(source):
    dest=OUT/'parts'/('%05d'%source)
    note=dest/'supplement.json'
    if note.exists() and json.loads(note.read_text()).get('complete'):
        return json.loads(note.read_text())
    if not (dest/'.native_complete').exists():
        return None
    expected=(OUT/'bin/supplement_light.sha256').read_text().strip()
    if (dest/'.native_complete').read_text().strip()!=expected:
        raise RuntimeError('Native output binary hash changed')
    dest,a,candidates=original(source)
    n=len(a['run'])
    job=json.loads((dest/'native_job.json').read_text())
    if job['original_sha256']!=w.ev.digest(dest/'original.npz') or job['requested']!=n:
        raise RuntimeError('Prepared event population changed')
    temporary=dest/'extra.partial.root'
    with uproot.open(str(temporary)) as f:
        extra=f['beta_extra'].arrays(library='np')
    if len(extra['run'])!=n or not np.array_equal(extra['run'],a['run']) or not np.array_equal(extra['event'],a['event']):
        raise RuntimeError('Incomplete or reordered native join for %d'%source)
    if not extra['found'].all():
        raise RuntimeError('Missing native events')
    valid=extra['selected']==1
    associations={k:int(np.sum(valid & (extra[k]!=a[k]))) for k in ('itrtrack','ibetah','irich')}
    if any(associations.values()):
        raise RuntimeError('Native/NEWTOF association mismatch: %s'%associations)
    joined=w.ev.sub(a,np.flatnonzero(valid))
    for key in EXTRA:
        joined[key]=extra[key][valid]
    w.set_edges(w.R_EDGES)
    joined=w.ev.sub(joined,np.flatnonzero(w.selection(joined)))
    with np.load(str(OUT/'training_runs.npz'),allow_pickle=False) as training:
        if np.intersect1d(joined['run'],training['run']).size:
            raise RuntimeError('Training run entered validation')
    np.savez_compressed(str(dest/'joined.npz'),**joined)
    temporary.replace(dest/'extra.root')
    record=dict(complete=True,source=source,requested=n,candidates=candidates,
                existing_reference_quality_preselection=True,
                extracted_rigidity_elements='all Z conventions at inner span only',
                extracted_error_elements='Kalman V6 inner span only; unused elements NaN',
                native_selected=int(valid.sum()),selected=len(joined['run']),
                producer_rejected=int((~valid).sum()),association_mismatches=associations,
                original_newtof_measurements_preserved=True,
                original_sha256=job['original_sha256'],supplement_sha256=w.ev.digest(dest/'extra.root'),
                worker_sha256=w.ev.digest(Path(__file__)),
                binary_sha256=w.ev.digest(OUT/'bin/supplement_light'),
                model_freeze_sha256=w.ev.digest(w.DEST/'freeze.json'))
    w.ev.dump(note,record)
    return record


def main():
    p=argparse.ArgumentParser()
    p.add_argument('stage',choices=['prepare','finalize'])
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--source',type=int)
    args=p.parse_args()
    w.check_freeze()
    proof=json.loads((OUT/'probe/batch_equivalence.json').read_text())
    digest=w.ev.digest(OUT/'bin/supplement_light')
    if not proof['required_measurements_equal'] or proof['light_binary_sha256']!=digest:
        raise RuntimeError('Native binary has not passed equivalence')
    indices=json.loads((OUT/'selection.json').read_text())['file_indices']
    tasks=[args.source] if args.source is not None else indices
    fn=prepare_job if args.stage=='prepare' else finalize_job
    records=[]
    with mp.Pool(args.workers) as pool:
        for i,record in enumerate(pool.imap_unordered(fn,tasks),1):
            if record is not None:records.append(record)
            if i%100==0 or i==len(tasks):
                print(args.stage,i,'/',len(tasks),'ready',len(records),flush=True)
    if args.stage=='prepare' and args.source is None:
        for worker in range(8):
            (OUT/'condor'/('tasks_%03d.txt'%worker)).write_text('\n'.join('%05d'%s for s in indices[worker::8])+'\n')
        (OUT/'bin/supplement_light.sha256').write_text(digest+'\n')
        w.ev.dump(OUT/'native_preparation.json',dict(complete=True,source_files=len(records),
                  candidates=sum(r['candidates'] for r in records),
                  reference_preselected=sum(r['requested'] for r in records),
                  binary_sha256=digest,records=records))
    elif args.stage=='finalize':
        w.ev.dump(OUT/'supplement_progress.json',dict(complete=len(records)==len(indices),
                  completed_files=len(records),expected_files=len(indices),
                  selected=sum(r['selected'] for r in records)))


if __name__=='__main__':main()
