#!/usr/bin/env python3
"""Evaluate the unchanged flight model on supplemented NEWTOF measurements."""
import argparse
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import flight_beta_wide as w
from prepare_tofcalib_validation import OUT

# Finer evaluation bins; all original R/beta boundaries remain included.
R_EDGES = np.unique(np.concatenate([np.geomspace(lo,hi,5) for lo,hi in zip(w.R_EDGES[:-1],w.R_EDGES[1:])]))
BETA_EDGES = np.array([.78,.80,.82,.84,.85,.87,.89,.90,.92,.94,.95,.955,.96,.965,.97,.975,.98,.985,.99,.9925,.995,.9965,.998,.999,1.,1.001,1.002,1.004,1.008,1.02])


MATERIAL = None

def initialize_predictor():
    global MATERIAL
    MATERIAL = w.rt.MaterialMap()

def predict_part(dest):
    with np.load(str(dest/'joined.npz'), allow_pickle=False) as f:
        a = {k:f[k] for k in f.files}
    w.set_edges(w.R_EDGES)
    d = w.rt.center_features(a)
    corrected = w.predict(d)
    geometry = MATERIAL.build_columns(a['tk_pos'],a['tk_dir'])
    delays=[]
    valid=geometry['valid'].copy()
    for mass in (11.17793,12.11221):
        nominal=np.sqrt(1+(mass/(6*d['R']))**2)
        physical=w.rt.propagate(nominal,d['R'],geometry,mass_override=mass)
        delays.append(1/physical['beta_planes'][:,-1]-nominal)
        valid &= physical['valid']
    delay=delays[0]
    raw, new = 1/(d['u']+delay), 1/(corrected+delay)
    reference = a['rich_beta'][:,0]
    mask = w.rt.rich_quality(a) & (a['rich_q'][:,1]>4) & (a['rich_q'][:,1]<8)
    mask &= (a['rich_npe'][:,2]>0) & valid
    mask &= np.isfinite(raw) & np.isfinite(new) & (raw>0) & (new>0)
    return dict(run=a['run'],event=a['event'],source=a['source'],
                   R=d['R'],rich_beta=reference,raw_u=d['u'],corrected_u=corrected,
                   raw_beta_rich=raw,new_beta_rich=new,reference_mask=mask,NaF=a['rich_NaF'],
                   delay_u_C12=delays[0],delay_u_C13=delays[1])


def predict(workers=8):
    w.check_freeze()
    config = json.loads((OUT/'selection.json').read_text())
    parts = []
    for source in config['file_indices']:
        dest = OUT/'parts'/('%05d'%source)
        if not (dest/'supplement.json').exists():
            raise RuntimeError('Supplement incomplete for source %d'%source)
        record = json.loads((dest/'supplement.json').read_text())
        if not record.get('complete'):
            raise RuntimeError('Incomplete supplement')
        if record['selected']:
            parts.append(dest)
    records = []
    with mp.Pool(workers, initializer=initialize_predictor) as pool:
        for i, record in enumerate(pool.imap(predict_part, parts), 1):
            records.append(record)
            if i%100==0 or i==len(parts):
                print('Predicted',i,'/',len(parts),flush=True)
    if not records:
        raise RuntimeError('No supplemented selected events')
    p = w.ev.concatenate(records)
    ids = np.rec.fromarrays([p['run'],p['event']],names='run,event')
    unique, first, inverse = np.unique(ids,return_index=True,return_inverse=True)
    duplicates = len(ids)-len(unique)
    if duplicates:
        # Duplicate events must have identical measurements/predictions; an
        # inconsistent event cannot be silently discarded from the evaluation.
        for key in p:
            if key=='source':continue
            if not np.allclose(p[key],p[key][first[inverse]],equal_nan=True,rtol=0,atol=0):
                raise RuntimeError('Conflicting duplicate event values: '+key)
        p = w.ev.sub(p,first)
    with np.load(str(OUT/'training_runs.npz'),allow_pickle=False) as f:
        if np.intersect1d(p['run'],f['run']).size:
            raise RuntimeError('Training run leakage')
    np.savez_compressed(str(OUT/'predictions.npz'),**p)
    w.ev.dump(OUT/'prediction_audit.json',dict(complete=True,selected=len(p['run']),
              reference_n=int(p['reference_mask'].sum()),identical_duplicates_removed=duplicates,
              actual_training_run_overlap=0,model_freeze_sha256=w.ev.digest(w.DEST/'freeze.json'),
              prediction_sha256=w.ev.digest(OUT/'predictions.npz')))


def evaluate():
    w.check_freeze()
    with np.load(str(OUT/'predictions.npz'),allow_pickle=False) as f:
        p={k:f[k] for k in f.files}
    mask=p['reference_mask']
    raw,new,reference=(p[k] for k in ('raw_beta_rich','new_beta_rich','rich_beta'))
    rows={}
    for axis,values,edges in [('rigidity',p['R'],R_EDGES),('beta',reference,BETA_EDGES)]:
        records=[]
        for low,high in zip(edges[:-1],edges[1:]):
            select=mask & (values>=low) & (values<high)
            n=int(select.sum());row=dict(low=low,high=high,n=n)
            if n>=2:
                row.update(w.ev.compare(raw[select],new[select],reference[select]))
                row['x']=float(np.median(values[select]))
                # Actual run is the statistical cluster: source files can
                # contain multiple runs and multiple files can share a run.
                row['bootstrap']=w.ev.bootstrap(raw[select],new[select],reference[select],p['run'][select])
            records.append(row)
            print(axis,low,high,'N',n,'gain',row.get('gain_beta_pct'),flush=True)
        if sum(r['n'] for r in records)!=int(mask.sum()):
            raise RuntimeError('Incomplete bin coverage')
        rows[axis]=records
    w.ev.dump(OUT/'results.json',dict(complete=True,selected=len(mask),reference_n=int(mask.sum()),
              rows=rows,rigidity_edges_GV=R_EDGES,beta_edges=BETA_EDGES,
              bootstrap_cluster='actual run',RICH_used_for_training=False,retrained=False,
              retained_raw_beta_above_one=int((mask & (raw>1)).sum()),
              retained_new_beta_above_one=int((mask & (new>1)).sum()),
              overall=w.ev.compare(raw[mask],new[mask],reference[mask]),
              source_selection=json.loads((OUT/'selection.json').read_text())))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['predict','evaluate'])
    parser.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    if args.stage=='predict':predict(args.workers)
    else:evaluate()


if __name__=='__main__':main()
