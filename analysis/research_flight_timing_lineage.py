#!/usr/bin/env python3
"""Audit stored timing branches and exact event identity, without RICH/truth."""
import argparse
import json
import hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import flight_timing_redundancy as m
import research_flight_timing_redundancy as r
ROOT=r.ROOT
OUT=ROOT/'results/flight_calibration_research/timing_lineage'


def one(record):
    source=int(record['source'])
    with np.load(str(ROOT/'results/time_local_beta/parts'/('%05d'%source)/'measured.npz'),allow_pickle=False) as f:
        wanted=set(zip(f['run'].astype(int),f['event'].astype(int)))
    with uproot.open(record['path']) as f:
        tree=f['amstreea'];a=tree.arrays(['run','event','tof_tl','tof_leng','tof_tl_NewCalib'],library='np')
    ix=np.array([i for i,k in enumerate(zip(a['run'].astype(int),a['event'].astype(int))) if k in wanted],int)
    a={k:v[ix] for k,v in a.items()};a['source']=np.full(len(ix),source,int)
    return a


def prepare(workers,full):
    OUT.mkdir(exist_ok=True)
    if full:
        sources=json.loads((r.OUT/'protocol.json').read_text())['sources']
        records=json.loads((ROOT/'results/tofcalib_validation/preparation.json').read_text())['records']
        records=[a for a in records if int(a['source']) in sources]
    else:records=json.loads((ROOT/'results/flight_calibration_research/input_inventory.json').read_text())['records']
    records=sorted(records,key=lambda a:int(a['source']))
    r.write(OUT/'protocol.json',dict(sources=[int(x['source']) for x in records],RICH_read=False,MC_truth_read=False,
        comparison='All stored NewCalib method/threshold slots; no outcome-based slot selection',
        baseline='tof_tl inherited from MIT input for bTofCalibDst=false',
        source_proof='TofCalib/src/isotope_selection.cxx CloneTree(0) lines229-231; tof_tl only added/replaced for bTofCalibDst=true at331-333 and1046; new times written separately at1249',
        producer_source_sha256=r.digest(ROOT.parent/'TofCalib/src/isotope_selection.cxx')))
    pieces=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i,a in enumerate(pool.map(one,records),1):
            pieces.append(a)
            if i%30==0 or i==len(records):print('TIMING_LINEAGE_READ',i,'/',len(records),flush=True)
    data={k:np.concatenate([a[k] for a in pieces]) for k in pieces[0]}
    ids=list(zip(data['run'],data['event']))
    if len(set(ids))!=len(ids):raise ValueError('Duplicate lineage events')
    np.savez_compressed(str(OUT/'measured.npz'),**data)


def evaluate():
    with np.load(str(OUT/'measured.npz'),allow_pickle=False) as f:data={k:f[k] for k in f.files}
    x=data['tof_leng'].astype(float)/29.9792458;t=data['tof_tl'].astype(float);n=len(x)
    p,w=m.projector(x);slots=data['tof_tl_NewCalib'].reshape(n,-1,4).astype(float)
    base_u=np.sum(w*t,axis=1);base_res=np.einsum('nij,nj->ni',p,t)
    rows=[]
    for j in range(slots.shape[1]):
        tt=slots[:,j];valid=np.isfinite(tt).all(axis=1)&(tt!=0).all(axis=1)&(np.abs(tt)<1000).all(axis=1)
        delta=tt[valid]-t[valid];delta-=delta.mean(axis=1)[:,None]
        u=np.sum(w[valid]*tt[valid],axis=1);res=np.einsum('nij,nj->ni',p[valid],tt[valid])
        rows.append(dict(slot=j,method=j//3,threshold=['LT','HT','SHT'][j%3],valid_events=int(valid.sum()),
            centered_time_change_ns=r.stats(delta.ravel()),u_change=r.stats(u-base_u[valid]),
            baseline_contrasts=[r.stats(base_res[valid,l]) for l in range(4)],new_contrasts=[r.stats(res[:,l]) for l in range(4)],
            beta_resolution_measured=False))
    duplicated=[]
    for j in range(3):
        if slots.shape[1]>j+3:
            duplicated.append(dict(threshold=j,equal_fraction=float(np.mean(slots[:,j]==slots[:,j+3]))))
    result=dict(complete=True,events=n,runs=len(np.unique(data['run'])),stored_shape=list(data['tof_tl_NewCalib'].shape),rows=rows,
        between_method_equality=duplicated,claim='Stored tof_tl is not the separately stored NewCalib time array; their common offsets alone do not explain all differences',
        no_RICH_or_truth=True,cache_sha256=r.digest(OUT/'measured.npz'))
    r.write(OUT/'audit.json',result)
    print(json.dumps(r.native(result),indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','evaluate','all']);p.add_argument('--full',action='store_true');p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    if a.stage in ['prepare','all']:prepare(a.workers,a.full)
    if a.stage in ['evaluate','all']:evaluate()
