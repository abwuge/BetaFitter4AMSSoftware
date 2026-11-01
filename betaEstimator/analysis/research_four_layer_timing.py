#!/usr/bin/env python3
"""Version 3: respect the native used-hit count before fitting four times.

This is an explicit four-used-layer domain, not an outcome-selected residual
cut. It inherits the original BetaH selection and is not an unbiased noise
sample without a selection model. V1 and V2 remain available as failure audits.
"""
import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
import research_flight_timing_redundancy as r
import research_timing_mc_closure as mc
ROOT=r.ROOT;BASE=ROOT/'results/flight_calibration_research';OUT=BASE/'used_layer_audit'


def read_quality(record):
    with uproot.open(record['path']) as f:
        a=f['amstreea'].arrays(['run','event','tof_hsumh','tof_hsumhu'],library='np')
    a['source']=np.full(len(a['run']),int(record['source']),int)
    return a


def audit_subset(data,quality):
    lookup={key:i for i,key in enumerate(zip(quality['run'],quality['event']))}
    if len(lookup)!=len(quality['run']):raise ValueError('Duplicate native quality identity')
    ix=np.array([lookup[key] for key in zip(data['run'],data['event'])],int)
    nused=quality['tof_hsumhu'][ix];nhit=quality['tof_hsumh'][ix]
    p,w=r.m.projector(data['x']);residual=np.einsum('nij,nj->ni',p,data['t']);norm=np.linalg.norm(residual,axis=1)
    counts=[]
    for used in sorted(np.unique(nused)):
        take=nused==used
        counts.append(dict(native_used_layers=int(used),events=int(take.sum()),contrast_norm_above_1ns=int((norm[take]>1).sum()),
            contrast_norm_above_10ns=int((norm[take]>10).sum()),layer_contrasts=[r.stats(residual[take,l]) for l in range(4)]))
    return nused,nhit,counts


def flight_prepare():
    OUT.mkdir(exist_ok=True)
    wanted=json.loads((BASE/'timing_redundancy/protocol.json').read_text())['sources']
    records=json.loads((ROOT/'results/tofcalib_validation/preparation.json').read_text())['records']
    records=sorted([a for a in records if int(a['source']) in wanted],key=lambda a:int(a['source']))
    pieces=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i,a in enumerate(pool.map(read_quality,records),1):
            pieces.append(a)
            if i%40==0 or i==len(records):print('USED_LAYER_AUDIT_READ',i,'/',len(records),flush=True)
    quality={k:np.concatenate([a[k] for a in pieces]) for k in pieces[0]}
    # Keep only the existing common detector-domain event identities.
    with np.load(str(BASE/'valid_timing_native/measured.npz'),allow_pickle=False) as f:data={k:f[k] for k in f.files}
    nused,nhit,counts=audit_subset(data,quality)
    np.savez_compressed(str(OUT/'flight_quality.npz'),run=data['run'],event=data['event'],tof_hsumh=nhit,tof_hsumhu=nused)
    keep=(nused==4)&(nhit==4)
    r.write(OUT/'flight_audit.json',dict(parent_events=len(keep),kept_four_used_events=int(keep.sum()),removed_events=int((~keep).sum()),counts_by_used_layers=counts,
        source_evidence='amsd69n/analysis_amsd69n.C stores GetSumHit/GetUseHit; include/root.h IsGoodBeta requires track and GetUseHit()==4; IsBetaUseHL is distinct from TestExistHL',
        rule='Require original native tof_hsumhu==4 and tof_hsumh==4 for a four-layer fit; no new time-residual, beta or truth cut',
        selection_caveat='BetaH used-layer choice depends on observed detector signals; this restricted domain does not independently identify unconditional timing-noise covariance',
        RICH_read=False,MC_read=False))
    for name in ['native','newcalib']:
        source=BASE/('valid_timing_'+name);dest=BASE/('four_layer_timing_'+name);dest.mkdir(exist_ok=True)
        with np.load(str(source/'measured.npz'),allow_pickle=False) as f:a={k:f[k] for k in f.files}
        if not np.array_equal(a['run'],data['run']) or not np.array_equal(a['event'],data['event']):raise ValueError('Baseline identity mismatch')
        a=r.subset(a,keep);np.savez_compressed(str(dest/'measured.npz'),**a)
        r.write(dest/'protocol.json',dict(version=3,baseline=name,parent_protocol_sha256=r.digest(source/'protocol.json'),
            quality_audit_sha256=r.digest(OUT/'flight_audit.json'),parent_events=len(keep),events=int(keep.sum()),
            criterion='Four original native used layers and all four existing finite non-sentinel times; shared native/NewCalib event set',
            model='Unchanged low-dimensional response and shared correction family; all parameters refit from measured values',
            selection_caveat='This is a restricted native-four-used-layer domain, not proof of unbiased noise selection',
            development_status='Previously inspected run cohort; not a fresh blind RICH test',RICH_read=False,MC_coefficients_used=False))
    print('FLIGHT_FOUR_USED_READY',int(keep.sum()),'of',len(keep),flush=True)


def mc_prepare():
    source=BASE/'timing_mc_closure';dest=BASE/'four_layer_timing_mc_closure';dest.mkdir(exist_ok=True)
    protocol=json.loads((source/'protocol.json').read_text())
    with np.load(str(source/'measured.npz'),allow_pickle=False) as f:data={k:f[k] for k in f.files}
    nused=np.full(len(data['run']),-1,int);nhit=nused.copy()
    for record in protocol['records']:
        ix=np.flatnonzero(data['source']==record['source'])
        if not len(ix):continue
        with uproot.open(record['path']) as f:a=f['amstreea'].arrays(['run','event','tof_hsumh','tof_hsumhu'],entry_stop=3000,library='np')
        entries=data['entry'][ix]
        if not np.array_equal(a['event'][entries],data['event'][ix]) or not np.array_equal(a['run'][entries],data['run'][ix]):raise ValueError('MC quality identity mismatch')
        nused[ix]=a['tof_hsumhu'][entries];nhit[ix]=a['tof_hsumh'][entries]
    keep=(nused==4)&(nhit==4)
    protocol.update(version=3,parent_events=len(keep),events=int(keep.sum()),removed_events=int((~keep).sum()),
        quality='Original native four fired and four used layers; no truth read for selection or training',
        counts_by_used_layers={str(int(i)):int((nused==i).sum()) for i in np.unique(nused)},
        development_status='Retrospective detector-validity repair on prior MC development cohort; not a new blind sample',
        quality_fields=['run','event','tof_hsumh','tof_hsumhu'])
    r.write(dest/'protocol.json',protocol)
    data=r.subset(data,keep);np.savez_compressed(str(dest/'measured.npz'),**data)
    print('MC_FOUR_USED_READY',int(keep.sum()),'of',len(keep),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','native','newcalib','mc']);p.add_argument('--replicates',type=int,default=32);a=p.parse_args()
    if a.stage=='prepare':flight_prepare();mc_prepare()
    elif a.stage in ['native','newcalib']:
        r.OUT=BASE/('four_layer_timing_'+a.stage);r.freeze();r.evaluate();r.bootstrap(a.replicates)
    elif a.stage=='mc':
        mc.OUT=BASE/'four_layer_timing_mc_closure';mc.fit();mc.evaluate();mc.bootstrap()
