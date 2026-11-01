#!/usr/bin/env python3
"""Measured flight reference cache; never opens MC or truth branches."""
import argparse, hashlib, json, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
import numpy as np
import uproot
from extract_inputs import CHARGES, BRANCHES, detector_selection, parse_sources
from timing_model import clean
REFERENCE_BRANCHES = ['irich','itrtrack','rich_itrtrack','rich_beta','rich_ebeta','rich_ebeta1',
 'rich_good','rich_clean','rich_NaF','rich_hit','rich_used','rich_pb',
 'rich_BetaConsistency','rich_pos','rich_tile','rich_distb','tk_rigidity']
ALL_BRANCHES = list(dict.fromkeys(BRANCHES+REFERENCE_BRANCHES))
def process(record):
    sid,(role,_,path)=record
    with uproot.open(path) as f:
        a=f['amstreea'].arrays(ALL_BRANCHES,library='np')
    charge,keep=detector_selection(a)
    out={k:a[k][keep] for k in ['run','event','tof_betah','tof_tl','tof_leng','tof_barid','tof_edep','tof_bpatt','tk_edep']+REFERENCE_BRANCHES}
    out['position']=a['tof_pos'][keep,:,:2]
    out['direction']=a['tof_dir'][keep,0]
    out['charge']=charge[keep]
    out['source_entry']=np.flatnonzero(keep)
    out['source_id']=np.full(int(keep.sum()),sid,dtype=np.int32)
    ok=clean(out)
    rng=np.random.RandomState(1981+sid)
    chunks={}
    stats={}
    for z in CHARGES:
        idx=np.flatnonzero(ok&(out['charge']==z))
        total=len(idx)
        # RICH availability only, not TOF-RICH agreement. Keep saturation anchors too.
        broad=(out['irich'][idx]>=0)&(out['rich_beta'][idx,0]>.6)
        saturation=np.isfinite(out['tk_rigidity'][idx])&(out['tk_rigidity'][idx]>100)
        idx=idx[broad|saturation]
        eligible=len(idx)
        if z<3 and len(idx)>3000:
            idx=np.sort(rng.choice(idx,3000,replace=False))
        chunks[z]={k:v[idx] for k,v in out.items()}
        stats[str(z)]={'clean_pool':total,'reference_pool':eligible,'stored':len(idx)}
    return sid,role,path,chunks,{'source_id':sid,'role':role,'path':path,'raw':len(charge),'charges':stats}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--sources',default='config/data_sources.conf')
    ap.add_argument('--output',required=True)
    ap.add_argument('--limit',type=int,default=0)
    args=ap.parse_args()
    dest=Path(args.output);dest.mkdir(parents=True,exist_ok=True)
    records=list(enumerate(parse_sources(args.sources),1))
    if args.limit: records=records[:args.limit]
    stores=defaultdict(lambda:defaultdict(list));audit=[];run_roles={}
    with ProcessPoolExecutor(max_workers=3) as pool:
        for sid,role,path,chunks,row in pool.map(process,records):
            for z,d in chunks.items():
                for run in np.unique(d['run']):
                    if run_roles.setdefault(int(run),role)!=role:
                        raise ValueError('Run overlap')
                for k,v in d.items():stores[(role,z)][k].append(v)
            audit.append(row)
            if len(audit)%20==0: print('reference extraction %d/%d'%(len(audit),len(records)),flush=True)
    counts={}
    for (role,z),s in sorted(stores.items()):
        d={k:np.concatenate(v) for k,v in s.items()}
        ids=np.rec.fromarrays([d['run'],d['event']])
        if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate events')
        name='%s_Z%d.npz'%(role,z)
        np.savez_compressed(str(dest/name),**d)
        counts[name]=len(d['run'])
        print(name,counts[name],flush=True)
    (dest/'reference_audit.json').write_text(json.dumps({'branches':ALL_BRANCHES,'truth_access':False,
        'manifest_sha256':hashlib.sha256(Path(args.sources).read_bytes()).hexdigest(),
        'counts':counts,'sources':audit,'unique_events':True,'disjoint_runs':True},indent=2))
if __name__=='__main__':main()
