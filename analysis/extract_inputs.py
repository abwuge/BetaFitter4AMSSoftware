#!/usr/bin/env python3
"""Extract aligned train/validation arrays using the existing system uproot environment."""
import argparse, json, os
from collections import defaultdict
import numpy as np
import uproot

CHARGES=(1,2,6,8)
BRANCHES=['run','event','mch','mmass','mevmom1','tof_betah','tof_tl','tof_leng','tof_barid','tof_pos','tof_dir','tof_edep','tof_bpatt','tk_edep','tk_pos','tk_q','tk_qin','tk_rigidity']

def parse_sources(path):
    out=[]
    for line in open(path):
        line=line.strip()
        if not line or line.startswith('#'): continue
        role,z,p=line.split(); out.append((role,int(z),p))
    return out

def load_raw(path):
    with uproot.open(path) as f:
        a=f['amstreea'].arrays(BRANCHES,library='np')
    charge=(np.where(a['tk_qin'][:,0,2]<2.5,a['tk_q'][:,1],a['tk_qin'][:,0,2])+0.5).astype(np.int32)
    keep=(a['mevmom1'][:,10]!=-1000)&np.isin(charge,CHARGES)&(a['mmass']>0)
    keep &= np.all(np.isfinite(a['tk_edep'][:,1:8])&(a['tk_edep'][:,1:8]>0),axis=1)
    keep &= np.all(np.isfinite(a['tk_pos'][:,1:8]),axis=(1,2))
    keep &= (a['tk_pos'][:,4,2]>0)&(a['tk_pos'][:,5,2]<0)
    p=a['mevmom1'][:,10].astype(np.float64);m=a['mmass'].astype(np.float64)
    mc=p/np.sqrt(p*p+m*m)
    out={}
    for k in ['run','event','mch','mmass','tof_betah','tof_tl','tof_leng','tof_barid','tof_edep','tof_bpatt','tk_edep','tk_rigidity']:
        out[k]=a[k][keep]
    out['position']=a['tof_pos'][keep,:,:2]
    out['direction']=a['tof_dir'][keep,0]
    out['charge']=charge[keep]
    out['mcBetaRaw']=mc[keep]
    return out

def append(store,d):
    for k,v in d.items(): store[k].append(v)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sources',default='config/sources.conf');ap.add_argument('--beta-fitter',default='/nas01/home/hxwu/AMSSoft/betaFitter');ap.add_argument('--output',default='results/comparison/cache');args=ap.parse_args()
    if not os.path.isdir(args.output): os.makedirs(args.output)
    data={(r,z):defaultdict(list) for r in ('train','validation') for z in CHARGES}
    validation_index=0; audit=[]
    records=parse_sources(args.sources)
    for i,(role,z,path) in enumerate(records,1):
        d=load_raw(path); d['source_id']=np.full(len(d['charge']),i,dtype=np.int16)
        if role=='validation':
            validation_index+=1
            part=os.path.join(args.beta_fitter,'results/hybrid_physical_20260913/validation/parts','part_%03d_Z%d.root'%(validation_index,z))
            with uproot.open(part) as f: b=f['betaTree'].arrays(library='np')
            assert len(b['Z'])==len(d['charge']), (part,len(b['Z']),len(d['charge']))
            assert np.array_equal(b['Z'],d['charge']), part
            assert np.array_equal(b['mcCharge'],d['mch'].astype(np.int32)), part
            assert np.allclose(b['linearBeta'],d['tof_betah'],rtol=0,atol=0,equal_nan=True), part
            assert np.allclose(b['innerRigidity'],d['tk_rigidity'],rtol=0,atol=0,equal_nan=True), part
            for k in ['mcBeta','linearBeta','nonlinearBeta','hybridBeta']: d[k]=b[k]
        append(data[(role,z)],d)
        audit.append(dict(role=role,Z=z,path=path,selected_rows=len(d['charge'])))
        if i%16==0: print('read %d/%d source files'%(i,len(records)),flush=True)
    for key,store in sorted(data.items()):
        role,z=key; out={k:np.concatenate(v) for k,v in store.items()}
        dest=os.path.join(args.output,'%s_Z%d.npz'%(role,z));np.savez_compressed(dest,**out)
        print('%s Z=%d rows=%d'%(role,z,len(out['charge'])),flush=True)
    with open(os.path.join(args.output,'extract_audit.json'),'w') as f: json.dump(audit,f,indent=2)
    print('all validation partitions aligned',flush=True)
if __name__=='__main__': main()
