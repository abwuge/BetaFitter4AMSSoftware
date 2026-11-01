#!/usr/bin/env python3
"""In-situ MC terminal audit. No production products are changed."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,json,hashlib,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,uproot
ROOT=Path('/nas01/home/hxwu/AMSSoft/betaEstimator')
OUT=ROOT/'results/terminal_ceiling_20261015'
C=29.9792458
BASE=['run','event','mmass','mpar','mch','mevmom1','mevcoo1','mevdir1','mtrmom','mtrpri','mtrpar','mtrcoo1','tof_tl','tof_leng','tof_pos','tof_barid','tof_edep','tk_edep','tk_ql','tk_q','tk_qin','tk_pos','tk_dir','tk_rigidityk','tof_bpatt','mtofdep','mtofpri','mtofbar','tof_betah','tof_betahmc']
OPTIONAL=['tof_Aadc','tof_adca','tof_rtime','tof_etl','tof_tl_HT','tof_tl_SHT','tk_chis1','tk_erigidity1','trd_dedx','trd_qk','ecal_en','rich_beta','rich_good','rich_clean','rich_NaF']

def read_one(record):
 role,z,path,sid=record
 with uproot.open(path) as file:
  tree=file['amstreea'];keys=set(tree.keys());fields=[k for k in BASE+OPTIONAL if k in keys]
  a=tree.arrays(fields,library='np')
 n=len(a['run']);t=a['tof_tl'].astype(float);x=a['tof_leng'].astype(float)/C
 xc=x-x.mean(1)[:,None];den=np.sum(xc*xc,axis=1)
 with np.errstate(all='ignore'):
  u=np.sum(xc*t,axis=1)/den
  q=np.where(a['tk_qin'][:,0,2]<2.5,a['tk_q'][:,1],a['tk_qin'][:,0,2])
  mm=a['mmass'].astype(float);pp=a['mtrmom'][:,4].astype(float)
  ut=np.hypot(1,mm/pp)
 good=np.isfinite(t).all(1)&np.isfinite(x).all(1)&(den>1)&np.isfinite(u)&(u>0)
 good&=(a['tof_bpatt']==4444)&np.isfinite(q)&(q>=z-.5)&(q<z+.5)
 good&=np.isfinite(a['tof_pos']).all((1,2))&np.isfinite(a['tk_pos'][:,1:8]).all((1,2))
 good&=(a['tk_pos'][:,4,2]>0)&(a['tk_pos'][:,5,2]<0)
 good&=np.isfinite(a['tof_edep']).all(1)&(a['tof_edep']>0).all(1)
 good&=np.isfinite(a['tk_edep'][:,1:8]).all(1)&(a['tk_edep'][:,1:8]>0).all(1)
 good&=((a['tof_barid']>=0)&(a['tof_barid']<10)).all(1)
 measured=int(good.sum())
 truth=(a['mch']==z)&(a['mtrpri'][:,4]==1)&(a['mtrpar'][:,4]==a['mpar'])
 truth&=np.isfinite(pp)&(pp>0)&np.isfinite(mm)&(mm>0)&np.isfinite(ut)&(ut>1)&(ut<=4)
 truth&=np.abs(a['mtrcoo1'][:,4,2]-1.698)<.04
 good&=truth
 b={k:v[good] for k,v in a.items()}
 b.update(u=u[good],u_true=ut[good],source=np.full(good.sum(),sid,np.int32),role=np.full(good.sum(),role,np.int8),Z=np.full(good.sum(),z,np.int8))
 note=dict(source=sid,role=role,Z=z,path=path,raw=n,measured=measured,accepted=int(good.sum()),missing=sorted(set(BASE+OPTIONAL)-keys),file_bytes=Path(path).stat().st_size)
 print('SOURCE',z,role,sid,n,measured,int(good.sum()),flush=True)
 return b,note

def run(z,workers):
 rows=[l.split() for l in (ROOT.parent/'betaFitter/config/hybrid_material_fine.sources').read_text().splitlines() if l and not l.startswith('#')]
 rec=[(0 if r[0]=='train' else 1,int(r[1]),r[2],i) for i,r in enumerate(rows) if int(r[1])==z]
 pieces=[];notes=[]
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for a,n in pool.map(read_one,rec):pieces.append(a);notes.append(n)
 keys=set.intersection(*(set(a) for a in pieces));a={k:np.concatenate([b[k] for b in pieces]) for k in sorted(keys)}
 ids=a['run'].astype(np.uint64)*np.uint64(2**32)+a['event'].astype(np.uint64)
 if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate identities')
 if np.intersect1d(a['run'][a['role']==0],a['run'][a['role']==1]).size:raise ValueError('Run overlap')
 cp=uproot.open(str(ROOT/('results/comparison/comparison_Z%d.root'%z)))['comparison'].arrays(library='np')
 ck=cp['run'].astype(np.uint64)*np.uint64(2**32)+cp['event'].astype(np.uint64)
 order=np.argsort(ck);j=np.searchsorted(ck[order],ids);jj=np.minimum(j,len(ck)-1);found=(j<len(ck))&(ck[order[jj]]==ids)
 a['previous_comparison']=found
 for k in ['mcBeta','hybridBeta','newBeta','linearBeta','nonlinearBeta']:
  a['previous_'+k]=np.where(found,cp[k][order[jj]],np.nan)
 np.savez_compressed(str(OUT/('broad_Z%d.npz'%z)),**a)
 stat=dict(Z=z,events=len(ids),train=int((a['role']==0).sum()),test=int((a['role']==1).sum()),previous_comparison=int(found.sum()),columns={k:list(v.shape) for k,v in a.items()},sources=notes,truth='Primary mtrmom[4], mass mmass; mtrcoo1[4] L5 plane; no tof_betahmc truth',created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
 (OUT/('broad_Z%d_manifest.json'%z)).write_text(json.dumps(stat,indent=2))
 print('COMPLETE',z,len(ids),stat['train'],stat['test'],stat['previous_comparison'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--z',type=int,required=True);p.add_argument('--workers',type=int,default=3);args=p.parse_args();run(args.z,args.workers)
