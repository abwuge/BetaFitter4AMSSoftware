#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,csv,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,uproot
from audit_common import HERE,write_json,C
R=HERE.parent
P=np.load(str(R/'c_joint_beta_rigidity_closure/expanded_validation/predictions.npz'))
MAN=list(csv.DictReader((R/'c_joint_beta_rigidity_closure/expanded_validation/manifest.tsv').open(),delimiter='\t'))
SPLIT={(r['isotope'],int(r['source'])):r['role'] for r in csv.DictReader((R/'c_mass_toy_oracle/split_manifest.tsv').open(),delimiter='\t')}
SHAPES={'tof_tl':(4,),'tof_leng':(4,),'tof_pos':(4,3),'tof_barid':(4,),'tof_edep':(4,),'tof_rtime':(4,2),'tof_Aadc':(4,2),'tof_adca':(4,2),'mtrmom':(9,),'mtrpri':(9,),'mtrpar':(9,),'mtrcoo1':(9,3),'mmass':(),'mpar':(),'mch':(),'mevmom1':(21,),'tk_edep':(9,),'tk_ql':(9,),'tk_rigidityk':(3,3),'tk_erigidity1':(3,3),'tk_chis1':(3,3,3,2),'tk_pos':(9,3),'mtofdep':(4,),'rich_beta':(2,)}

def one(r):
 sid=int(r['source']);label=0 if r['isotope']=='C12' else 1;ii=np.flatnonzero((P['source']==sid)&(P['label']==label));ids=P['run'][ii].astype(np.uint64)*np.uint64(2**32)+P['event'][ii].astype(np.uint64);order=np.argsort(ids);sid_sorted=ids[order]
 if not len(ii):return None,dict(source=sid,events=0)
 pieces=[];missing=[];shapes={}
 with uproot.open(r['path']) as file:
  tree=file['amstreea'];keys=set(tree.keys());fields=['run','event']+[k for k in SHAPES if k in keys]
  missing=sorted(set(SHAPES)-keys)
  for b in tree.iterate(fields,step_size=50000,library='np'):
   kk=b['run'].astype(np.uint64)*np.uint64(2**32)+b['event'].astype(np.uint64);j=np.searchsorted(sid_sorted,kk);jj=np.minimum(j,len(ids)-1);take=(j<len(ids))&(sid_sorted[jj]==kk)
   if not take.any():continue
   x={k:b[k][take] for k in b};x['index']=ii[order[j[take]]];pieces.append(x)
 if not pieces:raise RuntimeError('No targets '+r['path'])
 a={k:np.concatenate([b[k] for b in pieces]) for k in pieces[0]};n=len(a['index'])
 if n!=len(ii) or len(np.unique(a['index']))!=n:raise RuntimeError('Missing or duplicate target source %d'%sid)
 for k in SHAPES:
  if k not in a:a[k]=np.full((n,)+SHAPES[k],np.nan,np.float32)
 idx=a['index'];a['label']=P['label'][idx];a['source']=np.full(n,sid+256*label,np.int32);a['role']=np.full(n,{'train':0,'valid':2,'test':1}[SPLIT[(r['isotope'],sid)]],np.int8);a['Z']=np.full(n,6,np.int8)
 mass_fallback=not np.isfinite(a['mmass']).all()
 mass=np.where(np.isfinite(a['mmass'])&(a['mmass']>0),a['mmass'],np.where(a['label']==0,11.17793,12.11221))
 with np.errstate(all='ignore'):
  ut=np.hypot(1,mass/a['mtrmom'][:,4]);xc=a['tof_leng'].astype(float)/C;xc-=xc.mean(1)[:,None];u=np.sum(xc*a['tof_tl'],axis=1)/np.sum(xc*xc,axis=1)
 a['u_true']=ut;a['u']=u;a['true_mass']=mass
 a['primary_valid']=(a['mtrpri'][:,4]==1)&(a['mtrpar'][:,4]==a['mpar'])&(a['mtrmom'][:,4]>0)&np.isfinite(ut)&(ut>1)&(np.abs(a['mtrcoo1'][:,4,2]-1.698)<.04)
 for k in ['R','R_true','u_true','u_raw','u_beta_only_frozen','R_beta_only_frozen','u_joint_anchored','R_joint_anchored']:a['cached_'+k]=P[k][idx]
 a['cached_u']=P['u'][idx]
 note=dict(source=sid+256*label,original_source=sid,path=r['path'],isotope=r['isotope'],events=n,primary_valid=int(a['primary_valid'].sum()),missing=missing,mass_fallback=mass_fallback,max_u_difference=float(np.nanmax(np.abs(u-P['u_raw'][idx]))))
 print('MIXED_SOURCE',sid,n,note['primary_valid'],note['max_u_difference'],flush=True)
 return a,note

with ThreadPoolExecutor(max_workers=8) as pool:results=list(pool.map(one,MAN))
pieces=[a for a,n in results if a is not None];keys=set.intersection(*(set(a) for a in pieces));a={k:np.concatenate([b[k] for b in pieces]) for k in keys};ix=np.argsort(a['index']);a={k:v[ix] for k,v in a.items()}
if not np.array_equal(a['index'],np.arange(len(P['run']))):raise RuntimeError('Incomplete mixed extraction')
np.savez_compressed(str(HERE/'mixed.npz'),**a)
write_json(HERE/'mixed_manifest.json',dict(events=len(a['run']),primary_valid=int(a['primary_valid'].sum()),sources=[n for _,n in results],mass_caveat='C13 mmass is absent: retain the existing nominal MC mass 12.11221 GeV; no native C13 generator mass validation is claimed. C12 uses native-exported mmass.',columns={k:list(v.shape) for k,v in a.items()}))
print('MIXED_COMPLETE',len(a['run']),int(a['primary_valid'].sum()),flush=True)
