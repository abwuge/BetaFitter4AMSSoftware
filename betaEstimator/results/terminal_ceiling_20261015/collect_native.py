#!/usr/bin/env python3
"""Primary-identity native truth collection. No error-based event selection."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
from pathlib import Path
import numpy as np,uproot
from audit_common import HERE,C,write_json,metric

b=np.load(str(HERE/'broad_Z6.npz'));n=len(b['run'])
key=b['run'].astype(np.uint64)*np.uint64(2**32)+b['event'].astype(np.uint64);order=np.argsort(key)
a=dict(run=b['run'],event=b['event'],source=b['source'],role=b['role'],found=np.zeros(n,bool),association=np.zeros(n,bool),mass=np.full(n,np.nan),u_true=np.full(n,np.nan),acceptance_count=np.zeros((n,21),np.int16),acceptance_p=np.full((n,21),np.nan),acceptance_pos=np.full((n,21,3),np.nan),acceptance_dir=np.full((n,21,3),np.nan),tof_primary_count=np.zeros((n,4),np.int16),tof_count=np.zeros((n,4),np.int16),tof_time=np.full((n,4),np.nan),tof_first_any_time=np.full((n,4),np.nan),tof_beta=np.full((n,4),np.nan),tof_edep=np.full((n,4),np.nan),tof_total_edep=np.zeros((n,4)),tof_step=np.full((n,4),np.nan),tof_pos=np.full((n,4,3),np.nan),tof_dir_z=np.full((n,4),np.nan),native_tk_p=np.full(n,np.nan),native_tk_count=np.zeros(n,np.int16))
plan=json.loads((HERE/'native_plan.json').read_text())
start=time.monotonic()
for record in plan['records']:
 r=uproot.open(record['output'])['reference'].arrays(library='np')
 ids=r['run'].astype(np.uint64)*np.uint64(2**32)+r['event'].astype(np.uint64);ix=order[np.searchsorted(key[order],ids)]
 assert np.array_equal(key[ix],ids)
 for j,k in enumerate(ix):
  a['found'][k]=bool(r['found'][j]);a['association'][k]=bool(r['track_exists'][j] and r['beta_exists'][j] and r['beta_track'][j]==r['requested_track'][j])
  a['mass'][k]=r['gen_mass'][j]
  primary=(r['track'][j]==r['gen_track'][j])&(r['parent'][j]==r['gen_parent'][j])&(r['particle'][j]==r['gen_particle'][j])
  primary&=np.isfinite(r['momentum'][j])&(r['momentum'][j]>0)&(np.abs(r['mass'][j]-r['gen_mass'][j])<1e-4)&(r['charge'][j]==r['gen_charge'][j])
  pos=np.asarray(r['position'][j]).reshape(-1,3);dire=np.asarray(r['direction'][j]).reshape(-1,3)
  for plane in range(21):
   take=np.flatnonzero(primary&(r['plane'][j]==plane));a['acceptance_count'][k,plane]=len(take)
   if len(take)==1:
    t=take[0];a['acceptance_p'][k,plane]=r['momentum'][j][t];a['acceptance_pos'][k,plane]=pos[t];a['acceptance_dir'][k,plane]=dire[t]
  p=a['acceptance_p'][k,10]
  if p>0:a['u_true'][k]=np.hypot(1.,float(a['mass'][k])/p)
  tp=np.asarray(r['tof_position'][j]).reshape(-1,3)
  primary_tof=(r['tof_track'][j]==r['gen_track'][j])&(r['tof_parent'][j]==r['gen_parent'][j])&(r['tof_particle'][j]==r['gen_particle'][j])
  primary_tof&=np.isfinite(r['tof_time'][j])&(r['tof_time'][j]>0)&np.isfinite(r['tof_beta'][j])&(r['tof_beta'][j]>0)&(r['tof_beta'][j]<=1.00001)
  for layer in range(4):
   allix=np.flatnonzero(r['tof_layer'][j]==layer);take=np.flatnonzero(primary_tof&(r['tof_layer'][j]==layer))
   a['tof_count'][k,layer]=len(allix);a['tof_primary_count'][k,layer]=len(take)
   if len(allix):
    a['tof_first_any_time'][k,layer]=np.min(r['tof_time'][j][allix])*1e9;a['tof_total_edep'][k,layer]=np.sum(r['tof_edep'][j][allix])
   if not len(take):continue
   t=take[np.argmin(r['tof_time'][j][take])]
   for dest,src in [('tof_time','tof_time'),('tof_beta','tof_beta'),('tof_edep','tof_edep'),('tof_step','tof_step')]:a[dest][k,layer]=r[src][j][t]
   a['tof_time'][k,layer]*=1e9;a['tof_pos'][k,layer]=tp[t]
   valid=np.isfinite(a['acceptance_pos'][k,:,2])&np.isfinite(a['acceptance_dir'][k,:,2])
   if valid.sum()>=2:
    zz=a['acceptance_pos'][k,valid,2];vv=a['acceptance_dir'][k,valid,2];si=np.argsort(zz)
    a['tof_dir_z'][k,layer]=np.interp(tp[t,2],zz[si],vv[si])
  tk=(r['tk_primary'][j]==1)&(r['tk_track'][j]==r['gen_track'][j])&(r['tk_particle'][j]==r['gen_particle'][j])&(r['tk_momentum'][j]>0)
  a['native_tk_count'][k]=np.sum(tk)
  if np.sum(tk)==1:a['native_tk_p'][k]=r['tk_momentum'][j][tk][0]
 print('COLLECT',record['source'],len(ix),'seconds',round(time.monotonic()-start,1),flush=True)

zref=b['tof_pos'][:,:,2].astype(float);dz=zref-a['tof_pos'][:,:,2]
a['tof_time_plane']=a['tof_time']+dz/(C*a['tof_dir_z']*a['tof_beta'])
a['truth_valid']=a['found']&a['association']&np.isfinite(a['u_true'])&(a['u_true']>1)&(a['u_true']<4)&(np.abs(a['acceptance_pos'][:,10,2]-1.698)<.08)
a['primary_four']=a['truth_valid']&np.all(a['tof_primary_count']>0,axis=1)&np.isfinite(a['tof_time_plane']).all(1)&(a['tof_dir_z']<-.2).all(1)
a['primary_four_unique']=a['primary_four']&np.all(a['tof_primary_count']==1,axis=1)
a['primary_postS4']=a['primary_four']&(a['acceptance_count'][:,17]==1)&np.isfinite(a['acceptance_p'][:,17])
# All operations below are algebraic, with finite-data masks carried explicitly.
x=b['tof_leng'].astype(float)/C;xc=x-x.mean(1)[:,None];weights=xc/np.sum(xc*xc,axis=1)[:,None]
a['u_time_raw']=np.sum(weights*a['tof_time'],axis=1)
a['u_time_plane']=np.sum(weights*a['tof_time_plane'],axis=1)
a['delta_prop_raw']=a['u_time_raw']-a['u_true'];a['delta_prop']=a['u_time_plane']-a['u_true']
a['noise_vector']=b['tof_tl'].astype(float)-a['tof_time_plane'];a['noise_vector']-=a['noise_vector'].mean(1)[:,None]
a['noise_slope']=np.sum(weights*a['noise_vector'],axis=1)
a['noise_orthogonal']=a['noise_vector']-a['noise_slope'][:,None]*xc
closure=b['u']-a['u_true']-a['delta_prop']-a['noise_slope']
valid=a['primary_four'];maxclosure=float(np.max(np.abs(closure[valid])))
if maxclosure>1e-12:raise RuntimeError('Algebraic decomposition closure failed')
np.savez_compressed(str(HERE/'native_collected.npz'),**a)
summary=dict(total=n,found=int(a['found'].sum()),association=int(a['association'].sum()),native_L5=int(a['truth_valid'].sum()),primary_four=int(a['primary_four'].sum()),primary_four_unique=int(a['primary_four_unique'].sum()),primary_postS4=int(a['primary_postS4'].sum()),missing_primary_by_layer=np.sum(a['tof_primary_count']==0,axis=0),multiple_primary_by_layer=np.sum(a['tof_primary_count']>1,axis=0),max_algebraic_closure=maxclosure,reference_u_difference=metric((a['u_true']-b['u_true'])[a['truth_valid']]),native_mtr_p_difference=metric((a['native_tk_p']-b['mtrmom'][:,4])[a['native_tk_count']==1]),time_plane_shift_ns={str(l):metric((a['tof_time_plane']-a['tof_time'])[valid,l]) for l in range(4)},seconds=time.monotonic()-start,notes=['TOF times: earliest primary track+parent+particle record in the reconstructed bar.','Plane correction uses local native beta and interpolated primary direction; it is a local transport approximation, not a saved exact crossing time.','L5 acceptance and Tracker-cluster reference mismatch is quantified rather than ignored.','Missing/fragmented primary tracks are counted, never interpreted as a resolution improvement.'])
write_json(HERE/'native_collection_summary.json',summary)
print('NATIVE_COLLECTION_COMPLETE',json.dumps({k:v for k,v in summary.items() if isinstance(v,(int,float))}),flush=True)
