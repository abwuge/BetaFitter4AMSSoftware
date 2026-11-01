#!/usr/bin/env python3
"""Sparse native-trajectory quadrature and propagation component accounting."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import sys
COUNTER='--counter' in sys.argv
TAG='native_counter' if COUNTER else 'native'
from audit_common import HERE,C,metric,write_json
n=np.load(str(HERE/(TAG+'_collected.npz')));b=np.load(str(HERE/'broad_Z6.npz'))
N=len(n['run']);pos=n['acceptance_pos'];dire=n['acceptance_dir'];p=n['acceptance_p'];z=pos[:,:,2]
idx=np.argsort(np.where(np.isfinite(z),z,np.inf),axis=1);rr=np.arange(N)[:,None]
z=z[rr,idx];xy=pos[rr,idx,:2];di=dire[rr,idx];p=p[rr,idx]
with np.errstate(all='ignore'):slope=di[:,:,:2]/di[:,:,2:3]
length=np.zeros((N,4));chord_length=np.zeros((N,4));time_linear_p=np.zeros((N,4));time_linear_E=np.zeros((N,4));covered=np.zeros((N,4))
q,w=np.polynomial.legendre.leggauss(5);m=n['mass'];zref=n['z_counter'] if COUNTER else b['tof_pos'][:,:,2].astype(float)
for seg in range(20):
 za,zb=z[:,seg],z[:,seg+1];width=zb-za
 valid=np.isfinite(za)&np.isfinite(zb)&(width>1e-5)&np.isfinite(p[:,seg])&np.isfinite(p[:,seg+1])&(p[:,seg]>0)&(p[:,seg+1]>0)
 width=np.where(valid,width,1.);a0=np.where(valid[:,None],xy[:,seg],0);a1=np.where(valid[:,None],xy[:,seg+1],0);s0=np.where(valid[:,None],slope[:,seg],0);s1=np.where(valid[:,None],slope[:,seg+1],0)
 chord=np.sqrt(1+np.sum(((a1-a0)/width[:,None])**2,axis=1))
 ea=np.hypot(np.where(valid,p[:,seg],1),m);eb=np.hypot(np.where(valid,p[:,seg+1],1),m)
 for layer in range(1,4):
  lower=np.maximum(np.where(valid,za,0),zref[:,layer]);upper=np.minimum(np.where(valid,zb,0),zref[:,0]);span=np.where(valid,np.maximum(upper-lower,0),0)
  covered[:,layer]+=span;chord_length[:,layer]+=span*chord
  for node,weight in zip(q,w):
   zz=.5*(upper+lower)+.5*span*node;f=np.clip((zz-np.where(valid,za,0))/width,0,1)
   derivative=((6*f*f-6*f)[:,None]*a0+(-6*f*f+6*f)[:,None]*a1)/width[:,None]+(3*f*f-4*f+1)[:,None]*s0+(3*f*f-2*f)[:,None]*s1
   speed=np.sqrt(1+np.sum(derivative*derivative,axis=1))
   pp=(1-f)*np.where(valid,p[:,seg],1)+f*np.where(valid,p[:,seg+1],1);ee=(1-f)*ea+f*eb
   uu=np.hypot(1,m/pp);ue=ee/np.sqrt(np.maximum(ee*ee-m*m,1e-20))
   dlength=.5*span*weight*speed
   length[:,layer]+=dlength;time_linear_p[:,layer]+=dlength*uu/C;time_linear_E[:,layer]+=dlength*ue/C
 x=b['tof_leng'].astype(float)/C;xc=x-x.mean(1)[:,None];weights=xc/np.sum(xc*xc,axis=1)[:,None]
valid=n['primary_postS4']&np.isfinite(length).all(1)&np.all(np.abs(covered[:,1:]-(zref[:,0,None]-zref[:,1:]))<1e-4,axis=1)
valid&=(n['acceptance_count'][:,2]==1)&(n['acceptance_count'][:,17]==1)
ug=np.sum(weights*length/C,axis=1);ugch=np.sum(weights*chord_length/C,axis=1)
g=n['u_true']*(ug-1);loss=n['delta_prop']-g
res=dict(run=n['run'],event=n['event'],valid=valid,length=length,chord_length=chord_length,time_profile_p=time_linear_p,time_profile_E=time_linear_E,delta_geometry=g,delta_loss=loss,delta_geometry_chord=n['u_true']*(ugch-1),profile_p_minus_native=np.sum(weights*time_linear_p,axis=1)-n['u_time_plane'],profile_E_minus_native=np.sum(weights*time_linear_E,axis=1)-n['u_time_plane'])
np.savez_compressed(str(HERE/(TAG+'_geometry.npz')),**res)
summary={'valid':int(valid.sum()),'counts':{str(int(role)):int(np.sum(valid&(n['role']==role))) for role in [0,1]},'caveat':'Geometry is reconstructed from sparse primary truth checkpoints using cubic Hermite interpolation, compared to chord interpolation. True G4 step-by-step path length is not saved. Component split is approximate; total propagation/noise algebra is exact for the adopted TOF-time convention.'}
for key in ['delta_geometry','delta_loss','profile_p_minus_native','profile_E_minus_native']:
 summary[key]=metric(res[key][valid])
summary['geometry_Hermite_minus_chord']=metric((res['delta_geometry']-res['delta_geometry_chord'])[valid])
write_json(HERE/(TAG+'_geometry_summary.json'),summary);print('GEOMETRY_COMPLETE',summary,flush=True)
