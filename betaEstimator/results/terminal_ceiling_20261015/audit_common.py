import hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
C=29.9792458

def w68(x):
 x=np.asarray(x);x=x[np.isfinite(x)]
 return float(np.diff(np.percentile(x,[16,84]))[0]/2) if len(x) else float('nan')

def metric(e):
 e=np.asarray(e);e=e[np.isfinite(e)]
 if not len(e):return dict(n=0)
 return dict(n=len(e),mean=float(np.mean(e)),median=float(np.median(e)),w68=w68(e),std=float(np.std(e)),rmse=float(np.sqrt(np.mean(e*e))),tail02=float(np.mean(np.abs(e)>.02)),tail05=float(np.mean(np.abs(e)>.05)))

def primary_features(a):
 t=a['tof_tl'].astype(float);x=a['tof_leng'].astype(float)/C;xc=x-x.mean(1)[:,None];span=np.ptp(x,axis=1)
 u=np.sum(xc*t,axis=1)/np.sum(xc*xc,axis=1)
 r=(t-t.mean(1)[:,None]-u[:,None]*xc)/span[:,None]
 p=a['tof_pos'].astype(float);dz=p[:,3,2]-p[:,0,2]
 g=np.column_stack([p[:,:,0]/120,p[:,:,1]/120,(p[:,3,:2]-p[:,0,:2])/dz[:,None],xc,np.log(span/4),a['tof_barid']])
 z=a['Z'].astype(float)
 lt=np.log(np.maximum(a['tof_edep'].astype(float),1e-8)/(z[:,None]**2))
 lk=np.log(np.maximum(a['tk_edep'][:,1:8].astype(float),1e-8)/(z[:,None]**2))
 q=2*np.log(np.maximum(a['tk_ql'][:,1:8].astype(float),1e-8)/z[:,None])
 groups=dict(geometry=g,geometry_time=np.column_stack([g,r]),measured_tof=np.column_stack([g,r,lt]),measured_all=np.column_stack([g,r,lt,lk]))
 groups['beta_fed']=np.column_stack([groups['measured_all'],q])
 truth_dep=np.log(np.maximum(a['mtofdep'].astype(float),1e-8)/(z[:,None]**2))
 groups['truth_deposit']=np.column_stack([groups['measured_all'],truth_dep])
 adc_key='tof_Aadc' if 'tof_Aadc' in a else 'tof_adca' if 'tof_adca' in a else None
 raw=[]
 if adc_key:raw.append(np.log(np.maximum(a[adc_key].astype(float),1e-4)).reshape(len(t),-1))
 if 'tof_rtime' in a:
  rt=a['tof_rtime'].astype(float);cor=t+rt.mean(2);cor-=cor.mean(1)[:,None]
  raw += [cor/span[:,None],(rt[:,:,0]-rt[:,:,1])/10.]
 groups['raw_readout']=np.column_stack([groups['measured_all']]+raw)
 rr=a['tk_rigidityk'][:,1,1].astype(float)
 groups['unrestricted']=np.column_stack([groups['raw_readout'],u,np.arcsinh(rr/5)])
 if 'rich_beta' in a:
  rich=a['rich_beta'][:,0].astype(float);present=(rich>.5)&(rich<1.1)
  groups['rich_fusion']=np.column_stack([groups['unrestricted'],np.where(present,1/np.maximum(rich,.01),0),present])
 return groups

def partitions(a):
 train_sources=sorted(np.unique(a['source'][a['role']==0]).tolist());vsrc=np.array(train_sources[::5]);tune=(a['role']==0)&np.isin(a['source'],vsrc)
 train=(a['role']==0)&~tune;test=a['role']==1
 assert not np.intersect1d(a['run'][train],a['run'][tune]).size
 assert not np.intersect1d(a['run'][train|tune],a['run'][test]).size
 return train,tune,test

def strata(a,test):
 b=1/a['u_true'];r=a['tk_rigidityk'][:,1,1]
 ans={'all_beta_gt06':test,'R2p5_6':test&(r>=2.5)&(r<6),'beta06_09':test&(b<.9),'beta09_1':test&(b>=.9)}
 for lo,hi in zip([.6,.65,.7,.75,.8,.85,.9,.95],[.65,.7,.75,.8,.85,.9,.95,1.00000001]):ans['beta_%.2f_%.2f'%(lo,min(hi,1))]=test&(b>=lo)&(b<hi)
 return ans

def write_json(path,x):
 def default(v):
  if isinstance(v,np.ndarray):return v.tolist()
  if isinstance(v,(np.integer,np.floating,np.bool_)):return v.item()
  raise TypeError(str(type(v)))
 Path(path).write_text(json.dumps(x,indent=2,default=default,allow_nan=False))
