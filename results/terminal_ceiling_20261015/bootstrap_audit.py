#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from audit_common import HERE,w68,write_json
REPS=600

def width_sorted(values,weights,integer):
 take=weights>0;v=values[take];w=weights[take]
 if not len(v):return np.nan
 c=np.cumsum(w)
 if integer:
  rank=np.array([.16,.84])*(c[-1]-1);lo=np.floor(rank).astype(int);hi=np.ceil(rank).astype(int)
  vl=v[np.searchsorted(c,lo,side='right')];vh=v[np.searchsorted(c,hi,side='right')];q=vl+(vh-vl)*(rank-lo)
 else:q=np.interp([.16,.84],(c-.5*w)/c[-1],v)
 return float((q[1]-q[0])/2)

def bootstrap(errors,source,baseline,increment=None,weights=None,stratification=None,seed=102615):
 names=list(errors);n=len(source);uniq,inv=np.unique(source,return_inverse=True);K=len(uniq)
 integer=weights is None;basew=np.ones(n) if integer else np.asarray(weights,float)
 order={k:np.argsort(errors[k]) for k in names};sorted_values={k:errors[k][order[k]] for k in names}
 point={k:width_sorted(sorted_values[k],basew[order[k]],integer) for k in names}
 rng=np.random.RandomState(seed);draw=np.zeros((REPS,K),np.int16)
 if stratification is None:draw=rng.multinomial(K,np.ones(K)/K,size=REPS)
 else:
  ss=np.array([stratification[np.flatnonzero(inv==j)[0]] for j in range(K)])
  for group in np.unique(ss):
   ix=np.flatnonzero(ss==group);draw[:,ix]=rng.multinomial(len(ix),np.ones(len(ix))/len(ix),size=REPS)
 widths=np.empty((REPS,len(names)))
 for j in range(REPS):
  ww=basew*draw[j,inv]
  for k,name in enumerate(names):widths[j,k]=width_sorted(sorted_values[name],ww[order[name]],integer)
 ib=names.index(baseline);ii=names.index(increment) if increment else None;out={}
 for k,name in enumerate(names):
  gain=100*(1-widths[:,k]/widths[:,ib]);m=dict(w68=point[name],gain_pct=100*(1-point[name]/point[baseline]),gain95=np.percentile(gain,[2.5,97.5]),w68_95=np.percentile(widths[:,k],[2.5,97.5]))
  if increment:
   inc=100*(1-widths[:,k]/widths[:,ii]);m['additional_gain_pct']=100*(1-point[name]/point[increment]);m['additional_gain95']=np.percentile(inc,[2.5,97.5])
  out[name]=m
 return dict(n=n,sources=K,replicates=REPS,baseline=baseline,increment_reference=increment,weighted=not integer,methods=out,interval_scope='Paired whole-source bootstrap conditional on the fitted models; not model-class uncertainty or a certified information bound.')

def broad(z):
 a=np.load(str(HERE/('broad_Z%d.npz'%z)));take=(a['role']==1)&a['previous_comparison']&(a['previous_hybridBeta']>0)&(1/a['u_true']>.6)
 methods={'hybrid':1/a['previous_hybridBeta']-a['u_true']}
 paths=sorted(HERE.glob('fit_Z%d_*.json'%z))
 for p in paths:
  m=json.loads(p.read_text());q=np.load(str(HERE/('pred_'+m['tag']+'.npz')));corr=np.full(len(take),np.nan);corr[q['index']]=q['correction'][int(q['selected'])];methods[m['tag']]=a['u']-a['u_true']-corr
 out={}
 for st,mask in [('all',take),('R2p5_6',take&(a['tk_rigidityk'][:,1,1]>=2.5)&(a['tk_rigidityk'][:,1,1]<6))]:
  out[st]=bootstrap({k:v[mask] for k,v in methods.items()},a['source'][mask],'hybrid',seed=102615+z)
 write_json(HERE/('bootstrap_broad_Z%d.json'%z),out);print('BOOTSTRAP_BROAD',z,flush=True)

def native():
 a=np.load(str(HERE/'native_counter_analysis_errors.npz'));take=a['test'];names=['hybrid','measured_frozen','measured_refit','oracle_transport','transport_plus_measured','transport_plus_true_residual','transport_plus_true_residual_deposit','oracle_noise_removed','oracle_transport_halfstep_minus','oracle_transport_halfstep_plus'];out={}
 for st,mask in [('all',take),('R2p5_6',take&(a['R']>=2.5)&(a['R']<6)),('beta06_09',take&(1/a['u_true']<.9))]:
  out[st]=bootstrap({k:a[k][mask] for k in names},a['source'][mask],'hybrid',increment='measured_refit')
 write_json(HERE/'bootstrap_native.json',out);print('BOOTSTRAP_NATIVE',flush=True)

def mixed():
 a=np.load(str(HERE/'mixed_analysis_errors.npz'));take=a['role']==1;names=['ols','frozen_beta_only','frozen_joint','measured_tof','beta_fed_tracker','unrestricted','known_isotope','true_rigidity','shuffle'];out={}
 for lab in [-1,0,1]:
  mask=take if lab==-1 else take&(a['label']==lab)
  out['balanced' if lab==-1 else 'C%d'%(12+lab)]=bootstrap({k:a[k][mask] for k in names},a['source'][mask],'ols',increment='measured_tof',weights=a['weights'][mask],stratification=a['label'][mask])
 write_json(HERE/'bootstrap_mixed.json',out)
 a=np.load(str(HERE/'mixed_mass_errors.npz'));take=a['role']==1;names=['OLS_beta_measured_R','true_beta_measured_R','OLS_beta_true_R','measured_tof_beta_measured_R','unrestricted_beta_measured_R','known_isotope_beta_measured_R'];out={}
 for lab in [-1,0,1]:
  mask=take if lab==-1 else take&(a['label']==lab)
  out['balanced' if lab==-1 else 'C%d'%(12+lab)]=bootstrap({k:a[k][mask] for k in names},a['source'][mask],'OLS_beta_measured_R',weights=a['weights'][mask],stratification=a['label'][mask])
 write_json(HERE/'bootstrap_mass.json',out);print('BOOTSTRAP_MIXED',flush=True)
if __name__=='__main__':
 import sys
 tasks=[lambda:broad(1),lambda:broad(2),lambda:broad(6),lambda:broad(8),native]
 if '--mixed' in sys.argv:tasks=[mixed]
 with ThreadPoolExecutor(max_workers=5) as pool:list(pool.map(lambda fun:fun(),tasks))
