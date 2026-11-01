#!/usr/bin/env python3
"""Broad measured-feature stress test, not an intrinsic-TOF bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import histboost as hb
from audit_common import HERE,primary_features,metric,write_json
from bootstrap_audit import bootstrap
C=np.load(str(HERE.parent/'c_mass_toy_oracle/feature_cache.npz'));names=json.loads((HERE.parent/'c_mass_toy_oracle/feature_names.json').read_text())
raw=np.load(str(HERE/'mixed.npz'));old=np.load(str(HERE/'mixed_analysis_errors.npz'));ii=old['index'];key=C['run'].astype(np.uint64)*np.uint64(2**32)+C['event'].astype(np.uint64);order=np.argsort(key)
assert len(np.unique(key))==len(key)
ids=raw['run'][ii].astype(np.uint64)*np.uint64(2**32)+raw['event'][ii].astype(np.uint64);j=np.searchsorted(key[order],ids);jj=np.minimum(j,len(key)-1);found=(j<len(key))&(key[order[jj]]==ids)
index=ii[found];ci=order[j[found]];a={k:raw[k][index] for k in raw.files};assert np.array_equal(a['label'],C['label'][ci])
tr=a['role']==0;va=a['role']==2;te=a['role']==1;lab=a['label'];R=a['cached_R'];u=a['u'];ut=a['u_true'];e=u-ut
F=primary_features(a);common=F['measured_tof'];q=2*np.log(np.maximum(a['tk_ql'][:,1:8],1e-8)/6);small=np.column_stack([common,q,u,np.arcsinh(R/5)])
energy_keep=np.array([not name.startswith(('log_trd_qk','trd_lik','log_ecal_q','ecal_bdt')) for name in names['energy']])
noRich=np.column_stack([common,C['mass'][ci],C['timing'][ci],C['nuisance'][ci],C['energy'][ci][:,energy_keep]])
allRich=np.column_stack([noRich,C['rich'][ci]])
weights=np.ones(len(e));bins=np.digitize(R,[2.5,3,4,5,6])-1
for role in [0,1,2]:
 mask=a['role']==role
 for k in range(4):
  for label in [0,1]:
   s=mask&(bins==k)&(lab==label)
   if s.sum():weights[s]=mask.sum()/(8*s.sum())
def ww(x,w):
 s=np.argsort(x);x=x[s];w=w[s];return float(np.diff(np.interp([.16,.84],(np.cumsum(w)-.5*w)/w.sum(),x))[0]/2)
def strata(mask):return {'balanced':mask,'C12':mask&(lab==0),'C13':mask&(lab==1)}
plan=[('unrestricted_tracker',small,640,15,120),('all_noRich',noRich,640,15,120),('all_with_RICH',allRich,640,15,120),('all_noRich_capacity',noRich,1280,31,80),('all_noRich_known_isotope',np.column_stack([noRich,lab]),640,15,120)]
write_json(HERE/'extended_plan.json',dict(events=len(e),lost_cache_matches=int((~found).sum()),stages=[dict(name=n,features=X.shape[1],trees=t,leaves=l,minleaf=m) for n,X,t,l,m in plan],excluded_energy_features=[n for n,k in zip(names['energy'],energy_keep) if not k],lineage_caveat='The pre-existing measured-feature cache includes beta-fed Tracker charge, native measured TOF beta, alternative rigidity fits and detector quality. It contains no explicit MC truth target inputs. MC-special calibrated TRD charge/likelihood and ECAL charge/BDT blocks are excluded pending end-to-end producer validation. Remaining cached features are an extended measured-feature stress test, not a certified flight-domain model.'))

def run(item):
 name,X,trees,leaves,minleaf=item;start=time.monotonic();snaps=tuple(sorted(set([80,160,320,640,trees])))
 pred,model=hb.fit(X,np.clip(e,-.05,.05),tr,snapshots=snaps,leaves=leaves,minleaf=minleaf,weight=weights)
 score=[ww((e-p)[va],weights[va]) for p in pred];chosen=int(np.argmin(score));err=e-pred[chosen]
 out=dict(name=name,features=X.shape[1],selected_trees=snaps[chosen],tuning_weighted_w68=score,training_w68=[ww((e-p)[tr],weights[tr]) for p in pred],seconds=time.monotonic()-start,test={})
 for st,mask in strata(te).items():
  m=metric(err[mask]);m['weighted_w68']=ww(err[mask],weights[mask]);m['gain_pct']=100*(1-m['weighted_w68']/ww(e[mask],weights[mask]));out['test'][st]=m
 np.savez_compressed(str(HERE/('extended_model_'+name+'.npz')),**model)
 np.savez_compressed(str(HERE/('extended_prediction_'+name+'.npz')),index=index,prediction=pred,selected=np.array(chosen),snapshots=np.array(snaps))
 print('EXTENDED_FIT',name,out['test']['balanced']['gain_pct'],out['seconds'],flush=True)
 return name,err,out
with ThreadPoolExecutor(max_workers=4) as pool:fits=list(pool.map(run,plan))
errors={'ols':e,'measured_tof':old['measured_tof'][found],'simple_unrestricted':old['unrestricted'][found]};reports={}
for name,err,out in fits:errors[name]=err;reports[name]=out
stats={}
for st,mask in strata(te).items():stats[st]=bootstrap({k:v[mask] for k,v in errors.items()},a['source'][mask],'ols',increment='measured_tof',weights=weights[mask],stratification=lab[mask])
# Mass-squared performance of the selected beta estimator, keeping measured rigidity fixed.
mass={}
for name,err in errors.items():
 rel2=(6*R)**2*((ut+err)**2-1)/a['true_mass']**2-1
 mass[name]={st:dict(w68=ww(rel2[mask],weights[mask]),mean=float(np.average(rel2[mask],weights=weights[mask])),rmse=float(np.sqrt(np.average(rel2[mask]**2,weights=weights[mask])))) for st,mask in strata(te).items()}
write_json(HERE/'extended_summary.json',dict(events=len(e),counts={st:int(mask.sum()) for st,mask in strata(te).items()},models=reports,bootstrap=stats,mass=mass,claim_scope='Achieved MC prediction error in selected mixed-production samples. No end-to-end beta-field reprocessing/unit-response test or flight-domain validation is claimed for the high-dimensional cache models.'))
np.savez_compressed(str(HERE/'extended_errors.npz'),index=index,source=a['source'],role=a['role'],label=lab,R=R,u_true=ut,weights=weights,**errors)
print('EXTENDED_COMPLETE',len(e),int(te.sum()),flush=True)
