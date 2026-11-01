#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import time,json
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import histboost as hb
from audit_common import HERE,C,primary_features,metric,w68,write_json
raw=np.load(str(HERE/'mixed.npz'));good=raw['primary_valid']&np.isfinite(raw['u'])&(raw['cached_R']>=2.5)&(raw['cached_R']<6)&(1/raw['u_true']>.6)
indices=np.flatnonzero(good);a={k:raw[k][good] for k in raw.files};tr=a['role']==0;va=a['role']==2;te=a['role']==1
for m1,m2 in [(tr,va),(tr,te),(va,te)]:assert not np.intersect1d(a['run'][m1],a['run'][m2]).size
F=primary_features(a);e=a['u']-a['u_true'];R=a['cached_R'];truthR=a['mtrmom'][:,4]/6;labels=a['label'];u=a['u'];ut=a['u_true'];bins=np.digitize(R,[2.5,3,4,5,6])-1
weights=np.ones(len(e))
for role in [0,1,2]:
 mask=a['role']==role
 for k in range(4):
  for lab in [0,1]:
   take=mask&(bins==k)&(labels==lab);count=take.sum()
   if count:weights[take]=mask.sum()/(8*count)

def weighted_width(x,w):
 s=np.argsort(x);x=x[s];w=w[s];cdf=(np.cumsum(w)-.5*w)/np.sum(w)
 return float(np.diff(np.interp([.16,.84],cdf,x))[0]/2)

def strata(mask):
 out={'balanced_R2p5_6':mask,'C12_R2p5_6':mask&(labels==0),'C13_R2p5_6':mask&(labels==1)}
 for j,(lo,hi) in enumerate([(2.5,3),(3,4),(4,5),(5,6)]):
  for lab in [0,1]:out['C%d_R%g_%g'%(12+lab,lo,hi)]=mask&(labels==lab)&(bins==j)
 return out

def metrics(err):
 out={}
 for name,take in strata(te).items():
  if take.sum()<40:continue
  q=metric(err[take]);base=metric(e[take]);q['gain_vs_OLS_pct']=100*(1-q['w68']/base['w68']);q['baseline_OLS']=base
  q['weighted_w68']=weighted_width(err[take],weights[take]);q['weighted_gain_vs_OLS_pct']=100*(1-q['weighted_w68']/weighted_width(e[take],weights[take]));out[name]=q
 return out
common=F['measured_tof'];invar=F['geometry_time'];q=2*np.log(np.maximum(a['tk_ql'][:,1:8],1e-8)/6)
Xgroups={'geometry':F['geometry'],'time_residual':invar,'measured_tof':common,'beta_fed_tracker':np.column_stack([common,q]),'unrestricted':np.column_stack([common,u,np.arcsinh(R/5)]),'known_isotope':np.column_stack([common,u,np.arcsinh(R/5),labels]),'true_rigidity':np.column_stack([common,u,np.arcsinh(truthR/5)]),'shuffle':common}
(HERE/'mixed_model_plan.json').write_text(json.dumps(dict(stages=list(Xgroups),fit_weight='Equal isotope and measured-R-bin weights within each original whole-file role. Truth labels set only weights, except explicitly labelled isotope oracle.',truth_label='Primary mtrmom L5 plus native C12 mass or nominal C13 mass; no mixture likelihood supervision.',tuning='Weighted W68 on original validation-source files only.',test_counts={s:int(m.sum()) for s,m in strata(te).items()}),indent=2))

def fitone(item):
 name,X=item;target=np.clip(e,-.05,.05).copy()
 if name=='shuffle':
  ix=np.flatnonzero(tr);target[ix]=target[np.random.RandomState(1232016).permutation(ix)]
 snaps=(80,160,320,640);start=time.monotonic();pred,model=hb.fit(X,target,tr,snapshots=snaps,minleaf=120,weight=weights)
 scores=[weighted_width((e-p)[va],weights[va]) for p in pred];choice=int(np.argmin(scores));chosen=pred[choice];nt=snaps[choice]
 np.savez_compressed(str(HERE/('mixed_model_'+name+'.npz')),**model)
 np.savez_compressed(str(HERE/('mixed_prediction_'+name+'.npz')),index=indices,prediction=pred,selected=np.array(choice),snapshots=np.array(snaps))
 injection=[]
 sample=np.flatnonzero(te)[:1500];Xs=X[sample].copy();p0=hb.predict(Xs,model,trees=nt)
 if name in ['unrestricted','known_isotope','true_rigidity']:
  for delta in [-.01,.01]:
   Xa=Xs.copy();Xa[:,common.shape[1]]+=delta;pp=hb.predict(Xa,model,trees=nt);response=(delta-(pp-p0))/delta
   injection.append(dict(delta=delta,mean_response=float(np.mean(response)),median_response=float(np.median(response)),correction_change_rms=float(np.sqrt(np.mean((pp-p0)**2)))))
 out=dict(name=name,features=X.shape[1],train=int(tr.sum()),tune=int(va.sum()),test=int(te.sum()),selected_trees=nt,tuning_weighted_width=scores,training_width=[weighted_width((e-p)[tr],weights[tr]) for p in pred],test_metrics=metrics(e-chosen),injection=injection,seconds=time.monotonic()-start)
 print('MIXED_FIT',name,nt,out['test_metrics']['balanced_R2p5_6']['weighted_gain_vs_OLS_pct'],flush=True)
 return name,e-chosen,out
with ThreadPoolExecutor(max_workers=6) as pool:fits=list(pool.map(fitone,Xgroups.items()))
errors={'ols':e,'frozen_beta_only':a['cached_u_beta_only_frozen']-ut,'frozen_joint':a['cached_u_joint_anchored']-ut}
reports={name:dict(test_metrics=metrics(err)) for name,err in errors.items()}
for name,err,report in fits:errors[name]=err;reports[name]=report
# Exact mass-squared identity avoids discarding events whose reconstructed mass squared is negative.
ms_true=a['true_mass']**2
mass_methods={
 'OLS_beta_measured_R':(u,R),'true_beta_measured_R':(ut,R),'OLS_beta_true_R':(u,truthR),'all_truth':(ut,truthR),
 'measured_tof_beta_measured_R':(ut+errors['measured_tof'],R),'unrestricted_beta_measured_R':(ut+errors['unrestricted'],R),
 'known_isotope_beta_measured_R':(ut+errors['known_isotope'],R),'frozen_beta_only_measured_R':(a['cached_u_beta_only_frozen'],R),
}
mass_errors={};mass_metrics={}
for name,(ub,rr) in mass_methods.items():
 m2=(6*rr)**2*(ub*ub-1);rel2=m2/ms_true-1;mass_errors[name]=rel2;detail={}
 for s,take in strata(te).items():
  if take.sum()<40:continue
  mt=metric(rel2[take]);mt['negative_mass2_fraction']=float(np.mean(m2[take]<0));mt['weighted_w68']=weighted_width(rel2[take],weights[take])
  positive=take&(m2>0);mt['positive_mass_relative']=metric(np.sqrt(m2[positive]/ms_true[positive])-1)
  detail[s]=mt
 mass_metrics[name]=detail
cov={}
for name,take in strata(te).items():
 if take.sum()<40:continue
 er=R/truthR-1;du=(u-ut)*ut/(ut*ut-1)
 cov[name]=dict(n=int(take.sum()),linear_relative_mass_terms_covariance=np.cov(np.column_stack([er,du])[take].T,ddof=0),relative_rigidity=metric(er[take]),beta_mass_term=metric(du[take]))
prior_tests={}
for frac in [.01,.1,.5]:
 ww=np.where(labels==1,frac,1-frac).astype(float)
 for lab in [0,1]:ww[labels==lab]/=np.sum(te&(labels==lab))
 prior_tests[str(frac)]={name:dict(w68=weighted_width(err[te],ww[te]),mean=float(np.average(err[te],weights=ww[te])),rmse=float(np.sqrt(np.average(err[te]**2,weights=ww[te])))) for name,err in errors.items()}
np.savez_compressed(str(HERE/'mixed_analysis_errors.npz'),index=indices,source=a['source'],run=a['run'],role=a['role'],label=labels,R=R,u_true=ut,weights=weights,**errors)
np.savez_compressed(str(HERE/'mixed_mass_errors.npz'),index=indices,source=a['source'],role=a['role'],label=labels,R=R,weights=weights,**mass_errors)
write_json(HERE/'mixed_analysis_summary.json',dict(events=len(e),counts={s:int(m.sum()) for s,m in strata(te).items()},models=reports,mass=mass_metrics,mass_linear_covariance=cov,prior_weight_tests=prior_tests,reference_caveat='True rigidity from the same L5 momentum gives an exact identity with L5 beta. Measured rigidity is at nominal z=0, so true-beta/measured-R counterfactual includes reference-plane transport mismatch and Tracker error. C13 native mass is not saved; nominal 12.11221 GeV is used.',domain_caveat='C12 and C13 are different MC productions. Common TOF features exclude absent C13 tracker Edep/raw-time fields; beta-fed tracker and isotope features are separate diagnostic/oracle groups. Test results are not a flight-domain guarantee.'))
print('MIXED_ANALYSIS_COMPLETE',len(e),int(te.sum()),flush=True)
