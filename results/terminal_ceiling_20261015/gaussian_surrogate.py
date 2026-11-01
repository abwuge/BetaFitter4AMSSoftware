#!/usr/bin/env python3
"""Explicitly conditional Gaussian GLS/Fisher reference; never a universal bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import special
from types import SimpleNamespace
def standardized_moment(x,k):
 z=(x-np.mean(x))/np.std(x);return np.mean(z**k)
def ks_normal(x,unused=None):
 v=np.sort(x);cdf=special.ndtr(v);n=len(v);D=max(np.max(np.arange(1,n+1)/n-cdf),np.max(cdf-np.arange(n)/n))
 # Asymptotic Kolmogorov tail; diagnostic only (noise model estimated on other sources).
 k=np.arange(1,101);z=(np.sqrt(n)+.12+.11/np.sqrt(n))*D;p=np.clip(2*np.sum((-1.)**(k-1)*np.exp(-2*k*k*z*z)),0,1)
 return SimpleNamespace(statistic=D,pvalue=p)
def bracket_root(fun,lo,hi):
 for _ in range(70):
  mid=.5*(lo+hi)
  if fun(mid)>0:hi=mid
  else:lo=mid
 return .5*(lo+hi)
stats=SimpleNamespace(skew=lambda x:standardized_moment(x,3),kurtosis=lambda x:standardized_moment(x,4)-3,kstest=ks_normal)
optimize=SimpleNamespace(brentq=bracket_root)
import histboost as hb
from audit_common import HERE,C,primary_features,partitions,metric,w68,write_json
B=np.load(str(HERE/'broad_Z6.npz'));N=np.load(str(HERE/'native_counter_collected.npz'))
keep=N['primary_four']&(1/N['u_true']>.6);index=np.flatnonzero(keep);a={k:B[k][keep] for k in B.files};n={k:N[k][keep] for k in N.files}
tr,va,te=partitions(a);te&=a['previous_comparison']&(a['previous_hybridBeta']>0)
F=primary_features(a);X=np.column_stack([F['geometry'],F['measured_all'][:,23:34]])
x=a['tof_leng'].astype(float)/C;s=x[:,1:]-x[:,0,None];eps=n['noise_vector'][:,1:]-n['noise_vector'][:,0,None]

def meanfit(j):
 p,m=hb.fit(X,np.clip(eps[:,j],-.5,.5),tr,snapshots=(80,160,320),minleaf=160)
 score=[np.mean((eps[va,j]-v[va])**2) for v in p];best=int(np.argmin(score))
 np.savez_compressed(str(HERE/('gaussian_mean_model_%d.npz'%j)),**m)
 return p[best],dict(component=j,selected_trees=[80,160,320][best],tuning_mse=score)
with ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(meanfit,range(3)))
mu=np.column_stack([p for p,m in r]);noise=eps-mu;energy=np.mean(np.log(np.maximum(a['tof_edep'],1e-8)/36),axis=1);result={};stored={}
for label,scope,ncells in [('all_beta_gt06',np.ones(len(index),bool),8),('R2p5_6',(a['tk_rigidityk'][:,1,1]>=2.5)&(a['tk_rigidityk'][:,1,1]<6),4)]:
 calib=va&scope;test=te&scope;edges=np.percentile(energy[tr&scope],np.linspace(0,100,ncells+1)[1:-1]);cell=np.digitize(energy,edges)
 global_cov=np.cov(noise[calib].T);covs=[];means=[];counts=[]
 for k in range(ncells):
  sel=calib&(cell==k);counts.append(int(sel.sum()));m=np.mean(noise[sel],axis=0) if sel.sum() else np.zeros(3);cov=np.cov(noise[sel].T) if sel.sum()>=30 else global_cov
  cov=.9*cov+.1*global_cov+np.eye(3)*1e-10;covs.append(cov);means.append(m)
 covs=np.array(covs);means=np.array(means);inv=np.linalg.inv(covs);vec=np.einsum('nij,nj->ni',inv[cell],s);info=np.sum(s*vec,axis=1);variance=1/info
 estimate=np.sum(vec*(noise-means[cell]),axis=1)/info
 pull=estimate[test]/np.sqrt(variance[test]);pred_rms=float(np.sqrt(np.mean(variance[test])))
 fun=lambda value:np.mean(special.ndtr(value/np.sqrt(variance[test])))-.84
 pred_w68=float(optimize.brentq(fun,0,10*pred_rms));baseline=1/a['previous_hybridBeta'][test]-n['u_true'][test]
 ks=stats.kstest(pull,'norm');m=metric(estimate[test]);m['gain_vs_hybrid_w68_pct']=100*(1-m['w68']/w68(baseline));m['gain_vs_hybrid_rms_pct']=100*(1-m['rmse']/np.sqrt(np.mean(baseline**2)))
 result[label]=dict(n=int(test.sum()),calibration_events=int(calib.sum()),cells=ncells,cell_counts=counts,covariances=covs,gaussian_predicted_rms=pred_rms,gaussian_predicted_w68=pred_w68,conditional_gaussian_gain_rms_pct=100*(1-pred_rms/np.sqrt(np.mean(baseline**2))),conditional_gaussian_gain_w68_pct=100*(1-pred_w68/w68(baseline)),actual_GLS=m,pull=dict(mean=float(np.mean(pull)),std=float(np.std(pull)),skew=float(stats.skew(pull)),excess_kurtosis=float(stats.kurtosis(pull)),KS_stat=float(ks.statistic),KS_p=float(ks.pvalue),coverage_1=float(np.mean(np.abs(pull)<1)),coverage_2=float(np.mean(np.abs(pull)<2))),assumptions='Known deterministic transport; side information treated fixed; independent events; relative-time noise conditionally multivariate Gaussian with known covariance (here estimated from validation sources). No speed/rigidity prior, unbiased unit-response location estimator.')
 stored['error_'+label]=estimate;stored['variance_'+label]=variance
write_json(HERE/'gaussian_surrogate_summary.json',dict(results=result,mean_models=[m for p,m in r],interpretation='The inverse Fisher-information expression is exact only under the stated model. Estimated covariances and a failed Gaussian pull check do NOT establish a strict bound for the actual MC or for flight data. In particular an MSE lower bound is not an unconditional W68 lower bound.'))
np.savez_compressed(str(HERE/'gaussian_surrogate_predictions.npz'),index=index,test=te,**stored)
print('GAUSSIAN_COMPLETE',result,flush=True)
