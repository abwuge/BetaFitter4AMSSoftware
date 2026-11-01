#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,sys
TAG='native_counter' if '--counter' in sys.argv else 'native'
from pathlib import Path
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import histboost as hb
from audit_common import HERE,C,primary_features,partitions,strata,metric,w68,write_json
B=np.load(str(HERE/'broad_Z6.npz'));N=np.load(str(HERE/(TAG+'_collected.npz')));G=np.load(str(HERE/(TAG+'_geometry.npz')))
keep=N['primary_four']&(1/N['u_true']>.6)&(1/N['u_true']<=1)
index=np.flatnonzero(keep);a={k:B[k][keep] for k in B.files};n={k:N[k][keep] for k in N.files};g={k:G[k][keep] for k in G.files};a['u_true']=n['u_true']
tr,va,te=partitions(a);common=te&a['previous_comparison']&np.isfinite(a['previous_hybridBeta'])&(a['previous_hybridBeta']>0)
features=primary_features(a);base=features['measured_all'];geom=features['geometry'];energy=base[:,23:34];span=np.ptp(a['tof_leng']/C,axis=1)
residual=n['noise_orthogonal']/span[:,None]
x=a['tof_leng'].astype(float)/C;xc=x-x.mean(1)[:,None];weights=xc/np.sum(xc*xc,axis=1)[:,None]
e=a['u']-n['u_true'];noise=n['noise_slope'];prop=n['delta_prop']
methods={'ols':e,'hybrid':1/a['previous_hybridBeta']-n['u_true'],'previous_new':1/a['previous_newBeta']-n['u_true'],'oracle_transport':noise,'oracle_transport_uncorrected_time':a['u']-n['u_true']-n['delta_prop_raw'],'oracle_noise_removed':prop,'all_truth':np.zeros(len(e))}
for name,fn in [('measured_frozen','Z6_measured_all_t640_l15_m160_f100_core'),('unrestricted_frozen','Z6_unrestricted_t640_l15_m160_f100_core')]:
 p=np.load(str(HERE/('pred_'+fn+'.npz')));full=np.full(len(B['run']),np.nan);full[p['index']]=p['correction'][int(p['selected'])];methods[name]=e-full[index]
for sign in [-1,1]:
 dt=.5*n['tof_step']/(C*n['tof_beta'])*sign
 methods['oracle_transport_halfstep_%s'%('minus' if sign<0 else 'plus')]=noise-np.sum(weights*dt,axis=1)
# Identity and exact total covariance accounting, before optional geometry split.
results={};covariance={}
for name,mask in strata(a,common).items():
 if mask.sum()<40:continue
 E=np.column_stack([prop,noise]);cov=np.cov(E[mask].T,ddof=0)
 covariance[name]=dict(n=int(mask.sum()),prop=metric(prop[mask]),noise=metric(noise[mask]),covariance=cov,variance_total=float(np.var(e[mask])),variance_recombined=float(cov.sum()),prop_variance_share=float(cov[0,0]/np.var(e[mask])),noise_variance_share=float(cov[1,1]/np.var(e[mask])),cross_variance_share=float(2*cov[0,1]/np.var(e[mask])))

stages={
 'measured_refit':(base,e,e),
 'transport_plus_measured':(base,noise,noise),
 'transport_plus_true_residual':(np.column_stack([geom,energy,residual]),noise,noise),
 'transport_plus_true_residual_deposit':(np.column_stack([geom,energy,residual,np.log(np.maximum(n['tof_edep'],1e-8)),np.log(np.maximum(n['tof_total_edep'],1e-8))]),noise,noise),
 'transport_prediction_only':(base,prop,e),
}
(HERE/(TAG+'_model_plan.json')).write_text(json.dumps(dict(stages=list(stages),snapshots=[80,160,320,640],selection='Tune-source W68 only; oracle stages contain true transport or true orthogonal residual, never the true noise slope as an input. All output errors assessed against native primary L5 acceptance beta.'),indent=2))

def fitstage(item):
 name,(X,y,baseline)=item;started=time.monotonic();p,model=hb.fit(X,np.clip(y,-.05,.05),tr,snapshots=(80,160,320,640),minleaf=160)
 score=[metric((baseline-z)[va]) for z in p];chosen=int(np.argmin([s['w68'] for s in score]));nt=[80,160,320,640][chosen]
 # For propagation-only regression selection must assess prediction error, not accidentally exploit timing correlations.
 if name=='transport_prediction_only':
  score=[metric((y-z)[va]) for z in p];chosen=int(np.argmin([s['rmse'] for s in score]));nt=[80,160,320,640][chosen]
 np.savez_compressed(str(HERE/(TAG+'_model_'+name+'.npz')),**model)
 np.savez_compressed(str(HERE/(TAG+'_prediction_'+name+'.npz')),index=index,prediction=p,selected=np.array(chosen),snapshots=np.array([80,160,320,640]))
 out=dict(name=name,feature_count=X.shape[1],selected_trees=nt,tuning=score,train=int(tr.sum()),tune=int(va.sum()),test=int(common.sum()),seconds=time.monotonic()-started)
 print('NATIVE_FIT',name,nt,'gain',100*(1-w68((baseline-p[chosen])[common])/w68(methods['hybrid'][common])),flush=True)
 return name,baseline-p[chosen],out
with ThreadPoolExecutor(max_workers=5) as pool:fits=list(pool.map(fitstage,stages.items()))
fitnotes=[]
for name,err,note in fits:methods[name]=err;fitnotes.append(note)
# Also examine geometry-only and loss-only interventions on the identical sparse-trajectory-valid cohort.
methods['oracle_geometry_only']=e-g['delta_geometry'];methods['oracle_loss_only']=e-g['delta_loss']
for stratum,mask in strata(a,common).items():
 if mask.sum()<40:continue
 results[stratum]={}
 for name,err in methods.items():
  sel=mask&np.isfinite(err)
  if name in ['oracle_geometry_only','oracle_loss_only']:sel&=g['valid']
  if sel.sum()<40:continue
  mm=metric(err[sel]);mm['gain_vs_hybrid_pct']=100*(1-mm['w68']/w68(methods['hybrid'][sel]));mm['baseline_hybrid_same_events']=metric(methods['hybrid'][sel]);results[stratum][name]=mm
geom_results={}
for name,mask in strata(a,common&g['valid']).items():
 if mask.sum()<40:continue
 mat=np.column_stack([g['delta_geometry'],g['delta_loss'],noise]);cov=np.cov(mat[mask].T,ddof=0)
 geom_results[name]=dict(n=int(mask.sum()),terms={k:metric(mat[mask,j]) for j,k in enumerate(['geometry','energy_and_time_convention','effective_timing'])},covariance=cov,geometry_Hermite_minus_chord=metric((g['delta_geometry']-g['delta_geometry_chord'])[mask]),profile_p_closure=metric(g['profile_p_minus_native'][mask]),profile_E_closure=metric(g['profile_E_minus_native'][mask]))
np.savez_compressed(str(HERE/(TAG+'_analysis_errors.npz')),index=index,source=a['source'],run=a['run'],role=a['role'],test=common,train=tr,tune=va,u_true=n['u_true'],R=a['tk_rigidityk'][:,1,1],geometry_valid=g['valid'],**methods)
cohorts={}
for name,flag in [('all_L5',N['truth_valid']),('primary_four',N['primary_four']),('primary_postS4',N['primary_postS4'])]:
 sel=flag&(1/N['u_true']>.6)&(B['role']==1)&B['previous_comparison']&(B['previous_hybridBeta']>0)
 cohorts[name]=dict(n=int(sel.sum()),ols=metric((B['u']-N['u_true'])[sel]),hybrid=metric((1/B['previous_hybridBeta']-N['u_true'])[sel]))
write_json(HERE/(TAG+'_analysis_summary.json'),dict(cohorts=cohorts,results=results,covariance=covariance,geometry=geom_results,model_fits=fitnotes,reference='Native primary acceptance record at L5; nominal z=1.698cm; local MC TOF time-to-plane correction explicitly included.',limit_statement='These genie interventions and fitted conditional regressions provide measured noise budgets and empirical saturation, not a nonparametric certified upper bound on all estimators or on W68.'))
print('NATIVE_ANALYSIS_COMPLETE',len(index),int(common.sum()),flush=True)
