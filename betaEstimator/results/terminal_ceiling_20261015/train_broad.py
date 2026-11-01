#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import argparse,time,json,hashlib
from pathlib import Path
import numpy as np
import histboost as hb
from audit_common import HERE,primary_features,partitions,strata,metric,w68,write_json,C

def run(z,group,trees,leaves,minleaf,fraction,loss):
 raw=np.load(str(HERE/('broad_Z%d.npz'%z)));keep=(1/raw['u_true']>.6)&(1/raw['u_true']<=1)
 a={k:raw[k][keep] for k in raw.files};indices=np.flatnonzero(keep)
 tr,va,te=partitions(a);train_available=int(tr.sum())
 if fraction<1:
  rng=np.random.RandomState(102615);take=rng.rand(len(tr))<fraction;tr&=take
 groups=primary_features(a)
 source_group='measured_all' if group=='shuffle' else group
 X=groups[source_group];e=a['u']-a['u_true']
 target=e.copy() if loss=='mse' else np.clip(e,-{1:.15,2:.08,6:.05,8:.05}[z],{1:.15,2:.08,6:.05,8:.05}[z])
 if group=='shuffle':
  ix=np.flatnonzero(tr);target[ix]=target[np.random.RandomState(31516).permutation(ix)]
 snaps=tuple(sorted(set([min(80,trees),min(160,trees),min(320,trees),trees])))
 tag='Z%d_%s_t%d_l%d_m%d_f%d_%s'%(z,group,trees,leaves,minleaf,round(fraction*100),loss)
 if (HERE/('fit_'+tag+'.json')).exists():raise RuntimeError('Existing fit '+tag)
 start=time.monotonic();print('START',tag,int(tr.sum()),int(va.sum()),int(te.sum()),X.shape,flush=True)
 pred,model=hb.fit(X,target,tr,snapshots=snaps,leaves=leaves,minleaf=minleaf)
 score=[metric((e-p)[va]) for p in pred]
 selected=int(np.argmin([m['w68'] if loss=='core' else m['rmse'] for m in score]))
 chosen=pred[selected];nsel=snaps[selected]
 model['selected_trees']=np.array(nsel);np.savez_compressed(str(HERE/('model_'+tag+'.npz')),**model)
 np.savez_compressed(str(HERE/('pred_'+tag+'.npz')),index=indices,correction=pred,snapshots=np.array(snaps),selected=np.array(selected),train=tr,tune=va,test=te)
 tests={}
 # Final test requires the same available old-baseline comparison for every method.
 common=te&a['previous_comparison']&np.isfinite(a['previous_hybridBeta'])&(a['previous_hybridBeta']>0)
 for s,mask in strata(a,common).items():
  if mask.sum()<40:continue
  base=metric(e[mask]);m=metric((e-chosen)[mask]);m['gain_w68_pct']=100*(1-m['w68']/base['w68']);m['baseline']=base
  for name in ['linearBeta','hybridBeta','nonlinearBeta','newBeta']:
   m['previous_'+name]=metric(1/a['previous_'+name][mask]-a['u_true'][mask])
  m['gain_vs_hybrid_pct']=100*(1-m['w68']/m['previous_hybridBeta']['w68'])
  m['relative_beta_w68']=w68(a['u_true'][mask]/(a['u'][mask]-chosen[mask])-1)
  tests[s]=m
 injection=[]
 ix=np.flatnonzero(common)[:1000];small={k:v[ix].copy() for k,v in a.items()}
 original=hb.predict(primary_features(small)[source_group],model,trees=nsel)
 for delta in [-.03,.03]:
  altered={k:v.copy() for k,v in small.items()};shift=10+delta*small['tof_leng'].astype(float)/C
  altered['tof_tl']=small['tof_tl'].astype(float)+shift
  if 'tof_rtime' in small:altered['tof_rtime']=small['tof_rtime'].astype(float)-shift[:,:,None]
  p=hb.predict(primary_features(altered)[source_group],model,trees=nsel)
  injection.append(dict(delta_u=delta,correction_change_max=float(np.max(np.abs(p-original))),correction_change_rms=float(np.sqrt(np.mean((p-original)**2))),interpretation='Simultaneous time shift of corrected and raw timing; other fields fixed. Beta-fed observables are NOT reprocessed, so this is not an end-to-end invariance test for those groups.'))
 result=dict(tag=tag,Z=z,group=group,feature_count=X.shape[1],feature_lineage='audit_common.primary_features',train_available=train_available,train=int(tr.sum()),tune=int(va.sum()),test_events=int(te.sum()),paired_test=int(common.sum()),snapshots=snaps,selected_trees=nsel,selection='Minimum tuning-source W68' if loss=='core' else 'Minimum tuning-source RMSE',tuning_metrics=score,training_metrics=[metric((e-p)[tr]) for p in pred],test=tests,injection=injection,loss=loss,truth_used_only_as_training_target=group not in ['truth_deposit'],seconds=time.monotonic()-start,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 write_json(HERE/('fit_'+tag+'.json'),result)
 print('COMPLETE',tag,nsel,'seconds',result['seconds'],'all_gain',tests['all_beta_gt06']['gain_vs_hybrid_pct'],'lowR_gain',tests.get('R2p5_6',{}).get('gain_vs_hybrid_pct'),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--z',type=int,required=True);p.add_argument('--group',default='measured_all');p.add_argument('--trees',type=int,default=160);p.add_argument('--leaves',type=int,default=15);p.add_argument('--minleaf',type=int,default=160);p.add_argument('--fraction',type=float,default=1.);p.add_argument('--loss',choices=['core','mse'],default='core');x=p.parse_args();run(x.z,x.group,x.trees,x.leaves,x.minleaf,x.fraction,x.loss)
