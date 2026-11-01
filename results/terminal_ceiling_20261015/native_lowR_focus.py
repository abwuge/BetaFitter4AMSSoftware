#!/usr/bin/env python3
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import histboost as hb
from audit_common import HERE,C,primary_features,partitions,metric,write_json,w68
from bootstrap_audit import bootstrap
B=np.load(str(HERE/'broad_Z6.npz'));N=np.load(str(HERE/'native_counter_collected.npz'));keep=N['primary_four']&(1/N['u_true']>.6);index=np.flatnonzero(keep)
a={k:B[k][keep] for k in B.files};n={k:N[k][keep] for k in N.files};F=primary_features(a);tr,va,te=partitions(a);R=a['tk_rigidityk'][:,1,1];low=(R>=2.5)&(R<6);val=va&low;test=te&low&a['previous_comparison']&(a['previous_hybridBeta']>0)
e=a['u']-n['u_true'];noise=n['noise_slope'];span=np.ptp(a['tof_leng']/C,axis=1)
X0=F['measured_all'];X1=np.column_stack([F['geometry'],X0[:,23:34],n['noise_orthogonal']/span[:,None]])
X2=np.column_stack([X1,np.log(np.maximum(n['tof_edep'],1e-8)),np.log(np.maximum(n['tof_total_edep'],1e-8))])
groups={'measured':(X0,e),'transport_residual':(X1,noise),'transport_residual_deposit':(X2,noise)}
strategies=['lowR_only','lowR_weighted'];plan=[]
for strategy in strategies:
 for group,(X,target) in groups.items():plan.append((strategy,group,X,target))
write_json(HERE/'native_lowR_plan.json',dict(groups=list(groups),strategies=strategies,validation='Only R2.5-6 source-held-out tuning events choose checkpoint and strategy; no test error selection.',counts=dict(train_low=int((tr&low).sum()),train_all=int(tr.sum()),tune_low=int(val.sum()),test_low=int(test.sum()))))

def run(item):
 strategy,name,X,target=item;train=tr&low if strategy=='lowR_only' else tr.copy();weight=np.ones(len(e))
 if strategy=='lowR_weighted':
  weight[low]=.9*tr.sum()/max(1,(tr&low).sum());weight[~low]=.1*tr.sum()/max(1,(tr&~low).sum())
 p,m=hb.fit(X,np.clip(target,-.05,.05),train,snapshots=(80,160,320,640),leaves=15,minleaf=60,weight=weight)
 scores=[w68((target-v)[val]) for v in p];best=int(np.argmin(scores));nameall=strategy+'_'+name
 np.savez_compressed(str(HERE/('native_lowR_model_'+nameall+'.npz')),**m)
 np.savez_compressed(str(HERE/('native_lowR_prediction_'+nameall+'.npz')),index=index,prediction=p,selected=np.array(best))
 report=dict(strategy=strategy,group=name,selected_trees=[80,160,320,640][best],tuning_w68=scores,test=metric((target-p[best])[test]))
 print('LOWR_FOCUS',nameall,report,flush=True)
 return nameall,target-p[best],report
with ThreadPoolExecutor(max_workers=6) as pool:fits=list(pool.map(run,plan))
errors={'hybrid':1/a['previous_hybridBeta']-n['u_true']};reports={};selected={}
for name,err,report in fits:errors[name]=err;reports[name]=report
for group in groups:
 candidates=[s+'_'+group for s in strategies];best=min(candidates,key=lambda k:min(reports[k]['tuning_w68']));selected[group]=best;errors['selected_'+group]=errors[best]
boot=bootstrap({k:v[test] for k,v in errors.items()},a['source'][test],'hybrid',increment='selected_measured')
write_json(HERE/'native_lowR_focus_summary.json',dict(selected=selected,models=reports,bootstrap=boot))
np.savez_compressed(str(HERE/'native_lowR_focus_errors.npz'),index=index,test=test,source=a['source'],**errors)
print('LOWR_FOCUS_COMPLETE',selected,flush=True)
