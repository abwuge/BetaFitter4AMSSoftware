#!/usr/bin/env python3
"""Paired source bootstrap for the mass-score AUC of fixed estimators."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from audit_common import HERE,write_json
raw=np.load(str(HERE/'mixed.npz'));e=np.load(str(HERE/'extended_errors.npz'));take=e['role']==1
ix=e['index'][take];y=e['label'][take];r=e['R'][take];ut=e['u_true'][take];w=e['weights'][take];src=e['source'][take]
p=raw['mtrmom'][ix,4];tr=p/6
scores={name:((ut+e[name][take])**2-1)*(6*r)**2 for name in ['ols','measured_tof','simple_unrestricted','all_noRich','all_noRich_capacity','all_with_RICH']}
scores['true_beta_measured_R']=(ut*ut-1)*(6*r)**2
scores['ols_beta_true_R']=((ut+e['ols'][take])**2-1)*(6*tr)**2
names=list(scores);orders={n:np.argsort(s,kind='mergesort') for n,s in scores.items()};starts={n:np.r_[0,np.flatnonzero(np.diff(scores[n][orders[n]])!=0)+1] for n in names}
u,inv=np.unique(src,return_inverse=True);sl=np.array([y[np.flatnonzero(inv==k)[0]] for k in range(len(u))]);rng=np.random.RandomState(20261015);draw=np.zeros((600,len(u)),np.int16)
for l in [0,1]:
    j=np.flatnonzero(sl==l);draw[:,j]=rng.multinomial(len(j),np.ones(len(j))/len(j),size=600)
def calc(name,weight):
    o=orders[name];pos=weight[o]*(y[o]==1);neg=weight[o]*(y[o]==0);pg=np.add.reduceat(pos,starts[name]);ng=np.add.reduceat(neg,starts[name]);return float(np.sum(pg*(np.cumsum(ng)-.5*ng))/(pg.sum()*ng.sum()))
values=np.array([[calc(n,w*sample[inv]) for n in names] for sample in draw]);point=np.array([calc(n,w) for n in names]);bi=names.index('ols');si=names.index('measured_tof')
out={n:{'auc':point[k],'auc95':np.percentile(values[:,k],[2.5,97.5]),'delta_vs_OLS':point[k]-point[bi],'delta_vs_OLS95':np.percentile(values[:,k]-values[:,bi],[2.5,97.5]),'delta_vs_measured':point[k]-point[si],'delta_vs_measured95':np.percentile(values[:,k]-values[:,si],[2.5,97.5])} for k,n in enumerate(names)}
write_json(HERE/'particle_id_bootstrap.json',{'events':len(y),'sources':len(u),'replicates':600,'source_stratification':'isotope','scope':'Fixed-model, paired whole-source bootstrap. Source composition and production-domain systematic effects are not covered.','methods':out})
print('AUC_BOOTSTRAP_COMPLETE',len(y),len(u),flush=True)
for n,v in out.items():print(n,v,flush=True)
