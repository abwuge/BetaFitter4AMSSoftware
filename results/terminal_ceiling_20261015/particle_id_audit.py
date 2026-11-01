#!/usr/bin/env python3
"""Evaluate actual mass-squared ranking, never truth-normalized mass as a classifier."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from audit_common import HERE,write_json

raw=np.load(str(HERE/'mixed.npz'))
E=np.load(str(HERE/'extended_errors.npz'))
ix=E['index'];role=E['role'];lab=E['label'];R=E['R'];ut=E['u_true'];ww=E['weights'];te=role==1
assert np.array_equal(raw['label'][ix],lab)
truthR=raw['mtrmom'][ix,4]/6
mass=raw['true_mass'][ix]

def weighted_quantiles(x,w,qq):
    order=np.argsort(x);return np.interp(qq,(np.cumsum(w[order])-.5*w[order])/w.sum(),x[order])
def auc(score,y,w):
    order=np.argsort(score,kind='mergesort');s=score[order];p=w[order]*(y[order]==1);n=w[order]*(y[order]==0)
    starts=np.r_[0,np.flatnonzero(np.diff(s)!=0)+1]
    pg=np.add.reduceat(p,starts);ng=np.add.reduceat(n,starts)
    return float(np.sum(pg*(np.cumsum(ng)-.5*ng))/(p.sum()*n.sum()))
methods={name:((ut+E[name])**2-1)*(6*R)**2 for name in ['ols','measured_tof','simple_unrestricted','all_noRich','all_noRich_capacity','all_with_RICH','all_noRich_known_isotope']}
methods['true_beta_measured_R']=(ut*ut-1)*(6*R)**2
methods['ols_beta_true_R']=((ut+E['ols'])**2-1)*(6*truthR)**2
methods['all_truth']=(ut*ut-1)*(6*truthR)**2
scopes={'balanced_R2p5_6':te}
for lo,hi in [(2.5,3),(3,4),(4,5),(5,6)]:scopes['R%g_%g'%(lo,hi)]=te&(R>=lo)&(R<hi)
results={}
for scope,sel in scopes.items():
    r={}
    for name,score in methods.items():
        s=score[sel];y=lab[sel];w=ww[sel]
        c={}
        for l in [0,1]:
            v=s[y==l];wt=w[y==l];q=weighted_quantiles(v,wt,[.16,.5,.84])
            c['C%d'%(l+12)]={'n':len(v),'median_m2_GeV2':float(q[1]),'w68_m2_GeV2':float((q[2]-q[0])/2)}
        c['auc']=auc(s,y,w)
        c['separation_median_w68']=(c['C13']['median_m2_GeV2']-c['C12']['median_m2_GeV2'])/np.hypot(c['C12']['w68_m2_GeV2'],c['C13']['w68_m2_GeV2'])
        r[name]=c
    results[scope]=r
# Sensitivity of the C13 reference mass: fixed predictions, only reference changes.
reference={}
for dm in [-.0003,.0003]:
    mnew=mass*np.where(lab==1,1+dm,1);unew=np.hypot(1,mnew/(6*truthR))
    sel=te&(lab==1)
    reference[str(dm)]={name:{'w68_u':float(np.diff(weighted_quantiles((ut+E[name]-unew)[sel],ww[sel],[.16,.84]))[0]/2),'mean_u':float(np.average((ut+E[name]-unew)[sel],weights=ww[sel]))} for name in ['ols','measured_tof','all_noRich']}
write_json(HERE/'particle_id_audit.json',{'description':'Unretrained mass-squared score AUC; positive=C13; identical held-out events. No true-mass-normalized feature is used in AUC. C13 nominal generator mass caveat retained.','results':results,'C13_mass_reference_sensitivity':reference})
print('PID_COMPLETE',int(te.sum()),flush=True)
for name,r in results['balanced_R2p5_6'].items():print(name,'AUC',round(r['auc'],6),'D',round(r['separation_median_w68'],6),flush=True)
print('C13_MASS_REFERENCE',reference,flush=True)
