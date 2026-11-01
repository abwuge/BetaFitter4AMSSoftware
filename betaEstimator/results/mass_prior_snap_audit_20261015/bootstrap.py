#!/usr/bin/env python3
"""Fixed-model paired source bootstrap for the anchored energy correction."""
import json
from pathlib import Path
import numpy as np

ROOT=Path('/nas01/home/hxwu/AMSSoft/betaEstimator/results')
OUT=ROOT/'mass_prior_snap_audit_20261015'
p=np.load(ROOT/'c_joint_beta_rigidity_closure/expanded_validation/predictions.npz')
e=np.load(ROOT/'c_joint_beta_rigidity_closure/expanded_validation/energy_comparison/energy_predictions.npz')
m=np.load(ROOT/'terminal_ceiling_20261015/mixed.npz')
sel=p['truth_valid']&np.isfinite(e['u_anchored_energy'])
l=p['label'][sel];r=p['R'][sel];ut=p['u_true'][sel]
src=m['source'][sel]
uniq,inv=np.unique(src,return_inverse=True)
file_labels=np.array([l[np.flatnonzero(inv==j)[0]] for j in range(len(uniq))])
assert all(np.all(l[inv==j]==file_labels[j]) for j in range(len(uniq)))
score_raw=(6*r)**2*(p['u_raw'][sel]**2-1)
score_energy=(6*r)**2*(e['u_anchored_energy'][sel]**2-1)
err_raw=ut/p['u_raw'][sel]-1
err_energy=ut/e['u_anchored_energy'][sel]-1

def auc_prep(score):
    o=np.argsort(score,kind='stable')
    s=score[o]
    start=np.r_[0,np.flatnonzero(np.diff(s)!=0)+1]
    return o,start,l[o]==1
def auc(prep,w):
    o,start,pos=prep
    wp=w[o]*pos
    wn=w[o]*(~pos)
    pg=np.add.reduceat(wp,start)
    ng=np.add.reduceat(wn,start)
    return float(np.dot(pg,np.cumsum(ng)-.5*ng)/(pg.sum()*ng.sum()))

def w_prep(d,mask=None):
    idx=np.arange(len(l)) if mask is None else np.flatnonzero(mask)
    return idx[np.argsort(d[idx])]
def w68(prep,vec,w):
    vals=vec[prep]
    wc=w[prep]
    good=wc>0
    vals=vals[good];wc=wc[good]
    cw=(np.cumsum(wc)-wc*.5)/sum(wc)
    return float(np.diff(np.interp([.16,.84],cw,vals))[0]/2)

AR=auc_prep(score_raw);AE=auc_prep(score_energy)
o12r=w_prep(err_raw,l==0);o12e=w_prep(err_energy,l==0)
o13r=w_prep(err_raw,l==1);o13e=w_prep(err_energy,l==1)
omr=w_prep(err_raw);ome=w_prep(err_energy)
def evalw(w):
    a=auc(AR,w);b=auc(AE,w)
    c12r=w68(o12r,err_raw,w);c12e=w68(o12e,err_energy,w)
    c13r=w68(o13r,err_raw,w);c13e=w68(o13e,err_energy,w)
    adjusted=np.where(l==0,w*.9/(w[l==0].sum()),w*.1/(w[l==1].sum()))
    mixr=w68(omr,err_raw,adjusted);mixe=w68(ome,err_energy,adjusted)
    return [a,b,b-a,1-c12e/c12r,1-c13e/c13r,1-mixe/mixr]
point=evalw(np.ones(len(l)))
rng=np.random.RandomState(20261015)
replicates=600
val=[]
for _ in range(replicates):
    nw=np.zeros(len(uniq),dtype=np.int16)
    for k in [0,1]:
        fs=np.flatnonzero(file_labels==k)
        nw[fs]=rng.multinomial(len(fs),np.ones(len(fs))/len(fs))
    val.append(evalw(nw[inv].astype(float)))
val=np.asarray(val)
names=['AUC_raw','AUC_energy','delta_AUC','relative_W68_gain_C12',
       'relative_W68_gain_C13','relative_W68_gain_mix_pi13_0.1']
out={'N':len(l),'sources_total':len(uniq),'source_counts':{'C12':int((file_labels==0).sum()),'C13':int((file_labels==1).sum())},
    'replicates':replicates,'scheme':'Paired source-file resampling, stratified isotope. Models already frozen; excludes retraining and domain shift.',
    'metrics':{name:{'point':float(point[i]),'CI95':list(map(float,np.percentile(val[:,i],[2.5,97.5])))} for i,name in enumerate(names)}}
(OUT/'bootstrap.json').write_text(json.dumps(out,indent=2))
print('BOOTSTRAP',len(l),len(uniq),replicates,flush=True)
for n,v in out['metrics'].items():print(n,v,flush=True)
