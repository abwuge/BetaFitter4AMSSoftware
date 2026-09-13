#!/usr/bin/env python3
"""Train the standalone conditional estimator from extracted arrays.

The sdiat runtime uses a temporary Python 3.6-compatible binary environment;
no packages are installed into the user's system Python.
"""
import argparse, json, pickle
from pathlib import Path
import numpy as np
try:
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa: F401
except ImportError:
    pass
from sklearn.ensemble import HistGradientBoostingRegressor

C=29.9792458
CHARGES=(1,2,6,8)
CAP={1:.15,2:.08,6:.05,8:.05}
PARAMS=dict(max_iter=160,max_leaf_nodes=15,max_depth=6,min_samples_leaf=160,l2_regularization=20,learning_rate=.06,early_stopping=False,random_state=613)

def load(path):
    with np.load(path,allow_pickle=False) as f:return {k:f[k] for k in f.files}

def kin(d):
    x=d['tof_leng'].astype(float)/C;t=d['tof_tl'].astype(float);xc=x-x.mean(axis=1)[:,None]
    with np.errstate(invalid='ignore',divide='ignore'):
        u=(xc*t).sum(axis=1)/(xc*xc).sum(axis=1)
        r=t-t.mean(axis=1)[:,None]-xc*u[:,None]
    return x,u,r

def clean(d):
    x,u,_=kin(d);t=d['tof_tl'].astype(float)
    return (d['tof_bpatt']==4444)&np.all(np.isfinite(x)&np.isfinite(t),axis=1)&(np.ptp(x,axis=1)>0)&np.isfinite(u)&(u>0)

def features(d,z):
    x,_,r=kin(d);xc=x-x.mean(axis=1)[:,None]
    with np.errstate(invalid='ignore',divide='ignore'):
        rr=r/np.ptp(x,axis=1)[:,None]
        z2=float(z*z)
        dtof=np.log(np.where(d['tof_edep']>0,d['tof_edep'],np.nan)/z2)
        dtk=np.log(np.where(d['tk_edep'][:,1:8]>0,d['tk_edep'][:,1:8],np.nan)/z2)
    X=np.column_stack([d['position'].reshape(-1,8),d['direction'],d['tof_barid'],xc,rr,dtof,dtk]).astype(float)
    X[~np.isfinite(X)]=np.nan
    return X

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,default=Path('results/comparison/cache'));ap.add_argument('--output',type=Path,default=Path('results/comparison'));args=ap.parse_args()
    model_dir=args.output/'models';model_dir.mkdir(parents=True,exist_ok=True)
    summary={'parameters':PARAMS,'target_caps':CAP,'charges':{}}
    for z in CHARGES:
        tr=load(args.cache/f'train_Z{z}.npz');va=load(args.cache/f'validation_Z{z}.npz')
        _,ut,_=kin(tr)
        train_mask=(tr['charge']==z)&(tr['mch'].astype(int)==z)&(tr['mcBetaRaw']>=.25)&(tr['mcBetaRaw']<=1)&np.isfinite(tr['tof_betah'])&(tr['tof_betah']>0)&clean(tr)
        target=np.clip(ut[train_mask]-1/tr['mcBetaRaw'][train_mask],-CAP[z],CAP[z])
        model=HistGradientBoostingRegressor(**PARAMS);model.fit(features(tr,z)[train_mask],target)
        _,uv,_=kin(va);candidate=uv-model.predict(features(va,z))
        valid=clean(va)&np.isfinite(candidate)&(candidate>0)
        candidate=np.where(valid,candidate,1/va['tof_betah'].astype(float))
        np.save(args.output/f'new_inverse_beta_Z{z}.npy',candidate)
        with (model_dir/f'model_Z{z}.pkl').open('wb') as f:pickle.dump(model,f,pickle.HIGHEST_PROTOCOL)
        summary['charges'][str(z)]={'training_entries':int(train_mask.sum()),'validation_entries':int(len(candidate)),'fallback_entries':int((~valid).sum())}
        print('Z=%d train=%d validation=%d fallback=%d'%(z,train_mask.sum(),len(candidate),(~valid).sum()),flush=True)
    (args.output/'training_summary.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
