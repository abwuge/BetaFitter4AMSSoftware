#!/usr/bin/env python3
"""Test-only feature knockout; no retraining and no event-data transfer."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import json
import numpy as np
from audit_common import HERE,primary_features
import histboost as hb
C=np.load(str(HERE.parent/"c_mass_toy_oracle/feature_cache.npz"))
names=json.loads((HERE.parent/"c_mass_toy_oracle/feature_names.json").read_text())
a0=np.load(str(HERE/"mixed.npz"))
E=np.load(str(HERE/"extended_errors.npz"))
take=E["role"]==1
indices=E["index"][take]
a={k:a0[k][indices] for k in a0.files}
key=C["run"].astype(np.uint64)*np.uint64(2**32)+C["event"].astype(np.uint64)
order=np.argsort(key)
ids=a["run"].astype(np.uint64)*np.uint64(2**32)+a["event"].astype(np.uint64)
j=np.searchsorted(key[order],ids)
assert np.all(j<len(key)) and np.array_equal(key[order[j]],ids)
ci=order[j]
F=primary_features(a)
en_keep=np.array([not n.startswith(("log_trd_qk","trd_lik","log_ecal_q","ecal_bdt")) for n in names["energy"]])
base=np.column_stack([F["measured_tof"],C["mass"][ci],C["timing"][ci],C["nuisance"][ci],C["energy"][ci][:,en_keep]])
X=np.column_stack([base,C["rich"][ci]])
model=np.load(str(HERE/"extended_model_all_with_RICH.npz"))
pred=hb.predict(X,model,trees=640)
P=np.load(str(HERE/"extended_prediction_all_with_RICH.npz"))
assert np.array_equal(P["index"][take],indices)
original=P["prediction"][int(P["selected"])][take]
assert np.max(np.abs(pred-original))<1e-10
def ww(x,w):
 o=np.argsort(x)
 q=np.interp([.16,.84],(np.cumsum(w[o])-.5*w[o])/w.sum(),x[o])
 return (q[1]-q[0])/2
w=E["weights"][take]
err=E["ols"][take]
r=a["rich_beta"]
valid=np.isfinite(r[:,0])&(r[:,0]>.2)&(r[:,0]<1.2)
s0=ww(err-pred,w)
modified=X.copy()
modified[~valid,base.shape[1]:]=np.nan
pred_masked=hb.predict(modified,model,trees=640)
s1=ww(err-pred_masked,w)
modified2=X.copy()
for ix in [2,5,8]:
 modified2[~valid,base.shape[1]+ix]=np.nan
pred_masked2=hb.predict(modified2,model,trees=640)
s2=ww(err-pred_masked2,w)
print("TEST",len(err),"RICH_VALID",int(valid.sum()),"BAD_RICH",int((~valid).sum()),flush=True)
print("REPLAY_MAX",float(np.max(abs(pred-original))),flush=True)
print("WIDTH_ORIGINAL",s0,"WIDTH_MISSING_ALL_RICH",s1,"WIDTH_MISSING_STALE_BETA2_DERIVED",s2,flush=True)
print("PRED_CHANGE_ALL",float(np.mean(np.abs(pred-pred_masked))),"PRED_CHANGE_3",float(np.mean(np.abs(pred-pred_masked2))),flush=True)
print("ERROR_VALID",ww((err-pred)[valid],w[valid]),"ERROR_INVALID",ww((err-pred)[~valid],w[~valid]),flush=True)
print("MASKED_ERROR_VALID",ww((err-pred_masked)[valid],w[valid]),"MASKED_ERROR_INVALID",ww((err-pred_masked)[~valid],w[~valid]),flush=True)
(HERE/"rich_stale_probe.json").write_text(json.dumps({"n":len(err),"RICH_valid":int(valid.sum()),"W68_original":float(s0),"W68_mask_all_when_missing":float(s1),"W68_mask_three_stale_when_missing":float(s2),"replay_max":float(np.max(abs(pred-original)))},indent=2))
