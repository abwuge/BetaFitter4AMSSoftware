#!/usr/bin/env python3
from __future__ import print_function
import json
from pathlib import Path
import numpy as np
import uproot
HERE=Path(__file__).resolve().parent
MASS=HERE.parent/"c_mass_toy_oracle"/"predictions.npz"
STAGES=["strictgeom_c12train","strictgeom_c13train"]
def w68(x):
 x=np.asarray(x,float); q=np.percentile(x[np.isfinite(x)],[16,84]); return float((q[1]-q[0])/2)
def load(stage):
 chunks=[]
 for b in range(4):
  with uproot.open(str(HERE/"tmva_strict_transfer"/("%s_bin%d.root"%(stage,b)))) as f:
   key=next(k for k in f.keys(recursive=True) if "TestTree" in k)
   chunks.append(f[key].arrays(["event_index","BDTG"],library="np"))
 return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}
def main():
 with np.load(str(MASS),allow_pickle=False) as z:
  role=np.char.decode(z["role"],"ascii"); test=role=="test"; idx=np.flatnonzero(test)
  y=z["label"][test].astype(np.int8); src=z["source"][test].astype(np.int32)
  raw=np.log(z["R"][test]/z["R_true"][test])
 idxmap={int(e):i for i,e in enumerate(idx)}
 out={}
 errors={}
 for stage in STAGES:
  a=load(stage); pred=np.full(len(idx),np.nan)
  for ei,v in zip(a["event_index"].astype(np.int64),a["BDTG"]): pred[idxmap[int(ei)]]=float(v)
  err=raw-pred; errors[stage]=err
  d={}
  for lab,name in [(0,"C12"),(1,"C13")]:
   m=y==lab
   d[name]={"n":int(np.sum(m)),"w68":w68(err[m]),"gain_pct":100*(1-w68(err[m])/w68(raw[m])),"median":float(np.median(err[m]))}
  d["median_gap_C12_minus_C13"]=d["C12"]["median"]-d["C13"]["median"]
  out[stage]=d
 groups={}; keys={0:[],1:[]}
 for lab in (0,1):
  for s in sorted(int(x) for x in np.unique(src[y==lab])):
   k=(lab,s); groups[k]=np.flatnonzero((y==lab)&(src==s)); keys[lab].append(k)
 rng=np.random.RandomState(300933); reps=1000; boot={}
 for stage in STAGES:
  boot[stage]={}
  for lab,name in [(0,"C12"),(1,"C13")]:
   vals=[]
   for _ in range(reps):
    ii=np.concatenate([groups[keys[lab][j]] for j in rng.randint(0,len(keys[lab]),size=len(keys[lab]))])
    vals.append(100*(1-w68(errors[stage][ii])/w68(raw[ii])))
   a=np.asarray(vals)
   boot[stage][name]={"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))}
 result={"results":out,"bootstrap":boot,"note":"Strict geometry model trained on only one isotope's train+validation source files; both isotopes evaluated only on untouched test source files."}
 (HERE/"strict_transfer_results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("STRICT TRANSFER")
 for stage in STAGES:
  d=out[stage]
  print(stage,"C12",d["C12"],"C13",d["C13"],"gap",d["median_gap_C12_minus_C13"])
 print("BOOT",boot)
if __name__=="__main__":main()
