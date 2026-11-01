#!/usr/bin/env python3
from __future__ import print_function
import json
from pathlib import Path
import numpy as np
import uproot

HERE=Path(__file__).resolve().parent
MASS=HERE.parent/"c_mass_toy_oracle"/"predictions.npz"
STAGES=["morph_c12train","morph_c13train","charge_c12train","charge_c13train"]

def w68(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    q=np.percentile(x,[16,84]); return float((q[1]-q[0])/2)

def load_stage(stage):
    chunks=[]
    for b in range(4):
        p=HERE/"tmva_transfer"/("%s_bin%d.root"%(stage,b))
        with uproot.open(str(p)) as f:
            key=next(k for k in f.keys(recursive=True) if "TestTree" in k)
            chunks.append(f[key].arrays(["event_index","target","BDTG","label","source"],library="np"))
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}

def metric(error, raw, label, mask):
    return {
        "n":int(np.sum(mask)),
        "raw_w68":w68(raw[mask]),
        "corrected_w68":w68(error[mask]),
        "gain_pct":100*(1-w68(error[mask])/w68(raw[mask])),
        "raw_median":float(np.median(raw[mask])),
        "corrected_median":float(np.median(error[mask])),
        "rms":float(np.sqrt(np.mean(error[mask]**2))),
    }

def main():
    with np.load(str(MASS),allow_pickle=False) as z:
        roles=np.char.decode(z["role"],"ascii"); test=roles=="test"
        base_idx=np.flatnonzero(test)
        base={
            "label":z["label"][test].astype(np.int8),
            "source":z["source"][test].astype(np.int32),
            "bin":z["bin"][test].astype(np.int8),
            "raw":np.log(z["R"][test]/z["R_true"][test]),
        }
    idxmap={int(e):i for i,e in enumerate(base_idx)}
    results={}
    errors={}
    for stage in STAGES:
        a=load_stage(stage)
        pred=np.full(len(base_idx),np.nan)
        for ei,v in zip(a["event_index"].astype(np.int64),a["BDTG"]):
            pred[idxmap[int(ei)]]=float(v)
        if np.any(~np.isfinite(pred)): raise RuntimeError("missing "+stage)
        err=base["raw"]-pred; errors[stage]=err
        result={"overall":metric(err,base["raw"],base["label"],np.ones(len(err),bool)),"by_isotope":{},"bins":{}}
        for lab,name in [(0,"C12"),(1,"C13")]:
            result["by_isotope"][name]=metric(err,base["raw"],base["label"],base["label"]==lab)
        for b in range(4):
            result["bins"][str(b)]=metric(err,base["raw"],base["label"],base["bin"]==b)
        result["median_C12_minus_C13"]=result["by_isotope"]["C12"]["corrected_median"]-result["by_isotope"]["C13"]["corrected_median"]
        results[stage]=result

    # paired source-file bootstrap on transfer gain, separately for same/opposite isotope.
    groups={}; keys={0:[],1:[]}
    for lab in (0,1):
        for src in sorted(int(x) for x in np.unique(base["source"][base["label"]==lab])):
            k=(lab,src); groups[k]=np.flatnonzero((base["label"]==lab)&(base["source"]==src)); keys[lab].append(k)
    rng=np.random.RandomState(300931); reps=1000
    boot={}
    for stage in STAGES:
        vals={"C12":[],"C13":[]}
        for _ in range(reps):
            for lab,name in [(0,"C12"),(1,"C13")]:
                idx=np.concatenate([groups[keys[lab][j]] for j in rng.randint(0,len(keys[lab]),size=len(keys[lab]))])
                vals[name].append(100*(1-w68(errors[stage][idx])/w68(base["raw"][idx])))
        boot[stage]={}
        for name,a in vals.items():
            a=np.asarray(a)
            boot[stage][name]={"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))}
    out={"results":results,"bootstrap":boot,"note":"Each model is trained on train+validation files of only one isotope and evaluated on the untouched test files of both isotopes."}
    (HERE/"transfer_results.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("TRANSFER SUMMARY")
    for stage in STAGES:
        r=results[stage]
        print("%-16s overall gain %6.2f med %+0.5f | C12 gain %6.2f med %+0.5f | C13 gain %6.2f med %+0.5f | med-gap %+0.5f"%(
            stage,r["overall"]["gain_pct"],r["overall"]["corrected_median"],
            r["by_isotope"]["C12"]["gain_pct"],r["by_isotope"]["C12"]["corrected_median"],
            r["by_isotope"]["C13"]["gain_pct"],r["by_isotope"]["C13"]["corrected_median"],
            r["median_C12_minus_C13"]))
    print("BOOTSTRAP")
    for stage in STAGES: print(stage,boot[stage])

if __name__=="__main__": main()
