#!/usr/bin/env python3
from __future__ import print_function
import json
from pathlib import Path
import numpy as np
import uproot

HERE=Path(__file__).resolve().parent
MASS=HERE.parent/"c_mass_toy_oracle"/"predictions.npz"
STAGES=["strict_multifit","strict_fitdiag","strict_geometry","strict_morph"]
MASSES=np.array([11.17793,12.11221],dtype=float)
Z=6.0

def w68(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    q=np.percentile(x,[16,84]); return float((q[1]-q[0])/2)

def auc(score,label):
    score=np.asarray(score,float);label=np.asarray(label,np.int8)
    ok=np.isfinite(score);score=score[ok];label=label[ok]
    o=np.argsort(score,kind="mergesort");s=score[o];y=label[o]
    n1=int(np.sum(y==1));n0=int(np.sum(y==0))
    st=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1];en=np.r_[st[1:],len(s)]
    ranks=np.repeat((st+1.0+en)/2.0,en-st)
    return float((np.sum(ranks[y==1])-n1*(n1+1.0)/2.0)/(n0*n1))

def load_stage(stage):
    chunks=[]
    for b in range(4):
        p=HERE/"tmva_strict"/("%s_bin%d.root"%(stage,b))
        with uproot.open(str(p)) as f:
            key=next(k for k in f.keys(recursive=True) if "TestTree" in k)
            chunks.append(f[key].arrays(["event_index","target","BDTG","label","source"],library="np"))
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}

def metric(err,raw,mask):
    return {
        "n":int(np.sum(mask)),
        "median":float(np.median(err[mask])),
        "w68":w68(err[mask]),
        "rms":float(np.sqrt(np.mean(err[mask]**2))),
        "gain_pct":100*(1-w68(err[mask])/w68(raw[mask])),
    }

def main():
    with np.load(str(MASS),allow_pickle=False) as z:
        roles=np.char.decode(z["role"],"ascii"); test=roles=="test"; idx=np.flatnonzero(test)
        base={k:z[k][test] for k in ["label","source","bin","R","R_true","u","u_true"]}
    y=base["label"].astype(np.int8)
    raw=np.log(base["R"]/base["R_true"])
    idxmap={int(e):i for i,e in enumerate(idx)}
    pred={"raw":np.zeros(len(idx))}
    for stage in STAGES:
        a=load_stage(stage)
        p=np.full(len(idx),np.nan)
        for ei,v in zip(a["event_index"].astype(np.int64),a["BDTG"]): p[idxmap[int(ei)]]=float(v)
        if np.any(~np.isfinite(p)): raise RuntimeError("missing "+stage)
        pred[stage]=p
    methods=["raw"]+STAGES
    result={"overall":{},"by_isotope":{},"bins":{},"mass":{}}
    for m in methods:
        err=raw-pred[m]
        corrR=base["R"]*np.exp(-pred[m])
        result["overall"][m]=metric(err,raw,np.ones(len(raw),bool))
        result["overall"][m]["correction_label_auc"]=auc(pred[m],y) if m!="raw" else float("nan")
        result["by_isotope"][m]={}
        for lab,name in [(0,"C12"),(1,"C13")]:
            result["by_isotope"][m][name]=metric(err,raw,y==lab)
        result["overall"][m]["median_C12_minus_C13"]=result["by_isotope"][m]["C12"]["median"]-result["by_isotope"][m]["C13"]["median"]
        result["bins"][m]={}
        for b in range(4): result["bins"][m][str(b)]=metric(err,raw,base["bin"]==b)
        for beta_name,u in [("raw_beta",base["u"]),("truth_beta",base["u_true"])]:
            m2=(Z*corrR)**2*(u*u-1)
            e=m2-MASSES[y]**2
            result["mass"][m+"__"+beta_name]={
                "auc":auc(m2,y),"w68_mass2_error":w68(e),"median_mass2_error":float(np.median(e))
            }

    # Paired whole-source bootstrap.
    groups={};keys={0:[],1:[]}
    for lab in (0,1):
        for src in sorted(int(x) for x in np.unique(base["source"][y==lab])):
            k=(lab,src);groups[k]=np.flatnonzero((y==lab)&(base["source"]==src));keys[lab].append(k)
    rng=np.random.RandomState(300932); reps=1000
    gains={m:[] for m in STAGES}; bin_gains={m:[[] for _ in range(4)] for m in STAGES}
    errs={m:raw-pred[m] for m in STAGES}
    for _ in range(reps):
        draw=np.concatenate([groups[keys[lab][j]] for lab in (0,1) for j in rng.randint(0,len(keys[lab]),size=len(keys[lab]))])
        rw=w68(raw[draw])
        for m in STAGES:
            gains[m].append(100*(1-w68(errs[m][draw])/rw))
            for b in range(4):
                ib=draw[base["bin"][draw]==b]
                bin_gains[m][b].append(100*(1-w68(errs[m][ib])/w68(raw[ib])))
    result["bootstrap"]={}
    for m in STAGES:
        a=np.asarray(gains[m]); item={"overall":{"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))},"bins":[]}
        for b in range(4):
            q=np.asarray(bin_gains[m][b])
            item["bins"].append({"bin":b,"median":float(np.median(q)),"low95":float(np.percentile(q,2.5)),"high95":float(np.percentile(q,97.5))})
        result["bootstrap"][m]=item
    result["definition"]="Strict curvature-only excludes reconstructed R magnitude, GetAdvancedFitMass/Charge, tracker charge magnitude/RMS/hit charge, other-hit energy fraction, and all dE/dx-like cluster charge. The final strict_morph adds only cluster COG and external-hit distance to fit/residual/geometry information."
    (HERE/"strict_results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("STRICT SUMMARY")
    for m in methods:
        r=result["overall"][m]
        print("%-16s W68 %.6f gain %6.2f%% med %+0.5f gap %+0.5f labelAUC %s"%(
            m,r["w68"],r["gain_pct"],r["median"],r["median_C12_minus_C13"],
            ("%.3f"%r["correction_label_auc"] if np.isfinite(r["correction_label_auc"]) else "-")))
    print("MASS raw beta")
    for m in methods:
        q=result["mass"][m+"__raw_beta"]
        print("%-16s AUC %.6f W68m2 %.3f med %+0.3f"%(m,q["auc"],q["w68_mass2_error"],q["median_mass2_error"]))
    print("BOOT",result["bootstrap"])

if __name__=="__main__": main()
