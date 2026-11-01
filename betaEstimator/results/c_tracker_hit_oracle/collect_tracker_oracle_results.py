#!/usr/bin/env python3
from __future__ import print_function
import csv, json, math
from pathlib import Path
import numpy as np
import uproot

HERE=Path(__file__).resolve().parent
MASS_TOY=HERE.parent/"c_mass_toy_oracle"
PRED=MASS_TOY/"predictions.npz"
STAGES=["multifit","fitdiag","geometry","cluster"]
CHARGE=6.0
MASSES=np.array([11.17793,12.11221],dtype=float)
R_EDGES=np.array([2.5,3.0,4.0,5.0,6.0],dtype=float)

def w68(x):
    x=np.asarray(x,dtype=float); x=x[np.isfinite(x)]
    q=np.percentile(x,[16,84]); return float((q[1]-q[0])/2)

def rms(x):
    x=np.asarray(x,dtype=float); x=x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x*x)))

def auc(score,label):
    score=np.asarray(score,float); label=np.asarray(label,np.int8)
    ok=np.isfinite(score); score=score[ok]; label=label[ok]
    n1=int(np.sum(label==1)); n0=int(np.sum(label==0))
    order=np.argsort(score,kind="mergesort")
    s=score[order]; y=label[order]
    starts=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1]; ends=np.r_[starts[1:],len(s)]
    ranks=np.repeat((starts+1.0+ends)/2.0,ends-starts)
    return float((np.sum(ranks[y==1])-n1*(n1+1)/2.0)/(n0*n1))

def rejection50(score,label):
    score=np.asarray(score,float); label=np.asarray(label,np.int8)
    s1=np.sort(score[label==1])
    threshold=np.quantile(s1,0.5)
    acc=np.mean(score[label==0]>=threshold)
    return float(1.0/max(acc,1.0/max(1,np.sum(label==0))))

def write_tsv(path,rows):
    if not rows:return
    with open(str(path),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def load_tmva_predictions():
    pred={}
    for stage in STAGES:
        local={}
        for b in range(4):
            path=HERE/"tmva_full"/("%s_bin%d.root"%(stage,b))
            if not path.exists(): raise RuntimeError("missing "+str(path))
            with uproot.open(str(path)) as f:
                key=next(k for k in f.keys(recursive=True) if "TestTree" in k)
                a=f[key].arrays(["event_index","BDTG","target","label","source","R","R_true","u","u_true"],library="np")
            for i,e in enumerate(a["event_index"].astype(np.int64)):
                if e in local: raise RuntimeError("duplicate event "+str(e))
                local[int(e)]=float(a["BDTG"][i])
        pred[stage]=local
    return pred

def metrics_for(method,dlogr,label,binid,mask):
    rows=[]
    raw_width={}
    for b in [-1,0,1,2,3]:
        for labname,labmask in [("all",np.ones(len(mask),bool)),("C12",label==0),("C13",label==1)]:
            sel=mask & labmask & ((binid==b) if b>=0 else True)
            x=dlogr[sel]
            rows.append({"method":method,"bin":b,"sample":labname,"n":int(np.sum(sel)),
                         "median_dlogR":float(np.median(x)),"w68_dlogR":w68(x),"rms_dlogR":rms(x)})
    return rows

def source_bootstrap(method_errors,label,source,binid,test,reps=1000):
    groups={}
    keys={0:[],1:[]}
    for lab in (0,1):
        for src in sorted(int(x) for x in np.unique(source[test&(label==lab)])):
            key=(lab,src); groups[key]=np.flatnonzero(test&(label==lab)&(source==src)); keys[lab].append(key)
    rng=np.random.RandomState(300930)
    records={m:{"overall":[],"bins":[[] for _ in range(4)]} for m in method_errors if m!="raw"}
    for _ in range(reps):
        idx=np.concatenate([groups[keys[lab][j]] for lab in (0,1) for j in rng.randint(0,len(keys[lab]),size=len(keys[lab]))])
        raw=w68(method_errors["raw"][idx])
        for m in records:
            records[m]["overall"].append(100*(1-w68(method_errors[m][idx])/raw))
            for b in range(4):
                ib=idx[binid[idx]==b]
                rb=w68(method_errors["raw"][ib])
                records[m]["bins"][b].append(100*(1-w68(method_errors[m][ib])/rb))
    out={}
    for m,v in records.items():
        arr=np.asarray(v["overall"])
        out[m]={"overall_gain_pct":{"median":float(np.median(arr)),"low95":float(np.percentile(arr,2.5)),"high95":float(np.percentile(arr,97.5)),"positive_fraction":float(np.mean(arr>0))},"bins":[]}
        for b,a in enumerate(v["bins"]):
            a=np.asarray(a)
            out[m]["bins"].append({"bin":b,"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))})
    return out

def main():
    with np.load(str(PRED),allow_pickle=False) as z:
        d={k:z[k] for k in ["run","event","source","label","bin","R","R_true","u","u_true","role"]}
    role=np.char.decode(d["role"],"ascii"); test=role=="test"
    tmva=load_tmva_predictions()
    pred_error={m:np.full(len(test),np.nan,dtype=float) for m in STAGES}
    for m in STAGES:
        for idx,val in tmva[m].items(): pred_error[m][idx]=val
        if np.any(~np.isfinite(pred_error[m][test])): raise RuntimeError("missing test predictions "+m)
    target=np.log(d["R"]/d["R_true"])
    errors={"raw":target.copy()}
    rigidity={"raw":d["R"].astype(float).copy()}
    for m in STAGES:
        errors[m]=target-pred_error[m]
        rigidity[m]=d["R"]*np.exp(-pred_error[m])

    metric_rows=[]
    for m in ["raw"]+STAGES: metric_rows.extend(metrics_for(m,errors[m],d["label"],d["bin"],test))
    raw_lookup={(r["bin"],r["sample"]):r["w68_dlogR"] for r in metric_rows if r["method"]=="raw"}
    for r in metric_rows:
        r["improvement_vs_raw_pct"]=100*(1-r["w68_dlogR"]/raw_lookup[(r["bin"],r["sample"])])
    write_tsv(HERE/"rigidity_metrics.tsv",metric_rows)

    mass_rows=[]
    true_m2=MASSES[d["label"]]**2
    for beta_name,u in [("raw_beta",d["u"]),("truth_beta",d["u_true"])]:
        for m in ["raw"]+STAGES:
            m2=(CHARGE*rigidity[m])**2*(u*u-1)
            for b in [-1,0,1,2,3]:
                sel=test & ((d["bin"]==b) if b>=0 else True)
                err=m2[sel]-true_m2[sel]
                mass_rows.append({"beta":beta_name,"method":m,"bin":b,"n":int(np.sum(sel)),
                                  "mass2_auc":auc(m2[sel],d["label"][sel]),
                                  "rejection50":rejection50(m2[sel],d["label"][sel]),
                                  "median_mass2_error":float(np.median(err)),
                                  "w68_mass2_error":w68(err)})
    write_tsv(HERE/"mass_metrics.tsv",mass_rows)

    bootstrap=source_bootstrap(errors,d["label"],d["source"],d["bin"],test,1000)
    (HERE/"bootstrap.json").write_text(json.dumps(bootstrap,indent=2,sort_keys=True)+"\n")

    payload={"event_index":np.flatnonzero(test),"run":d["run"][test],"event":d["event"][test],
             "source":d["source"][test],"label":d["label"][test],"bin":d["bin"][test],
             "R_raw":d["R"][test],"R_true":d["R_true"][test],"u":d["u"][test],"u_true":d["u_true"][test],
             "dlogR_raw":errors["raw"][test]}
    for m in STAGES:
        payload["pred_dlogR_"+m]=pred_error[m][test]
        payload["R_"+m]=rigidity[m][test]
        payload["dlogR_"+m]=errors[m][test]
    np.savez_compressed(str(HERE/"predictions.npz"),**payload)

    overall={r["method"]:r for r in metric_rows if r["bin"]==-1 and r["sample"]=="all"}
    byiso={m:{s:next(r for r in metric_rows if r["method"]==m and r["bin"]==-1 and r["sample"]==s) for s in ["C12","C13"]} for m in ["raw"]+STAGES}
    mass_overall=[r for r in mass_rows if r["bin"]==-1]
    result={"test_events":int(np.sum(test)),"test_c12":int(np.sum(test&(d["label"]==0))),"test_c13":int(np.sum(test&(d["label"]==1))),
            "overall_rigidity":overall,"by_isotope":byiso,"mass_overall":mass_overall,
            "bootstrap":bootstrap,
            "training_note":"TMVA BDTG, per reconstructed-R bin, trained on train+validation whole files; test source files were never used for fitting. Huber loss, 300 trees, randomized feature subset.",
            "limitation":"C12 B1308 and C13 B1315 remain different productions; this is an MC oracle information-bound diagnostic, not a flight-transfer result."}
    (HERE/"results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("OVERALL")
    for m in ["raw"]+STAGES:
        r=overall[m]
        print("%-10s W68 %.6f gain %6.2f%% median %+0.5f RMS %.6f"%(m,r["w68_dlogR"],r["improvement_vs_raw_pct"],r["median_dlogR"],r["rms_dlogR"]))
    print("MASS raw-beta AUC")
    for m in ["raw"]+STAGES:
        r=next(x for x in mass_overall if x["beta"]=="raw_beta" and x["method"]==m)
        print("%-10s AUC %.6f W68m2 %.3f rejection50 %.3f"%(m,r["mass2_auc"],r["w68_mass2_error"],r["rejection50"]))

if __name__=="__main__": main()
