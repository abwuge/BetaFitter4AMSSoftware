#!/usr/bin/env python3
from __future__ import print_function

import csv
import json
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
INPUT = HERE / "tkadc_shape_features.root"
TMVA = HERE / "tmva_tkadc"
STAGES = ["tkadc_base","tkadc_adcshape","tkadc_adcstrip","tkadc_morph","tkadc_shape_only"]
Z = 6.0
MASS_C12 = 11.17793


def w68(values):
    values=np.asarray(values,dtype=float)
    values=values[np.isfinite(values)]
    q=np.percentile(values,[16.0,84.0])
    return float((q[1]-q[0])/2.0)


def rms(values):
    values=np.asarray(values,dtype=float)
    values=values[np.isfinite(values)]
    return float(np.sqrt(np.mean(values*values)))


def corr(x,y):
    x=np.asarray(x,dtype=float); y=np.asarray(y,dtype=float)
    ok=np.isfinite(x)&np.isfinite(y)
    x=x[ok]; y=y[ok]
    if len(x)<3: return float("nan")
    x=x-np.mean(x); y=y-np.mean(y)
    den=np.sqrt(np.sum(x*x)*np.sum(y*y))
    return float(np.sum(x*y)/den) if den>0 else float("nan")


def within_bin_corr(x,y,bins):
    xs=[]; ys=[]
    for b in range(4):
        take=(bins==b)&np.isfinite(x)&np.isfinite(y)
        if np.sum(take)<3: continue
        xx=np.asarray(x[take],float); yy=np.asarray(y[take],float)
        xs.append(xx-np.mean(xx)); ys.append(yy-np.mean(yy))
    if not xs: return float("nan")
    return corr(np.concatenate(xs),np.concatenate(ys))


def metric(error, raw_error, mask):
    e=np.asarray(error[mask],dtype=float)
    r=np.asarray(raw_error[mask],dtype=float)
    return {
        "n":int(np.sum(mask)),
        "median":float(np.median(e)),
        "w68":w68(e),
        "rms":rms(e),
        "gain_vs_raw_pct":100.0*(1.0-w68(e)/w68(r)),
    }


def load_tmva(stage, tree_token):
    chunks=[]
    for b in range(4):
        path=TMVA/("%s_bin%d.root"%(stage,b))
        with uproot.open(str(path)) as f:
            key=next(k for k in f.keys(recursive=True) if tree_token in k)
            chunks.append(f[key].arrays(["event_index","target","BDTG","source"],library="np"))
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}


def write_tsv(path, rows):
    if not rows: return
    with open(str(path),"w",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0].keys()),delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def bootstrap(test_data, predictions, replicates=1000):
    source=test_data["source"].astype(np.int32)
    bins=test_data["bin"].astype(np.int8)
    raw=test_data["target"].astype(float)
    keys=sorted(int(x) for x in np.unique(source))
    groups={s:np.flatnonzero(source==s) for s in keys}
    rng=np.random.RandomState(300933)
    direct={stage:[] for stage in STAGES}
    incremental={stage:[] for stage in STAGES if stage!="tkadc_base"}
    direct_bins={stage:[[] for _ in range(4)] for stage in STAGES}
    increment_bins={stage:[[] for _ in range(4)] for stage in STAGES if stage!="tkadc_base"}
    base_error=raw-predictions["tkadc_base"]
    for _ in range(replicates):
        draw=np.concatenate([groups[keys[j]] for j in rng.randint(0,len(keys),size=len(keys))])
        raw_w=w68(raw[draw]); base_w=w68(base_error[draw])
        for stage in STAGES:
            err=raw-predictions[stage]
            direct[stage].append(100.0*(1.0-w68(err[draw])/raw_w))
            if stage!="tkadc_base":
                incremental[stage].append(100.0*(1.0-w68(err[draw])/base_w))
            for b in range(4):
                ib=draw[bins[draw]==b]
                if len(ib)<10: continue
                direct_bins[stage][b].append(100.0*(1.0-w68(err[ib])/w68(raw[ib])))
                if stage!="tkadc_base":
                    increment_bins[stage][b].append(100.0*(1.0-w68(err[ib])/w68(base_error[ib])))
    out={"replicates":replicates,"direct_vs_raw":{},"incremental_vs_base":{}}
    for stage,values in direct.items():
        a=np.asarray(values)
        out["direct_vs_raw"][stage]={
            "median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),
            "high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0)),
            "bins":[]
        }
        for b in range(4):
            q=np.asarray(direct_bins[stage][b])
            out["direct_vs_raw"][stage]["bins"].append({
                "bin":b,"median":float(np.median(q)),"low95":float(np.percentile(q,2.5)),
                "high95":float(np.percentile(q,97.5))
            })
    for stage,values in incremental.items():
        a=np.asarray(values)
        out["incremental_vs_base"][stage]={
            "median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),
            "high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0),
            ),"bins":[]
        }
        for b in range(4):
            q=np.asarray(increment_bins[stage][b])
            out["incremental_vs_base"][stage]["bins"].append({
                "bin":b,"median":float(np.median(q)),"low95":float(np.percentile(q,2.5)),
                "high95":float(np.percentile(q,97.5))
            })
    return out


def bootstrap_mass(test_data, predictions, replicates=1000):
    source=test_data["source"].astype(np.int32)
    keys=sorted(int(s) for s in np.unique(source))
    groups={s:np.flatnonzero(source==s) for s in keys}
    errors={}
    draws={}
    for beta_name, branch in [("raw_beta", "u"), ("truth_beta", "u_true")]:
        u=test_data[branch].astype(float)
        errors[beta_name]={}
        draws[beta_name]={stage:[] for stage in STAGES if stage!="tkadc_base"}
        for stage in STAGES:
            rigidity=test_data["R"].astype(float)*np.exp(-predictions[stage])
            errors[beta_name][stage]=(Z*rigidity)**2*(u*u-1)-MASS_C12**2
    rng=np.random.RandomState(300933)
    for _ in range(replicates):
        ix=np.concatenate([groups[keys[j]] for j in rng.randint(0,len(keys),size=len(keys))])
        for beta_name in errors:
            denominator=w68(errors[beta_name]["tkadc_base"][ix])
            for stage in draws[beta_name]:
                draws[beta_name][stage].append(100*(1-w68(errors[beta_name][stage][ix])/denominator))
    report={"replicates":replicates,"source_files":len(keys),"seed":300933,"incremental_vs_base":{}}
    for beta_name in draws:
        report["incremental_vs_base"][beta_name]={}
        for stage, values in draws[beta_name].items():
            a=np.asarray(values)
            report["incremental_vs_base"][beta_name][stage]={
                "direct":100*(1-w68(errors[beta_name][stage])/w68(errors[beta_name]["tkadc_base"])),
                "median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),
                "high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0)),
            }
    return report


def main():
    with uproot.open(str(INPUT)) as f:
        tree=f["events"]
        names=["event_index","role","source","bin","R","R_true","u","u_true","target","adc_log_scale","qlqs_log_scale"]
        full=tree.arrays(names,library="np")
    role=full["role"].astype(np.int8)
    test_mask=role==2
    train_mask=role<2
    test={k:v[test_mask] for k,v in full.items()}
    train={k:v[train_mask] for k,v in full.items()}
    if len(np.unique(full["event_index"])) != len(full["event_index"]):
        raise RuntimeError("duplicate feature event identities")
    source_roles = [set(full["source"][role == r].tolist()) for r in range(3)]
    if any(source_roles[i] & source_roles[j] for i in range(3) for j in range(i)):
        raise RuntimeError("source file overlap across train/valid/test")
    for r, name in enumerate(["train", "valid", "test"]):
        with np.load(str(HERE.parent/"c_mass_toy_oracle"/"predictions.npz")) as original:
            expected = np.flatnonzero((original["label"] == 0) & (original["role"] == name.encode()))
        if not np.array_equal(np.sort(full["event_index"][role == r]), expected):
            raise RuntimeError("changed original sample for " + name)

    predictions_test={}
    predictions_train={}
    for stage in STAGES:
        for token,destination,reference in [
            ("TestTree",predictions_test,test),
            ("TrainTree",predictions_train,train),
        ]:
            a=load_tmva(stage,token)
            if len(a["event_index"]) != len(reference["event_index"]) or not np.array_equal(
                np.sort(a["event_index"]), np.sort(reference["event_index"])
            ):
                raise RuntimeError("prediction coverage mismatch %s %s" % (stage, token))
            if len(np.unique(a["event_index"])) != len(a["event_index"]):
                raise RuntimeError("duplicate model predictions")
            index_map={int(e):i for i,e in enumerate(reference["event_index"])}
            p=np.full(len(reference["event_index"]),np.nan)
            for event,value in zip(a["event_index"].astype(np.int64),a["BDTG"]):
                p[index_map[int(event)]]=float(value)
            if np.any(~np.isfinite(p)):
                raise RuntimeError("missing predictions %s %s"%(stage,token))
            destination[stage]=p

    rows=[]
    result={"test":{},"train":{},"bins":{},"diagnostics":{},"mass":{}}
    result["sample_audit"] = {
        "role_counts": {name: int(np.sum(role == r)) for r, name in enumerate(["train", "valid", "test"])},
        "source_counts": {name: len(source_roles[r]) for r, name in enumerate(["train", "valid", "test"])},
        "source_overlap": False, "original_event_coverage": 1.0,
        "training_roles": "train+valid", "preprocessor_fit_role": "train",
        "test_used_for_training": False,
    }
    for split, data in [("test", test), ("train", train)]:
        result[split]["raw"] = metric(data["target"], data["target"], np.ones(len(data["target"]), dtype=bool))
    for beta_name, u in [("raw_beta", test["u"]), ("truth_beta", test["u_true"])]:
        error = (Z * test["R"])**2 * (u*u - 1.0) - MASS_C12**2
        result["mass"]["raw__" + beta_name] = {"median_mass2_error": float(np.median(error)), "w68_mass2_error": w68(error)}
    for split_name,data,predictions in [("test",test,predictions_test),("train",train,predictions_train)]:
        raw=data["target"].astype(float)
        base_error=raw-predictions["tkadc_base"]
        for stage in STAGES:
            error=raw-predictions[stage]
            m=metric(error,raw,np.ones(len(raw),dtype=bool))
            m["incremental_gain_vs_base_pct"]=100.0*(1.0-w68(error)/w68(base_error))
            result[split_name][stage]=m
            if split_name=="test":
                result["bins"][stage]={}
                for b in range(4):
                    bm=metric(error,raw,data["bin"]==b)
                    bm["incremental_gain_vs_base_pct"]=100.0*(
                        1.0-w68(error[data["bin"]==b])/w68(base_error[data["bin"]==b])
                    )
                    result["bins"][stage][str(b)]=bm
                pred=predictions[stage]
                delta=pred-predictions["tkadc_base"]
                result["diagnostics"][stage]={
                    "prediction_corr_target_within_Rbin":within_bin_corr(pred,raw,data["bin"]),
                    "increment_corr_base_residual_within_Rbin":within_bin_corr(delta,base_error,data["bin"]),
                    "increment_corr_adc_scale_within_Rbin":within_bin_corr(delta,data["adc_log_scale"],data["bin"]),
                    "increment_corr_qlqs_scale_within_Rbin":within_bin_corr(delta,data["qlqs_log_scale"],data["bin"]),
                    "increment_corr_u_true_within_Rbin":within_bin_corr(delta,data["u_true"],data["bin"]),
                }
                corrected_R=data["R"]*np.exp(-pred)
                for beta_name,u in [("raw_beta",data["u"]),("truth_beta",data["u_true"])]:
                    mass2=(Z*corrected_R)**2*(u*u-1.0)
                    e=mass2-MASS_C12**2
                    result["mass"][stage+"__"+beta_name]={
                        "median_mass2_error":float(np.median(e)),
                        "w68_mass2_error":w68(e),
                    }
    for stage in STAGES:
        for split in ["test","train"]:
            row={"stage":stage,"split":split}
            row.update(result[split][stage]); rows.append(row)
    write_tsv(HERE/"tkadc_metrics.tsv",rows)

    boot=bootstrap(test,predictions_test,replicates=1000)
    result["bootstrap"]=boot
    mass_boot=bootstrap_mass(test,predictions_test,replicates=1000)
    result["mass_bootstrap"]=mass_boot
    (HERE/"tkadc_mass_bootstrap.json").write_text(json.dumps(mass_boot,indent=2,sort_keys=True)+"\n")
    plot_rows=[]
    for comparison, key in [("raw", "direct_vs_raw"), ("base", "incremental_vs_base")]:
        for stage, item in boot[key].items():
            for b, interval in enumerate(item["bins"]):
                metric_key = "gain_vs_raw_pct" if comparison == "raw" else "incremental_gain_vs_base_pct"
                plot_rows.append({"comparison": comparison, "stage": stage, "bin": b,
                                 "direct": result["bins"][stage][str(b)][metric_key],
                                 "median": interval["median"], "low95": interval["low95"], "high95": interval["high95"]})
    write_tsv(HERE/"tkadc_gain_plot.tsv", plot_rows)
    payload = {k: v for k, v in test.items()}
    for stage, prediction in predictions_test.items():
        payload[stage] = prediction.astype(np.float32)
    with uproot.recreate(str(HERE/"tkadc_test_predictions.root")) as output:
        output["events"] = payload
    result["feature_counts"]={}
    for stage,filename in [
        ("tkadc_base","features_tkadc_base.txt"),
        ("tkadc_adcshape","features_tkadc_adcshape.txt"),
        ("tkadc_adcstrip","features_tkadc_adcstrip.txt"),
        ("tkadc_morph","features_tkadc_morph.txt"),
        ("tkadc_shape_only","features_tkadc_shape_only.txt"),
    ]:
        result["feature_counts"][stage]=len([x for x in (HERE/filename).read_text().splitlines() if x.strip()])
    result["definition"]=(
        "C12-only, whole-source-file train/valid/test. The primary new-feature tests remove absolute "
        "tracker charge magnitude: ADC is normalized within layer/side, qlqs uses only centered log ratios, "
        "and stripx supplies local morphology/address metadata. Omitted absolute ADC/qlqs scales are spectators "
        "used only to diagnose kinematic leakage."
    )
    (HERE/"tkadc_results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    (HERE/"tkadc_bootstrap.json").write_text(json.dumps(boot,indent=2,sort_keys=True)+"\n")

    print("TKADC TEST")
    for stage in STAGES:
        m=result["test"][stage]
        d=result["diagnostics"][stage]
        print("%-18s W68 %.6f rawgain %6.2f%% incbase %+6.2f%% med %+0.5f  dCorrADC %+0.3f dCorrU %+0.3f"%(
            stage,m["w68"],m["gain_vs_raw_pct"],m["incremental_gain_vs_base_pct"],m["median"],
            d["increment_corr_adc_scale_within_Rbin"],d["increment_corr_u_true_within_Rbin"]))
    print("BOOT incremental")
    for stage,item in boot["incremental_vs_base"].items():
        print(stage,item["median"],item["low95"],item["high95"],item["positive_fraction"])


if __name__=="__main__":
    main()
