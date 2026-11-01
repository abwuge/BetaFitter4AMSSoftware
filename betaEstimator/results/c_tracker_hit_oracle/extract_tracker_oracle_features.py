#!/usr/bin/env python3
"""Extract progressively deeper tracker features for a truth-supervised rigidity upper-bound scan.

This is deliberately an MC-oracle diagnostic.  MC truth is used only for the
regression target and final evaluation.  Input variables are reconstructed
tracker quantities common to the C12 B1308 and C13 B1315 ntuples.
"""
from __future__ import print_function

import csv
import json
import math
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
MASS_TOY = HERE.parent / "c_mass_toy_oracle"
PRED = MASS_TOY / "predictions.npz"
MANIFEST = MASS_TOY / "split_manifest.tsv"
SEED = 300930

GROUPS = ["multifit", "fitdiag", "geometry", "cluster"]
BRANCHES = [
    "run","event",
    "tk_rigidity1","tk_rigidityk","tk_erigidity1","tk_rigidityi","tk_chisi",
    "tk_chis1","tk_bcor","tk_cdif","tk_fzm","tk_qrms","tk_qhit","tk_hitb",
    "tk_res","tk_res1","tk_hitc","tk_pos","tk_dir","tk_oel","tk_ohitl",
    "tk_rigiditys","tk_chiss","tk_pos2s","tk_dir2s","tk_hitbs",
    "tk_q","tk_ql","tk_qln","tk_lid","tk_cad","tk_cof","tk_iso",
    "tk_exql","tk_exqln","tk_exqls","tk_exdis","tk_exlid","tk_excad",
    "tk_qis","tk_qirmss","tk_qlss",
]


def log(msg):
    print(msg, flush=True)


def write_tsv(path, rows):
    if not rows:
        return
    with open(str(path), "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def flatten_named(array, prefix):
    array = np.asarray(array, dtype=np.float64)
    flat = array.reshape(len(array), -1)
    shape = array.shape[1:]
    names = []
    for idx in np.ndindex(shape):
        names.append(prefix + "_" + "_".join(str(i) for i in idx))
    if not shape:
        names = [prefix]
        flat = array.reshape(len(array), 1)
    return flat, names


def add(cols, names, array, prefix):
    flat, local_names = flatten_named(array, prefix)
    cols.append(flat)
    names.extend(local_names)


def clean_finite(array, limit=1000.0):
    x = np.asarray(array, dtype=np.float64).copy()
    x[~np.isfinite(x)] = np.nan
    x[np.abs(x) >= limit] = np.nan
    return x


def clean_positive_log(array, upper=1000.0):
    x = np.asarray(array, dtype=np.float64).copy()
    good = np.isfinite(x) & (x > 0.0) & (x < upper)
    out = np.full(x.shape, np.nan, dtype=np.float64)
    out[good] = np.log(x[good])
    return out


def log_ratio(array, baseline):
    x = np.asarray(array, dtype=np.float64)
    base = np.asarray(baseline, dtype=np.float64)
    expand = (len(base),) + (1,) * (x.ndim - 1)
    b = base.reshape(expand)
    good = np.isfinite(x) & np.isfinite(b) & (np.abs(x) > 1e-8) & (b > 0.0)
    out = np.full(x.shape, np.nan, dtype=np.float64)
    out[good] = np.log(np.abs(x[good] / np.broadcast_to(b, x.shape)[good]))
    out[np.abs(out) > 1.5] = np.nan
    return out


def bit_columns(values, bits=9):
    x = np.asarray(values, dtype=np.int64)
    out = []
    for axis in range(x.shape[1]):
        for bit in range(bits):
            out.append(((x[:, axis] >> bit) & 1).astype(np.float64))
    return np.column_stack(out)


def cad_features(cad, prefix):
    cad = np.asarray(cad, dtype=np.int64)
    n, layers = cad.shape
    vals = []
    names = []
    for layer in range(layers):
        c = cad[:, layer]
        valid = c >= 0
        head = np.where(valid, c // 100000000, -1)
        xaddr = np.where(valid, c % 10000, -1)
        yaddr = np.where(valid, (c // 10000) % 10000, -1)
        k7 = np.where(valid, head // 10, -1)
        mult = np.where(valid, head % 10, -1)
        for suffix, v in [("k7",k7),("mult",mult),("xaddr",xaddr),("yaddr",yaddr)]:
            a = v.astype(np.float64)
            a[~valid] = np.nan
            vals.append(a)
            names.append("%s_L%d_%s" % (prefix, layer+1, suffix))
    return np.column_stack(vals), names


def excad_features(cad):
    cad = np.asarray(cad, dtype=np.int64)
    vals = []
    names = []
    for side in range(cad.shape[1]):
        c = cad[:, side]
        valid = c >= 0
        xaddr = np.where(valid, c % 10000, -1).astype(np.float64)
        yaddr = np.where(valid, (c // 10000) % 10000, -1).astype(np.float64)
        xaddr[~valid] = np.nan
        yaddr[~valid] = np.nan
        vals.extend([xaddr,yaddr])
        names.extend(["excad_%d_xaddr"%side, "excad_%d_yaddr"%side])
    return np.column_stack(vals), names


def make_groups(a, baseline_r):
    n = len(baseline_r)
    out = {}
    out_names = {}

    # Stage 1: close to the previous multi-fit oracle information boundary.
    cols=[]; names=[]
    add(cols,names,log_ratio(a["tk_rigidity1"],baseline_r),"logratio_r1")
    add(cols,names,log_ratio(a["tk_rigidityk"],baseline_r),"logratio_rk")
    add(cols,names,log_ratio(a["tk_rigidityi"],baseline_r),"logratio_leave1out")
    relerr = np.abs(np.asarray(a["tk_erigidity1"],dtype=np.float64) * baseline_r[:,None,None])
    add(cols,names,clean_positive_log(relerr,upper=10.0),"log_relerr")
    add(cols,names,clean_positive_log(a["tk_chisi"],upper=1000.0),"log_chisi")
    add(cols,names,clean_finite(a["tk_bcor"],limit=10.0),"bcor")
    add(cols,names,clean_finite(a["tk_cdif"],limit=100.0),"pg_ciemat_diff")
    add(cols,names,clean_finite(a["tk_fzm"],limit=1e3),"advanced_fit")
    add(cols,names,clean_finite(a["tk_qrms"],limit=1e3),"qrms")
    add(cols,names,clean_finite(a["tk_qhit"],limit=1e3),"qhit")
    add(cols,names,np.log(np.asarray(baseline_r,dtype=np.float64)),"logR")
    out["multifit"]=np.column_stack(cols); out_names["multifit"]=names

    # Stage 2: per-fit chi2, residual/influence and second reconstructed fit.
    cols=[]; names=[]
    chis = np.asarray(a["tk_chis1"],dtype=np.float64)[...,:2]
    add(cols,names,clean_positive_log(chis,upper=1000.0),"log_chis1_xy")
    res=clean_finite(a["tk_res"],limit=50.0)
    res1=clean_finite(a["tk_res1"],limit=50.0)
    add(cols,names,res,"res_full")
    add(cols,names,res1,"res_leave1out")
    add(cols,names,res1-res,"res_influence")
    add(cols,names,log_ratio(a["tk_rigiditys"],baseline_r),"logratio_particle_refit")
    add(cols,names,clean_positive_log(a["tk_chiss"],upper=1000.0),"log_chiss")
    hb=bit_columns(a["tk_hitb"],9)
    add(cols,names,hb,"hitbit_primary")
    hbs=bit_columns(a["tk_hitbs"],9)
    add(cols,names,hbs,"hitbit_secondary")
    add(cols,names,clean_finite(a["tk_oel"],limit=100.0),"other_hit_energy_fraction")
    add(cols,names,clean_finite(a["tk_ohitl"],limit=1e3),"other_hit_count")
    out["fitdiag"]=np.column_stack(cols); out_names["fitdiag"]=names

    # Stage 3: actual hit and fitted-track geometry.
    cols=[]; names=[]
    hit=clean_finite(a["tk_hitc"],limit=1000.0)
    pos=clean_finite(a["tk_pos"],limit=1000.0)
    direction=clean_finite(a["tk_dir"],limit=10.0)
    add(cols,names,hit-pos,"hit_minus_track_xyz")
    add(cols,names,hit,"hit_xyz")
    add(cols,names,pos,"track_pos_xyz")
    add(cols,names,direction,"track_dir_xyz")
    add(cols,names,clean_finite(a["tk_pos2s"],limit=1000.0),"particle_track_pos_xyz")
    add(cols,names,clean_finite(a["tk_dir2s"],limit=10.0),"particle_track_dir_xyz")
    out["geometry"]=np.column_stack(cols); out_names["geometry"]=names

    # Stage 4: charge/cluster morphology and local-address information.
    cols=[]; names=[]
    add(cols,names,clean_positive_log(a["tk_q"],upper=100.0),"log_q")
    add(cols,names,clean_positive_log(a["tk_ql"],upper=100.0),"log_ql")
    add(cols,names,clean_positive_log(a["tk_qln"],upper=100.0),"log_qln")
    add(cols,names,clean_finite(a["tk_cof"],limit=10.0),"cluster_cog")
    add(cols,names,clean_finite(a["tk_iso"],limit=100.0),"cluster_isolation")
    add(cols,names,clean_positive_log(a["tk_exql"],upper=100.0),"log_external_q")
    add(cols,names,clean_positive_log(a["tk_exqln"],upper=100.0),"log_external_qn")
    add(cols,names,clean_finite(a["tk_exdis"],limit=100.0),"external_hit_distance")
    add(cols,names,clean_positive_log(a["tk_qis"],upper=100.0),"log_particle_q")
    add(cols,names,clean_positive_log(a["tk_qirmss"],upper=100.0),"log_particle_qrms")
    add(cols,names,clean_positive_log(a["tk_qlss"],upper=100.0),"log_particle_layer_q")
    add(cols,names,clean_finite(a["tk_lid"],limit=1e9),"ladder_id")
    cad,cadn=cad_features(a["tk_cad"],"cad"); cols.append(cad); names.extend(cadn)
    exc,exn=excad_features(a["tk_excad"]); cols.append(exc); names.extend(exn)
    add(cols,names,clean_finite(a["tk_exlid"],limit=1e9),"external_ladder_id")
    add(cols,names,clean_finite(a["tk_exqls"],limit=1e9),"external_q_status")
    out["cluster"]=np.column_stack(cols); out_names["cluster"]=names
    return out,out_names


def extract_one(task):
    label,source,path,rows,target_run,target_event,target_r = task
    started=time.time()
    with uproot.open(path) as f:
        tree=f["amstreea"]
        missing=[b for b in BRANCHES if b not in tree.keys()]
        if missing:
            raise RuntimeError("%s missing branches %s"%(path,missing))
        a=tree.arrays(BRANCHES,library="np")
    ids=np.rec.fromarrays([a["run"].astype(np.int64),a["event"].astype(np.int64)],names="run,event")
    tids=np.rec.fromarrays([target_run.astype(np.int64),target_event.astype(np.int64)],names="run,event")
    common,ie,it=np.intersect1d(ids,tids,assume_unique=False,return_indices=True)
    if len(common)!=len(rows):
        raise RuntimeError("match failure label=%d source=%d got=%d need=%d"%(label,source,len(common),len(rows)))
    order=np.argsort(it); ie=ie[order]; rows_out=rows[it[order]]
    local={k:v[ie] for k,v in a.items()}
    branch_r=np.asarray(local["tk_rigidityk"][:,1,1],dtype=np.float64)
    r_replay=float(np.max(np.abs(branch_r-target_r[it[order]]))) if len(rows) else 0.0
    groups,names=make_groups(local,branch_r)
    audit={"label":label,"source":source,"path":path,"events":len(rows),"seconds":time.time()-started,"R_replay_max":r_replay}
    return rows_out,groups,names,audit


def preprocess(raw,names,train_mask):
    train=np.asarray(raw[train_mask],dtype=np.float64)
    finite_fraction=np.mean(np.isfinite(train),axis=0)
    keep=finite_fraction>=0.005
    indices=np.flatnonzero(keep)
    train=train[:,keep]
    names=[n for n,k in zip(names,keep) if k]
    finite_fraction=finite_fraction[keep]
    low=np.nanpercentile(train,0.5,axis=0)
    high=np.nanpercentile(train,99.5,axis=0)
    clipped=np.minimum(np.maximum(train,low),high)
    center=np.nanmedian(clipped,axis=0)
    q25,q75=np.nanpercentile(clipped,[25,75],axis=0)
    robust=0.741301109*(q75-q25)
    std=np.nanstd(clipped,axis=0)
    scale=np.maximum(robust,0.1*std)
    good=np.isfinite(center)&np.isfinite(scale)&(scale>1e-7)
    indices=indices[good]; low=low[good]; high=high[good]; center=center[good]; scale=scale[good]
    names=[n for n,g in zip(names,good) if g]
    finite_fraction=finite_fraction[good]
    miss=(1.0-finite_fraction)>=0.005
    values=np.asarray(raw[:,indices],dtype=np.float64)
    missing=~np.isfinite(values)
    values=np.minimum(np.maximum(values,low),high)
    values=np.where(missing,center,values)
    z=np.clip((values-center)/scale,-6,6)
    columns=[z.astype(np.float32)]
    out_names=list(names)
    if np.any(miss):
        columns.append(missing[:,miss].astype(np.float32))
        out_names += ["missing_"+names[i] for i in np.flatnonzero(miss)]
    state={"indices":indices.tolist(),"low":low.tolist(),"high":high.tolist(),"center":center.tolist(),"scale":scale.tolist(),"missing_indicator":miss.tolist(),"input_names":names,"output_names":out_names}
    return np.column_stack(columns),state


def main():
    started=time.time()
    with np.load(str(PRED),allow_pickle=False) as z:
        data={k:z[k] for k in ["run","event","source","label","bin","R","R_true","u","u_true","role"]}
    roles=np.char.decode(data["role"],"ascii")
    if not np.all(np.isfinite(data["R_true"])) or not np.all(data["R_true"]>0):
        raise RuntimeError("prediction input is not truth-valid")
    train_mask=roles=="train"
    rows_manifest=list(csv.DictReader(open(str(MANIFEST)),delimiter="\t"))
    iso_map={"C12":0,"C13":1}
    paths={(iso_map[r["isotope"]],int(r["source"])):r["path"] for r in rows_manifest}
    tasks=[]
    for label in (0,1):
        for source in sorted(int(x) for x in np.unique(data["source"][data["label"]==label])):
            rows=np.flatnonzero((data["label"]==label)&(data["source"]==source))
            tasks.append((label,source,paths[(label,source)],rows,data["run"][rows],data["event"][rows],data["R"][rows]))
    raw_groups=None; group_names=None; filled=np.zeros(len(roles),dtype=bool); audits=[]
    log("extracting %d files / %d events"%(len(tasks),len(roles)))
    with mp.Pool(processes=6) as pool:
        for done,(rows,local,names,audit) in enumerate(pool.imap_unordered(extract_one,tasks,chunksize=1),1):
            if raw_groups is None:
                group_names=names
                raw_groups={g:np.empty((len(roles),local[g].shape[1]),dtype=np.float32) for g in GROUPS}
            for g in GROUPS:
                if names[g]!=group_names[g]:
                    raise RuntimeError("schema mismatch "+g)
                raw_groups[g][rows]=local[g]
            filled[rows]=True; audits.append(audit)
            if done%24==0 or done==len(tasks):
                log("files %d/%d elapsed %.1fs"%(done,len(tasks),time.time()-started))
    if not np.all(filled):
        raise RuntimeError("unfilled events")
    write_tsv(HERE/"extraction_audit.tsv",sorted(audits,key=lambda x:(x["label"],x["source"])))
    max_replay=max(x["R_replay_max"] for x in audits)
    if max_replay>2e-6:
        raise RuntimeError("R replay mismatch %g"%max_replay)

    transformed=[]; transformed_names=[]; states={}; stage_counts={}; group_counts={}
    for group in GROUPS:
        x,state=preprocess(raw_groups[group],group_names[group],train_mask)
        states[group]=state
        transformed.append(x)
        transformed_names.extend(["%s__%s"%(group,n) for n in state["output_names"]])
        stage_counts[group]=sum(v.shape[1] for v in transformed)
        group_counts[group]={"raw":int(raw_groups[group].shape[1]),"retained":len(state["input_names"]),"transformed":int(x.shape[1])}
        log("%-10s raw=%d retained=%d transformed=%d cumulative=%d"%(group,raw_groups[group].shape[1],len(state["input_names"]),x.shape[1],stage_counts[group]))
    X=np.column_stack(transformed).astype(np.float32)
    del raw_groups, transformed
    target=np.log(data["R"]/data["R_true"]).astype(np.float32)
    event_index=np.arange(len(target),dtype=np.int32)
    role_code=np.where(roles=="train",0,np.where(roles=="valid",1,2)).astype(np.int8)
    payload={
        "event_index":event_index,"role":role_code,"label":data["label"].astype(np.int8),
        "bin":data["bin"].astype(np.int8),"source":data["source"].astype(np.int32),
        "R":data["R"].astype(np.float32),"R_true":data["R_true"].astype(np.float32),
        "u":data["u"].astype(np.float32),"u_true":data["u_true"].astype(np.float32),
        "target":target,
    }
    for i in range(X.shape[1]):
        payload["f%04d"%i]=X[:,i]
    with uproot.recreate(str(HERE/"tracker_features.root"),compression=uproot.ZLIB(4)) as f:
        f["events"]=payload
    write_tsv(HERE/"feature_names.tsv",[{"index":i,"branch":"f%04d"%i,"feature":name} for i,name in enumerate(transformed_names)])
    meta={
        "groups":GROUPS,"stage_counts":stage_counts,"group_counts":group_counts,
        "preprocessing":states,"events":len(target),
        "train_events":int(np.sum(role_code==0)),"valid_events":int(np.sum(role_code==1)),"test_events":int(np.sum(role_code==2)),
        "target":"log(R_reco/R_true)","max_R_replay":max_replay,
        "input_exclusions":["MC truth as feature","isotope label as feature","run/event/source as feature","production version"],
        "elapsed_seconds":time.time()-started,
    }
    (HERE/"feature_meta.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    log("WROTE %s with %d features in %.1fs"%(HERE/"tracker_features.root",X.shape[1],time.time()-started))

if __name__=="__main__":
    main()
