#!/usr/bin/env python3
from __future__ import print_function
import csv, json
from pathlib import Path
import numpy as np
import uproot

HERE=Path(__file__).resolve().parent
MASS=HERE.parent/"c_mass_toy_oracle"/"predictions.npz"
MASSES=np.array([11.17793,12.11221],dtype=float)
Z=6.0
STAGES=[
 ("raw",None,None),
 ("multifit","tmva_full","multifit"),
 ("otherhit","tmva_ablation","otherhit"),
 ("noother","tmva_ablation","noother"),
 ("fitdiag","tmva_full","fitdiag"),
 ("geometry","tmva_full","geometry"),
 ("morph","tmva_ablation","morph"),
 ("cluster","tmva_full","cluster"),
 ("chargeprobe","tmva_ablation","chargeprobe"),
]

def w68(x):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 q=np.percentile(x,[16,84]); return float((q[1]-q[0])/2)

def rms(x):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 return float(np.sqrt(np.mean(x*x)))

def auc(score,label):
 score=np.asarray(score,float);label=np.asarray(label,np.int8)
 ok=np.isfinite(score);score=score[ok];label=label[ok]
 o=np.argsort(score,kind="mergesort");s=score[o];y=label[o]
 n1=int(np.sum(y==1));n0=int(np.sum(y==0))
 st=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1];en=np.r_[st[1:],len(s)]
 ranks=np.repeat((st+1.0+en)/2.0,en-st)
 return float((np.sum(ranks[y==1])-n1*(n1+1.0)/2.0)/(n0*n1))

def write_tsv(path,rows):
 if not rows:return
 with open(str(path),"w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0].keys()),delimiter="\t",lineterminator="\n")
  w.writeheader();w.writerows(rows)

def load_stage(directory,stage,which="TestTree"):
 chunks=[]
 for b in range(4):
  path=HERE/directory/("%s_bin%d.root"%(stage,b))
  with uproot.open(str(path)) as f:
   key=next(k for k in f.keys(recursive=True) if which in k)
   a=f[key].arrays(["event_index","target","BDTG","label","source","R","R_true","u","u_true"],library="np")
  chunks.append(a)
 out={}
 for key in chunks[0]: out[key]=np.concatenate([x[key] for x in chunks])
 return out

def main():
 with np.load(str(MASS),allow_pickle=False) as z:
  role=np.char.decode(z["role"],"ascii")
  test=role=="test"
  base={k:z[k][test] for k in ["run","event","source","label","bin","R","R_true","u","u_true"]}
  base_idx=np.flatnonzero(test)
 index_to_local={int(v):i for i,v in enumerate(base_idx)}
 n=len(base_idx)
 raw=np.log(base["R"]/base["R_true"])
 predictions={"raw":np.zeros(n)}
 stage_arrays={}
 for name,directory,stage in STAGES[1:]:
  a=load_stage(directory,stage)
  pred=np.full(n,np.nan)
  for ei,v in zip(a["event_index"].astype(np.int64),a["BDTG"]):
   pred[index_to_local[int(ei)]]=float(v)
  if np.any(~np.isfinite(pred)):raise RuntimeError("missing "+name)
  predictions[name]=pred
  stage_arrays[name]=a

 rows=[]
 massrows=[]
 for name,_,_ in STAGES:
  err=raw-predictions[name]
  corrR=base["R"]*np.exp(-predictions[name])
  for b in [-1,0,1,2,3]:
   for sample,sm in [("all",np.ones(n,bool)),("C12",base["label"]==0),("C13",base["label"]==1)]:
    sel=sm & ((base["bin"]==b) if b>=0 else True)
    rows.append({"method":name,"bin":b,"sample":sample,"n":int(np.sum(sel)),
                 "median_dlogR":float(np.median(err[sel])),"w68_dlogR":w68(err[sel]),"rms_dlogR":rms(err[sel]),
                 "correction_label_auc":auc(predictions[name],base["label"]) if b==-1 and sample=="all" and name!="raw" else float("nan")})
  for beta_name,u in [("raw_beta",base["u"]),("truth_beta",base["u_true"])]:
   m2=(Z*corrR)**2*(u*u-1)
   e=m2-MASSES[base["label"]]**2
   massrows.append({"method":name,"beta":beta_name,"n":n,"mass2_auc":auc(m2,base["label"]),"w68_mass2_error":w68(e),"median_mass2_error":float(np.median(e))})
 raw_lookup={(r["bin"],r["sample"]):r["w68_dlogR"] for r in rows if r["method"]=="raw"}
 for r in rows:r["gain_vs_raw_pct"]=100*(1-r["w68_dlogR"]/raw_lookup[(r["bin"],r["sample"])])
 write_tsv(HERE/"ablation_metrics.tsv",rows)
 write_tsv(HERE/"ablation_mass_metrics.tsv",massrows)

 # Train/test gap from TMVA trees.
 gap=[]
 for name,directory,stage in STAGES[1:]:
  for which in ["TrainTree","TestTree"]:
   a=load_stage(directory,stage,which)
   e=a["target"]-a["BDTG"]
   gap.append({"method":name,"sample":which,"n":len(e),"raw_w68":w68(a["target"]),"corrected_w68":w68(e),"gain_pct":100*(1-w68(e)/w68(a["target"])),
               "median_dlogR":float(np.median(e))})
 write_tsv(HERE/"train_test_gap.tsv",gap)

 # Whole-source bootstrap, isotope-stratified and paired across methods.
 groups={};keys={0:[],1:[]}
 for lab in (0,1):
  for src in sorted(int(x) for x in np.unique(base["source"][base["label"]==lab])):
   k=(lab,src);groups[k]=np.flatnonzero((base["label"]==lab)&(base["source"]==src));keys[lab].append(k)
 rng=np.random.RandomState(300930)
 methods=[x[0] for x in STAGES[1:]]
 gains={m:{"overall":[],"bins":[[] for _ in range(4)]} for m in methods}
 reps=1000
 errors={m:raw-predictions[m] for m in methods}
 for _ in range(reps):
  idx=np.concatenate([groups[keys[lab][j]] for lab in (0,1) for j in rng.randint(0,len(keys[lab]),size=len(keys[lab]))])
  rw=w68(raw[idx])
  for m in methods:
   gains[m]["overall"].append(100*(1-w68(errors[m][idx])/rw))
   for b in range(4):
    ib=idx[base["bin"][idx]==b]
    gains[m]["bins"][b].append(100*(1-w68(errors[m][ib])/w68(raw[ib])))
 boot={"replicates":reps,"methods":{}}
 for m in methods:
  a=np.asarray(gains[m]["overall"])
  boot["methods"][m]={"overall":{"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))},"bins":[]}
  for b,x in enumerate(gains[m]["bins"]):
   a=np.asarray(x)
   boot["methods"][m]["bins"].append({"bin":b,"median":float(np.median(a)),"low95":float(np.percentile(a,2.5)),"high95":float(np.percentile(a,97.5)),"positive_fraction":float(np.mean(a>0))})
 (HERE/"ablation_bootstrap.json").write_text(json.dumps(boot,indent=2,sort_keys=True)+"\n")

 overall={r["method"]:r for r in rows if r["bin"]==-1 and r["sample"]=="all"}
 byiso={m:{s:next(r for r in rows if r["method"]==m and r["bin"]==-1 and r["sample"]==s) for s in ["C12","C13"]} for m,_,_ in STAGES}
 mass={r["method"]+"__"+r["beta"]:r for r in massrows}
 result={"overall":overall,"by_isotope":byiso,"mass":mass,"bootstrap":boot,
         "interpretation":{"detector_specific_candidate":"morph","unsafe_charge_shortcut":"chargeprobe",
         "charge_shortcut_mechanism":"Tracker charge/dE/dx plus reconstructed R estimates beta-like kinematics. Without isotope identity, the regression chooses an intermediate mass hypothesis, producing opposite C12/C13 rigidity-scale biases while artificially narrowing each class."}}
 (HERE/"ablation_results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("SUMMARY")
 for m,_,_ in STAGES:
  r=overall[m]
  print("%-12s W68 %.6f gain %6.2f med %+0.5f labelAUC %s"%(m,r["w68_dlogR"],r["gain_vs_raw_pct"],r["median_dlogR"],("%.3f"%r["correction_label_auc"] if np.isfinite(r["correction_label_auc"]) else "-")))
 print("BOOTSTRAP morph",boot["methods"]["morph"]["overall"])
 print("BOOTSTRAP noother",boot["methods"]["noother"]["overall"])

if __name__=="__main__":main()
