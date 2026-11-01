#!/usr/bin/env python3
import numpy as np,uproot,json
from pathlib import Path
H=Path(__file__).resolve().parent
stages=['rev_noeloss','rev_eloss','base_fair','base_eloss']
def load(stage):
 out=[]
 for b in range(4):
  with uproot.open(str(H/'tmva_reverse'/('%s_bin%d.root'%(stage,b)))) as f:
   k=next(k for k in f.keys(recursive=True) if 'TestTree' in k)
   out.append(f[k].arrays(['event_index','role','source','R','R_true','target','BDTG'],library='np'))
 return {k:np.concatenate([x[k] for x in out]) for k in out[0]}
D={s:load(s) for s in stages}
# Align all stages by event_index.
ids=np.sort(D[stages[0]]['event_index'].astype(np.int64))
P={}
for s in stages:
 d=D[s]; order=np.argsort(d['event_index']); assert np.array_equal(d['event_index'][order],ids)
 if s==stages[0]:
  role=d['role'][order].astype(np.int8); source=d['source'][order].astype(np.int32); raw=d['target'][order].astype(float)
 P[s]=np.asarray(d['BDTG'][order]).reshape(-1).astype(float)
valid=role==1; test=role==2
def metric(x,mask):
 v=x[mask];q=np.percentile(v,[2.5,16,50,84,97.5]);return {'n':int(mask.sum()),'median':float(q[2]),'w68':float((q[3]-q[1])/2),'w95':float((q[4]-q[0])/2),'std':float(np.std(v))}
alphas=np.arange(0,1.51,0.1)
sel={}
errors={}
for s in stages:
 scan=[]
 for a in alphas:
  e=raw-a*P[s]; scan.append((a,metric(e,valid)['w68']))
 a=min(scan,key=lambda x:x[1])[0]; sel[s]=float(a); errors[s]=raw-a*P[s]
 print(s,'alpha',a,'valid',metric(errors[s],valid),'test',metric(errors[s],test))
# source bootstrap on test
src=source[test]; keys=np.unique(src); groups=[np.flatnonzero(src==k) for k in keys]; rng=np.random.RandomState(20261002)
def w68(v):
 q=np.percentile(v,[16,84]);return (q[1]-q[0])/2
boot={s:[] for s in stages}; inc=[]
rt=raw[test]
for _ in range(1000):
 ix=np.concatenate([groups[i] for i in rng.randint(len(groups),size=len(groups))])
 rw=w68(rt[ix])
 for s in stages: boot[s].append(100*(1-w68(errors[s][test][ix])/rw))
 inc.append(100*(1-w68(errors['base_eloss'][test][ix])/w68(errors['base_fair'][test][ix])))
def ci(v): return [float(x) for x in np.percentile(v,[2.5,50,97.5])]
# Truth transport diagnostic (not model input)
A=np.load(str(H/'aligned.npz')); Z=np.load(str(H/'decomposition.npz')); lookup={int(v):i for i,v in enumerate(A['ref_event_index'])}; ix=np.array([lookup[int(v)] for v in ids])
diag={}
for s in stages:
 q=test & Z['validtruth'][ix]
 diag[s]={
  'corr_pred_raw_error':float(np.corrcoef(P[s][q],raw[q])[0,1]),
  'corr_pred_truth_transport_projection':float(np.corrcoef(P[s][q],Z['transport'][ix[q]])[0,1]),
  'corr_corrected_truth_transport_projection':float(np.corrcoef(errors[s][q],Z['transport'][ix[q]])[0,1])
 }
# Measured-visible scatter subsets: use preprocessed energy-loss kink_absmax feature's raw precursor reconstructed from reverse refs.
with uproot.open(str(H/'reverse.root')) as f:E=f['reverse'].arrays(['event_index','yref'],library='np')
assert np.array_equal(E['event_index'],A['ref_event_index'])
y=A['tk_hitc'][:,1:8,1].astype(float); z=A['tk_hitc'][:,1:8,2].astype(float); res=y-E['yref']; slopes=np.diff(res,axis=1)/np.diff(z,axis=1); kval=np.max(np.abs(np.diff(slopes,axis=1)),axis=1)
kv=kval[ix]
cut=np.percentile(kv[valid],[50,75,90])
sub={}
for label,m in [('all',test),('top50',test&(kv>=cut[0])),('top25',test&(kv>=cut[1])),('top10',test&(kv>=cut[2]))]:
 sub[label]={'raw':metric(raw,m)}
 for s in stages:sub[label][s]=metric(errors[s],m)
# bin metrics
bins=np.empty_like(role)
# recover from source file root
with uproot.open(str(H/'reverse_features.root')) as f:
 T=f['events'].arrays(['event_index','bin'],library='np')
 bm={int(e):int(b) for e,b in zip(T['event_index'],T['bin'])}
bins=np.array([bm[int(e)] for e in ids],dtype=np.int8)
bybin=[]
for b in range(4):
 m=test&(bins==b); row={'bin':b,'raw':metric(raw,m)}
 for s in stages:row[s]=metric(errors[s],m)
 bybin.append(row)
out={'alpha_selected_on_validation':sel,'raw_test':metric(raw,test),'test':{s:metric(errors[s],test) for s in stages},
     'gain_vs_raw_bootstrap_CI_pct':{s:ci(boot[s]) for s in stages},
     'base_eloss_incremental_vs_base_fair_CI_pct':ci(inc),
     'diagnostics_truth_only':diag,'visible_kink_subsets':sub,'by_bin':bybin,
     'definition':'BDTG trained on original train sources only; one global correction scale selected on original validation; final metrics on source-file-disjoint test. Reverse features use reconstructed R-conditioned native-field trajectory; eloss variant includes train-derived mean dlogp profile. MC truth enters only targets/evaluation.'}
(H/'reverse_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
