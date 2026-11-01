#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import run_closure as rc
import run_energy_residualized as er

base=Path(__file__).resolve().parent
train,valid=er.load_cache(base/'controlled_cache.npz')
first=np.load(str(base/'predictions.npz'))
state=rc.make_feature_state(train,'no_edep')
xb_tr=rc.apply_feature_state(train,state)
xb_va=rc.apply_feature_state(valid,state)
response=er.response_state(train)
xtr,xva=er.append_energy(xb_tr,xb_va,er.energy_raw(train,'fitted',response),er.energy_raw(valid,'fitted',response))
control=json.load(open(str(base/'energy_best_control_results.json')))
theta=np.asarray(control['models']['nominal']['theta'],dtype=float)
arrays={
 'raw':valid['u'],
 'no_edep':first['u_weak_no_edep'],
 'fitted':valid['u']-xva.dot(theta),
}
truth=valid['truth_valid']
y=valid['label']
b=valid['bin']
u_true=valid['u_true']

def width(values):
 q=np.percentile(values,[16.,84.]); return float((q[1]-q[0])/2.)

def one(indices, group):
 if group=='C12': indices=indices[y[indices]==0]
 elif group=='C13': indices=indices[y[indices]==1]
 return {k:width(v[indices]-u_true[indices]) for k,v in arrays.items()}

rng=np.random.RandomState(1301315)
strata=[np.flatnonzero(truth&(y==label)&(b==ibin)) for label in (0,1) for ibin in range(rc.N_BINS)]
replicates=[]
for i in range(5000):
 idx=np.concatenate([s[rng.randint(0,len(s),len(s))] for s in strata])
 item={}
 for group in ('all','C12','C13'):
  w=one(idx,group)
  item[group]={
   'raw':w['raw'],'no_edep':w['no_edep'],'fitted':w['fitted'],
   'no_edep_vs_raw_percent':100.*(1.-w['no_edep']/w['raw']),
   'fitted_vs_raw_percent':100.*(1.-w['fitted']/w['raw']),
   'fitted_vs_no_edep_percent':100.*(1.-w['fitted']/w['no_edep']),
   'incremental_percentage_points':100.*(w['no_edep']-w['fitted'])/w['raw'],
  }
 replicates.append(item)

summary={}
for group in ('all','C12','C13'):
 summary[group]={}
 for metric in replicates[0][group]:
  values=np.asarray([r[group][metric] for r in replicates])
  summary[group][metric]={
   'median':float(np.median(values)),
   'p2p5':float(np.percentile(values,2.5)),
   'p97p5':float(np.percentile(values,97.5)),
   'probability_positive':float(np.mean(values>0.)) if ('percent' in metric or 'points' in metric) else None,
  }
observed={group:one(np.flatnonzero(truth),group) for group in ('all','C12','C13')}
out={'observed_widths':observed,'paired_stratified_event_bootstrap_5000':summary,'note':'Fixed trained models; resamples validation events within isotope x rigidity-bin strata. Does not include training or production-version uncertainty.'}
(base/'bootstrap_fixed_validation.json').write_text(json.dumps(out,indent=2,sort_keys=True))
print(json.dumps(out,indent=2,sort_keys=True))
