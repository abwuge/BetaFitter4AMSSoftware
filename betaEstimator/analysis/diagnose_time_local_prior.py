import sys,json,numpy as np
from pathlib import Path
sys.path.insert(0,'/nas01/home/hxwu/AMSSoft/betaEstimator/analysis')
import time_local_beta as t
out=t.OUT
plan=json.loads((out/'plan.json').read_text()); p=t.load(out/'predictions.npz'); c=t.load(t.OLD/'diagnostic_covariates.npz',['continuous','bars'])
rng=np.random.RandomState(20261011); shuffled=np.full(len(p['run']),np.nan); slope=np.full(len(p['run']),np.nan)
for f in plan['folds']:
 rows=np.flatnonzero(p['fold']==f['fold'])
 for bi,edges in enumerate(t.w.BLOCKS):
  ix=rows[(p['R'][rows]>=edges[0])&(p['R'][rows]<edges[-1])]
  if not len(ix):continue
  t.w.set_edges(edges); dest=out/'folds'/('%04d'%f['fold'])
  state=t.w.rt.arrays(json.loads((dest/('block_%d_model.json'%bi)).read_text()))
  d=dict(u=p['raw_u'][ix],R=p['R'][ix],continuous=c['continuous'][ix].astype(float),bars=c['bars'][ix])
  d['bin']=np.searchsorted(edges,d['R'],side='right')-1
  coef=state['parameters']; x=t.tr.feature_matrix(d,state); correction=x.dot(coef[:x.shape[1]])
  assert np.allclose(d['u']-correction,p['local_corrected_u'][ix],atol=1e-12,rtol=0)
  shifted=dict(d,u=d['u']+1e-4); xs=t.tr.feature_matrix(shifted,state)
  slope[ix]=1-(xs.dot(coef[:xs.shape[1]])-correction)/1e-4
  mixed=d['continuous'].copy()
  for b in range(4):
   loc=np.flatnonzero(d['bin']==b); mixed[loc,18:29]=d['continuous'][rng.permutation(loc),18:29]
  ds=dict(d,continuous=mixed); xs=t.tr.feature_matrix(ds,state)
  shuffled[ix]=1/(d['u']-xs.dot(coef[:xs.shape[1]])+p['delay_u_C12'][ix])
assert np.isfinite(shuffled).all() and np.isfinite(slope).all()
def stats(x):
 q=np.percentile(x,[16,84]);return dict(w68=float((q[1]-q[0])/2),mean=float(x.mean()),rmse=float(np.sqrt(np.mean(x*x))))
records=[]
for lo,hi in [(2.5,6),(6,20),(20,200)]:
 m=p['reference_mask']&(p['R']>=lo)&(p['R']<hi);ref=p['rich_beta'][m]
 data={name:stats(beta[m]-ref) for name,beta in [('raw',p['raw_beta_rich']),('local',p['local_beta_rich']),('energy_shuffled',shuffled)]}
 row=dict(low=lo,high=hi,n=int(m.sum()),metrics=data,energy_shuffled_gain=100*(1-data['energy_shuffled']['w68']/data['raw']['w68']),direct_u_path_slope_quantiles=np.percentile(slope[m],[5,25,50,75,95]).tolist())
 records.append(row);print(json.dumps(row),flush=True)
t.w.ev.dump(out/'prior_shortcut_diagnostic.json',dict(posthoc_diagnostic=True,models_retrained=False,predictions_sha256=t.w.ev.digest(out/'predictions.npz'),energy_shuffle='joint TOF/Tracker energy columns 18:29 within target fold and training rigidity bin; preserve geometry, TOF timing residuals and raw u',slope='finite difference corrected_u versus raw_u through explicit feature u path only; fixed other measured features, epsilon=1e-4; not total physical response',rows=records))
