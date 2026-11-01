#!/usr/bin/env python3
"""Controlled hit/transport decomposition and measured-only GLS closure.
The pure-field truth trajectory is diagnostic ONLY. The deployable GLS uses
only reconstructed anchor state, hit positions and train-learned covariance.
"""
import json,csv
from pathlib import Path
import numpy as np
import uproot
H=Path(__file__).resolve().parent

def metric(x):
 x=np.asarray(x); x=x[np.isfinite(x)]
 if not len(x):return {'n':0}
 q=np.percentile(x,[2.5,16,50,84,97.5])
 return dict(n=len(x),median=float(q[2]),w68=float((q[3]-q[1])/2),w95=float((q[4]-q[0])/2),std=float(np.std(x)))
def boot(a,b,src,n=600):
 rng=np.random.RandomState(20261002);groups=[np.flatnonzero(src==s) for s in np.unique(src)];v=[]
 for _ in range(n):
  k=np.concatenate([groups[i] for i in rng.randint(len(groups),size=len(groups))]);v.append(100*(1-metric(b[k])['w68']/metric(a[k])['w68']))
 return [float(v) for v in np.percentile(v,[2.5,50,97.5])]
def main():
 a=np.load(str(H/'aligned.npz'))
 with uproot.open(str(H/'native.root')) as f:p=f['probe'].arrays(library='np')
 assert np.array_equal(p['event_index'],a['ref_event_index'])
 n=len(a['ref_R']);ids=np.array([7,8,9,10,11,12,13]);h=a['tk_hitc'][:,1:8,:].astype(float)
 c=a['mevcoo1'][:,ids,:].astype(float);d=a['mevdir1'][:,ids,:].astype(float)
 truth=c+d/d[:,:,2,None]*(h[:,:,2]-c[:,:,2])[:,:,None]
 validtruth=(a['mevmom1'][:,ids]>0).all(axis=1)&(a['mtrpri'][:,1:8]==1).all(axis=1)
 validtruth&=((a['tk_hitb'][:,1]&254)==254)
 roles=a['ref_role'];tr=(roles==b'train')&validtruth;va=(roles==b'valid')&validtruth;te=(roles==b'test')&validtruth
 err=np.log(a['ref_R']/a['ref_R_true']); bins=a['ref_bin'];
 ide=p['ideal'][:,1:8,1]; iderec=p['idealrec'][:,1:8,1]
 measurement=h[:,:,1]-truth[:,:,1];transport=truth[:,:,1]-ide;total=h[:,:,1]-ide
 # y0, dy/dz, logR; x state is held fixed, not fitted from y alone.
 J=p['jac'][:,1:8,1,:][:,:,[1,3,4]]
 Jr=p['jacrec'][:,1:8,1,:][:,:,[1,3,4]]
 terms={k:np.full(n,np.nan) for k in ['transport','measurement','total','gls_delta','OLS_transport','OLS_measurement']}
 weights=np.zeros((n,7));covs=[];means=[]
 for b in range(4):
  train=tr&(bins==b);rows=np.flatnonzero(bins==b)
  X=total[train]; lo,hi=np.percentile(X,[1,99],axis=0);xc=np.clip(X,lo,hi)
  C=np.cov(xc,rowvar=False);C+=np.eye(7)*2.5e-9
  mu=np.median(X,axis=0);W=np.linalg.inv(C)
  covs.append(C.tolist());means.append(mu.tolist())
  for jj,kind in [(J,'truth'),(Jr,'reco')]:
   mat=jj[rows];mtw=np.einsum('nji,jk->nik',mat,W)
   A=np.matmul(np.linalg.inv(np.matmul(mtw,mat)),mtw)[:,2,:]
   if kind=='truth':
    for key,x in [('transport',transport),('measurement',measurement),('total',total)]:terms[key][rows]=np.einsum('ni,ni->n',A,x[rows])
   else:
    weights[rows]=A;terms['gls_delta'][rows]=np.einsum('ni,ni->n',A,h[rows,:,1]-iderec[rows]-mu)
  mat=J[rows];mt=mat.transpose(0,2,1);A=np.matmul(np.linalg.inv(np.matmul(mt,mat)),mt)[:,2,:]
  for key,x in [('OLS_transport',transport),('OLS_measurement',measurement)]:terms[key][rows]=np.einsum('ni,ni->n',A,x[rows])
 # Select just one interpolation strength on validation; zero is included.
 scan=[]
 for alpha in [0,.25,.5,.75,1.]:
  den=1-alpha*terms['gls_delta'];r=np.where(den>0,err-np.log(np.maximum(den,1e-12)),np.nan)
  scan.append({'alpha':alpha,'valid':metric(r[va])})
 alpha=min(scan,key=lambda x:x['valid']['w68'])['alpha']
 den=1-alpha*terms['gls_delta'];ecorr=np.where(den>0,err-np.log(np.maximum(den,1e-12)),err)
 measurements=[]
 for l in range(7):
  good=(roles==b'test')&(a['mtrpri'][:,l+1]==1)&(a['mevmom1'][:,ids[l]]>0)
  measurements.append({'layer':l+2,'error_y_um':metric(measurement[good,l]*1e4),'transport_from_L5_y_um':metric(transport[good,l]*1e4),'naive_y_difference_um':metric((h[good,l,1]-c[good,l,1])*1e4),'energy_log_ratio_L5':metric(np.log(a['mevmom1'][good,ids[l]]/(a['ref_R_true'][good]*6)))})
 result={'definition':'W68 of ln(R/Rtruth at L5); native controlled fits use constant errors, NOT exact baseline replay','counts':{s:int(np.sum((roles==s.encode())&validtruth)) for s in ['train','valid','test']},'original_counts':{s:int(np.sum(roles==s.encode())) for s in ['train','valid','test']},'per_layer':measurements,'raw_test_full':metric(err[roles==b'test']),'raw_test_common':metric(err[te]),'reference_R_truth_max_abs_diff':float(np.max(abs(a['mevmom1'][:,10]/6-a['ref_R_true']))),'gls':{'validation_scan':scan,'selected_alpha':alpha,'test':metric(ecorr[te]),'test_gain_bootstrap_pct':boot(err[te],ecorr[te],a['ref_source'][te]),'covariances_cm2':covs,'means_cm':means},'native_fits':{},'projection':{k:metric(v[te]) for k,v in terms.items() if k!='gls_delta'}}
 fits=p['fits'];support=te&(fits>0).all(axis=1)&np.isfinite(fits).all(axis=1)
 for j,name in enumerate(['measured_noMS','measured_MS','truth_hits_noMS','truth_hits_MS']):
  result['native_fits'][name]=metric(np.log(fits[support,j]/a['ref_R_true'][support]))
 for name,first,last in [('ideal_hits_noMS',0,2),('ideal_hits_MS',1,3)]:
  e1=np.log(fits[support,first]/a['ref_R_true'][support]);e2=np.log(fits[support,last]/a['ref_R_true'][support]);result['native_fits'][name+'_gain_CI_pct']=boot(e1,e2,a['ref_source'][support]);result['native_fits'][name+'_difference']=metric(np.log(fits[support,first]/fits[support,last]))
 result['native_fits']['raw_same_support']=metric(err[support])
 # Variance accounting is additive with covariance. W68 values are not variances.
 q=te&np.isfinite(terms['total']);x=terms['measurement'][q];y=terms['transport'][q]
 result['projection']['variance_accounting']={'measurement_var':float(np.var(x)),'transport_var':float(np.var(y)),'twice_cov':float(2*np.mean((x-x.mean())*(y-y.mean()))),'total_var':float(np.var(x+y)),'correlation_transport_rawR':float(np.corrcoef(y,err[q])[0,1]),'correlation_measurement_rawR':float(np.corrcoef(x,err[q])[0,1]),'closure_max':float(np.max(abs(terms['total'][q]-x-y)))}
 # Transport-cancelled and hit-error-cancelled counterfactuals using native linear response.
 for key in ['measurement','transport']:
  er=err-terms[key];result['projection']['raw_minus_'+key]=metric(er[te])
 result['by_bin']=[]
 for b in range(4):
  q=te&(bins==b);s=support&(bins==b)
  result['by_bin'].append({'bin':b,'raw':metric(err[q]),'gls':metric(ecorr[q]),'measurement_projection':metric(terms['measurement'][q]),'transport_projection':metric(terms['transport'][q]),'native_meas_MS':metric(np.log(fits[s,1]/a['ref_R_true'][s])),'native_truth_MS':metric(np.log(fits[s,3]/a['ref_R_true'][s]))})
 (H/'results.json').write_text(json.dumps(result,indent=2))
 np.savez_compressed(str(H/'decomposition.npz'),measurement_hit_cm=measurement,transport_hit_cm=transport,truth=truth,validtruth=validtruth,weights=weights,ecorr=ecorr,**terms)
 print(json.dumps({k:v for k,v in result.items() if k not in ['gls','per_layer']},indent=2),flush=True)
 print('GLS',json.dumps({k:v for k,v in result['gls'].items() if k not in ['covariances_cm2','means_cm']}),flush=True)
if __name__=='__main__':main()
