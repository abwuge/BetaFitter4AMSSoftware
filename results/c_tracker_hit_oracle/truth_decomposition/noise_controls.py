import json
from pathlib import Path
import numpy as np,uproot
from analyze import metric,boot
H=Path(__file__).resolve().parent;a=np.load(str(H/'aligned.npz'));z=np.load(str(H/'decomposition.npz'))
with uproot.open(str(H/'native.root')) as f:p=f['probe'].arrays(library='np')
roles=a['ref_role'];tr=(roles==b'train')&z['validtruth'];va=(roles==b'valid')&z['validtruth'];te=(roles==b'test')&z['validtruth'];n=len(tr)
err=np.log(a['ref_R']/a['ref_R_true']);bins=a['ref_bin'];j=p['jacrec'][:,1:8,1,:][:,:,[1,3,4]];ir=p['idealrec'][:,1:8,1]
result=[]
for scale in [1.,.5,0.]:
 trainerr=z['transport_hit_cm']+scale*z['measurement_hit_cm'];points=z['truth'][:,:,1]+scale*z['measurement_hit_cm'];delta=np.zeros(n)
 for b in range(4):
  q=tr&(bins==b);rows=np.flatnonzero(bins==b);xx=trainerr[q];lo,hi=np.percentile(xx,[1,99],axis=0);C=np.cov(np.clip(xx,lo,hi),rowvar=False)+np.eye(7)*1e-12;W=np.linalg.inv(C);mu=np.median(xx,axis=0)
  mat=j[rows];mw=np.einsum('nji,jk->nik',mat,W);A=np.matmul(np.linalg.inv(np.matmul(mw,mat)),mw)[:,2,:]
  delta[rows]=np.einsum('ni,ni->n',A,points[rows]-ir[rows]-mu)
 scan=[];ee=[]
 for alpha in [0,.25,.5,.75,1.]:
  den=1-alpha*delta;e=err-np.log(np.maximum(den,1e-12));e[den<=0]=err[den<=0];scan.append({'alpha':alpha,'valid':metric(e[va])});ee.append(e)
 choice=min(range(5),key=lambda k:scan[k]['valid']['w68']);e=ee[choice]
 result.append({'measurement_noise_scale':scale,'full_GLS_test':metric(ee[-1][te]),'validation_scan':scan,'selected_alpha':scan[choice]['alpha'],'selected_test':metric(e[te]),'gain_CI_percent_vs_raw':boot(err[te],e[te],a['ref_source'][te])})
 if scale==1:ref=ee[-1][te]
 else:result[-1]['gain_CI_full_GLS_vs_scale1']=boot(ref,ee[-1][te],a['ref_source'][te])
(H/'noise_controls.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
