import numpy as np, uproot, json
A=np.load('aligned.npz'); Z=np.load('decomposition.npz')
with uproot.open('native.root') as f: P=f['probe'].arrays(['jacrec','idealrec'],library='np')
with uproot.open('intervals.root') as f: K=f['intervals'].arrays(['kick','offset'],library='np')
roles=A['ref_role']; valid=Z['validtruth']; tr=(roles==b'train')&valid; va=(roles==b'valid')&valid; te=(roles==b'test')&valid
R=A['ref_R'].astype(float); err=np.log(R/A['ref_R_true']); src=A['ref_source']
y=A['tk_hitc'][:,1:8,1].astype(float); z=A['tk_hitc'][:,1:8,2].astype(float)
ir=P['idealrec'][:,1:8,1].astype(float); J=P['jacrec'][:,1:8,1,:][:,:,[1,3,4]].astype(float)
Merr=Z['measurement_hit_cm']; lo,hi=np.percentile(Merr[tr],[1,99],axis=0); M=np.cov(np.clip(Merr[tr],lo,hi),rowvar=False)+np.eye(7)*1e-12
off=K['offset'][:,1:7,1].astype(float); kick=K['kick'][:,1:7,1].astype(float)
U=np.concatenate([off*R[:,None],kick*R[:,None]],axis=1); lo,hi=np.percentile(U[tr],[1,99],axis=0); Q=np.cov(np.clip(U[tr],lo,hi),rowvar=False)+np.eye(12)*1e-16
res=y-ir; mu=np.median(res[tr],axis=0); zmed=np.median(z[tr],axis=0)
G=np.zeros((7,12))
for j in range(7):
 for k in range(6):
  if k<j:
   G[j,k]=1.
   if j>k+1: G[j,6+k]=zmed[j]-zmed[k+1]
weak=[1,3,5]
def metric(x):
 x=x[np.isfinite(x)]; q=np.percentile(x,[16,50,84]); return {'n':len(x),'median':float(q[1]),'w68':float((q[2]-q[0])/2),'std':float(np.std(x))}
def delta_for(lam, mask):
 S=np.ones(12)
 for k in weak: S[k]=S[6+k]=np.sqrt(lam)
 Qm=Q*S[:,None]*S[None,:]; Cp=G@Qm@G.T; idx=np.flatnonzero(mask); out=np.full(len(R),np.nan)
 for c0 in range(0,len(idx),1500):
  ii=idx[c0:c0+1500]; C=M[None,:,:]+Cp[None,:,:]/(R[ii,None,None]**2); W=np.linalg.inv(C); Ji=J[ii]
  JW=np.einsum('npi,npq->niq',Ji,W); info=JW@Ji; Aop=(np.linalg.inv(info)@JW)[:,2,:]
  out[ii]=np.einsum('ni,ni->n',Aop,res[ii]-mu)
 return out
lams=[0.25,0.5,1,2,4,8,16]; alphas=[0,.25,.5,.75,1,1.25]; allmask=tr|va|te; rows=[]; ds={}; es={}
for lam in lams:
 d=delta_for(lam,allmask); ds[lam]=d; best=None
 for a in alphas:
  den=1-a*d; e=err-np.log(np.maximum(den,1e-9)); e[den<=0]=err[den<=0]; m=metric(e[va])
  if best is None or m['w68']<best[1]['w68']: best=(a,m)
 a=best[0]; den=1-a*d; e=err-np.log(np.maximum(den,1e-9)); e[den<=0]=err[den<=0]; es[lam]=e
 rows.append({'lambda_weak':lam,'alpha':a,'valid':best[1],'test':metric(e[te])}); print(rows[-1],flush=True)
choice=min(rows,key=lambda x:x['valid']['w68']); lam=choice['lambda_weak']; ec=es[lam]; base=[x for x in rows if x['lambda_weak']==1][0]; eb=es[1]
rng=np.random.RandomState(20261002); st=src[te]; ids=np.unique(st); groups=[np.flatnonzero(st==s) for s in ids]
def boot(x0,x1,n=600):
 vals=[]; x0=x0[te];x1=x1[te]
 for _ in range(n):
  kk=np.concatenate([groups[i] for i in rng.randint(len(groups),size=len(groups))]); vals.append(100*(1-metric(x1[kk])['w68']/metric(x0[kk])['w68']))
 return np.percentile(vals,[2.5,50,97.5]).tolist()
bend=[];sig=[]
with uproot.open('intervals.root') as f: KK=f['intervals'].arrays(['bend','kick'],library='np')
for k in range(1,7):
 b=np.median(np.abs(KK['bend'][te,k,1])); x=KK['kick'][te,k,1]; q=np.percentile(x,[16,84]); bend.append(float(b));sig.append(float((q[1]-q[0])/2))
snr=np.array(bend)/np.array(sig); info=snr**2
each_ci={str(x['lambda_weak']):boot(err,es[x['lambda_weak']]) for x in rows}
out={'rows':rows,'selected':choice,'raw_test':metric(err[te]),'lambda1_test':base['test'],'each_gain_vs_raw_CI_pct':each_ci,'selected_gain_vs_lambda1_CI_pct':boot(eb,ec),'selected_gain_vs_raw_CI_pct':boot(err,ec),'snr':snr.tolist(),'info_fraction':(info/info.sum()).tolist(),'strong_info_fraction':float(info[[0,2,4]].sum()/info.sum()),'zmed':zmed.tolist()}
open('interval_weight_scan.json','w').write(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
