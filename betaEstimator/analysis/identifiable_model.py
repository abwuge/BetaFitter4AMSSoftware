"""Keep only the slope correction identifiable under independent layer noise."""
import numpy as np
from timing_model import clean,kinematics,contrasts,predict as timing_predict
def predict(d,model):
    x,t,u0=kinematics(d)
    q=contrasts(x)
    _,mu,v=timing_predict(d,model)
    a0=x/np.sum(x*x,axis=1)[:,None]
    n=np.ones_like(x)
    for i in range(4):
        for j in range(4):
            if i!=j:n[:,i]*=x[:,i]-x[:,j]
    nn=np.linalg.norm(n,axis=1)
    safe=nn>1e-14
    n=n/np.maximum(nn,1e-14)[:,None]
    gn=np.einsum('nri,ni->nr',q,n*a0)
    g=np.einsum('nri,ni->nr',q,v*a0)
    cov=np.einsum('nri,ni,nsi->nrs',q,v,q)
    h=np.stack([-gn[:,1],gn[:,0]],axis=1)
    hnorm=np.linalg.norm(h,axis=1)
    safe&=hnorm>1e-12
    h=h/np.maximum(hnorm,1e-12)[:,None]
    den=np.einsum('nr,nrs,ns->n',h,cov,h)
    coef=np.sum(h*g,axis=1)/np.maximum(den,1e-20)
    r=np.einsum('nri,ni->nr',q,t-mu)
    correction=coef*np.sum(h*r,axis=1)
    correction[~safe]=0
    return u0-correction,correction
def reconstruct(d,model):
    valid=clean(d)
    with np.errstate(invalid='ignore',divide='ignore'):u=1/d['tof_betah'].astype(float)
    idx=np.flatnonzero(valid)
    for start in range(0,len(idx),50000):
        ii=idx[start:start+50000]
        part={k:v[ii] for k,v in d.items()}
        cand,_=predict(part,model)
        ok=np.isfinite(cand)&(cand>0)
        u[ii[ok]]=cand[ok]
        valid[ii[~ok]]=False
    return u,valid
