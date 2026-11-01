"""Auditable minimal histogram boosting wrapper, with replay and null tests."""
import ctypes as c,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def quantize(X,train,edges=None):
 X=np.asarray(X,float)
 if edges is None:
  edges=np.empty((X.shape[1],255));med=[]
  for j in range(X.shape[1]):
   v=X[train,j];v=v[np.isfinite(v)];v=v if len(v) else np.array([0.])
   edges[j]=np.percentile(v,np.linspace(0,100,257)[1:-1]);med.append(float(np.median(v)))
  state=dict(edges=edges,median=np.array(med))
 else:state=edges
 out=np.empty(X.shape,np.uint8,order='F')
 for j in range(X.shape[1]):out[:,j]=np.searchsorted(state['edges'][j],np.where(np.isfinite(X[:,j]),X[:,j],state['median'][j]),side='right')
 return out,state

def fit(X,y,train,snapshots=(160,320,640),leaves=15,depth=6,minleaf=160,l2=20.,rate=.06,weight=None):
 Xq,state=quantize(X,train);y=np.ascontiguousarray(y,dtype=np.float64);n,d=Xq.shape
 ix=np.ascontiguousarray(np.flatnonzero(train),dtype=np.int32)
 w=np.ascontiguousarray(np.ones(n) if weight is None else weight,dtype=np.float64)
 snaps=np.array(snapshots,dtype=np.int32);trees=int(max(snaps));nodes=2*leaves-1
 out=np.full((len(snaps),n),np.nan);fi=np.full((trees,nodes),-1,np.int32);cu=fi.copy();le=fi.copy();ri=fi.copy();va=np.zeros((trees,nodes));init=np.zeros(1)
 lib=c.CDLL(str(HERE/'libhistboost.so'));fun=lib.fit_hgb
 fun.argtypes=[c.c_void_p,c.c_int,c.c_int,c.c_void_p,c.c_void_p,c.c_void_p,c.c_int,c.c_int,c.c_int,c.c_int,c.c_int,c.c_double,c.c_double,c.c_void_p,c.c_int,c.c_void_p,c.c_void_p,c.c_void_p,c.c_void_p,c.c_void_p,c.c_void_p,c.c_void_p]
 fun.restype=c.c_int
 ptr=lambda a:a.ctypes.data_as(c.c_void_p)
 rc=fun(ptr(Xq),n,d,ptr(y),ptr(w),ptr(ix),len(ix),trees,leaves,depth,minleaf,l2,rate,ptr(snaps),len(snaps),ptr(out),ptr(fi),ptr(cu),ptr(le),ptr(ri),ptr(va),ptr(init))
 if rc or not np.isfinite(out).all():raise RuntimeError('Histogram booster failed: %d'%rc)
 model=dict(feature=fi,cut=cu,left=le,right=ri,value=va,initial=init,rate=np.array(rate),edges=state['edges'],median=state['median'],snapshots=snaps)
 return out,model

def predict(X,model,trees=None):
 Xq,_=quantize(X,np.ones(len(X),bool),edges=model);out=np.full(len(Xq),float(model['initial'][0]));rows=np.arange(len(Xq))
 ntrees=len(model['feature']) if trees is None else trees
 for t in range(ntrees):
  node=np.zeros(len(Xq),np.int32)
  while True:
   fs=model['feature'][t,node];active=fs>=0
   if not active.any():break
   r=rows[active];nn=node[active];f=fs[active]
   node[active]=np.where(Xq[r,f]<=model['cut'][t,nn],model['left'][t,nn],model['right'][t,nn])
  out+=float(model['rate'])*model['value'][t,node]
 return out

def self_test():
 rng=np.random.RandomState(141026);X=rng.normal(size=(6000,4));train=np.arange(len(X))<4000
 signal=.7*np.sin(2*X[:,0])+.3*(X[:,1]>0);y=signal+rng.normal(0,.1,len(X))
 p,m=fit(X,y,train,snapshots=(80,160),leaves=12,minleaf=30,l2=1.)
 mse=float(np.mean((p[-1,~train]-y[~train])**2));replay=float(np.max(np.abs(p[-1]-predict(X,m))))
 const,mc=fit(np.zeros((1000,2)),np.arange(1000.)/1000,np.arange(1000)<600,snapshots=(5,),minleaf=30)
 constant_error=float(np.max(np.abs(const-np.mean(np.arange(600.)/1000))))
 yn=rng.normal(size=len(X));pn,mn=fit(X,yn,train,snapshots=(160,),leaves=8,minleaf=100)
 null_r2=float(1-np.mean((yn[~train]-pn[-1,~train])**2)/np.var(yn[~train]))
 assert .007<mse<.04,(mse,'bad signal fit')
 assert replay<1e-12 and constant_error<1e-12,(replay,constant_error)
 assert null_r2<.03,null_r2
 result=dict(passed=True,known_noise_variance=.01,heldout_mse=mse,replay_error=replay,constant_error=constant_error,null_heldout_r2=null_r2)
 (HERE/'learner_self_test.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':self_test()
