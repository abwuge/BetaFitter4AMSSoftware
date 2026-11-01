"""Flight-signal-invariant readout contrasts and shared timing corrections.

The estimator keeps the event time slope unrestricted. A smaller residual
contrast is a diagnostic, not proof that the reconstructed slope is better.
"""
import numpy as np


def projector(x):
    x=np.asarray(x,float)
    if x.ndim!=2 or x.shape[1]!=4 or not np.isfinite(x).all():
        raise ValueError('Four finite path/c coordinates are required')
    z=x-x.mean(axis=1)[:,None];den=np.sum(z*z,axis=1)
    if np.any(den<=0):raise ValueError('Zero path span')
    w=z/den[:,None]
    p=np.eye(4)[None,:,:]-0.25-z[:,:,None]*w[:,None,:]
    return p,w


def correction(h,channel,theta):
    h=np.asarray(h,float);channel=np.asarray(channel,int)
    nf=h.shape[2]
    index=channel[:,:,None]*nf+np.arange(nf)[None,None,:]
    return np.sum(h*theta[index],axis=2)


def fit_correction(x,t,h,channel,event_weight=None,robust_iterations=8):
    """Fit sparse per-channel coefficients with both local parameters removed."""
    x,t,h=[np.asarray(v,float) for v in (x,t,h)]
    channel=np.asarray(channel,int)
    if t.shape!=x.shape or h.shape[:2]!=x.shape or channel.shape!=x.shape:
        raise ValueError('Array shape mismatch')
    if not (np.isfinite(t).all() and np.isfinite(h).all()):raise ValueError('Nonfinite input')
    if np.any((channel<0)|(channel>=40)):raise ValueError('Channel outside 4x10 layout')
    p,w=projector(x);pt=np.einsum('nij,nj->ni',p,t)
    nf=h.shape[2];k=40*nf
    index=channel[:,:,None]*nf+np.arange(nf)[None,None,:]
    base=np.ones(len(x)) if event_weight is None else np.asarray(event_weight,float)
    if np.any(base<0) or base.sum()<=0:raise ValueError('Invalid event weights')
    weight=base.copy();theta=np.zeros(k)
    for iteration in range(robust_iterations+1):
        fisher=np.zeros((k,k));rhs=np.zeros(k)
        for l in range(4):
            for j in range(nf):
                rhs+=np.bincount(index[:,l,j],weights=weight*h[:,l,j]*pt[:,l],minlength=k)
                for m in range(4):
                    for q in range(nf):
                        pairs=index[:,l,j]*k+index[:,m,q]
                        fisher+=np.bincount(pairs,weights=weight*p[:,l,m]*h[:,l,j]*h[:,m,q],minlength=k*k).reshape(k,k)
        active=np.diag(fisher)>1.e-12
        scales=np.sqrt(np.diag(fisher)[active]);scaled=fisher[np.ix_(active,active)]/scales[:,None]/scales[None,:]
        eigen=np.linalg.eigvalsh(scaled)
        if len(eigen)==0 or eigen.min()<1.e-8:raise ValueError('Timing correction is not identified')
        revised=np.zeros(k);revised[active]=np.linalg.solve(scaled,rhs[active]/scales)/scales
        residual=pt-np.einsum('nij,nj->ni',p,correction(h,channel,revised))
        radius=np.sqrt(np.sum(residual*residual,axis=1))
        sample=radius[base>0]
        scale=max(1.e-8,np.median(sample)/np.sqrt(2*np.log(2)))
        if iteration==robust_iterations or np.max(np.abs(revised-theta))<1.e-8:
            theta=revised;break
        theta=revised;weight=base*np.minimum(1.,3*scale/np.maximum(radius,1.e-12))
    # Run-level resampling, rather than inverse Hessian, supplies calibration errors.
    return dict(theta=theta,active=active,eigenvalues=eigen,fisher=fisher,scale=scale,
                weighted_residual_sum_squares=float(np.sum(weight[:,None]*residual**2)),iterations=iteration+1)


def reconstruct(x,t,h,channel,theta):
    p,w=projector(x);c=correction(h,channel,theta);new=t-c
    return dict(u=np.sum(w*new,axis=1),correction=c,
                residual=np.einsum('nij,nj->ni',p,new),slope_weights=w)


def robust_fit(x,y,weight=None,iterations=12):
    weight=np.ones(len(y)) if weight is None else np.asarray(weight,float)
    use=weight>0
    x=x[use];y=y[use];weight=weight[use]
    if len(y)<x.shape[1]+5:raise ValueError('Insufficient response observations')
    b=np.linalg.lstsq(x*np.sqrt(weight)[:,None],y*np.sqrt(weight),rcond=1.e-10)[0]
    for _ in range(iterations):
        r=y-x.dot(b);s=max(1.e-7,1.4826*np.median(np.abs(r-np.median(r))))
        ww=weight*np.minimum(1.,2.5*s/np.maximum(np.abs(r),1.e-12))
        new=np.linalg.lstsq(x*np.sqrt(ww)[:,None],y*np.sqrt(ww),rcond=1.e-10)[0]
        if np.max(np.abs(new-b))<1.e-8:b=new;break
        b=new
    return b


def position_basis(xy,adc,center,axis):
    z=np.sum((xy-center)*axis,axis=1)/100.
    adc=np.maximum(np.asarray(adc,float),1.e-8)
    al=np.sqrt(100./adc[:,0]);ar=np.sqrt(100./adc[:,1])
    return np.column_stack([np.ones(len(z)),z,z*z,z*z*z,al,ar,al*al,ar*ar,z*(al-ar)])


def raw_contrasts(raw):
    """raw shape: event,layer,side,threshold(LT,HT,SHT), in ns."""
    raw=np.asarray(raw,float)
    valid=np.isfinite(raw)&(np.abs(raw)<1000)&(raw!=0)
    values=np.zeros((len(raw),4,3));ok=np.zeros(values.shape,bool)
    values[:,:,0]=raw[:,:,0,0]-raw[:,:,1,0];ok[:,:,0]=np.all(valid[:,:,:,0],axis=2)
    for j in [1,2]:
        values[:,:,j]=np.mean(raw[:,:,:,j]-raw[:,:,:,0],axis=2)
        ok[:,:,j]=np.all(valid[:,:,:,j]&valid[:,:,:,0],axis=2)
    return np.where(ok,values,0.),ok


def fit_response(data,train,weight=None):
    """Fit response from a separate run cohort; return serializable parameters."""
    y,valid=raw_contrasts(data['raw']);models=[]
    weight=np.ones(len(train)) if weight is None else np.asarray(weight,float)
    for layer in range(4):
        for bar in range(10):
            group=(data['bar'][:,layer]==bar)&data['good'][:,layer]
            select=group&train&(weight>0)
            if select.sum()<100:continue
            center=np.mean(data['xy'][select,layer],axis=0)
            _,_,vt=np.linalg.svd(data['xy'][select,layer]-center,full_matrices=False)
            axis=vt[0]
            basis=position_basis(data['xy'][:,layer],data['adc'][:,layer],center,axis)
            for feature in range(3):
                take=select&valid[:,layer,feature]
                supported=[int(e) for e in np.unique(data['epoch'][take]) if (take&(data['epoch']==e)).sum()>=30]
                take&=np.isin(data['epoch'],supported)
                if take.sum()<100 or not supported:continue
                ix=np.flatnonzero(take);epochs=data['epoch'][ix]
                design=np.column_stack([basis[ix]]+[(epochs==e).astype(float) for e in supported[1:]])
                b=robust_fit(design,y[ix,layer,feature],weight[ix])
                pooled=robust_fit(basis[ix],y[ix,layer,feature],weight[ix])
                res=y[ix,layer,feature]-design.dot(b)
                scale=max(1.e-4,1.4826*np.median(np.abs(res-np.median(res))))
                models.append(dict(layer=layer,bar=bar,feature=feature,epochs=supported,center=center.tolist(),axis=axis.tolist(),
                    coef=b.tolist(),pooled=pooled.tolist(),scale=float(scale),n=int(take.sum())))
    return models


def apply_response(data,models):
    y,valid=raw_contrasts(data['raw']);h=np.zeros(y.shape);pooled=np.full(y.shape,np.nan);raw_residual=pooled.copy()
    supported=np.zeros(y.shape,bool)
    for m in models:
        layer=m['layer'];bar=m['bar'];feature=m['feature'];epochs=m['epochs']
        take=(data['bar'][:,layer]==bar)&data['good'][:,layer]&valid[:,layer,feature]&np.isin(data['epoch'],epochs)
        ix=np.flatnonzero(take)
        if not len(ix):continue
        basis=position_basis(data['xy'][ix,layer],data['adc'][ix,layer],np.asarray(m['center']),np.asarray(m['axis']))
        design=np.column_stack([basis]+[(data['epoch'][ix]==e).astype(float) for e in epochs[1:]])
        r=y[ix,layer,feature]-design.dot(m['coef'])
        # A bounded input contrast preserves the coherent-flight response.
        h[ix,layer,feature]=np.clip(r,-5*m['scale'],5*m['scale'])
        raw_residual[ix,layer,feature]=r
        pooled[ix,layer,feature]=y[ix,layer,feature]-basis.dot(m['pooled'])
        supported[ix,layer,feature]=True
    return h,supported,raw_residual,pooled
