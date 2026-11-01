"""Selected-sample joint timing / direct end-ADC working response model.

The continuous event estimator does not use an isotope, rigidity or source
posterior. Combining ADC with time changes the response to time-only injection;
that response is reported explicitly, not mislabeled as pure-TOF equivariance.
"""
import ctypes as ct
from pathlib import Path
import numpy as np
from scipy.special import gammaln,expit
import flight_conditional_likelihood as base

ROOT=base.ROOT
OUT=ROOT/'results/flight_calibration_research/continuation_20261013/joint_adc'
UREF=1.15
TRI=[(i,j) for i in range(7) for j in range(i+1)]
GEOMETRY_NAMES=['L5_x_over50','L5_y_over50','L5_x2_over2500','L5_y2_over2500','L5_xy_over2500']


def response_features(data,train,normalization=None):
    p=np.asarray(data['tk_pos'],float)[:,4,:2]/50.;xx=p[:,0];yy=p[:,1]
    raw=np.column_stack([xx,yy,xx*xx,yy*yy,xx*yy])
    if normalization is None:
        center=raw[train].mean(0);scale=raw[train].std(0)
        if np.any(scale<1e-9):raise ValueError('Degenerate response geometry')
        normalization=dict(names=GEOMETRY_NAMES,center=center.tolist(),scale=scale.tolist(),clip=3.)
    if normalization['names']!=GEOMETRY_NAMES:raise ValueError('Geometry feature order differs')
    geo=np.clip((raw-np.array(normalization['center']))/np.array(normalization['scale']),-3.,3.)
    bar=np.asarray(data['bar'],int)
    if bar.shape!=(len(raw),4) or np.any((bar<0)|(bar>=10)):raise ValueError('Invalid four-layer bar association')
    features=np.concatenate([np.ones((len(raw),4,1)),np.broadcast_to(geo[:,None,:],(len(raw),4,5)),
        (bar[:,:,None]==np.arange(10)[None,None,:]).astype(float)],axis=2)
    if not np.isfinite(features).all():raise ValueError('Nonfinite response features')
    return features,normalization


def precision(theta):
    theta=np.asarray(theta,float)
    if len(theta)<39 or not np.isfinite(theta).all():raise ValueError('Invalid joint response parameters')
    l=np.zeros((7,7))
    for k,(i,j) in enumerate(TRI):l[i,j]=10*(np.exp(theta[2+k]) if i==j else theta[2+k])
    return l,l.dot(l.T)


def parameters_initial(features,adc,train,lambda0=2.,timing_theta=None):
    nf=features.shape[2];theta=np.zeros(39+4*nf)
    old=base.initial(1,.5) if timing_theta is None else np.asarray(timing_theta,float)
    theta[:2]=old[:2];theta[30:35]=old[8:13]
    old_l,_=base.precision(old)
    for k,(i,j) in enumerate(TRI):
        if i<3:theta[2+k]=np.log(old_l[i,j]/10.) if i==j else old_l[i,j]/10.
        elif i==j:theta[2+k]=np.log(.5)
    theta[35:39]=np.log(lambda0)
    signal=.5*np.sum(np.log(adc[train]),axis=2)
    alpha=theta[39:].reshape(4,nf)
    alpha[:,0]=np.median(signal,axis=0)
    return theta


def bounds(nfeatures):
    b=[(0.,2.)]*2
    b += [(-3.,5.) if i==j else (-20.,20.) for i,j in TRI]
    b += base.BOUNDS[8:13]+[(-3.,3.)]*4
    b += [(0.,12.) if j==0 else (-2.,2.) for layer in range(4) for j in range(nfeatures)]
    return b


def observation(data,adc,time_field='t_native'):
    adc=np.asarray(adc,float)
    if adc.shape!=(len(data[time_field]),4,2) or not np.isfinite(adc).all() or np.any(adc<=0):raise ValueError('Finite positive end ADC required')
    return np.column_stack([np.asarray(data[time_field],float).dot(base.D.T),.5*np.sum(np.log(adc),axis=2)])


class JointADC(base.ConditionalModel):
    def __init__(self,data,geometry,select,features,adc,penalty=.01,**kwargs):
        super().__init__(data,geometry,select,**kwargs)
        choose=np.flatnonzero(select) if np.asarray(select).dtype==bool else np.asarray(select,int)
        self.features=np.ascontiguousarray(np.asarray(features,float)[choose]);self.nfeatures=self.features.shape[2]
        if self.features.shape[:2]!=(len(choose),4) or not np.isfinite(self.features).all():raise ValueError('Invalid joint response features')
        ad=np.asarray(adc,float)[choose]
        if not np.isfinite(ad).all() or np.any(ad<=0):raise ValueError('Invalid selected ADC')
        self.y=np.ascontiguousarray(np.column_stack([self.y,.5*np.sum(np.log(ad),axis=2)]))
        self.arrays[0]=self.y;self.pointers=[v.ctypes.data_as(self.ptr) for v in self.arrays]
        self.nparameters=39+4*self.nfeatures;self.penalty=float(penalty)
        self.penalty_mask=np.zeros(self.nparameters);self.penalty_mask[39:]=np.tile(np.r_[0.,np.ones(self.nfeatures-1)],4)
        self.library=ct.CDLL(str(OUT/'libjointadc.so'));self.function=self.library.joint_adc_loss
        self.function.argtypes=[self.ptr,ct.c_int,ct.c_int,ct.c_int]+[self.ptr]*15+[ct.c_int,ct.c_double]
        self.function.restype=ct.c_double
        self.feature_pointer=self.features.ctypes.data_as(self.ptr)

    def evaluate(self,theta,events=False,regularized=True):
        theta=np.ascontiguousarray(theta,dtype=np.float64)
        if theta.shape!=(self.nparameters,) or not np.isfinite(theta).all():raise ValueError('Joint parameter shape/values')
        grad=np.zeros(self.nparameters);per=np.zeros(len(self.y))
        value=self.function(theta.ctypes.data_as(self.ptr),len(self.y),self.nodes,self.nfeatures,self.feature_pointer,*self.pointers,
            grad.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),int(self.threads),self.rigidity_nu)
        if regularized:
            value+=.5*self.penalty*np.dot(self.penalty_mask,theta*theta);grad+=self.penalty*self.penalty_mask*theta
        if not np.isfinite(value) or not np.isfinite(grad).all():raise FloatingPointError('Nonfinite joint response likelihood')
        return (value,grad,per) if events else (value,grad)

    def numpy_event(self,theta,index,y_override=None):
        y=self.y[index] if y_override is None else np.asarray(y_override,float);l,Q=precision(theta);z=self.z
        hmix=theta[32]+theta[33]*z
        mix=np.stack([-np.logaddexp(0,hmix),-np.logaddexp(0,-hmix)])
        rr=(self.r[index]-self.q)/(self.sigma[index]*np.exp(theta[34]))
        rlog=-.5*rr*rr if self.rigidity_nu<=0 else -.5*(self.rigidity_nu+1)*np.log1p(rr*rr/self.rigidity_nu)
        h=self.logw+theta[30]*z+theta[31]*z*z+rlog+mix
        mean=np.zeros((2,self.nodes,7));mean[:,:,:3]=self.x[index][None,None,:]*self.u[:,:,None]
        for g in range(2):mean[:,:,:3]+=theta[g]*(self.moment0[index,:,g][None,None,:]*self.a[:,:,None]+self.moment1[index,:,g][None,None,:]*self.b[:,:,None])
        offset=np.sum(self.features[index]*theta[39:].reshape(4,self.nfeatures),axis=1)
        mean[:,:,3:]=offset[None,None,:]+(self.u[:,:,None]-UREF)*np.exp(theta[35:39])[None,None,:]
        residual=y-mean;radius=np.sum(np.matmul(residual,l)**2,axis=-1)
        logtime=gammaln(6)-gammaln(2.5)-3.5*np.log(5*np.pi)+np.log(np.diag(l)).sum()-6*np.log1p(radius/5)
        norm=base.logsumexp(h);loglike=base.logsumexp(h+logtime)-norm
        return dict(loglike=float(loglike),means=mean,posterior=np.exp(h+logtime-base.logsumexp(h+logtime)),conditional_rigidity_weights=np.exp(h-norm))


def reconstruct(theta,x,y,features):
    """Continuous seven-observable GLS; isotope/source/R parameters are unused."""
    x=np.asarray(x,float);y=np.asarray(y,float);features=np.asarray(features,float)
    if x.ndim!=2 or x.shape[1]!=4 or y.shape!=(len(x),7) or features.shape[:2]!=(len(x),4):raise ValueError('Joint reconstruction shapes')
    if np.asarray(theta).shape!=(39+4*features.shape[2],):raise ValueError('Joint reconstruction parameter dimension')
    l,Q=precision(theta);cov=np.linalg.inv(Q);lam=np.exp(theta[35:39]);at=x.dot(base.D.T)
    aa=np.broadcast_to(lam,(len(x),4));direction=np.column_stack([at,aa])
    offset=np.sum(features*theta[39:].reshape(4,features.shape[2])[None,:,:],axis=2)-UREF*lam[None,:]
    centered=y-np.column_stack([np.zeros((len(x),3)),offset])
    weights=direction.dot(Q);weights/=np.sum(weights*direction,axis=1)[:,None]
    joint=np.sum(weights*centered,axis=1)
    qt=np.linalg.inv(cov[:3,:3]);wt=at.dot(qt);wt/=np.sum(wt*at,axis=1)[:,None]
    qa=np.linalg.inv(cov[3:,3:]);wa=aa.dot(qa);wa/=np.sum(wa*aa,axis=1)[:,None]
    time=np.sum(wt*y[:,:3],axis=1);amplitude=np.sum(wa*centered[:,3:],axis=1)
    information=np.einsum('ni,ij,nj->n',direction,Q,direction);timeinfo=np.einsum('ni,ij,nj->n',at,qt,at)
    condition=cov[3:,3:]-cov[3:,:3].dot(qt).dot(cov[:3,3:]);contrast=aa-at.dot(qt).dot(cov[:3,3:])
    gain=np.einsum('ni,ij,nj->n',contrast,np.linalg.inv(condition),contrast)
    return dict(joint_u=joint,time_u=time,adc_u=amplitude,weights=weights,time_weights=wt.dot(base.D),adc_weights=wa,
        time_injection_gain=np.sum(weights[:,:3]*at,axis=1),adc_injection_gain=np.sum(weights[:,3:]*aa,axis=1),
        information=information,time_information=timeinfo,additional_adc_information=gain,covariance=cov,direction=direction,centered=centered)
