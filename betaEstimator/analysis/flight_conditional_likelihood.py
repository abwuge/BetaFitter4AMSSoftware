"""Conditional raw-time working model with mass and true-rigidity integration.

The density is normalized over the three orthonormal time contrasts given
measured inverse rigidity, its stored error and geometry. It models the
selected sample as a conditional population. Native hit-selection effects,
time-calibration covariance, deposit response and source-period variation
still need separate treatment before use as a flight calibration.
"""
import ctypes as ct
import hashlib
import math
from pathlib import Path
import numpy as np
from scipy.special import gammaln,expit
import flight_mean_loss as physics

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/flight_calibration_research/mean_loss_likelihood/conditional'
D=np.array([[1,-1,0,0],[1,1,-2,0],[1,1,1,-3]],float)/np.sqrt(np.array([2,6,12]))[:,None]
MASSES=np.array([12.,13.00335483507])*.9314940954-6*.0005109989461
PARAMETERS=['kappa_TOF','kappa_inner','precision_log00','precision_log11','precision_log22',
            'precision_off10','precision_off20','precision_off21','spectrum_logq','spectrum_logq2',
            'C13_logodds','C13_logodds_logq','log_rigidity_error_scale']
BOUNDS=[(0.,2.),(0.,2.)]+[(-2.,5.)]*3+[(-20.,20.)]*3+[(-10.,10.)]*2+[(-8.,8.)]*2+[(math.log(.25),math.log(4.))]

def logsumexp(value,axis=None):
    value=np.asarray(value,float);maximum=np.max(value,axis=axis,keepdims=True)
    result=maximum+np.log(np.sum(np.exp(value-maximum),axis=axis,keepdims=True))
    return result.reshape(()).item() if axis is None else np.squeeze(result,axis=axis)

def kernels(q):
    q=np.asarray(q,float);u=np.sqrt(1+(MASSES[:,None]*q[None,:]/6)**2)
    beta=1/u;b2=beta*beta;g2=1/(1-b2);me=.511;mass=MASSES[:,None]*1000
    tmax=2*me*b2*g2/(1+2*np.sqrt(g2)*me/mass+(me/mass)**2)
    prefactor=(MASSES[:,None]**2*q[None,:]**3/6**3)*(.000307075*36/b2)
    aa=prefactor*(.5*np.log(2*me*b2*g2*tmax)-b2)
    bb=-prefactor
    return u,aa,bb

def initial(kappa=1.,fraction=.1):
    t=np.zeros(13);t[:2]=kappa;t[2:5]=np.log(2.5);t[8]=1.;t[10]=np.log(fraction/(1-fraction));return t

def precision(theta):
    l=np.zeros((3,3));l[np.diag_indices(3)]=10*np.exp(theta[2:5]);l[1,0]=10*theta[5];l[2,0]=10*theta[6];l[2,1]=10*theta[7]
    return l,l.dot(l.T)

class ConditionalModel:
    def __init__(self,data,geometry,select,time_field='t_native',nodes=128,qrange=(.02,.7),rigidity_nu=0.,threads=4):
        choose=np.flatnonzero(select) if np.asarray(select).dtype==bool else np.asarray(select,int)
        if not len(choose):raise ValueError('Empty conditional sample')
        self.identities={k:np.asarray(data[k])[choose].copy() for k in ['run','event']}
        self.y=np.asarray(data[time_field][choose],float).dot(D.T)
        self.x=np.asarray(geometry['x'][choose],float).dot(D.T)
        mm=np.einsum('il,nlgk->nigk',D,geometry['moment'][choose])
        self.moment0=np.sum(mm*physics.ZA,axis=-1)
        self.moment1=np.sum(mm*physics.ZA*np.log(physics.I_MEV),axis=-1)
        self.r=1/np.asarray(data['R'][choose],float)
        self.sigma=np.asarray(data['relR'][choose],float)*self.r
        leg,w=np.polynomial.legendre.leggauss(nodes);self.q=qrange[0]+(leg+1)*.5*(qrange[1]-qrange[0])
        self.logw=np.log(w*.5*(qrange[1]-qrange[0]));self.z=np.log(self.q/.25)
        self.u,self.a,self.b=kernels(self.q)
        names=['y','x','moment0','moment1','r','sigma','q','logw','z','u','a','b']
        self.arrays=[]
        for name in names:
            value=np.ascontiguousarray(getattr(self,name),dtype=np.float64)
            if not np.isfinite(value).all():raise ValueError('Nonfinite '+name)
            setattr(self,name,value);self.arrays.append(value)
        if np.any(self.sigma<=0):raise ValueError('Invalid rigidity uncertainty')
        self.nodes=nodes;self.qrange=qrange;self.rigidity_nu=float(rigidity_nu);self.threads=threads
        self.library=ct.CDLL(str(OUT/'libconditional.so'))
        self.function=self.library.conditional_loss;ptr=ct.POINTER(ct.c_double)
        self.function.argtypes=[ptr,ct.c_int,ct.c_int]+[ptr]*14+[ct.c_int,ct.c_double]
        self.function.restype=ct.c_double
        self.pointers=[v.ctypes.data_as(ptr) for v in self.arrays];self.ptr=ptr

    def evaluate(self,theta,events=False):
        theta=np.ascontiguousarray(theta,dtype=np.float64)
        if theta.shape!=(13,) or not np.isfinite(theta).all():raise ValueError('Thirteen finite parameters required')
        grad=np.zeros(13);per=np.zeros(len(self.y))
        loss=self.function(theta.ctypes.data_as(self.ptr),len(self.y),self.nodes,*self.pointers,
             grad.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),int(self.threads),self.rigidity_nu)
        if not np.isfinite(loss) or not np.isfinite(grad).all():raise FloatingPointError('Nonfinite likelihood')
        return (loss,grad,per) if events else (loss,grad)

    def numpy_event(self,theta,index,y_override=None):
        """Independent density implementation for numerical verification."""
        y=self.y[index] if y_override is None else np.asarray(y_override,float)
        l,p=precision(theta);z=self.z
        hmix=theta[10]+theta[11]*z;pi=expit(hmix)
        mix=np.stack([1-pi,pi]);rr=(self.r[index]-self.q)/(self.sigma[index]*np.exp(theta[12]))
        rlog=-.5*rr*rr if self.rigidity_nu<=0 else -.5*(self.rigidity_nu+1)*np.log1p(rr*rr/self.rigidity_nu)
        h=self.logw+theta[8]*z+theta[9]*z*z+rlog+np.log(mix)
        base=self.x[index][None,None,:]*self.u[:,:,None]
        response=np.zeros_like(base)
        for g in range(2):response+=theta[g]*(self.moment0[index,:,g][None,None,:]*self.a[:,:,None]+self.moment1[index,:,g][None,None,:]*self.b[:,:,None])
        residual=y-base-response;radius=np.sum(np.matmul(residual,l)**2,axis=-1)
        logtime=gammaln(4)-gammaln(2.5)-1.5*np.log(5*np.pi)+np.log(np.diag(l)).sum()-4*np.log1p(radius/5)
        norm=logsumexp(h);loglike=logsumexp(h+logtime)-norm
        posterior=np.exp(h+logtime-logsumexp(h+logtime));prior=np.exp(h-norm)
        return dict(loglike=float(loglike),posterior=posterior,conditional_rigidity_weights=prior,means=base+response)
