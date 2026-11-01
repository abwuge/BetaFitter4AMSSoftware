"""Geometry/era-conditioned time noise in the normalized seven-dimensional model."""
import ctypes as ct
import numpy as np
from scipy.special import gammaln
import joint_adc_likelihood as joint
import periodic_joint_adc as period
import flight_conditional_likelihood as base

OUT=period.OUT/'conditional_noise'
VARIANTS=['constant','geometry','geometry_period']


def features(data,train,dataset,variant,normalization=None):
    if variant not in VARIANTS:raise ValueError('Unknown conditional-noise variant')
    pos=np.asarray(data['tk_pos'],float)[:,4,:2]/50.;x,y=pos.T
    geometry=np.column_stack([x,y,x*x,y*y,x*y]);n=len(x)
    era=period.period_basis(data['run'],dataset,'step')
    raw=np.empty((n,0)) if variant=='constant' else (geometry if variant=='geometry' else np.column_stack([geometry,era]))
    if normalization is None:
        center=raw[train].mean(0);scale=raw[train].std(0)
        if np.any(scale<1e-10):raise ValueError('Degenerate conditional-noise feature')
        normalization=dict(dataset=dataset,variant=variant,center=center.tolist(),scale=scale.tolist(),clip=3.,
            geometry_names=joint.GEOMETRY_NAMES,period_features='Four step indicators with the first UTC era as reference; MC eras use the separately declared deterministic run hash')
    if normalization['dataset']!=dataset or normalization['variant']!=variant:raise ValueError('Conditional-noise normalization mismatch')
    result=np.clip((raw-np.asarray(normalization['center']))/np.asarray(normalization['scale']),-3,3)
    if not np.isfinite(result).all():raise ValueError('Nonfinite conditional-noise features')
    return result,normalization


def matrices(theta,nresponse,noise_features):
    theta=np.asarray(theta,float);f=np.asarray(noise_features,float);nbase=39+4*nresponse
    if f.ndim!=2 or len(theta)!=nbase+3*f.shape[1]:raise ValueError('Conditional-noise parameter dimensions')
    l,_=joint.precision(theta[:nbase]);eta=f.dot(theta[nbase:].reshape(3,f.shape[1]).T)
    scales=np.column_stack([np.exp(eta),np.ones((len(f),4))])
    ls=l[None,:,:]*scales[:,None,:];q=np.matmul(ls,np.swapaxes(ls,1,2))
    return ls,q


class ConditionalNoise(joint.JointADC):
    def __init__(self,data,geometry,select,response_features,adc,noise_features,noise_penalty=.1,**kwargs):
        super().__init__(data,geometry,select,response_features,adc,**kwargs)
        chosen=np.flatnonzero(select) if np.asarray(select).dtype==bool else np.asarray(select,int)
        self.noise_features=np.ascontiguousarray(np.asarray(noise_features,float)[chosen]);self.nnoise=self.noise_features.shape[1]
        self.nbase=self.nparameters;self.nparameters+=3*self.nnoise
        self.regularization=np.r_[self.penalty*self.penalty_mask,np.full(3*self.nnoise,noise_penalty)]
        self.noise_pointer=self.noise_features.ctypes.data_as(self.ptr)
        self.library=ct.CDLL(str(OUT/'libconditionalnoise.so'));self.function=self.library.conditional_time_noise_loss
        self.function.argtypes=[self.ptr,ct.c_int,ct.c_int,ct.c_int,ct.c_int]+[self.ptr]*16+[ct.c_int,ct.c_double]
        self.function.restype=ct.c_double

    def evaluate(self,theta,events=False,regularized=True):
        theta=np.ascontiguousarray(theta,dtype=np.float64)
        if theta.shape!=(self.nparameters,) or not np.isfinite(theta).all():raise ValueError('Invalid conditional-noise parameters')
        gradient=np.zeros(self.nparameters);per=np.zeros(len(self.y))
        value=self.function(theta.ctypes.data_as(self.ptr),len(self.y),self.nodes,self.nfeatures,self.nnoise,
            self.feature_pointer,self.noise_pointer,*self.pointers,gradient.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),int(self.threads),self.rigidity_nu)
        if regularized:
            value+=.5*np.dot(self.regularization,theta*theta);gradient+=self.regularization*theta
        if not np.isfinite(value) or not np.isfinite(gradient).all():raise FloatingPointError('Nonfinite conditional-noise objective')
        return (value,gradient,per) if events else (value,gradient)

    def numpy_event(self,theta,index,y_override=None):
        prior=super().numpy_event(np.asarray(theta)[:self.nbase],index,y_override)
        y=self.y[index] if y_override is None else np.asarray(y_override,float)
        l,_=matrices(theta,self.nfeatures,self.noise_features[index:index+1]);l=l[0]
        radius=np.sum(np.matmul(y-prior['means'],l)**2,axis=-1)
        density=gammaln(6)-gammaln(2.5)-3.5*np.log(5*np.pi)+np.log(np.diag(l)).sum()-6*np.log1p(radius/5)
        with np.errstate(divide='ignore'):logweights=np.log(prior['conditional_rigidity_weights'])
        return float(base.logsumexp(logweights+density))


def reconstruct(theta,x,y,response_features,noise_features):
    theta=np.asarray(theta,float);x=np.asarray(x,float);y=np.asarray(y,float);f=np.asarray(response_features,float)
    nbase=39+4*f.shape[2];original=joint.reconstruct(theta[:nbase],x,y,f)
    _,q=matrices(theta,f.shape[2],noise_features);cov=np.linalg.inv(q);d=original['direction'];centered=original['centered']
    w=np.einsum('nij,nj->ni',q,d);info=np.sum(w*d,axis=1);w/=info[:,None]
    qt=np.linalg.inv(cov[:,:3,:3]);wt=np.einsum('nij,nj->ni',qt,d[:,:3]);it=np.sum(wt*d[:,:3],axis=1);wt/=it[:,None]
    result=dict(joint_u=np.sum(w*centered,axis=1),weights=w,time_u=np.sum(wt*y[:,:3],axis=1),
        time_weights=wt.dot(base.D),information=info,time_information=it,time_injection_gain=np.sum(w[:,:3]*d[:,:3],axis=1),
        direction=d,centered=centered,precision=q,covariance=cov)
    for variant in ['speed_null','speed_gain_null']:
        columns=[np.column_stack([d[:,:3],np.zeros((len(d),4))]),np.column_stack([np.zeros((len(d),3)),d[:,3:]])]
        if variant=='speed_gain_null':columns.append(np.column_stack([np.zeros((len(d),3)),np.ones((len(d),4))]))
        c=np.stack(columns,axis=2);qc=np.einsum('nij,njk->nik',q,c);gram=np.einsum('nik,nil->nkl',c,qc)
        target=np.zeros((len(d),len(columns)));target[:,0]=1
        wc=np.einsum('nik,nk->ni',qc,np.linalg.solve(gram,target));response=np.einsum('nik,ni->nk',c,wc)
        if np.max(abs(response-target))>2e-10:raise ValueError('Conditional-noise output constraints failed')
        result[variant]=dict(u=np.sum(wc*centered,axis=1),weights=wc,time_weights=wc[:,:3].dot(base.D),response=response)
    return result
