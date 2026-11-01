"""Layer-contamination mixture with reference-free equivariant event estimates.

Calibration integrates the original physical conditional density. Prediction
instead marginalizes unrestricted continuous response coordinates analytically:
its component weights depend only on redundant residuals, never on the event's
isotope or source posterior. The two-coordinate branch has exact unit response
to time-slope shifts and exact zero response to the declared ADC-speed shift.
"""
import ctypes as ct
import numpy as np
from scipy.special import gammaln
import flight_conditional_likelihood as timing
import joint_adc_likelihood as joint
import free_timing_transport as material
import periodic_joint_adc as period

OUT=period.OUT/'layer_mixture'
BOUNDS=[(-10.,2.)]*4+[(np.log(.01),np.log(3.)),(np.log(.5),np.log(2.))]
VARIANTS=['joint','speed_null','speed_gain_null']


def probabilities(theta):
    theta=np.asarray(theta,float)
    if theta.shape!=(6,) or not np.isfinite(theta).all():raise ValueError('Six finite mixture parameters required')
    logits=np.r_[0.,theta[:4]];logits-=np.max(logits);weights=np.exp(logits);return weights/weights.sum()


def matrices(base_theta,theta):
    weights=probabilities(theta);l,_=joint.precision(base_theta);l[:,:3]*=np.exp(theta[5])
    core=np.linalg.inv(l.dot(l.T));v=np.column_stack([timing.D.T,np.zeros((4,4))])
    cov=np.broadcast_to(core,(5,7,7)).copy();cov[1:]+=np.exp(2*theta[4])*v[:,:,None]*v[:,None,:]
    q=np.linalg.inv(cov);sign,ld=np.linalg.slogdet(cov)
    if not (sign>0).all():raise ValueError('Nonpositive mixture covariance')
    return weights,cov,q,ld


def make_cache(model,base_theta,chunk=128):
    theta=np.asarray(base_theta,float);l,_=joint.precision(theta);n=len(model.y);s=2*model.nodes
    cache=np.empty((n,s,6));prior=np.empty((n,s));normalization=np.empty(n)
    directions=timing.D.T.dot(l[:3,:3]);h2=np.sum(directions*directions,axis=1)
    q=model.q;z=model.z;hmix=theta[32]+theta[33]*z
    logmix=np.stack([-np.logaddexp(0,hmix),-np.logaddexp(0,-hmix)])
    constant=float(gammaln(6)-gammaln(2.5)-3.5*np.log(5*np.pi)+np.log(np.diag(l)).sum())
    for low in range(0,n,chunk):
        high=min(n,low+chunk);rr=(model.r[low:high,None]-q)/(model.sigma[low:high,None]*np.exp(theta[34]))
        rlog=-.5*rr*rr if model.rigidity_nu<=0 else -.5*(model.rigidity_nu+1)*np.log1p(rr*rr/model.rigidity_nu)
        hp=model.logw+theta[30]*z+theta[31]*z*z+rlog[:,None,:]+logmix[None,:,:]
        prior[low:high]=hp.reshape(high-low,s);normalization[low:high]=timing.logsumexp(prior[low:high],axis=1)
        mean=np.empty((high-low,2,model.nodes,7));mean[:,:,:,:3]=model.x[low:high,None,None,:]*model.u[None,:,:,None]
        for group in range(2):
            mean[:,:,:,:3]+=theta[group]*(model.moment0[low:high,:,group][:,None,None,:]*model.a[None,:,:,None]+model.moment1[low:high,:,group][:,None,None,:]*model.b[None,:,:,None])
        offset=np.sum(model.features[low:high]*theta[39:].reshape(4,model.nfeatures),axis=2)
        mean[:,:,:,3:]=offset[:,None,None,:]+(model.u[None,:,:,None]-joint.UREF)*np.exp(theta[35:39])
        whiten=np.matmul(model.y[low:high,None,None,:]-mean,l).reshape(high-low,s,7)
        cache[low:high,:,0]=np.sum(whiten[:,:,:3]**2,axis=2);cache[low:high,:,1]=np.sum(whiten[:,:,3:]**2,axis=2)
        cache[low:high,:,2:]=np.matmul(whiten[:,:,:3],directions.T)
    if not all(np.isfinite(value).all() for value in [cache,prior,normalization,h2]):raise ValueError('Nonfinite mixture cache')
    return dict(cache=cache,logprior=prior,norm=normalization,h2=h2,logconstant=np.array(constant))


class LayerMixture:
    def __init__(self,cache,threads=4):
        self.values={key:np.ascontiguousarray(value,dtype=np.float64) for key,value in cache.items()}
        self.n,self.states,self.columns=self.values['cache'].shape
        if self.columns!=6 or self.values['logprior'].shape!=(self.n,self.states):raise ValueError('Mixture cache shape differs')
        self.threads=threads;self.ptr=ct.POINTER(ct.c_double);self.library=ct.CDLL(str(OUT/'liblayermixture.so'))
        self.function=self.library.layer_noise_mixture_loss
        self.function.argtypes=[self.ptr,ct.c_int,ct.c_int]+[self.ptr]*4+[ct.c_double,self.ptr,self.ptr,ct.c_int]
        self.function.restype=ct.c_double
        self.pointers=[self.values[key].ctypes.data_as(self.ptr) for key in ['cache','norm','logprior','h2']]

    def evaluate(self,theta,events=False):
        theta=np.ascontiguousarray(theta,dtype=np.float64);probabilities(theta);gradient=np.zeros(6);per=np.zeros(self.n)
        value=self.function(theta.ctypes.data_as(self.ptr),self.n,self.states,*self.pointers,float(self.values['logconstant']),
            gradient.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),self.threads)
        if not np.isfinite(value) or not np.isfinite(gradient).all():raise FloatingPointError('Nonfinite mixture likelihood')
        return (value,gradient,per) if events else (value,gradient)


def constraint_matrix(direction,variant):
    d=np.asarray(direction,float);n=len(d)
    if d.shape!=(n,7) or variant not in VARIANTS:raise ValueError('Mixture response dimensions')
    if variant=='joint':return d[:,:,None]
    c=np.zeros((n,7,2 if variant=='speed_null' else 3));c[:,:3,0]=d[:,:3];c[:,3:,1]=d[:,3:]
    if variant=='speed_gain_null':c[:,3:,2]=1.
    return c


def estimate_centered(centered,direction,base_theta,theta,variant='speed_null'):
    y=np.asarray(centered,float);d=np.asarray(direction,float)
    if y.shape!=d.shape or y.ndim!=2 or y.shape[1]!=7:raise ValueError('Mixture estimator dimensions')
    prior,cov,precision,logdet=matrices(base_theta,theta);c=constraint_matrix(d,variant);rank=c.shape[2];n=len(y)
    locations=np.empty((n,5));logweights=np.empty((n,5));weights=np.empty((n,5,7))
    for k,q in enumerate(precision):
        qc=np.einsum('ij,njk->nik',q,c);gram=np.einsum('nik,nil->nkl',c,qc)
        sign,ld=np.linalg.slogdet(gram)
        if not (sign>0).all():raise ValueError('Dependent free response coordinates')
        fitted=np.linalg.solve(gram,np.einsum('nik,ni->nk',qc,y))
        residual=y-np.einsum('nik,nk->ni',c,fitted);radius=np.einsum('ni,ij,nj->n',residual,q,residual)
        locations[:,k]=fitted[:,0]
        target=np.zeros((n,rank));target[:,0]=1.;weights[:,k]=np.einsum('nik,nk->ni',qc,np.linalg.solve(gram,target))
        logweights[:,k]=np.log(prior[k])-.5*logdet[k]-.5*ld-.5*(12-rank)*np.log1p(np.maximum(radius,0)/5)
    logweights-=timing.logsumexp(logweights,axis=1)[:,None];posterior=np.exp(logweights)
    effective=np.sum(posterior[:,:,None]*weights,axis=1);u=np.sum(posterior*locations,axis=1)
    response=np.einsum('nik,ni->nk',c,effective);target=np.zeros_like(response);target[:,0]=1
    if np.max(abs(response-target))>2e-9 or not np.isfinite(u).all():raise ValueError('Mixture estimator response failed')
    return dict(u=u,weights=effective,time_weights=effective[:,:3].dot(timing.D),component_posterior=posterior,
        component_u=locations,response=response,free_coordinate_count=rank)


def reconstruct(base_theta,theta,x,y,features,variant='speed_null'):
    original=joint.reconstruct(base_theta,x,y,features)
    result=estimate_centered(original['centered'],original['direction'],base_theta,theta,variant)
    result.update(centered=original['centered'],direction=original['direction'])
    return result


def reconstruct_transport(base_theta,theta,x,y,features,rigidity,moment,kappa=None,variant='speed_null',tolerance=2e-10,max_iterations=60):
    rec=reconstruct(base_theta,theta,x,y,features,variant);k=np.asarray(base_theta[:2] if kappa is None else kappa,float)
    initial=material.solve_center(rec['u'],rigidity,rec['time_weights'],moment,k)
    if not initial['converged'].all():raise ValueError('Initial mixture mean-loss inverse failed')
    u=initial['u'].copy();m0,m1=material.compress(moment);delta=np.full(len(u),np.inf);last=rec
    for iteration in range(max_iterations):
        loss=np.einsum('nlg,g->nl',material.continuous_kernel(u,rigidity,m0,m1),k).dot(timing.D.T)
        shifted=rec['centered'].copy();shifted[:,:3]-=loss
        last=estimate_centered(shifted,rec['direction'],base_theta,theta,variant);delta=last['u']-u
        if np.max(abs(delta))<=tolerance:u=last['u'];break
        u+=.7*delta
    converged=abs(delta)<=tolerance
    return dict(u=u,converged=converged,fixed_point_error=delta,iterations=iteration+1,
        component_posterior=last['component_posterior'],time_weights=last['time_weights'],raw_u=rec['u'])
