"""Geometry-dependent timing precision with a time-equivariant free output."""
import ctypes as ct
import numpy as np
import flight_conditional_likelihood as base

OUT=base.ROOT/'results/flight_calibration_research/feature_noise'
FEATURE_NAMES=['L5_x','L5_y','L5_x_squared','L5_y_squared','direction_x_over_z','direction_y_over_z']


def raw_features(data):
    pos=np.asarray(data['tk_pos'],float)[:,4];direction=np.asarray(data['tk_dir'],float)[:,4]
    if not np.isfinite(pos).all() or not np.isfinite(direction).all() or np.any(abs(direction[:,2])<.1):raise ValueError('Invalid measured Tracker geometry')
    xx=pos[:,0]/50.;yy=pos[:,1]/50.
    return np.column_stack([xx,yy,xx*xx,yy*yy,direction[:,0]/direction[:,2],direction[:,1]/direction[:,2]])


def prepare_features(data,train,normalization=None):
    raw=raw_features(data)
    if normalization is None:
        center=raw[train].mean(0);scale=raw[train].std(0)
        if np.any(scale<=1e-10):raise ValueError('Constant feature')
        normalization=dict(names=FEATURE_NAMES,center=center.tolist(),scale=scale.tolist(),clip=3.,source='Training measurements only')
    center=np.asarray(normalization['center']);scale=np.asarray(normalization['scale'])
    if normalization['names']!=FEATURE_NAMES:raise ValueError('Feature definition differs')
    transformed=np.clip((raw-center)/scale,-normalization['clip'],normalization['clip'])
    return transformed,normalization


def local_parameters(theta,features):
    theta=np.asarray(theta,float);features=np.asarray(features,float)
    if features.ndim!=2 or theta.shape!=(13+6*features.shape[1],):raise ValueError('Feature coefficient shape differs')
    local=np.broadcast_to(theta[:13],(len(features),13)).copy()
    local[:,2:8]+=features.dot(theta[13:].reshape(6,-1).T)
    return local


def precision_weights(theta,x,features):
    local=local_parameters(theta,features);n=len(x);ll=np.zeros((n,3,3))
    for j in range(3):ll[:,j,j]=10*np.exp(local[:,2+j])
    ll[:,1,0]=10*local[:,5];ll[:,2,0]=10*local[:,6];ll[:,2,1]=10*local[:,7]
    precision=np.matmul(ll,np.swapaxes(ll,1,2));a=np.asarray(x).dot(base.D.T)
    qa=np.einsum('nij,nj->ni',precision,a);den=np.einsum('ni,ni->n',a,qa)
    if np.any(den<=0) or not np.isfinite(den).all():raise FloatingPointError('Invalid output precision')
    weights=qa.dot(base.D)/den[:,None]
    return weights,precision


class FeatureConditional(base.ConditionalModel):
    def __init__(self,data,geometry,select,features,penalty=0.,**kwargs):
        super().__init__(data,geometry,select,**kwargs)
        choose=np.flatnonzero(select) if np.asarray(select).dtype==bool else np.asarray(select,int)
        self.features=np.ascontiguousarray(np.asarray(features)[choose],dtype=np.float64)
        if self.features.ndim!=2 or not np.isfinite(self.features).all():raise ValueError('Finite features required')
        self.nf=self.features.shape[1];self.penalty=float(penalty)
        if self.penalty<0:raise ValueError('Penalty must be nonnegative')
        self.libfeature=ct.CDLL(str(OUT/'libfeature.so'));self.feature_function=self.libfeature.conditional_feature_loss
        self.feature_function.argtypes=[self.ptr,ct.c_int,ct.c_int,ct.c_int]+[self.ptr]*15+[ct.c_int,ct.c_double]
        self.feature_function.restype=ct.c_double

    def evaluate(self,theta,events=False,regularized=True):
        theta=np.ascontiguousarray(theta,dtype=np.float64)
        if theta.shape!=(13+6*self.nf,) or not np.isfinite(theta).all():raise ValueError('Invalid feature-model parameter vector')
        grad=np.zeros_like(theta);per=np.zeros(len(self.y))
        loss=self.feature_function(theta.ctypes.data_as(self.ptr),len(self.y),self.nodes,self.nf,self.features.ctypes.data_as(self.ptr),*self.pointers,
            grad.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),self.threads,self.rigidity_nu)
        if regularized:
            loss+=.5*self.penalty*np.dot(theta[13:],theta[13:]);grad[13:]+=self.penalty*theta[13:]
        if not np.isfinite(loss) or not np.isfinite(grad).all():raise FloatingPointError('Nonfinite feature likelihood')
        return (loss,grad,per) if events else (loss,grad)

    def numpy_event(self,theta,index,y_override=None):
        local=local_parameters(theta,self.features[index:index+1])[0]
        return super().numpy_event(local,index,y_override)
