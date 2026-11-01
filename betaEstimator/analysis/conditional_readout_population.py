"""Normalized conditional readout-population model; no posterior in speed output."""
import ctypes as ct
from pathlib import Path
import numpy as np
import flight_conditional_likelihood as base
import feature_conditional as feature

ROOT=base.ROOT
OUT=ROOT/'results/flight_calibration_research/continuation_20261013/readout_population'


def local_parameters(theta,features):
    theta=np.asarray(theta,float);features=np.asarray(features,float)
    if features.ndim!=2 or theta.shape!=(13+10*features.shape[1],):raise ValueError('Population parameter/feature shape mismatch')
    value=np.broadcast_to(theta[:13],(len(features),13)).copy()
    value[:,2:12]+=features.dot(theta[13:].reshape(10,features.shape[1]).T)
    return value


def precision_weights(theta,x,features):
    # Deliberately strip all latent-population coefficients before calling output.
    nf=np.asarray(features).shape[1];theta=np.asarray(theta,float)
    if theta.shape!=(13+10*nf,):raise ValueError('Population parameter dimension')
    return feature.precision_weights(theta[:13+6*nf],x,features)


class PopulationConditional(base.ConditionalModel):
    def __init__(self,data,geometry,select,features,penalty=.03,**kwargs):
        super().__init__(data,geometry,select,**kwargs)
        choose=np.flatnonzero(select) if np.asarray(select).dtype==bool else np.asarray(select,int)
        self.features=np.ascontiguousarray(np.asarray(features,float)[choose]);self.nfeatures=self.features.shape[1]
        if self.features.ndim!=2 or not np.isfinite(self.features).all():raise ValueError('Finite 2-D features required')
        self.nparameters=13+10*self.nfeatures;self.penalty=float(penalty)
        self.library=ct.CDLL(str(OUT/'libpopulation.so'));self.function=self.library.conditional_population_loss
        self.function.argtypes=[self.ptr,ct.c_int,ct.c_int,ct.c_int]+[self.ptr]*15+[ct.c_int,ct.c_double]
        self.function.restype=ct.c_double
        self.feature_pointer=self.features.ctypes.data_as(self.ptr)

    def evaluate(self,theta,events=False,regularized=True):
        theta=np.ascontiguousarray(theta,dtype=np.float64)
        if theta.shape!=(self.nparameters,) or not np.isfinite(theta).all():raise ValueError('Finite population parameter vector required')
        gradient=np.zeros(self.nparameters);per=np.zeros(len(self.y))
        loss=self.function(theta.ctypes.data_as(self.ptr),len(self.y),self.nodes,self.nfeatures,self.feature_pointer,*self.pointers,
             gradient.ctypes.data_as(self.ptr),per.ctypes.data_as(self.ptr),int(self.threads),self.rigidity_nu)
        if regularized:
            loss+=.5*self.penalty*np.dot(theta[13:],theta[13:]);gradient[13:]+=self.penalty*theta[13:]
        if not np.isfinite(loss) or not np.isfinite(gradient).all():raise FloatingPointError('Nonfinite population likelihood')
        return (loss,gradient,per) if events else (loss,gradient)

    def numpy_event(self,theta,index,y_override=None):
        local=local_parameters(theta,self.features[index:index+1])[0]
        return base.ConditionalModel.numpy_event(self,local,index,y_override=y_override)
