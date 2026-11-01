#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
import flight_conditional_likelihood as base
import feature_conditional as feature
from research_conditional_output import weights as base_weights


class FeatureTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.RandomState(913);n=8
        x=np.tile([-2.2,-1.8,1.7,2.15],(n,1))*rng.uniform(1,1.3,(n,1))
        pos=rng.normal(0,20,(n,9,3));direction=rng.normal(0,.1,(n,9,3));direction[:,:,2]=-.95
        self.a=dict(run=np.arange(n),event=np.arange(n),R=rng.uniform(3,5,n),relR=np.full(n,.1),
            t_native=x*1.13+rng.normal(0,.025,(n,4)),tk_pos=pos,tk_dir=direction)
        self.g=dict(x=x,moment=rng.uniform(.01,.1,(n,4,2,9)))
        self.choose=np.ones(n,bool);self.features,self.norm=feature.prepare_features(self.a,self.choose)
        self.theta=np.r_[base.initial(),np.zeros(36)]

    def test_base_model_and_gradient(self):
        for nu in [0,5]:
            original=base.ConditionalModel(self.a,self.g,self.choose,nodes=64,rigidity_nu=nu,threads=1)
            new=feature.FeatureConditional(self.a,self.g,self.choose,self.features,nodes=64,rigidity_nu=nu,threads=1)
            first=original.evaluate(self.theta[:13],events=True);second=new.evaluate(self.theta,events=True)
            self.assertAlmostEqual(first[0],second[0],places=12)
            np.testing.assert_allclose(first[1],second[1][:13],atol=1e-12,rtol=1e-12)
            np.testing.assert_allclose(first[2],second[2],atol=1e-12,rtol=1e-12)

    def test_independent_density_and_all_gradients(self):
        rng=np.random.RandomState(41);theta=self.theta.copy();theta[13:]=rng.normal(0,.03,36)
        for nu in [0,5]:
            obj=feature.FeatureConditional(self.a,self.g,self.choose,self.features,nodes=64,rigidity_nu=nu,threads=1,penalty=.01)
            loss,grad,per=obj.evaluate(theta,events=True)
            independent=np.array([-obj.numpy_event(theta,i)['loglike'] for i in range(len(self.choose))])
            np.testing.assert_allclose(independent,per,atol=1e-12,rtol=1e-12)
            numeric=[]
            for j in range(len(theta)):
                up=theta.copy();dn=theta.copy();up[j]+=1e-5;dn[j]-=1e-5
                numeric.append((obj.evaluate(up)[0]-obj.evaluate(dn)[0])/2e-5)
            np.testing.assert_allclose(grad,numeric,atol=2e-7,rtol=2e-6)

    def test_output_baseline_and_coherent_response(self):
        rng=np.random.RandomState(10);x=self.g['x'];t=self.a['t_native']
        weights,q=feature.precision_weights(self.theta,x,self.features)
        np.testing.assert_allclose(weights,base_weights(self.theta[:13],x),atol=1e-13,rtol=1e-13)
        theta=self.theta.copy();theta[13:]=rng.normal(0,.05,36)
        weights,q=feature.precision_weights(theta,x,self.features)
        common=rng.normal(size=len(x));slope=rng.normal(size=len(x));injected=t+common[:,None]+slope[:,None]*x
        np.testing.assert_allclose(np.sum(weights*(injected-t),axis=1),slope,atol=2e-13,rtol=0)
        self.assertTrue(np.all(np.linalg.eigvalsh(q)>0))

    def test_variance_gap_identity(self):
        rng=np.random.RandomState(5);x=self.g['x'];theta=self.theta.copy();theta[13:]=rng.normal(0,.05,36)
        w,q=feature.precision_weights(theta,x,self.features);covariance=np.linalg.inv(q)
        a=x.dot(base.D.T);v=w.dot(base.D.T)
        delta=rng.normal(size=a.shape);delta-=a*np.sum(delta*a,axis=1)[:,None]/np.sum(a*a,axis=1)[:,None]
        proposed=v+delta
        first=np.einsum('ni,nij,nj->n',v,covariance,v);second=np.einsum('ni,nij,nj->n',proposed,covariance,proposed)
        gap=np.einsum('ni,nij,nj->n',delta,covariance,delta)
        np.testing.assert_allclose(second-first,gap,atol=1e-13,rtol=1e-12)
        self.assertTrue(np.all(gap>0))

    def test_normalization_frozen_and_invalid_shapes(self):
        changed={k:v.copy() for k,v in self.a.items()};changed['t_native']+=100
        ff,norm=feature.prepare_features(changed,self.choose,self.norm)
        np.testing.assert_array_equal(ff,self.features)
        with self.assertRaises(ValueError):feature.local_parameters(self.theta[:-1],self.features)


if __name__=='__main__':unittest.main()
