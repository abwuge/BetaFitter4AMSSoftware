#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
from joint_adc_scaled_optimizer import solve_scaled,projected_gradient

class Quadratic:
    def __init__(self,precision,target):self.precision=np.asarray(precision,float);self.target=np.asarray(target,float)
    def evaluate(self,x,regularized=True):
        d=np.asarray(x)-self.target
        return float(.5*d.dot(self.precision).dot(d)),self.precision.dot(d)

class ScaledOptimizerTests(unittest.TestCase):
    def test_ill_conditioned_quadratic(self):
        weights=np.logspace(-6,6,12);target=np.linspace(-.8,.8,12);obj=Quadratic(np.diag(weights),target)
        result=solve_scaled(obj,np.zeros(12),[(-1,1)]*12)
        self.assertTrue(result['success']);self.assertLess(result['penalized_nll'],1e-10)
        np.testing.assert_allclose(result['theta'],target,atol=2e-5,rtol=2e-5)

    def test_bound_projection(self):
        target=np.array([-3.,2.,.4]);obj=Quadratic(np.diag([.01,100.,1]),target)
        result=solve_scaled(obj,np.zeros(3),[(-1.,1.)]*3)
        np.testing.assert_allclose(result['theta'],[-1,1,.4],atol=2e-8)
        self.assertLess(result['projected_gradient_max'],1e-7)
        self.assertGreater(result['gradient_max'],10)

    def test_positive_definite_coupled(self):
        rng=np.random.RandomState(382);matrix=rng.normal(size=(8,8));precision=matrix.T.dot(matrix)+.1*np.eye(8);target=rng.uniform(-.5,.5,8)
        result=solve_scaled(Quadratic(precision,target),np.zeros(8),[(-1.,1.)]*8)
        self.assertTrue(result['success']);np.testing.assert_allclose(result['theta'],target,atol=1e-6,rtol=1e-6)

if __name__=='__main__':unittest.main()
