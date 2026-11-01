#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import free_timing_transport as m
import flight_conditional_likelihood as f
import flight_mean_loss as p


class FreeTransportTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.RandomState(27181);self.n=32
        self.moment=rng.uniform(.1,2,(self.n,4,2,9));self.R=np.linspace(2.5,6,self.n)
        self.x=np.tile([-2.2,-2.05,2.1,2.25],(self.n,1))
        self.w=m.weights(f.initial(),self.x)

    def test_scalar_kernel(self):
        u=np.linspace(1.001,1.7,self.n);a,b=m.compress(self.moment)
        vec=m.continuous_kernel(u,self.R,a,b)
        for i in range(self.n):
            mass=6*self.R[i]*np.sqrt(u[i]*u[i]-1)
            np.testing.assert_allclose(vec[i],p.first_order(u[i],mass,self.moment[i]),rtol=1e-11,atol=1e-12)

    def test_zero_loss_exact(self):
        s=np.linspace(.8,1.6,self.n)
        result=m.solve_center(s,self.R,self.w,self.moment,[0,0])
        np.testing.assert_array_equal(s,result['u']);self.assertTrue(result['converged'].all())

    def test_beta_above_one_preserved(self):
        s=np.linspace(.8,1,self.n)
        result=m.solve_center(s,self.R,self.w,self.moment,[2,2])
        np.testing.assert_array_equal(result['u'],s)

    def test_forward_inverse(self):
        u=np.linspace(.95,1.5,self.n);a,b=m.compress(self.moment);k=np.array([1.2,.8])
        s=u+np.einsum('nl,nlg,g->n',self.w,m.continuous_kernel(u,self.R,a,b),k)
        result=m.solve_center(s,self.R,self.w,self.moment,k)
        np.testing.assert_allclose(result['u'],u,atol=3e-12,rtol=0)

    def test_common_time_shift(self):
        theta=f.initial();t=1.2*self.x
        a=m.center_from_times(theta,self.x,t,self.R,self.moment)
        b=m.center_from_times(theta,self.x,t+np.arange(self.n)[:,None],self.R,self.moment)
        np.testing.assert_allclose(a['u'],b['u'],atol=1e-12,rtol=0)

    def test_free_weight_constraints(self):
        np.testing.assert_allclose(np.sum(self.w,axis=1),0,atol=1e-15)
        np.testing.assert_allclose(np.sum(self.w*self.x,axis=1),1,atol=1e-15)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):m.solve_center(np.ones(self.n),self.R,self.w,self.moment,[-1,1])
        a,b=m.compress(self.moment)
        with self.assertRaises(ValueError):m.continuous_kernel(np.ones(self.n),-self.R,a,b)


if __name__=='__main__':unittest.main()
