"""Algebraic checks for the flight-calibration research primitive."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import flight_calibration_core as c


class GlobalCalibrationTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.RandomState(917)
        self.x=np.array([-2.2,-2.1,2.1,2.2])[None,:]*rng.uniform(1,1.4,(25,1))
        self.b=rng.normal(size=(25,4,3))
        self.sigma=rng.uniform(.03,.07,(25,4))
        self.t=rng.normal(size=(25,4))

    def test_profile_matches_joint_least_squares(self):
        fit=c.solve_global(self.x,self.t,self.b,self.sigma)
        n=len(self.x)
        design=np.zeros((n,4,3+2*n))
        design[:,:,:3]=self.b
        for e in range(n):
            design[e,:,3+2*e]=1
            design[e,:,3+2*e+1]=self.x[e]
        whitened=(design/self.sigma[:,:,None]).reshape(-1,3+2*n)
        target=(self.t/self.sigma).ravel()
        solution=np.linalg.lstsq(whitened,target,rcond=None)[0]
        np.testing.assert_allclose(solution[:3],fit['theta'],atol=2e-12)
        covariance=np.linalg.inv(whitened.T.dot(whitened))[:3,:3]
        np.testing.assert_allclose(covariance,fit['covariance_on_estimable_subspace'],atol=2e-12)

    def test_affine_time_response(self):
        theta=np.array([.02,-.01,.03])
        raw=c.reconstruct(self.x,self.t,self.b,self.sigma,theta)
        shifted=c.reconstruct(self.x,self.t+.004*self.x+7,self.b,self.sigma,theta)
        np.testing.assert_allclose(shifted['u']-raw['u'],.004,atol=2e-14)
        np.testing.assert_allclose(shifted['t0']-raw['t0'],7,atol=2e-14)

    def test_unidentified_slope_rejected(self):
        b=np.concatenate([self.b,self.x[:,:,None]],axis=-1)
        with self.assertRaises(ValueError):c.solve_global(self.x,self.t,b,self.sigma)
        fit=c.solve_global(self.x,self.t,b,self.sigma,allow_rank_deficient=True,penalty=np.eye(4))
        self.assertEqual(fit['data_rank'],3)
        self.assertFalse(fit['covariance_valid_all_parameters'])

    def test_penalty_covariance_is_sandwich(self):
        lam=100*np.eye(3)
        fit=c.solve_global(self.x,self.t,self.b,self.sigma,penalty=lam)
        inv=np.linalg.inv(fit['fisher']+lam)
        np.testing.assert_allclose(inv,fit['penalized_precision_inverse'],rtol=1e-10)
        np.testing.assert_allclose(inv.dot(fit['fisher']).dot(inv),fit['covariance_on_estimable_subspace'],rtol=1e-10)

    def test_projection_and_reference_origin(self):
        x=self.x
        basis=np.stack([np.ones_like(x),x],axis=-1)
        np.testing.assert_allclose(c.project(x,basis,self.sigma),0,atol=5e-14)
        b0=c.impulse_kernel(x,0,center=0)
        b1=c.impulse_kernel(x,0,center=-1)
        np.testing.assert_allclose(c.project(x,b0,self.sigma),c.project(x,b1,self.sigma),atol=5e-14)

    def test_invalid_noise_and_penalty(self):
        with self.assertRaises(ValueError):c.solve_global(self.x,self.t,self.b,0.)
        with self.assertRaises(ValueError):c.solve_global(self.x,self.t,self.b,1.,penalty=-np.eye(3))


if __name__=='__main__':unittest.main()
