import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import free_curve_timing as curve
import free_timing_transport as material


class FreeCurveTiming(unittest.TestCase):
    def fixture(self):
        rng=np.random.RandomState(273)
        n=20;d=rng.normal(size=(n,7));d[:,:3]*=3;d[:,3:]=.5
        b=rng.normal(size=(7,7));q=b.dot(b.T)+np.eye(7)
        r=np.linspace(2.5,6,n);moment=rng.uniform(0,.015,(n,4,2,9))
        moment[:,:2]*=-1
        return rng,d,q,r,moment

    def test_zero_loss_recovers_exact_linear_solution(self):
        rng,d,q,r,m=self.fixture();y=rng.normal(size=d.shape)+d*1.2
        exact=np.einsum('ni,ij,nj->n',d,q,y)/np.einsum('ni,ij,nj->n',d,q,d)
        result=curve.refine(exact+.04,y,d,q,r,m,np.zeros(2))
        np.testing.assert_allclose(result['u'],exact,rtol=0,atol=2e-10)
        self.assertTrue(result['converged'].all())

    def test_noiseless_material_signal_is_recovered(self):
        _,d,q,r,m=self.fixture();u=np.linspace(1.01,1.5,len(r));k=np.array([1.,2.])
        y=curve.mean_curve(u,r,d,*material.compress(m),k)
        result=curve.refine(u+.03,y,d,q,r,m,k)
        np.testing.assert_allclose(result['u'],u,rtol=0,atol=3e-10)
        self.assertTrue(result['converged'].all())

    def test_objective_never_increases(self):
        rng,d,q,r,m=self.fixture();u=np.linspace(1.05,1.5,len(r));k=np.ones(2)
        y=curve.mean_curve(u,r,d,*material.compress(m),k)+rng.normal(0,.07,d.shape)
        result=curve.refine(u,y,d,q,r,m,k)
        self.assertLessEqual(np.max(result['objective']-result['initial_objective']),1e-10)
        self.assertTrue(result['converged'].all())

    def test_beta_above_one_is_not_clipped(self):
        _,d,q,r,m=self.fixture();u=np.linspace(.7,.98,len(r));y=u[:,None]*d
        result=curve.refine(u-.01,y,d,q,r,m,np.ones(2))
        np.testing.assert_allclose(result['u'],u,rtol=0,atol=3e-10)
        self.assertTrue((1/result['u']>1).all())

    def test_bad_shapes_rejected(self):
        _,d,q,r,m=self.fixture()
        with self.assertRaises(ValueError):curve.refine(np.ones(20),d,d,q,r,m,np.ones(3))
        with self.assertRaises(ValueError):curve.refine(np.ones(20),d,d,-q,r,m,np.ones(2))


if __name__=='__main__':unittest.main()
