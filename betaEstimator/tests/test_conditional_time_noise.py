import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import conditional_time_noise as noise
import joint_adc_likelihood as joint
import periodic_joint_adc as period


class ConditionalTimeNoise(unittest.TestCase):
    def fixture(self):
        rng=np.random.RandomState(315);n=100
        data=dict(run=np.linspace(period.EDGES[0]+10,period.EDGES[-1]-10,n).astype(np.int64),
            tk_pos=rng.normal(0,15,(n,9,3)),bar=rng.randint(0,10,(n,4)))
        train=np.arange(n)%5!=0;f,_=joint.response_features(data,train)
        theta=joint.parameters_initial(f,np.full((n,4,2),200.),train)
        theta[35:39]=np.log([.4,.5,.65,.6])
        for k,(i,j) in enumerate(joint.TRI):
            if i>=3 and j<3:theta[2+k]=rng.normal(0,.1)
        x=np.tile([-2.3,-1.8,1.7,2.4],(n,1));y=rng.normal(0,.1,(n,7))
        h,_=noise.features(data,train,'flight','geometry_period')
        return rng,data,train,f,theta,x,y,h

    def test_zero_feature_coefficients_reproduce_base(self):
        rng,data,train,f,theta,x,y,h=self.fixture()
        actual=noise.reconstruct(np.r_[theta,np.zeros(3*h.shape[1])],x,y,f,h)
        expected=joint.reconstruct(theta,x,y,f)
        for key in ['joint_u','time_u','weights','time_weights']:
            np.testing.assert_allclose(actual[key],expected[key],rtol=0,atol=2e-12)

    def test_marginal_adc_and_regression_are_preserved(self):
        rng,data,train,f,theta,x,y,h=self.fixture();par=np.r_[theta,rng.normal(0,.08,3*h.shape[1])]
        ls,q=noise.matrices(par,f.shape[2],h);cov=np.linalg.inv(q)
        reference=np.linalg.inv(joint.precision(theta)[1]);regression=reference[:3,3:].dot(np.linalg.inv(reference[3:,3:]))
        np.testing.assert_allclose(cov[:,3:,3:],np.broadcast_to(reference[3:,3:],cov[:,3:,3:].shape),rtol=1e-11,atol=1e-12)
        actual=np.matmul(cov[:,:3,3:],np.linalg.inv(cov[:,3:,3:]))
        np.testing.assert_allclose(actual,np.broadcast_to(regression,actual.shape),rtol=1e-11,atol=1e-12)
        self.assertTrue((np.linalg.eigvalsh(cov)>0).all())

    def test_exact_continuous_signal_response(self):
        rng,data,train,f,theta,x,y,h=self.fixture();par=np.r_[theta,rng.normal(0,.05,3*h.shape[1])]
        first=noise.reconstruct(par,x,y,f,h);delta=np.linspace(-.03,.03,len(x))
        second=noise.reconstruct(par,x,y+delta[:,None]*first['direction'],f,h)
        np.testing.assert_allclose(second['joint_u']-first['joint_u'],delta,rtol=0,atol=2e-12)
        for name in ['speed_null','speed_gain_null']:
            target=np.zeros_like(first[name]['response']);target[:,0]=1
            np.testing.assert_allclose(first[name]['response'],target,rtol=0,atol=2e-10)

    def test_heldout_positions_do_not_refit_normalization(self):
        rng,data,train,f,theta,x,y,h=self.fixture();_,norm=noise.features(data,train,'flight','geometry_period')
        changed={k:v.copy() for k,v in data.items()};changed['tk_pos'][~train]*=8
        _,again=noise.features(changed,train,'flight','geometry_period');self.assertEqual(norm,again)

    def test_bad_dimension_is_rejected(self):
        rng,data,train,f,theta,x,y,h=self.fixture()
        with self.assertRaises(ValueError):noise.matrices(theta,f.shape[2],h)
        with self.assertRaises(ValueError):noise.features(data,train,'flight','absent')


if __name__=='__main__':unittest.main()
