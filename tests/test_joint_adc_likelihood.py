#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
from scipy.special import gammaln
from scipy.integrate import quad
import joint_adc_likelihood as joint
import test_feature_conditional as fixture

class JointADCTests(unittest.TestCase):
    def setUp(self):
        f=fixture.FeatureTests();f.setUp();self.a=f.a;self.g=f.g;self.choose=f.choose
        rng=np.random.RandomState(818);n=len(self.choose)
        self.a['bar']=rng.randint(0,8,(n,4));self.adc=np.exp(rng.normal(6,.15,(n,4,2)))
        self.ff,self.norm=joint.response_features(self.a,self.choose)
        self.theta=joint.parameters_initial(self.ff,self.adc,self.choose,lambda0=1.6)
        for k,(i,j) in enumerate(joint.TRI):
            if i!=j:self.theta[2+k]=rng.normal(0,.05)
        self.theta[39:].reshape(4,16)[:,1:]=rng.normal(0,.02,(4,15))
        self.y=joint.observation(self.a,self.adc)

    def test_numpy_density_all_103_gradients(self):
        self.assertEqual(len(self.theta),103)
        for nu in [0.,5.]:
            obj=joint.JointADC(self.a,self.g,self.choose,self.ff,self.adc,nodes=64,threads=1,rigidity_nu=nu,penalty=.01)
            loss,grad,per=obj.evaluate(self.theta,events=True)
            other=np.array([-obj.numpy_event(self.theta,i)['loglike'] for i in range(len(self.choose))])
            np.testing.assert_allclose(per,other,atol=2e-12,rtol=2e-12)
            numeric=[]
            for j in range(len(self.theta)):
                up=self.theta.copy();dn=self.theta.copy();up[j]+=1e-5;dn[j]-=1e-5
                numeric.append((obj.evaluate(up)[0]-obj.evaluate(dn)[0])/2e-5)
            np.testing.assert_allclose(grad,numeric,atol=2e-7,rtol=5e-6)

    def test_joint_and_separate_injections(self):
        first=joint.reconstruct(self.theta,self.g['x'],self.y,self.ff)
        delta=np.linspace(-.1,.1,len(self.y));second=joint.reconstruct(self.theta,self.g['x'],self.y+delta[:,None]*first['direction'],self.ff)
        np.testing.assert_allclose(second['joint_u']-first['joint_u'],delta,atol=2e-13)
        injected=self.y.copy();injected[:,:3]+=delta[:,None]*first['direction'][:,:3]
        time=joint.reconstruct(self.theta,self.g['x'],injected,self.ff)
        np.testing.assert_allclose(time['joint_u']-first['joint_u'],delta*first['time_injection_gain'],atol=2e-13)
        np.testing.assert_allclose(time['time_u']-first['time_u'],delta,atol=2e-13)
        np.testing.assert_allclose(first['time_injection_gain']+first['adc_injection_gain'],1,atol=2e-13)
        above=joint.reconstruct(self.theta,self.g['x'],self.y-.6*first['direction'],self.ff)
        self.assertTrue((above['joint_u']<1).any())

    def test_information_decomposition_and_optimality_gap(self):
        a=joint.reconstruct(self.theta,self.g['x'],self.y,self.ff)
        np.testing.assert_allclose(a['information'],a['time_information']+a['additional_adc_information'],atol=1e-9,rtol=1e-12)
        self.assertTrue((a['additional_adc_information']>=0).all())
        rng=np.random.RandomState(71);delta=rng.normal(size=a['weights'].shape);v=a['direction']
        delta-=v*np.sum(delta*v,axis=1)[:,None]/np.sum(v*v,axis=1)[:,None]
        old=np.einsum('ni,ij,nj->n',a['weights'],a['covariance'],a['weights'])
        new=np.einsum('ni,ij,nj->n',a['weights']+delta,a['covariance'],a['weights']+delta)
        gap=np.einsum('ni,ij,nj->n',delta,a['covariance'],delta)
        np.testing.assert_allclose(new-old,gap,atol=1e-12,rtol=2e-12)
        np.testing.assert_allclose(old,1/a['information'],atol=1e-13,rtol=1e-12)

    def test_no_mass_or_source_posterior_in_output(self):
        theta=self.theta.copy();theta[:2]=[0,0];theta[30:35]=[-1,3,4,-5,.7]
        a=joint.reconstruct(self.theta,self.g['x'],self.y,self.ff);b=joint.reconstruct(theta,self.g['x'],self.y,self.ff)
        for k in ['joint_u','time_u','adc_u','weights']:np.testing.assert_array_equal(a[k],b[k])

    def test_seven_dimensional_radial_normalization(self):
        d=7.;nu=5.
        logconstant=gammaln((nu+d)/2)-gammaln(nu/2)-d/2*np.log(nu*np.pi)
        surface=2*np.pi**(d/2)/np.exp(gammaln(d/2))
        value,error=quad(lambda r:surface*np.exp(logconstant)*r**(d-1)*(1+r*r/nu)**(-(nu+d)/2),0,np.inf,epsabs=1e-11,epsrel=1e-12)
        self.assertAlmostEqual(value,1,places=10)
        self.assertLess(error,1e-9)

if __name__=='__main__':unittest.main()
