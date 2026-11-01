import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import layer_noise_mixture as mixture
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import flight_conditional_likelihood as timing


class LayerNoiseMixture(unittest.TestCase):
    def fixture(self,n=40):
        rng=np.random.RandomState(612);x=np.tile([-2.3,-1.8,1.7,2.4],(n,1))
        features=np.ones((n,4,1));theta=joint.parameters_initial(features,np.full((n,4,2),100.),np.ones(n,bool))
        theta[35:39]=np.log([.4,.5,.65,.7]);theta[39:]=4.6
        for k,(i,j) in enumerate(joint.TRI):
            if i>=3 and j<3:theta[2+k]=rng.normal(0,.07)
        par=np.r_[np.full(4,np.log(.05/.8)),np.log(.3),0.]
        y=np.column_stack([x.dot(timing.D.T)*1.17,np.ones((n,4))*4.6])+rng.normal(0,.05,(n,7))
        return rng,x,features,theta,par,y

    def test_zero_contamination_recovers_linear_outputs(self):
        rng,x,f,theta,par,y=self.fixture();par[4]=-50
        for variant in mixture.VARIANTS:
            actual=mixture.reconstruct(theta,par,x,y,f,variant)
            expected=joint.reconstruct(theta,x,y,f)['joint_u'] if variant=='joint' else constrained.reconstruct(theta,x,y,f,variant)['u']
            np.testing.assert_allclose(actual['u'],expected,rtol=0,atol=2e-11)
            np.testing.assert_allclose(actual['component_posterior'],np.broadcast_to(mixture.probabilities(par),actual['component_posterior'].shape),rtol=0,atol=2e-11)

    def test_response_and_responsibilities_are_shift_invariant(self):
        rng,x,f,theta,par,y=self.fixture()
        for variant in mixture.VARIANTS:
            first=mixture.reconstruct(theta,par,x,y,f,variant);c=mixture.constraint_matrix(first['direction'],variant)
            coefficients=rng.normal(0,.07,(len(y),c.shape[2]));shift=np.einsum('nik,nk->ni',c,coefficients)
            second=mixture.reconstruct(theta,par,x,y+shift,f,variant)
            np.testing.assert_allclose(second['u']-first['u'],coefficients[:,0],rtol=0,atol=2e-10)
            np.testing.assert_allclose(second['component_posterior'],first['component_posterior'],rtol=0,atol=2e-10)

    def test_integrated_density_matches_joint_component_weights(self):
        rng,x,f,theta,par,y=self.fixture(1);rec=mixture.reconstruct(theta,par,x,y,f,'joint')
        prior,cov,q,ld=mixture.matrices(theta,par);d=rec['direction'][0];z=rec['centered'][0]
        values=[];constant=gammaln(6)-gammaln(2.5)-3.5*np.log(5*np.pi)
        for k in range(5):
            center=d.dot(q[k]).dot(z)/d.dot(q[k]).dot(d)
            def density(offset):
                residual=z-d*(center+offset)
                return np.exp(constant-.5*ld[k]-6*np.log1p(residual.dot(q[k]).dot(residual)/5))
            integral,error=quad(density,-np.inf,np.inf,epsabs=1e-8,epsrel=2e-10,limit=200)
            values.append(prior[k]*integral)
        values=np.asarray(values);values/=values.sum()
        np.testing.assert_allclose(rec['component_posterior'][0],values,rtol=1e-8,atol=1e-10)

    def test_positive_definite_normalized_components(self):
        rng,x,f,theta,par,y=self.fixture();prior,cov,q,ld=mixture.matrices(theta,par)
        self.assertAlmostEqual(float(prior.sum()),1,places=14)
        self.assertTrue((prior>0).all());self.assertTrue((np.linalg.eigvalsh(cov)>0).all())
        np.testing.assert_allclose(np.matmul(cov,q),np.broadcast_to(np.eye(7),(5,7,7)),rtol=0,atol=3e-12)
        for k in range(1,5):self.assertLess(np.linalg.matrix_rank(cov[k]-cov[0],tol=1e-10),2)

    def test_material_inverse_with_zero_loss_is_identical(self):
        rng,x,f,theta,par,y=self.fixture();mom=np.zeros((len(x),4,2,9));r=np.full(len(x),4.)
        for variant in mixture.VARIANTS:
            direct=mixture.reconstruct(theta,par,x,y,f,variant)
            actual=mixture.reconstruct_transport(theta,par,x,y,f,r,mom,kappa=np.zeros(2),variant=variant)
            self.assertTrue(actual['converged'].all())
            np.testing.assert_allclose(actual['u'],direct['u'],rtol=0,atol=2e-11)

    def test_beta_above_one_is_not_clipped(self):
        rng,x,f,theta,par,y=self.fixture();base=joint.reconstruct(theta,x,y,f)
        shift=(.9-base['time_u'])[:,None]*np.column_stack([base['direction'][:,:3],np.zeros((len(x),4))])
        actual=mixture.reconstruct(theta,par,x,y+shift,f,'speed_null')
        self.assertTrue((actual['u']<1).all())

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):mixture.probabilities(np.zeros(5))
        with self.assertRaises(ValueError):mixture.constraint_matrix(np.ones((1,7)),'isotope_posterior')


if __name__=='__main__':unittest.main()
