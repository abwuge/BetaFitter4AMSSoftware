#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
import feature_conditional as feature
import conditional_readout_population as population
import test_amplitude_noise_20261013 as fixture
import bootstrap_amplitude_noise_20261013 as boot

class PopulationTests(unittest.TestCase):
    def setUp(self):
        f=fixture.AmplitudeTests();f.setUp()
        self.a=f.a;self.g=f.g;self.ff=f.ff;self.theta=np.r_[f.theta,np.zeros(4*self.ff.shape[1])];self.domain=f.domain;self.norm=f.norm

    def test_reduces_to_noise_only(self):
        for nu in [0.,5.]:
            old=feature.FeatureConditional(self.a,self.g,self.domain,self.ff,nodes=64,rigidity_nu=nu,threads=1,penalty=.03)
            new=population.PopulationConditional(self.a,self.g,self.domain,self.ff,nodes=64,rigidity_nu=nu,threads=1,penalty=.03)
            before=old.evaluate(self.theta[:121],events=True);after=new.evaluate(self.theta,events=True)
            self.assertAlmostEqual(before[0],after[0],places=12)
            np.testing.assert_allclose(before[1],after[1][:121],atol=1e-12,rtol=1e-12)
            np.testing.assert_allclose(before[2],after[2],atol=1e-12,rtol=1e-12)

    def test_independent_density_all_193_gradients(self):
        rng=np.random.RandomState(1729);theta=self.theta.copy();theta[121:]=rng.normal(0,.08,72)
        for nu in [0.,5.]:
            obj=population.PopulationConditional(self.a,self.g,self.domain,self.ff,nodes=64,rigidity_nu=nu,threads=1,penalty=.03)
            loss,grad,per=obj.evaluate(theta,events=True)
            other=np.array([-obj.numpy_event(theta,i)['loglike'] for i in range(len(self.ff))])
            np.testing.assert_allclose(per,other,atol=2e-12,rtol=1e-12)
            numeric=[]
            for j in range(len(theta)):
                up=theta.copy();dn=theta.copy();up[j]+=1e-5;dn[j]-=1e-5
                numeric.append((obj.evaluate(up)[0]-obj.evaluate(dn)[0])/2e-5)
            np.testing.assert_allclose(grad,numeric,atol=2e-7,rtol=4e-6)

    def test_population_never_enters_event_output(self):
        rng=np.random.RandomState(20);other=self.theta.copy();other[8:13]=rng.normal(0,10,5);other[121:]=rng.normal(0,10,72)
        first=population.precision_weights(self.theta,self.g['x'],self.ff)
        second=population.precision_weights(other,self.g['x'],self.ff)
        for a,b in zip(first,second):np.testing.assert_array_equal(a,b)
        np.testing.assert_allclose(first[0].sum(1),0,atol=1e-13)
        np.testing.assert_allclose(np.sum(first[0]*self.g['x'],axis=1),1,atol=1e-13)

    def test_population_node_weights_normalize(self):
        theta=self.theta.copy();theta[121:]=np.linspace(-.1,.1,72)
        obj=population.PopulationConditional(self.a,self.g,self.domain,self.ff,nodes=64,threads=1)
        for i in range(len(self.ff)):
            result=obj.numpy_event(theta,i)
            self.assertAlmostEqual(result['conditional_rigidity_weights'].sum(),1,places=12)
            self.assertAlmostEqual(result['posterior'].sum(),1,places=12)

    def test_bootstrap_rebase_without_clipping(self):
        theta=self.theta[:121];old=self.norm;new=dict(old);new['center']=(np.array(old['center'])+.03).tolist();new['scale']=(np.array(old['scale'])*1.1).tolist()
        rebased=boot.rebase(theta,old,new)
        rng=np.random.RandomState(22);raw=np.array(old['center'])+rng.uniform(-.1,.1,(4,18))*np.array(old['scale'])
        fo=(raw-np.array(old['center']))/np.array(old['scale']);fn=(raw-np.array(new['center']))/np.array(new['scale'])
        np.testing.assert_allclose(feature.local_parameters(theta,fo),feature.local_parameters(rebased,fn),atol=1e-13,rtol=1e-13)

if __name__=='__main__':unittest.main()
