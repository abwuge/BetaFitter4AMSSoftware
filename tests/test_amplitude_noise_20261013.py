#!/usr/bin/env python3
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
import flight_conditional_likelihood as base
import feature_conditional as feature
import research_amplitude_noise_20261013 as amplitude
from test_feature_conditional import FeatureTests

class AmplitudeTests(unittest.TestCase):
    def setUp(self):
        obj=FeatureTests();obj.setUp();self.a=obj.a;self.g=obj.g
        rng=np.random.RandomState(814);n=len(self.a['R'])
        self.a['edep']=np.exp(rng.normal(4,.2,(n,4)))
        self.adc=np.exp(rng.normal(6,.4,(n,4,2)))
        self.domain=np.ones(n,bool);self.train=np.arange(n)<6
        self.ff,self.norm=amplitude.feature_arrays(self.a,self.adc,self.domain,self.train,'readout')
        self.theta=np.r_[base.initial(),rng.normal(0,.02,6*self.ff.shape[1])]

    def test_full_121_gradient_and_numpy_density(self):
        self.assertEqual(self.ff.shape[1],18);self.assertEqual(len(self.theta),121)
        for nu in [0.,5.]:
            obj=feature.FeatureConditional(self.a,self.g,self.domain,self.ff,nodes=64,threads=1,rigidity_nu=nu,penalty=.03)
            loss,grad,per=obj.evaluate(self.theta,events=True)
            independently=np.array([-obj.numpy_event(self.theta,i)['loglike'] for i in range(len(self.ff))])
            np.testing.assert_allclose(per,independently,rtol=1e-12,atol=1e-12)
            numeric=[]
            for j in range(len(self.theta)):
                up=self.theta.copy();dn=self.theta.copy();up[j]+=1e-5;dn[j]-=1e-5
                numeric.append((obj.evaluate(up)[0]-obj.evaluate(dn)[0])/2e-5)
            np.testing.assert_allclose(grad,numeric,atol=2e-7,rtol=3e-6)

    def test_holdout_not_used_in_normalization(self):
        a={k:v.copy() for k,v in self.a.items()};adc=self.adc.copy()
        a['edep'][~self.train]*=1000;adc[~self.train]*=10000
        ff,norm=amplitude.feature_arrays(a,adc,self.domain,self.train,'readout')
        self.assertEqual(norm,self.norm)
        np.testing.assert_array_equal(ff[self.train],self.ff[self.train])

    def test_time_shift_invariance(self):
        a={k:v.copy() for k,v in self.a.items()};a['t_native']+=17
        ff,norm=amplitude.feature_arrays(a,self.adc,self.domain,self.train,'readout',self.norm)
        np.testing.assert_array_equal(ff,self.ff)
        w,q=feature.precision_weights(self.theta,self.g['x'],ff)
        np.testing.assert_allclose(w.sum(1),0,atol=2e-13)
        np.testing.assert_allclose(np.sum(w*self.g['x'],axis=1),1,atol=2e-13)
        self.assertTrue((np.linalg.eigvalsh(q)>0).all())

    def test_missing_ADC_does_not_create_selected_features(self):
        adc=self.adc.copy();adc[-1]=np.nan;domain=self.domain.copy();domain[-1]=False
        ff,norm=amplitude.feature_arrays(self.a,adc,domain,self.train,'readout',self.norm)
        np.testing.assert_array_equal(ff[-1],0)
        np.testing.assert_array_equal(ff[:-1],self.ff[:-1])

if __name__=='__main__':unittest.main()
