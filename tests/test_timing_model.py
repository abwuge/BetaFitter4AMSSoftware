"""Tests of the truth boundary, objective gradient, and detector response."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'analysis'))
from extract_inputs import ARRAY_KEYS, BRANCHES, load_arrays
from timing_model import (C, clean, contrasts, fit, kinematics, likelihood,
                          predict, prepare_features, reconstruct)


def sample(n=5000, seed=45):
    rng = np.random.RandomState(seed)
    x = np.tile([-2.3, -2.1, 2.1, 2.3], (n, 1)) * rng.uniform(1, 1.3, (n, 1))
    edep = np.exp(rng.normal(1, .55, (n,4)))
    pos = rng.uniform(-50, 50, (n,4,2))
    # Latent event speed is used only to synthesize times and score held-out
    # reconstruction. It is never a field passed to training.
    u = rng.uniform(.9, 1.8, n)
    sigma = .08 * np.sqrt(3/edep)
    t = rng.normal(0, 10, (n,1)) + x*u[:,None] + sigma*rng.normal(size=(n,4))
    xc = x-x.mean(axis=1)[:,None]
    measured_u = np.sum(xc*t,axis=1)/np.sum(xc*xc,axis=1)
    d = dict(tof_leng=x*C, tof_tl=t, tof_edep=edep, position=pos,
             tof_barid=rng.randint(0,8,(n,4)), tof_bpatt=np.full(n,4444),
             tof_betah=1/measured_u)
    return d,u


class TimingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.train, _ = sample()
        cls.model = fit(cls.train, 1)

    def test_gradient(self):
        d,_ = sample(24)
        x,t,_ = kinematics(d)
        q = contrasts(x)
        f,_ = prepare_features(d,1)
        y = np.einsum('nri,ni->nr',q,t)
        p = np.random.RandomState(19).normal(0,.07,8*f.shape[2])
        _,g = likelihood(p,f,q,y,.08)
        for i in range(0,len(p),7):
            step=np.zeros_like(p);step[i]=1e-6
            numeric=(likelihood(p+step,f,q,y,.08)[0]-likelihood(p-step,f,q,y,.08)[0])/2e-6
            self.assertAlmostEqual(g[i],numeric,places=6)

    def test_slope_and_offset_response(self):
        d,_=sample(200,91)
        original,mean,_=predict(d,self.model)
        shifted=dict(d);shifted['tof_tl']=d['tof_tl']+10+.07*d['tof_leng']/C
        result,_,_=predict(shifted,self.model)
        np.testing.assert_allclose(result,original+.07,rtol=0,atol=2e-13)
        x,_,_=kinematics(d)
        np.testing.assert_allclose(mean.sum(axis=1),0,atol=1e-13)
        np.testing.assert_allclose((mean*x).sum(axis=1),0,atol=1e-13)

    def test_training_loss_ignores_event_speed_and_time_zero(self):
        d,_=sample(100)
        x,t,_=kinematics(d)
        q=contrasts(x)
        f,_=prepare_features(d,1)
        y=np.einsum('nri,ni->nr',q,t)
        altered=t+np.arange(100)[:,None]+np.linspace(-.2,.4,100)[:,None]*x
        shifted=np.einsum('nri,ni->nr',q,altered)
        p=np.asarray(self.model['parameters'])
        original,gradient=likelihood(p,f,q,y,.1)
        changed,new_gradient=likelihood(p,f,q,shifted,.1)
        self.assertAlmostEqual(original,changed,places=10)
        np.testing.assert_allclose(gradient,new_gradient,atol=1e-10,rtol=0)

    def test_reconstruction_on_independent_synthetic_sample(self):
        d,truth=sample(4000,102)
        u,_,_=predict(d,self.model)
        _,_,linear=kinematics(d)
        self.assertLess(np.std(u-truth),.98*np.std(linear-truth))

    def test_fallback(self):
        d,_=sample(50)
        d['tof_bpatt'][0]=2444
        d['tof_edep'][1,2]=0
        u,valid=reconstruct(d,self.model)
        self.assertFalse(valid[0]);self.assertFalse(valid[1])
        self.assertEqual(u[0],1/d['tof_betah'][0])
        self.assertEqual(u[1],1/d['tof_betah'][1])

    def test_inference_across_batch_boundary(self):
        d,_=sample(100003, 92)
        expected=predict(d,self.model)[0]
        result,valid=reconstruct(d,self.model)
        self.assertTrue(np.all(valid))
        np.testing.assert_allclose(result,expected,atol=1e-13,rtol=0)

    def test_truth_never_enters_cache(self):
        self.assertFalse(any(k.startswith('mc') or k in ('mmass','mch','mevmom1') for k in BRANCHES))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'cache.npz'
            d={k:np.zeros(2) for k in ARRAY_KEYS}
            np.savez(str(p),**d)
            self.assertEqual(set(load_arrays(p)),set(ARRAY_KEYS))
            d['mcBeta']=np.ones(2)
            np.savez(str(p),**d)
            with self.assertRaises(ValueError):load_arrays(p)

    def test_model_serialization(self):
        d,_=sample(20)
        restored=json.loads(json.dumps(self.model))
        np.testing.assert_array_equal(predict(d,self.model)[0],predict(d,restored)[0])


if __name__=='__main__':unittest.main()
