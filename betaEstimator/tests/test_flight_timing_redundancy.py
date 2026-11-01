#!/usr/bin/env python3
"""Algebra, independent-end closure and observational-equivalence tests."""
import sys
import json
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import flight_timing_redundancy as m

class RedundancyTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.RandomState(823)

    def geometry(self,n):
        return np.array([-2.,-1.9,1.9,2.])[None,:]*(1+self.rng.uniform(0,.2,(n,1)))

    def test_projector(self):
        x=self.geometry(100);p,w=m.projector(x)
        np.testing.assert_allclose(p,p.transpose(0,2,1),atol=1.e-14)
        np.testing.assert_allclose(np.einsum('nij,njk->nik',p,p),p,atol=1.e-14)
        np.testing.assert_allclose(p.sum(axis=2),0,atol=1.e-14)
        np.testing.assert_allclose(np.einsum('nij,nj->ni',p,x),0,atol=1.e-14)
        np.testing.assert_allclose(np.sum(w*x,axis=1),1,atol=1.e-14)

    def test_sparse_schur_equals_dense(self):
        n=800;x=self.geometry(n);h=self.rng.normal(size=(n,4,2))
        ch=np.tile(np.array([0,10,20,30]),(n,1))+self.rng.randint(0,2,(n,4))
        t=self.rng.normal(size=(n,4));fit=m.fit_correction(x,t,h,ch,robust_iterations=0)
        b=np.zeros((n,4,80))
        for e in range(n):
            for l in range(4):
                b[e,l,2*ch[e,l]:2*ch[e,l]+2]=h[e,l]
        p,_=m.projector(x);bp=np.einsum('nij,njk->nik',p,b);pt=np.einsum('nij,nj->ni',p,t)
        fisher=np.einsum('nmi,nmj->ij',bp,bp)
        rhs=np.einsum('nmi,nm->i',bp,pt)
        np.testing.assert_allclose(fit['fisher'],fisher,atol=1.e-10)
        np.testing.assert_allclose(fit['theta'],np.linalg.pinv(fisher).dot(rhs),atol=1.e-10)

    def test_independent_end_closure(self):
        n=24000;x=self.geometry(n);u=1+self.rng.uniform(0,.4,n)
        sl=np.array([.10,.04,.08,.05]);sr=np.array([.04,.10,.05,.08])
        left=self.rng.normal(size=(n,4))*sl;right=self.rng.normal(size=(n,4))*sr
        h=(left-right)[:,:,None];t=x*u[:,None]+(left+right)/2
        ch=np.tile(np.array([0,10,20,30]),(n,1))
        fit=m.fit_correction(x[:16000],t[:16000],h[:16000],ch[:16000],robust_iterations=0)
        oracle=(sl*sl-sr*sr)/(2*(sl*sl+sr*sr))
        np.testing.assert_allclose(fit['theta'][[0,10,20,30]],oracle,atol=.025)
        reco=m.reconstruct(x[16000:],t[16000:],h[16000:],ch[16000:],fit['theta'])
        _,w=m.projector(x[16000:]);raw=np.sum(w*t[16000:],axis=1)-u[16000:]
        improved=reco['u']-u[16000:]
        self.assertLess(np.std(improved),.87*np.std(raw))
        self.assertLess(abs(np.mean(improved)),.0006)

    def test_full_readout_coherent_response(self):
        n=500;x=self.geometry(n);raw=50+self.rng.normal(size=(n,4,2,3))
        raw[::11,0,0,1]=0;raw[::13,2,1,2]=-9999
        h,good=m.raw_contrasts(raw);ch=np.tile(np.array([0,10,20,30]),(n,1))
        theta=self.rng.normal(scale=.1,size=120);t=x+self.rng.normal(scale=.05,size=(n,4))
        a=self.rng.normal(size=n);b=self.rng.uniform(-.2,.2,n);shift=a[:,None]+b[:,None]*x
        # Raw times are FT-LT; a physical later hit shifts all thresholds down.
        valid=np.isfinite(raw)&(np.abs(raw)<1000)&(raw!=0)
        changed=np.where(valid,raw-shift[:,:,None,None],raw)
        hh,gg=m.raw_contrasts(changed)
        np.testing.assert_allclose(h,hh,atol=2.e-14);np.testing.assert_array_equal(good,gg)
        r0=m.reconstruct(x,t,h,ch,theta);r1=m.reconstruct(x,t+shift,hh,ch,theta)
        np.testing.assert_allclose(r1['u']-r0['u'],b,atol=1.e-13)

    def test_equivalent_slope_worlds(self):
        n=1000;x=self.geometry(n);h=self.rng.normal(size=(n,4,1));u=1+self.rng.uniform(0,.1,n)
        noise=self.rng.normal(scale=.03,size=(n,4));delta=.025*h[:,0,0]
        time1=x*u[:,None]+noise
        time2=x*(u+delta)[:,None]+(noise-x*delta[:,None])
        np.testing.assert_allclose(time1,time2,atol=2.e-15)
        p,_=m.projector(x)
        np.testing.assert_allclose(np.einsum('nij,nj->ni',p,noise),np.einsum('nij,nj->ni',p,noise-x*delta[:,None]),atol=2.e-15)
        self.assertGreater(np.std(delta),.02)

    def test_invalid_readouts(self):
        raw=np.ones((3,4,2,3))*50
        raw[0,0,0,1]=0;raw[1,0,1,0]=-9999;raw[2,0,1,2]=np.nan
        h,ok=m.raw_contrasts(raw)
        self.assertFalse(ok[0,0,1]);self.assertFalse(ok[1,0,0]);self.assertFalse(ok[2,0,2])
        self.assertTrue(np.isfinite(h).all())
        with self.assertRaises(ValueError):m.projector(np.zeros((3,4)))

if __name__=='__main__':unittest.main(verbosity=2)
