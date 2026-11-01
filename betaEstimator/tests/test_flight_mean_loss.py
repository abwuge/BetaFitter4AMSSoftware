#!/usr/bin/env python3
"""Algebra and numerical limits of the coefficient-free mean-loss propagator."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import flight_mean_loss as m

class MeanLossTests(unittest.TestCase):
    def setUp(self):
        high=np.array([10.,5.,0.,-5.]);low=high-5.;length=np.full(4,5.)
        col=np.zeros((4,9));col[:,1]=.5
        self.args=(high,low,length,col,np.array([0,1,1,0]),np.array([0,1,-1,2]),np.array([1,-1,2,3]),0.,np.array([10.,5.,-5.,-10.]))
        self.x=np.array([-10.,-5.,5.,10.])/m.C
        nodes=np.array([-10.,-5.,0.,5.,10.])/m.C;middle=.5*(nodes[:-1]+nodes[1:])
        h=np.maximum(self.x[:,None]-middle,0)-np.maximum(-middle,0)-self.x[:,None]*(middle<=0)
        self.moment=np.zeros((4,2,9))
        for g in range(2):self.moment[:,g]=h[:,self.args[4]==g].dot(col[self.args[4]==g])

    def test_zero_loss(self):
        t,l=m.propagate_cells(1.2,12.,np.zeros(2),*self.args)
        np.testing.assert_allclose(t,1.2*self.x,rtol=0,atol=1e-14);np.testing.assert_array_equal(l,0)

    def test_charge_squared(self):
        np.testing.assert_allclose(m.stopping_power(.8,12.,6),4*m.stopping_power(.8,12.,3),rtol=1e-14)

    def test_small_loss_derivative(self):
        b=m.first_order(1.2,12.,self.moment);eps=1e-5
        for g in range(2):
            scale=np.zeros(2);scale[g]=eps;t,_=m.propagate_cells(1.2,12.,scale,*self.args,substeps=8)
            np.testing.assert_allclose((t-1.2*self.x)/eps,b[:,g],atol=2e-9,rtol=1e-5)

    def test_refinement(self):
        one,_=m.propagate_cells(1.2,12.,np.ones(2),*self.args,substeps=1)
        four,_=m.propagate_cells(1.2,12.,np.ones(2),*self.args,substeps=4)
        sixteen,_=m.propagate_cells(1.2,12.,np.ones(2),*self.args,substeps=16)
        self.assertLess(np.max(abs(four-sixteen)),np.max(abs(one-sixteen))/8)

    def test_unphysical_states_rejected(self):
        for u in [.9,1.,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):m.propagate_cells(u,12.,np.ones(2),*self.args)
        for scales in [np.array([-1.,1.]),np.array([np.nan,1.])]:
            with self.assertRaises(ValueError):m.propagate_cells(1.2,12.,scales,*self.args)

    def test_stopped_particle_flag(self):
        t,l=m.propagate_cells(1.2,12.,np.full(2,1e5),*self.args)
        self.assertTrue(np.isnan(t).all());self.assertTrue(np.isnan(l).all())

    def test_energy_derivative(self):
        energy=22.;mass=12.;eps=1e-4;p=np.sqrt(energy*energy-mass*mass)
        derivative=(m.inverse_beta(energy+eps,mass)-m.inverse_beta(energy-eps,mass))/(2*eps)
        self.assertAlmostEqual(derivative,-mass*mass/p**3,places=10)

if __name__=='__main__':unittest.main()
