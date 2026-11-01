"""Physics and frozen-model interface checks for center-to-RICH propagation."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import rich_transport as rt


class TransportTests(unittest.TestCase):
    def geometry(self, value=.02, n=3):
        return {'column':np.full((n,7,9),value,np.float32),'segment':np.arange(7),
                'valid':np.ones(n,bool)}

    def test_zero_material_identity(self):
        u=np.array([1.02,1.10,1.25])
        p=rt.propagate(u,np.array([5.,4.,3.]),self.geometry(0))
        np.testing.assert_allclose(p['beta_planes'],np.repeat((1/u)[:,None],8,axis=1),rtol=1e-13)
        np.testing.assert_array_equal(p['loss_GeV'],0)

    def test_energy_conservation_and_monotonic_speed(self):
        p=rt.propagate(np.array([1.02,1.1,1.25]),np.array([5.,4.,3.]),self.geometry())
        np.testing.assert_allclose(p['initial_energy_GeV']-p['final_energy_GeV'],p['loss_GeV'].sum(axis=1),atol=3e-13)
        self.assertTrue((p['loss_GeV']>=0).all())
        self.assertTrue((np.diff(p['beta_planes'],axis=1)<=1e-14).all())
        energy=p['mass_GeV']/np.sqrt(1-p['beta_planes'][:,-1]**2)
        np.testing.assert_allclose(energy,p['final_energy_GeV'],rtol=1e-13)

    def test_unphysical_input_not_snapped(self):
        p=rt.propagate(np.array([.99,1.,1.1]),np.array([4.,4.,4.]),self.geometry())
        np.testing.assert_array_equal(p['valid'],[False,False,True])
        self.assertTrue(np.isnan(p['beta_planes'][:2]).all())

    def test_stopping_charge_scaling(self):
        b=np.array([.8,.9]);m=np.array([12.,12.])
        np.testing.assert_allclose(rt.stopping_power(b,m,6),36*rt.stopping_power(b,m,1))

    def test_straight_track_interpolation(self):
        z=np.array([158.92,53.06,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882])
        pos=np.zeros((3,9,3));pos[:,:,2]=z;pos[:,:,0]=2+.2*z;pos[:,:,1]=3-.1*z
        direction=np.broadcast_to(np.array([-.2,.1,-1.]),pos.shape).copy()
        xy,metric=rt.track_at_z(pos,direction,np.array([-1.,-20.,-75.]))
        np.testing.assert_allclose(xy,np.column_stack([2+.2*np.array([-1.,-20.,-75.]),3-.1*np.array([-1.,-20.,-75.])]),atol=1e-12)
        np.testing.assert_allclose(metric,np.sqrt(1.05))

    def test_feature_state_replay_and_forbidden_keys(self):
        train,valid=rt.er.load_cache(rt.BASE/'controlled_cache.npz')
        d={k:valid[k][:20] for k in ['u','R','bin','continuous','bars']}
        model=rt.CenterModel();p=model.predict(d)
        d.update({'u_true':np.full(20,99.),'rich_beta':np.full(20,-1.),'isotope':np.arange(20)})
        np.testing.assert_array_equal(p,model.predict(d))
        self.assertTrue(np.isfinite(p).all())
        np.testing.assert_array_equal(p,rt.CenterModel().predict(d))

    def test_domain_rejection(self):
        train,valid=rt.er.load_cache(rt.BASE/'controlled_cache.npz')
        d={k:valid[k][:2] for k in ['u','R','bin','continuous','bars']}
        d['bin']=np.array([-1,4])
        with self.assertRaises(ValueError):rt.CenterModel().predict(d)


if __name__=='__main__':unittest.main()
