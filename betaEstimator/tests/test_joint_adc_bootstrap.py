#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
from bootstrap_joint_adc_20261013 import rebase,prediction
import joint_adc_likelihood as joint
import test_joint_adc_likelihood as fixture

class JointBootstrapTests(unittest.TestCase):
    def setUp(self):
        f=fixture.JointADCTests();f.setUp();self.f=f

    def test_geometry_rebase_exact_away_from_clip(self):
        f=self.f;old=f.norm;new=dict(old);new['center']=(np.asarray(old['center'])+.003).tolist();new['scale']=(np.asarray(old['scale'])*1.04).tolist()
        theta=rebase(f.theta,old,new);raw=np.asarray(old['center'])[None,:]+np.linspace(-.1,.1,len(f.ff))[:,None]*np.asarray(old['scale'])[None,:]
        ff0=f.ff.copy();ff1=f.ff.copy();ff0[:,:,1:6]=((raw-np.asarray(old['center']))/np.asarray(old['scale']))[:,None,:]
        ff1[:,:,1:6]=((raw-np.asarray(new['center']))/np.asarray(new['scale']))[:,None,:]
        alpha=f.theta[39:].reshape(4,16);alpha1=theta[39:].reshape(4,16)
        np.testing.assert_allclose(np.sum(ff0*alpha[None,:,:],axis=2),np.sum(ff1*alpha1[None,:,:],axis=2),atol=1e-13,rtol=1e-13)
        np.testing.assert_array_equal(theta[:39],f.theta[:39])

    def test_prediction_keeps_ids_and_original_direction(self):
        f=self.f;n=len(f.choose);f.a['run']=np.repeat([100,101],n//2);f.a['event']=np.arange(n);f.g['moment']=np.zeros((n,4,2,9));theta=f.theta.copy();theta[:2]=0
        expected=joint.reconstruct(theta,f.g['x'],f.y,f.ff);actual=prediction(theta,f.a,f.g,f.adc,f.ff,f.choose,'t_native')
        np.testing.assert_array_equal(actual['run'],f.a['run'])
        np.testing.assert_array_equal(actual['event'],f.a['event'])
        for key in ['joint_u','time_u','adc_u']:np.testing.assert_allclose(actual[key],expected[key],atol=1e-14,rtol=1e-14)
        np.testing.assert_allclose(actual['joint_transport_u'],actual['joint_u'],atol=1e-14,rtol=1e-14)

    def test_normalization_counts_duplicate_sampled_runs(self):
        f=self.f;index=np.array([0,1,2,0,1,2,4,5]);features,norm=joint.response_features(f.a,index)
        pos=f.a['tk_pos'][:,4,:2]/50.;x=pos[:,0];y=pos[:,1];raw=np.column_stack([x,y,x*x,y*y,x*y])
        np.testing.assert_allclose(norm['center'],raw[index].mean(0),atol=1e-14,rtol=1e-14)
        np.testing.assert_allclose(norm['scale'],raw[index].std(0),atol=1e-14,rtol=1e-14)

if __name__=='__main__':unittest.main()
