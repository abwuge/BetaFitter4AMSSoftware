#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import numpy as np
from evaluate_additional_mc_20261013 import matrix_metrics,matrix_difference,METRICS,DELTA
from evaluate_conditional_mc import metrics,pair_delta
from evaluate_amplitude_noise_20261013 import join

class AdditionalMetricsTests(unittest.TestCase):
    def test_all_metrics_equal_scalar(self):
        rng=np.random.RandomState(235);e=rng.standard_t(5,(671,21))*.014+rng.normal(0,.001,(1,21))
        vector=matrix_metrics(e)
        for i in range(e.shape[1]):
            expected=metrics(e[:,i]);np.testing.assert_allclose(vector[i],[expected[k] for k in METRICS],atol=1e-15,rtol=1e-13)
        vector_delta=matrix_difference(vector[0:1],vector)
        for i in range(e.shape[1]):
            expected=pair_delta(e[:,0],e[:,i]);np.testing.assert_allclose(vector_delta[i],[expected[k] for k in DELTA],atol=1e-14,rtol=1e-12)

    def test_nested_bootstrap_broadcast(self):
        rng=np.random.RandomState(84);e=rng.normal(0,.02,(219,5));draws=[rng.randint(0,len(e),len(e)) for _ in range(9)]
        mm=np.array([matrix_metrics(e[j]) for j in draws]);dd=matrix_difference(mm[:,1:2],mm)
        for b,index in enumerate(draws):
            for i in range(5):
                expected=pair_delta(e[index,1],e[index,i]);np.testing.assert_allclose(dd[b,i],[expected[k] for k in DELTA],atol=1e-14,rtol=1e-12)

    def test_composite_identity_exact_join_and_failures(self):
        left=dict(run=np.array([10,11,10],np.int64),event=np.array([100,100,200],np.int64))
        order=np.array([2,0,1]);right={k:v[order] for k,v in left.items()}
        np.testing.assert_array_equal(join(left,right),[1,2,0])
        bad={k:v.copy() for k,v in right.items()};bad['event'][0]=999
        with self.assertRaises(ValueError):join(left,bad)
        duplicate={k:np.r_[v,v[:1]] for k,v in right.items()}
        with self.assertRaises(ValueError):join(left,duplicate)

    def test_nonfinite_errors_never_silently_filtered(self):
        for value in [np.nan,np.inf,-np.inf]:
            e=np.zeros((8,2));e[1,1]=value
            with self.assertRaises(ValueError):matrix_metrics(e)

if __name__=='__main__':unittest.main()
