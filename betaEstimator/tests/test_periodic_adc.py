"""Period response features and padding, using synthetic measurements only."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import periodic_joint_adc as period
import joint_adc_likelihood as joint


class PeriodicADC(unittest.TestCase):
    def data(self):
        rng = np.random.RandomState(417)
        n = 120
        return dict(run=np.linspace(period.EDGES[0]+100, period.EDGES[-1]-100, n).astype(np.int64),
                    tk_pos=rng.normal(size=(n, 9, 3))*15,
                    bar=rng.randint(0, 10, size=(n, 4)))

    def test_step_partition_and_edges(self):
        run = period.EDGES[:-1] + 10
        actual = period.period_basis(run, 'flight', 'step')
        np.testing.assert_array_equal(actual, np.eye(5)[:, 1:])
        self.assertEqual(period.period_basis(run, 'flight', 'static').shape, (5, 0))

    def test_spline_is_continuous(self):
        middle = period.KNOTS[1:-1]
        left = period.period_basis(middle-1e-3, 'flight', 'spline')
        right = period.period_basis(middle+1e-3, 'flight', 'spline')
        self.assertLess(float(np.max(abs(left-right))), 1e-9)
        np.testing.assert_allclose(period.period_basis(period.KNOTS, 'flight', 'spline'), np.eye(5)[:,1:])

    def test_mc_run_hash_is_order_independent(self):
        run = np.array([12, 55, 12, 200, 55, 91])
        coordinate = period.period_coordinates(run, 'mc')
        np.testing.assert_array_equal(coordinate[::-1], period.period_coordinates(run[::-1], 'mc'))
        self.assertEqual(coordinate[0], coordinate[2])
        self.assertEqual(coordinate[1], coordinate[4])

    def test_frozen_normalization_and_reordered_prediction(self):
        data = self.data(); train = np.arange(120) % 5 != 0
        for variant in period.VARIANTS:
            value, norm = period.features(data, train, 'flight', variant)
            permutation = np.arange(120)[::-1]
            changed = {key: val[permutation] for key, val in data.items()}
            replay, same = period.features(changed, np.ones(120, bool), 'flight', variant, norm)
            np.testing.assert_allclose(value[permutation], replay)
            self.assertEqual(norm, same)

    def test_heldout_features_do_not_retrain_normalization(self):
        data = self.data(); train = np.arange(120) % 5 != 0
        _, before = period.features(data, train, 'flight', 'spline')
        changed = {key: val.copy() for key, val in data.items()}
        changed['run'][~train] = int(period.KNOTS[0])
        changed['tk_pos'][~train] *= 10
        _, after = period.features(changed, train, 'flight', 'spline')
        self.assertEqual(before, after)

    def test_zero_period_coefficients_recover_original_response(self):
        data = self.data(); train = np.arange(120) % 5 != 0
        old_features, _ = joint.response_features(data, train)
        new_features, _ = period.features(data, train, 'flight', 'spline')
        theta = joint.parameters_initial(old_features, np.full((120,4,2),100.), train)
        theta[39:] += np.random.RandomState(913).normal(0, .01, len(theta)-39)
        padded = period.lift(theta, new_features.shape[2])
        old_mean = np.einsum('nlf,lf->nl', old_features, theta[39:].reshape(4,-1))
        new_mean = np.einsum('nlf,lf->nl', new_features, padded[39:].reshape(4,-1))
        np.testing.assert_allclose(old_mean, new_mean)
        np.testing.assert_array_equal(theta[:39], padded[:39])

    def test_invalid_domains_and_norms_are_rejected(self):
        with self.assertRaises(ValueError):
            period.period_coordinates(np.array([period.EDGES[0]-1]), 'flight')
        with self.assertRaises(ValueError):
            period.period_basis(np.array([period.KNOTS[0]]), 'flight', 'unknown')
        data = self.data(); train = np.ones(120,bool)
        _, norm = period.features(data, train, 'flight', 'step')
        with self.assertRaises(ValueError):
            period.features(data, train, 'flight', 'spline', norm)
        with self.assertRaises(ValueError):
            period.lift(np.zeros(103), 15)


if __name__ == '__main__':
    unittest.main()
