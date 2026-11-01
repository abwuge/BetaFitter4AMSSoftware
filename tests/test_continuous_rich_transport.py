"""Measurement extension and physical one-way mean-loss transport checks."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import continuous_rich_transport as transfer
import rich_transport as inherited


def geometry(n, carbon=.02):
    # Seven material cells are assigned to the existing seven output segments.
    column = np.zeros((n, 7, 9)); column[:, :, 1] = carbon
    return dict(column=column, segment=np.arange(7), valid=np.ones(n, bool),
                target_z=np.full(n, -75.552), NaF_predicted=np.ones(n, bool))


class ContinuousRichTest(unittest.TestCase):
    def test_zero_loss_exact_and_unrestricted_extension(self):
        u = np.array([-.8, .85, 1., 1.01, 1.2, 25.])
        R = np.full(len(u), 4.)
        result = transfer.propagate(u, R, geometry(len(u)), material_scale=0.)
        np.testing.assert_array_equal(result['u'], u)
        self.assertTrue(result['valid'].all())
        result = transfer.propagate(u[:5], R[:5], geometry(5))
        np.testing.assert_array_equal(result['u'][:3], u[:3])
        self.assertGreater(result['beta'][1], 1.)
        self.assertTrue(result['valid'].all())
        self.assertTrue(np.all(result['u'][3:] > u[3:5]))

    def test_inherited_physical_result_exact(self):
        u = np.array([1.03, 1.1, 1.25]); R = np.array([2.6, 3.5, 5.8]); g = geometry(3)
        old = inherited.propagate(u, R, g)
        new = transfer.propagate(u, R, g)
        np.testing.assert_array_equal(new['beta'], old['beta_planes'][:, -1])
        np.testing.assert_array_equal(new['loss_GeV'], old['loss_GeV'])
        np.testing.assert_array_equal(new['valid'], old['valid'])

    def test_subdivision_and_positive_loss(self):
        u = np.array([1.03, 1.1, 1.25]); R = np.array([2.6, 3.5, 5.8]); g = geometry(3, carbon=.1)
        one = transfer.propagate(u, R, g)
        two = transfer.propagate(u, R, transfer.subdivide(g, 2))
        four = transfer.propagate(u, R, transfer.subdivide(g, 4))
        self.assertTrue(one['valid'].all() and two['valid'].all() and four['valid'].all())
        self.assertTrue(np.all(one['loss_GeV'] > 0))
        self.assertLess(np.max(abs(four['u'] - two['u'])), np.max(abs(two['u'] - one['u'])))
        self.assertLess(np.max(abs(four['u'] - one['u'])), 1e-7)

    def test_limit_at_light_speed(self):
        u = 1 + np.array([1e-4, 1e-6, 1e-8]); R = np.full(3, 4.)
        result = transfer.propagate(u, R, geometry(3))
        shift = result['u'] - u
        self.assertTrue(result['valid'].all())
        self.assertTrue(np.all(shift > 0))
        self.assertTrue(np.all(np.diff(shift) < 0))
        self.assertLess(shift[-1], 1e-8)

    def test_invalid_states_are_flagged_not_removed(self):
        u = np.array([1.2, 25., .95]); R = np.array([4., 4., 4.])
        result = transfer.propagate(u, R, geometry(3))
        self.assertEqual(len(result['u']), 3)
        np.testing.assert_array_equal(result['valid'], [True, False, True])
        self.assertTrue(np.isnan(result['u'][1]))
        for bad in [np.array([0., 1., 2.]), np.array([np.nan, 1., 2.])]:
            with self.assertRaises(ValueError):
                transfer.propagate(bad, R, geometry(3))
        with self.assertRaises(ValueError):
            transfer.propagate(u, -R, geometry(3))
        with self.assertRaises(ValueError):
            transfer.propagate(u, R, geometry(3), material_scale=-1)


if __name__ == '__main__':
    unittest.main()
