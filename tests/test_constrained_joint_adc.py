"""Response and minimum-variance checks for ADC-assisted free timing slopes."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from constrained_joint_adc import constrained_weights


class ConstrainedADCTest(unittest.TestCase):
    def setUp(self):
        rng = np.random.RandomState(815)
        a = rng.normal(size=(7, 7))
        self.cov = a.dot(a.T) + np.eye(7)
        self.q = np.linalg.inv(self.cov)
        self.t = rng.normal(size=(12, 3))
        self.adc = np.tile([.45, .50, .61, .70], (12, 1))
        self.rng = rng

    def test_constraints_and_injections(self):
        for gain in [False, True]:
            out = constrained_weights(self.q, self.t, self.adc, gain)
            w = out['weights']; c = out['constraint_matrix']
            target = np.zeros((len(w), c.shape[-1])); target[:, 0] = 1.
            np.testing.assert_allclose(np.einsum('ni,nik->nk', w, c), target, atol=1e-12)
            signals = self.rng.normal(size=(len(w), c.shape[-1]))
            injected = np.einsum('nik,nk->ni', c, signals)
            np.testing.assert_allclose(np.sum(w * injected, axis=1), signals[:, 0], atol=1e-12)

    def test_optimality_identity(self):
        for gain in [False, True]:
            out = constrained_weights(self.q, self.t, self.adc, gain)
            w = out['weights']; c = out['constraint_matrix']
            h = self.rng.normal(size=w.shape)
            ctc = np.einsum('nik,nil->nkl', c, c)
            projection = np.linalg.solve(ctc, np.einsum('nik,ni->nk', c, h))
            h -= np.einsum('nik,nk->ni', c, projection)
            np.testing.assert_allclose(np.einsum('nik,ni->nk', c, h), 0, atol=1e-12)
            delta = np.einsum('ni,ij,nj->n', w + h, self.cov, w + h) - out['scale_variance']
            expected = np.einsum('ni,ij,nj->n', h, self.cov, h)
            np.testing.assert_allclose(delta, expected, rtol=1e-12, atol=1e-12)
            self.assertTrue(np.all(delta >= 0))

    def test_uncorrelated_adc_reduces_to_time_gls(self):
        cov = self.cov.copy(); cov[:3, 3:] = 0; cov[3:, :3] = 0
        q = np.linalg.inv(cov)
        wt = self.t.dot(np.linalg.inv(cov[:3, :3]))
        wt /= np.sum(wt * self.t, axis=1)[:, None]
        for gain in [False, True]:
            w = constrained_weights(q, self.t, self.adc, gain)['weights']
            np.testing.assert_allclose(w[:, :3], wt, atol=1e-12)
            np.testing.assert_allclose(w[:, 3:], 0, atol=1e-12)

    def test_variance_ordering(self):
        speed = constrained_weights(self.q, self.t, self.adc)
        gain = constrained_weights(self.q, self.t, self.adc, True)
        qt = np.linalg.inv(self.cov[:3, :3])
        time_variance = 1 / np.einsum('ni,ij,nj->n', self.t, qt, self.t)
        direction = np.column_stack([self.t, self.adc])
        joint_variance = 1 / np.einsum('ni,ij,nj->n', direction, self.q, direction)
        self.assertTrue(np.all(joint_variance <= speed['scale_variance'] + 1e-12))
        self.assertTrue(np.all(speed['scale_variance'] <= gain['scale_variance'] + 1e-12))
        self.assertTrue(np.all(gain['scale_variance'] <= time_variance + 1e-12))

    def test_domain_checks(self):
        with self.assertRaises(ValueError):
            constrained_weights(np.eye(6), self.t, self.adc)
        bad = self.q.copy(); bad[0, 0] = -10
        with self.assertRaises(ValueError):
            constrained_weights(bad, self.t, self.adc)
        badt = self.t.copy(); badt[0, 0] = np.nan
        with self.assertRaises(ValueError):
            constrained_weights(self.q, badt, self.adc)


if __name__ == '__main__':
    unittest.main()
