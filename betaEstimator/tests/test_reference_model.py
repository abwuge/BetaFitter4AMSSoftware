"""Meaningful synthetic checks; no experiment files or MC fields are loaded.

The predictive test fixes its algorithm before exposing held-out latent speed:
training uses only raw TOF inputs and a statistically independent noisy
reference.  Latent speed appears only in the final, frozen-model score.
Run: OPENBLAS_NUM_THREADS=1 python -m unittest -v test_reference_model
"""
from __future__ import division

import copy
import json
import unittest

import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import reference_model as rm


def synthetic_sample(n, seed):
    rng = np.random.RandomState(seed)
    z = 8.0
    latent_u = 1.0 + rng.uniform(0.0, 0.45, n)
    pathscale = rng.uniform(0.92, 1.18, n)
    lengths = np.array([0.0, 45.0, 96.0, 141.0])[None, :] * pathscale[:, None]
    lengths += rng.uniform(-20.0, 20.0, n)[:, None]
    x = lengths / rm.C_CM_PER_NS
    xc = x - x.mean(axis=1)[:, None]
    span = np.ptp(x, axis=1)
    impact = rng.uniform(-35.0, 35.0, (n, 2))
    tilt = rng.normal(0.0, 0.075, (n, 2))
    position = impact[:, None, :] + (lengths - lengths[:, :1])[:, :, None] * tilt[:, None, :]
    direction = np.column_stack((tilt, np.ones(n)))
    bars = rng.randint(0, 8, (n, 4))
    tof = z * z * rng.lognormal(0.9, 0.25, (n, 4))
    tk = z * z * rng.lognormal(-0.5, 0.3, (n, 9))
    systematic = (0.009 * impact[:, 0] / 35.0
                  + 0.004 * impact[:, 1] / 35.0
                  + 0.006 * (np.log1p(tof[:, 0] / z ** 2) - 1.25)
                  + 0.0008 * (bars[:, 0] - 3.5))
    # Use exactly two timing residual directions, both perpendicular to 1,x.
    w = np.dot(xc / np.sqrt(np.sum(xc * xc, axis=1))[:, None], rm._H)
    h = w.copy()
    h[:, 0] += np.where(w[:, 0] >= 0.0, 1.0, -1.0)
    h /= np.sqrt(np.sum(h * h, axis=1))[:, None]
    hh = np.dot(h, rm._H.T)
    q0 = rm._H[:, 1][None, :] - 2.0 * hh * h[:, 1:2]
    q1 = rm._H[:, 2][None, :] - 2.0 * hh * h[:, 2:3]
    r0, r1 = rng.normal(0.0, 0.009, (2, n))
    slope_noise = ((0.55 + 0.2 * impact[:, 0] / 35.0) * r0
                   - 0.35 * r1 + rng.normal(0.0, 0.002, n))
    t0 = rng.normal(12.0, 1.0, n)
    t = (t0[:, None] + (latent_u + systematic + slope_noise)[:, None] * x
         + (span * r0)[:, None] * q0 + (span * r1)[:, None] * q1)
    d = {
        "tof_tl": t, "tof_leng": lengths, "position": position,
        "direction": direction, "tof_barid": bars, "tof_edep": tof,
        "tk_edep": tk, "tof_bpatt": np.full(n, 4444),
        "charge": np.full(n, z),
    }
    # Reference uncertainty is independent of TOF errors and input features.
    reference_u = latent_u + rng.normal(0.0, 0.004, n)
    outlier = rng.uniform(size=n) < 0.015
    reference_u[outlier] += rng.normal(0.0, 0.1, np.sum(outlier))
    return d, reference_u, latent_u


class ForbiddenInput(object):
    def __array__(self, *args, **kwargs):
        raise AssertionError("a forbidden beta/reference/truth input was read")


class ReferenceModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.small, cls.reference, unused_latent = synthetic_sample(500, 11)
        cls.target = rm.baseline(cls.small) - cls.reference
        cls.model = rm.fit(cls.small, 8, cls.target, 0.1,
                           mode="full", batch_size=113, max_iter=4)

    def test_affine_time_injection_equivariance(self):
        changed = dict(self.small)
        n = len(self.target)
        rng = np.random.RandomState(31)
        delta = rng.uniform(-0.75, 1.25, n)
        t0 = rng.uniform(-100.0, 100.0, n)
        changed["tof_tl"] = (self.small["tof_tl"] + t0[:, None]
                              + delta[:, None] * self.small["tof_leng"] / rm.C_CM_PER_NS)
        u, correction = rm.predict(self.small, self.model)
        shifted, shifted_correction = rm.predict(changed, self.model)
        np.testing.assert_allclose(shifted, u + delta, rtol=0, atol=5e-12)
        np.testing.assert_allclose(shifted_correction, correction, rtol=0, atol=5e-12)
        # Negative reconstructed u is also retained; no physical-range clamp.
        large = dict(self.small)
        large["tof_tl"] = self.small["tof_tl"] - 3.0 * self.small["tof_leng"] / rm.C_CM_PER_NS
        large_u, unused = rm.predict(large, self.model)
        self.assertTrue(np.all(large_u < 0.0))

    def test_forbidden_fields_are_never_read(self):
        poisoned = dict(self.small)
        for key in ("tof_betah", "rich_beta", "rigidity", "beta_true", "mc_beta", "target"):
            poisoned[key] = ForbiddenInput()
        f, names = rm.features(self.small, self.model)
        fp, namesp = rm.features(poisoned, self.model)
        np.testing.assert_array_equal(f, fp)
        self.assertEqual(names, namesp)
        out = rm.predict(poisoned, self.model)
        np.testing.assert_array_equal(out[0], rm.predict(self.small, self.model)[0])
        refit = rm.fit(poisoned, 8, self.target[::-1], 0.1,
                       mode="full", batch_size=113, max_iter=1)
        self.assertEqual(refit["feature_state"], self.model["feature_state"])
        self.assertFalse(np.allclose(refit["coef"], self.model["coef"]))

    def test_geometry_ablation_has_no_timing_features(self):
        model = rm.fit(self.small, 8, self.target, 0.1, mode="geometry", max_iter=1)
        changed = dict(self.small)
        rng = np.random.RandomState(83)
        changed["tof_tl"] = rng.normal(0.0, 100.0, self.small["tof_tl"].shape)
        f, names = rm.features(self.small, model)
        fp, unused = rm.features(changed, model)
        np.testing.assert_array_equal(f, fp)
        self.assertFalse(any("residual" in name for name in names))
        self.assertLess(len(names), 200)
        self.assertGreater(len(self.model["coef"]), 200)
        self.assertLess(len(self.model["coef"]), 400)

    def test_json_roundtrip_and_batch_invariance(self):
        loaded = json.loads(json.dumps(self.model, allow_nan=False))
        u, c = rm.predict(self.small, self.model, batch_size=500)
        ul, cl = rm.predict(self.small, loaded, batch_size=37)
        np.testing.assert_allclose(ul, u, rtol=0, atol=1e-13)
        np.testing.assert_allclose(cl, c, rtol=0, atol=1e-13)
        self.assertEqual(self.model["training"]["n_fit"], len(self.target))

    def test_partial_layer_and_model_fallbacks(self):
        d = {key: np.array(value[:5], copy=True) for key, value in self.small.items()}
        d["tof_tl"][0, 2:] = np.nan
        d["tof_tl"][1, 1:] = np.nan
        d["tof_leng"][2] = 10.0
        d["tof_bpatt"][3] = 3444
        u, correction = rm.predict(d, self.model)
        raw = rm.baseline(d)
        self.assertTrue(np.isfinite(u[0]))
        self.assertTrue(np.isnan(u[1]))
        self.assertTrue(np.isnan(u[2]))
        self.assertTrue(np.isfinite(u[3]))
        np.testing.assert_array_equal(correction[:4], np.zeros(4))
        np.testing.assert_allclose(u[:4], raw[:4], equal_nan=True)
        fallback, details = rm.reconstruct(d, {"version": -1}, return_details=True)
        np.testing.assert_allclose(fallback, raw, equal_nan=True)
        self.assertIn("unsupported", details["fallback_reason"])
        no_model, details = rm.reconstruct(d, return_details=True)
        np.testing.assert_allclose(no_model, raw, equal_nan=True)
        self.assertEqual(details["fallback_reason"], "no_model")

        d["tof_betah"] = np.array([0.8, 0.75, 0.9, 0.95, 0.85])
        producer, details = rm.reconstruct(d, self.model, return_details=True)
        np.testing.assert_array_equal(details["producer_fallback"],
                                      [True, True, True, True, False])
        np.testing.assert_array_equal(producer[:4], 1.0 / d["tof_betah"][:4])
        self.assertEqual(producer[4], u[4])
        # Neither numerical digits nor binary bits are interpreted as masks.
        pattern_data = rm.select_data(self.small, np.arange(4))
        pattern_data["tof_bpatt"] = np.array([4444, 3444, 15, 4443])
        kin = rm.kinematics(pattern_data)
        np.testing.assert_array_equal(kin["valid_model"], [True, False, False, False])
        self.assertTrue(np.all(kin["layer_valid"]))

    def test_event_charge_model_and_unseen_bars(self):
        event_model = rm.fit(self.small, self.small["charge"], self.target, 0.1,
                             mode="geometry", max_iter=1)
        d = dict(self.small)
        d["tof_barid"] = np.full_like(d["tof_barid"], 999)
        u, c = rm.predict(d, event_model)
        self.assertTrue(np.all(np.isfinite(u)))
        self.assertTrue(np.all(np.isfinite(c)))
        del d["charge"]
        with self.assertRaises(ValueError):
            rm.predict(d, event_model)
        out, details = rm.reconstruct(d, event_model, return_details=True)
        np.testing.assert_array_equal(out, rm.baseline(d))
        self.assertIn("charge", details["fallback_reason"])

    def test_independent_reference_improves_frozen_latent_score(self):
        train, reference, unused_train_latent = synthetic_sample(12000, 415)
        test, unused_test_reference, final_latent = synthetic_sample(6000, 892)
        target = rm.baseline(train) - reference
        # Prespecified models: no alpha/mode choice is made using latent truth.
        geometry = rm.fit(train, 8, target, 0.1, mode="geometry",
                          batch_size=2048, max_iter=4)
        full = rm.fit(train, 8, target, 0.1, mode="full",
                      batch_size=2048, max_iter=4)
        frozen = [json.loads(json.dumps(model)) for model in (geometry, full)]
        # Only here does the test inspect held-out latent speed.
        rms = []
        for u in (rm.baseline(test), rm.predict(test, frozen[0])[0],
                  rm.predict(test, frozen[1])[0]):
            rms.append(float(np.sqrt(np.mean((u - final_latent) ** 2))))
        print("\nSynthetic frozen latent RMSE: baseline={:.6f}, geometry={:.6f}, full={:.6f}".format(*rms))
        self.assertLess(rms[1], 0.90 * rms[0])
        self.assertLess(rms[2], 0.60 * rms[0])
        self.assertLess(rms[2], 0.80 * rms[1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
