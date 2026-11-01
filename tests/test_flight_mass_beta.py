"""Numerical and supervision-boundary checks for flight mass training."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import train_flight_mass_beta as tr


class FlightMassTrainingTests(unittest.TestCase):
    def test_fraction_and_width_gradient(self):
        rng = np.random.RandomState(20261009)
        rigidity = rng.uniform(2.5, 6., 400)
        data = {'R': rigidity, 'bin': np.searchsorted(tr.rc.R_EDGES, rigidity) - 1,
                'u': np.sqrt(1 + (11.17793 / (6 * rigidity))**2) + rng.normal(0, .02, 400),
                'relR': rng.uniform(.02, .1, 400)}
        density = tr.RawMixture(data)
        p = np.r_[np.log([.01, .015, .02, .03]), [.02, .1, .4, .8]]
        _, gradient = density.value_gradient(p)
        for j in range(len(p)):
            hi, lo = p.copy(), p.copy()
            hi[j] += 1e-6
            lo[j] -= 1e-6
            numerical = (density.value_gradient(hi)[0] - density.value_gradient(lo)[0]) / 2e-6
            self.assertAlmostEqual(gradient[j], numerical, delta=1e-6 * max(1., abs(numerical)))

    def test_training_rejects_answer_bearing_cache(self):
        for answer in ('rich_beta', 'u_true', 'isotope_label'):
            with tempfile.TemporaryDirectory(prefix='flight-beta-boundary-') as folder:
                folder = Path(folder)
                arrays = {key: np.zeros(1) for key in tr.TRAIN_KEYS}
                arrays[answer] = np.zeros(1)
                np.savez(str(folder / 'training.npz'), **arrays)
                with patch.object(tr, 'CACHE', folder), patch.object(tr, 'DEST', folder / 'results'):
                    with self.assertRaisesRegex(ValueError, 'Unexpected training fields'):
                        tr.train()
                    self.assertFalse((tr.DEST / 'freeze.json').exists())

    def test_extractor_rejects_holdout_source_before_open(self):
        with patch.object(tr.uproot, 'open', side_effect=AssertionError('must not open')):
            with self.assertRaisesRegex(ValueError, 'validation source'):
                tr.extract_part((0, 'validation', '/unused/holdout.root'))


if __name__ == '__main__':
    unittest.main()
