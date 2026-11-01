"""Time exclusion and supervision-boundary checks for adjacent-period fits."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import time_local_beta as t


class LocalTimeTests(unittest.TestCase):
    def test_window_boundaries_and_guard(self):
        low = 1500000000
        high = low+t.PERIOD
        width = 180*t.DAY
        run = np.array([low-t.GUARD-width-1, low-t.GUARD-width,
                        low-t.GUARD-1, low-t.GUARD, low, high-1,
                        high, high+t.GUARD-1, high+t.GUARD,
                        high+t.GUARD+width-1, high+t.GUARD+width])
        self.assertEqual(t.training_mask(run,low,high,180).tolist(),
                         [False,True,True,False,False,False,False,False,True,True,False])

    def test_training_rejects_reference_fields_before_fit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'folds').mkdir()
            run = np.array([1500000000])
            np.savez(str(root/'training.npz'), run=run, rich_beta=np.ones(1))
            t.w.ev.dump(root/'plan.json',dict(dependencies={},folds=[dict(fold=0)],training_sha256='irrelevant'))
            with patch.object(t,'OUT',root), patch.object(t.tr,'fit_raw',side_effect=AssertionError('fit reached')):
                with self.assertRaisesRegex(RuntimeError,'Unexpected training fields'):
                    t.fit_one(0)

    def test_all_models_must_be_frozen_before_reference_open(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'folds').mkdir()
            t.w.ev.dump(root/'plan.json',dict(folds=[dict(fold=0)]))
            with patch.object(t,'OUT',root), patch.object(t,'load',side_effect=AssertionError('reference opened')):
                with self.assertRaises(FileNotFoundError):
                    t.evaluate()


if __name__ == '__main__': unittest.main()
