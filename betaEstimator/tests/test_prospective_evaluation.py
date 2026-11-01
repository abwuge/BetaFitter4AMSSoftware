"""Cohort-preservation and fixed metric checks without opening real references."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from evaluate_prospective_flight import prediction_coverage
from evaluate_additional_mc_20261013 import matrix_metrics, matrix_difference
from rich_transport import rich_quality


class ProspectiveEvaluationTest(unittest.TestCase):
    def test_unrestricted_predictions_are_preserved(self):
        predictions = {'candidate': np.array([.95, 1.2, -.8, 1.])}
        result = prediction_coverage(predictions, 'candidate', np.ones(4, bool))
        self.assertTrue(result['complete'])
        self.assertEqual(result['events'], 4)
        self.assertEqual(result['invalid_events'], 0)

    def test_failure_is_not_a_cohort_filter(self):
        predictions = {'candidate': np.array([1.1, np.nan, 1.2, 0.]),
                       'candidate_valid': np.array([True, False, True, False])}
        mask = np.array([True, True, True, False])
        result = prediction_coverage(predictions, 'candidate', mask)
        self.assertFalse(result['complete'])
        self.assertEqual(result['events'], 3)
        self.assertEqual(result['invalid_events'], 1)
        np.testing.assert_array_equal(mask, [True, True, True, False])
        predictions['candidate'] = np.array([1.1, 1.3, 1.2, np.nan])
        predictions['candidate_valid'] = np.array([True, True, True, False])
        self.assertTrue(prediction_coverage(predictions, 'candidate', mask)['complete'])

    def test_metric_scaling_is_exact(self):
        error = np.linspace(-.05, .05, 201)
        values = np.column_stack([error, .8 * error])
        metrics = matrix_metrics(values)
        change = matrix_difference(metrics[0], metrics[1])
        self.assertAlmostEqual(change[0], .2, places=13)
        self.assertAlmostEqual(change[1], .2, places=13)
        self.assertAlmostEqual(change[2], 0., places=13)
        self.assertLessEqual(change[3], 0.)

    def test_existing_reference_quality_and_npe_columns(self):
        n = 3
        a = dict(rich_beta=np.tile([.88, .88, .88], (n, 1)), rich_NaF=np.ones(n, bool),
            rich_pos=np.zeros((n, 3)), rich_npe=np.tile([50., 2., 100.], (n, 1)),
            rich_BetaConsistency=np.full(n, .001), rich_PMTChargeConsistency=np.ones(n),
            irich=np.arange(n), itrtrack=np.zeros(n, int), rich_itrtrack=np.zeros(n, int),
            rich_good=np.ones(n, bool), rich_clean=np.ones(n, bool), rich_pb=np.full(n, .2),
            rich_pmt=np.full(n, 5), rich_used=np.full(n, 8), rich_ebeta1=np.full(n, .001), rich_tile=np.zeros(n, int))
        np.testing.assert_array_equal(rich_quality(a, strict=True), [True, True, True])
        a['rich_npe'][1, 0] = 30.  # Ring accounts for less than 40 percent of the event light.
        a['rich_npe'][2, 1] = .5   # Expected ring photoelectron count is below one.
        np.testing.assert_array_equal(rich_quality(a, strict=True), [True, False, False])


if __name__ == '__main__':
    unittest.main()
