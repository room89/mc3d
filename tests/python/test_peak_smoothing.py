import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
try:
    import numpy as np
    from experiment_peak_smoothing import METHODS, evaluate, summarize, training_profile
    from experiment_template_training import evaluate as evaluate_training
    from experiment_crossfit_gradient import fixture
    from shock_alignment import GRID
except ImportError:
    np = None


@unittest.skipIf(np is None, 'Requires numpy')
class PeakSmoothingTest(unittest.TestCase):
    def data(self):
        return np.tile(1/(1+np.exp(-GRID)), (3, 2, 1))

    def test_raw_preserves_profile(self):
        y = np.sin(GRID)
        np.testing.assert_array_equal(training_profile(y, 'raw'), y)

    def test_monotone_and_constant(self):
        y = np.sin(GRID)
        self.assertTrue(np.all(np.diff(training_profile(y, 'monotone')) >= 0))
        for method in METHODS:
            np.testing.assert_allclose(training_profile(np.ones_like(GRID), method), 1.)

    def test_heldout_does_not_select_own_peak(self):
        data = self.data()
        before = evaluate(data, 'logistic', 4.)
        data[0] += .02*np.exp(-((GRID-1)/.2)**2)
        after = evaluate(data, 'logistic', 4.)
        for method in METHODS:
            self.assertEqual(before[method]['peak_errors'][0], after[method]['peak_errors'][0])

    def test_clean_symmetric_profile(self):
        result = evaluate(self.data(), 'logistic', 4.)
        for method in METHODS:
            self.assertLess(abs(result[method]['normalized_slope']-1), .001)
            self.assertLess(max(abs(x) for x in result[method]['peak_errors']), .051)
        summary = summarize([result, result])
        for value in summary['paired_thickness_contrasts'].values():
            self.assertAlmostEqual(value['mean'], 0.)
            self.assertAlmostEqual(value['mc_se'], 0.)

    def test_training_control_matches_fixed_front_experiment(self):
        trials, oracle = fixture('gompertz', 3600000, .8, .2)
        first = evaluate(oracle, 'gompertz', np.e)
        second = evaluate_training(trials, oracle)
        self.assertEqual(first['smooth_monotone'], second['exact_training'])

    def test_training_template_excludes_heldout(self):
        trials, oracle = fixture('gompertz', 3600000, 0., 0., noise=0.)
        before = evaluate_training(trials, oracle)
        trials[0][0]['density_normalized'] = [0.] * 281
        # The other folds would use the altered seed; inspect seed zero's template directly.
        from shock_alignment import template
        from shock_crossfit_gradient import peak_position
        x, y = template(trials, 0)
        error = peak_position(x, y) - np.log(np.log(2))
        self.assertEqual(error, before['crossing_training_recentered']['peak_errors'][0])

    def test_invalid_input(self):
        with self.assertRaises(ValueError):training_profile(GRID, 'unknown')
        with self.assertRaises(ValueError):evaluate(self.data()[:2], 'logistic', 4.)
        data = self.data();data[0,0,0] = np.nan
        with self.assertRaises(ValueError):evaluate(data, 'logistic', 4.)


if __name__ == '__main__':
    unittest.main()
