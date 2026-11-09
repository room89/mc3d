import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
try:
    import numpy as np
    from decompose_gradient_bias import analytic_control, decompose, slopes_for
    from experiment_crossfit_gradient import fixture
except ImportError:
    np = None


@unittest.skipIf(np is None, 'Requires numpy')
class GradientDecompositionTest(unittest.TestCase):
    def test_exact_additive_identity(self):
        d = decompose([.8, 1., 1.2])
        self.assertAlmostEqual(d['mean_relative_slope'], 0.)
        self.assertGreater(d['reciprocal_jensen_gap'], 0.)
        self.assertAlmostEqual(d['thickness_bias'],
            d['inverse_mean_slope_bias'] + d['reciprocal_jensen_gap'])

    def test_constant_bias_has_no_jensen_gap(self):
        d = decompose([.9]*8)
        self.assertAlmostEqual(d['reciprocal_jensen_gap'], 0.)
        self.assertAlmostEqual(d['inverse_mean_slope_bias'], 1/.9-1)

    def test_invalid_slopes(self):
        for values in ([1.], [0., 1.], [-1., 1.], [float('nan'), 1.]):
            with self.assertRaises(ValueError):
                decompose(values)

    def test_noiseless_controls(self):
        for kind, exact in [('logistic', 4.), ('gompertz', np.e), ('flat_peak', 4.)]:
            self.assertLess(abs(analytic_control(kind, exact)['thickness_bias']), .005)
            trials, oracle = fixture(kind, 2400000, 0., 0., noise=0.)
            slopes = slopes_for(trials, oracle, kind, exact)
            for value in slopes.values():
                self.assertLess(abs(value-1), .005)

    def test_default_fixture_unchanged(self):
        a, x = fixture('logistic', 99, .8, .2)
        b, y = fixture('logistic', 99, .8, .2, noise=.03)
        self.assertEqual(a, b)
        np.testing.assert_array_equal(x, y)


if __name__ == '__main__':
    unittest.main()
