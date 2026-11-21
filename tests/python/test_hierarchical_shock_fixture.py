from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from shock_hierarchical_fixture import fixture
from validate_shock_estimator import shape


class HierarchicalFixtureTest(unittest.TestCase):
    def test_noiseless_profiles_follow_true_motion(self):
        for kind in ('logistic', 'gompertz'):
            runs, info = fixture(kind, 1234, windows=32, run_sigma=0, window_sigma=0)
            for run, source in zip(runs, info['sources']):
                for window, shift in zip(run, source['true_shifts']):
                    np.testing.assert_array_equal(window['density_normalized'], shape(np.array(window['x'])-shift, kind))

    def test_persistent_effect_is_not_removed_or_resampled_per_window(self):
        runs, _ = fixture('logistic', 1234, windows=32, window_sigma=0, motion_amplitude=0)
        for run in runs:
            np.testing.assert_array_equal(run[0]['density_normalized'], run[-1]['density_normalized'])
        self.assertNotEqual(runs[0][0]['density_normalized'], runs[1][0]['density_normalized'])
        residual = np.mean([np.array(r[0]['density_normalized'])-shape(np.array(r[0]['x']), 'logistic') for r in runs], axis=0)
        self.assertGreater(float(max(abs(residual))), 0)

    def test_reproducible_prefix_and_separate_motion_stream(self):
        a, _ = fixture('logistic', 1234, windows=32)
        b, _ = fixture('logistic', 1234, windows=64)
        c, _ = fixture('logistic', 1234, windows=32, motion_amplitude=0)
        self.assertEqual(a, [run[:32] for run in b])
        _, ia = fixture('logistic', 1234, windows=32)
        for moving, fixed, source in zip(a, c, ia['sources']):
            for m, f, shift in zip(moving, fixed, source['true_shifts']):
                x = np.array(m['x'])
                np.testing.assert_allclose(np.array(m['density_normalized'])-shape(x-shift, 'logistic'), np.array(f['density_normalized'])-shape(x, 'logistic'), rtol=0, atol=8*np.finfo(float).eps)
        different, _ = fixture('logistic', 1235, windows=32)
        self.assertNotEqual(a[0][0]['density_normalized'], different[0][0]['density_normalized'])

    def test_invalid_parameters_rejected(self):
        for options in ({'temporal_phi':1}, {'window_sigma':-1}, {'run_sigma':float('nan')}, {'windows':16}, {'downstream_ratio':0}, {'motion_amplitude':2}):
            with self.assertRaises(ValueError):fixture('logistic', 1234, **options)


if __name__ == '__main__':
    unittest.main()
