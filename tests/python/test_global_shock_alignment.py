from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from shock_global_alignment import global_shift
from shock_alignment import GRID, shift
from shock_ensemble_estimator import PROFILE_GRID
from shock_hierarchical_fixture import fixture
from validate_shock_estimator import shape


class GlobalAlignmentTest(unittest.TestCase):
    def test_known_translations(self):
        x = np.linspace(-13.975, 13.975, 560)
        for kind in ('logistic', 'gompertz'):
            for truth in (-.7, 0., .4):
                window = dict(x=x.tolist(), density_normalized=shape(x-truth, kind).tolist())
                front, info = global_shift(window, GRID, shape(GRID, kind), PROFILE_GRID)
                self.assertLess(abs(front-truth), 2e-4)
                self.assertEqual(info['method'], 'support_bounded_cellwise_huber_v1')

    def test_replayed_boundary_cases(self):
        # Fixed regression seeds from the selected failure replay, not validation.
        for ci, rep, group, source, index in ((3,9,8,2,13),(5,1,1,0,13),(7,9,8,1,19),(7,15,3,2,30)):
            kind = 'logistic' if ci == 3 else 'gompertz'
            motion = 0. if ci == 5 else .05
            runs, _ = fixture(kind,192000000+ci*100000+rep*100+group, windows=32,run_sigma=.01,motion_amplitude=motion)
            window = runs[source][index]
            with self.assertRaisesRegex(ValueError, 'Boundary minimum'):
                shift(window, GRID, shape(GRID,kind))
            front, info = global_shift(window, GRID, shape(GRID,kind), PROFILE_GRID)
            self.assertTrue(np.isfinite(front))
            self.assertLess(info['support_bounds'][0], front)
            self.assertGreater(info['support_bounds'][1], front)

    def test_replayed_ambiguity_is_not_silently_accepted(self):
        runs, _ = fixture('logistic',192000000+3*100000+8*100,windows=32,run_sigma=.01,motion_amplitude=.05)
        with self.assertRaisesRegex(ValueError, 'Ambiguous minimum'):
            global_shift(runs[1][13], GRID, shape(GRID,'logistic'), PROFILE_GRID)

    def test_actual_support_boundary_and_invalid_data(self):
        x=np.linspace(-13.975,13.975,560)
        with self.assertRaisesRegex(ValueError, 'Boundary minimum'):
            global_shift(dict(x=x,density_normalized=shape(x-3.5,'logistic')),GRID,shape(GRID,'logistic'),PROFILE_GRID)
        with self.assertRaises(ValueError):
            global_shift(dict(x=x,density_normalized=np.full(len(x),np.nan)),GRID,shape(GRID,'logistic'),PROFILE_GRID)
        with self.assertRaises(ValueError):
            global_shift(dict(x=x[::-1],density_normalized=shape(x,'logistic')),GRID,shape(GRID,'logistic'),PROFILE_GRID)


if __name__ == '__main__':
    unittest.main()
