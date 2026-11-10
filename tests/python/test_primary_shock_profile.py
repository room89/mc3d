import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from compare_primary_shock_profile import compare, GRID


class PrimaryProfileTest(unittest.TestCase):
    def inputs(self):
        primary = dict(reference_length_m=1., stated_relative_length_accuracy=.1,
                       stated_normalized_density_accuracy=.01,
                       points=[dict(pixel_row=1, x_paper_interval=[-2., -1.], rho_interval=[.4, .5])])
        summary = dict(center_profile=True, density=((GRID+10)/20).tolist(), density_simultaneous_radius=.02)
        return primary, summary

    def test_negative_coordinate_scale_bounds_and_no_acceptance(self):
        p, s = self.inputs(); r = compare(p, s, 1.)
        self.assertEqual(r['points'][0]['model_x_interval'], [-2.2, -.9])
        self.assertEqual(r['nonoverlap_count'], 0)
        self.assertEqual(r['physical_acceptance'], 'incomplete')
        p['points'][0]['rho_interval'] = [.8, .9]
        self.assertGreater(compare(p, s, 1.)['maximum_nonoverlap_gap'], .3)

    def test_interior_extrema_are_retained(self):
        p, s = self.inputs(); s['density'][np.argmin(abs(GRID+1.5))] = .9
        self.assertAlmostEqual(compare(p, s, 1.)['points'][0]['model_density_interval'][1], .92)

    def test_wrong_grid_and_extrapolation_are_rejected(self):
        p, s = self.inputs(); wrong = copy.deepcopy(s); wrong['density'] = [0.] * 441
        with self.assertRaisesRegex(ValueError, '401-point'): compare(p, wrong, 1.)
        p['points'][0]['x_paper_interval'] = [-11., -10.]
        with self.assertRaisesRegex(ValueError, 'extrapolation'): compare(p, s, 1.)


if __name__ == '__main__': unittest.main()
