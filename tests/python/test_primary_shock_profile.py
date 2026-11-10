import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from compare_primary_shock_profile import compare, GRID, nominal_axes, monotone_width_bounds, compare_width


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

    def test_shared_nominal_axes_use_pixels_without_changing_original(self):
        p, _ = self.inputs()
        p.update(axis_origin_x_interval=[9., 11.], axis_pixels_per_lambda_interval=[1., 3.],
                 axis_y0_interval=[99., 101.], axis_y1_interval=[-1., 1.], pixel_extent_padding=.5)
        p['points'] = [dict(pixel_row=25., pixel_black_run=[12., 14.]),
                       dict(pixel_row=75., pixel_black_run=[6., 8.])]
        before = copy.deepcopy(p); fixed = nominal_axes(p)
        self.assertEqual(p, before)
        self.assertEqual(fixed['points'][0]['x_paper_interval'], [.75, 2.25])
        self.assertEqual(fixed['points'][1]['x_paper_interval'], [-2.25, -.75])
        self.assertEqual(fixed['points'][0]['rho_interval'], [.745, .755])
        self.assertEqual(fixed['points'][1]['rho_interval'], [.245, .255])
        self.assertEqual(fixed['stated_relative_length_accuracy'], 0.)

    def test_monotone_width_bounds_and_accuracy_widening(self):
        p, s = self.inputs()
        p['points'] = [dict(pixel_row=i, x_paper_interval=[i/2, i/2],
                            rho_interval=[i/20, i/20]) for i in range(21)]
        narrow = monotone_width_bounds(p, False)['width_reference_length_interval']
        self.assertEqual(narrow, [7., 9.])
        p['stated_normalized_density_accuracy'] = .05
        broad = monotone_width_bounds(p, True)['width_reference_length_interval']
        self.assertLess(broad[0], narrow[0]); self.assertGreater(broad[1], narrow[1])
        s['width_10_90'] = dict(mean=8., statistical_radius=.1, processing_sensitivity=.2)
        result = compare_width(p, s, 1.)
        self.assertEqual(result['nonoverlap_gap_m'], 0.)
        self.assertEqual(result['model_width_interval_m'], [7.7, 8.3])
        p['points'] = p['points'][1:-1]
        with self.assertRaisesRegex(ValueError, 'extrapolation'): monotone_width_bounds(p, True)

    def test_wrong_grid_and_extrapolation_are_rejected(self):
        p, s = self.inputs(); wrong = copy.deepcopy(s); wrong['density'] = [0.] * 441
        with self.assertRaisesRegex(ValueError, '401-point'): compare(p, wrong, 1.)
        p['points'][0]['x_paper_interval'] = [-11., -10.]
        with self.assertRaisesRegex(ValueError, 'extrapolation'): compare(p, s, 1.)


if __name__ == '__main__': unittest.main()
