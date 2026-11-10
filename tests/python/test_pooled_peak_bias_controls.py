from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from diagnose_pooled_peak_bias import slopes, summarize
from shock_alignment import GRID
from shock_ensemble_estimator import PROFILE_GRID, centered_density
from shock_group_pooling import crossfit_pool_points
from validate_shock_estimator import shape


class PooledPeakControlsTest(unittest.TestCase):
    def test_clean_logistic_known_peak_and_order(self):
        y=np.tile(shape(GRID,'logistic'),(6,1))
        cubic,peaks=slopes(y,position=0.)
        quintic,_=slopes(y,position=0.,quintic=True)
        self.assertEqual(peaks,[0.]*6)
        self.assertLess(abs(quintic-.25),abs(cubic-.25))
        self.assertAlmostEqual(slopes(y)[0],cubic,places=13)

    def test_control_reproduces_production_with_distinct_subgroups(self):
        points=[]
        for shift in np.linspace(-.1,.1,6):
            y=shape(PROFILE_GRID-shift,'logistic')
            points.append(dict(density=y.tolist(),density_cubic=y.tolist(),slopes=[.25]*3,slope_cubic=.25))
        expected=crossfit_pool_points(points)['slopes'][1]
        actual,_=slopes(centered_density([p['density'] for p in points]))
        self.assertAlmostEqual(actual,expected,places=14)

    def test_quintic_derivative_preserves_known_polynomial_for_each_width(self):
        y=np.tile(.5+.1*GRID+1e-6*GRID**5,(6,1))
        for width in (1.5,2.,2.5):
            value,_=slopes(y,position=0.,quintic=True,width=width)
            self.assertAlmostEqual(value,.1,places=12)

    def test_paired_contrasts_and_errors_remain_visible(self):
        keys=('oracle_known_cubic','oracle_known_quintic','oracle_selected_cubic',
              'refined_uncentered_cubic','refined_centered_cubic',
              'refined_centered_quintic_same_peaks','refined_cubic_interpolation')
        rows=[dict(relative_errors=dict(zip(keys,[i,i,i+1,i+2,i+3,i+4,i+5]))) for i in (0.,1.)]
        rows.append(dict(repetition=2,error='Ambiguous minimum'))
        result=summarize(rows)
        self.assertEqual(result['valid'],2);self.assertEqual(len(result['errors']),1)
        self.assertEqual(result['summary']['centering']['mean'],1.)
        self.assertEqual(result['summary']['centering']['two_mc_se'],0.)


if __name__=='__main__':unittest.main()
