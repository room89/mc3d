import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
    from shock_alignment import GRID
    from shock_ensemble_estimator import estimate,GROUPS,PROFILE_GRID,centered_density,width_density_interval,cubic_profile
except ImportError:
    np=None

@unittest.skipIf(np is None,'Requires numpy')
class EnsembleTest(unittest.TestCase):
    def test_cubic_interpolation_reproduces_cubic_without_extrapolation(self):
        x=np.linspace(-3,3,61);points=np.linspace(-3,3,123)
        f=lambda x:.3+x+.1*x*x-.01*x**3
        np.testing.assert_allclose(cubic_profile(x,f(x),points),f(points),atol=1e-12)
        with self.assertRaises(ValueError):cubic_profile(x,f(x),np.array([3.1]))

    def test_interpolation_budget_is_separate_from_sampling_radius(self):
        groups=self.groups()
        for g in groups:
            g['density_cubic']=(np.asarray(g['density'])+.0002).tolist()
            g['slope_cubic']=.251
        baseline=estimate(groups,repeats=128)
        result=estimate(groups,repeats=128,width_band=True,interpolation_sensitivity=True)
        self.assertAlmostEqual(result['density_processing_sensitivity'],.0002)
        self.assertAlmostEqual(result['density_statistical_radius'],baseline['density_simultaneous_radius'])
        self.assertGreater(result['gradient_thickness']['processing_sensitivity'],0)
        self.assertGreater(result['width_10_90']['processing_sensitivity'],0)

    def test_width_bounds_contain_truth_when_density_band_does(self):
        from shock_mean_profile import measure
        rng=np.random.default_rng(8080)
        for truth in (1/(1+np.exp(-GRID)),np.exp(-np.exp(-GRID))):
            exact=measure(truth)[0,0]
            for _ in range(20):
                profile=truth+rng.uniform(-.002,.002,len(GRID))
                lo,hi=width_density_interval(profile,.002)
                self.assertLessEqual(lo,exact)
                self.assertGreaterEqual(hi,exact)
        with self.assertRaises(ValueError):width_density_interval(truth,.1)

    def groups(self):
        return [dict(density=(1/(1+np.exp(-GRID))).tolist(),slopes=[.25]*3) for _ in range(GROUPS)]

    def test_identical_groups(self):
        r=estimate(self.groups(),repeats=128)
        self.assertAlmostEqual(r['gradient_thickness']['mean'],4.)
        self.assertLess(r['gradient_thickness']['statistical_radius'],1e-12)
        self.assertLess(r['density_simultaneous_radius'],1e-12)

    def test_slopes_are_combined_before_inversion(self):
        groups=self.groups()
        for i,g in enumerate(groups):g['slopes']=[(.20,.25,.30)[i%3]]*3
        r=estimate(groups,repeats=128)
        self.assertAlmostEqual(r['gradient_thickness']['mean'],4.)
        self.assertGreater(r['gradient_thickness']['statistical_radius'],0.)

    def test_final_density_is_centered_without_scaling(self):
        y=1/(1+np.exp(-(PROFILE_GRID-.37)))
        centered=centered_density(y)[0]
        self.assertAlmostEqual(centered[len(GRID)//2],.5)
        np.testing.assert_allclose(centered,1/(1+np.exp(-GRID)),atol=4e-5)
        np.testing.assert_allclose(centered_density(y,interpolation='cubic')[0],1/(1+np.exp(-GRID)),atol=1e-7)
        groups=[dict(density=y.tolist(),slopes=[.25]*3) for _ in range(GROUPS)]
        r=estimate(groups,repeats=128,center_profile=True)
        self.assertAlmostEqual(r['density'][len(GRID)//2],.5)
        self.assertAlmostEqual(r['gradient_thickness']['mean'],4.)

    def test_centering_does_not_extrapolate(self):
        with self.assertRaises(ValueError):centered_density(1/(1+np.exp(-(PROFILE_GRID-2))))

    def test_invalid_groups(self):
        with self.assertRaises(ValueError):estimate(self.groups()[:3],repeats=128)
        groups=self.groups();groups[0]['slopes'][1]=float('nan')
        with self.assertRaises(ValueError):estimate(groups,repeats=128)

if __name__=='__main__':unittest.main()
