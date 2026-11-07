import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
except ImportError:
    np=None
if np is not None:
    from shock_mean_profile import GRID,measure,aligned,bootstrap_profiles,estimate
    from validate_shock_estimator import shape,truth,synthetic


@unittest.skipIf(np is None,'Experimental reanalysis requires numpy')
class MeanProfileTest(unittest.TestCase):
    def test_known_symmetric_and_asymmetric_profiles(self):
        for kind in ('logistic','gompertz'):
            for scale in (1.,2.):
                actual=measure(shape(GRID-.73,kind,scale))[0]
                self.assertTrue(np.all(abs(actual/truth(kind,scale)-1)<.001))

    def test_unresolved_narrow_profile_is_not_qualified(self):
        actual=measure(shape(GRID,'gompertz',.25))[0]
        self.assertGreater(abs(actual[1]/truth('gompertz',.25)[1]-1),.05)

    def test_derivative_preserves_linear_slope(self):
        result=measure((GRID+10)/20)[0]
        self.assertAlmostEqual(result[0],16)
        self.assertAlmostEqual(result[1],20)

    def test_constant_trials_have_zero_resampling_spread(self):
        data=np.tile(shape(GRID,'logistic'),(3,32,1))
        result=estimate(data,repeats=32)
        for key in ('width_10_90','gradient_thickness'):
            self.assertLess(result[key]['statistical_radius'],1e-11)

    def test_seed_variation_is_not_treated_as_96_independent_windows(self):
        data=np.array([np.tile(shape(GRID,'logistic',s),(32,1)) for s in (.9,1.,1.1)])
        result=estimate(data,repeats=32)
        self.assertGreater(result['width_10_90']['seed_t_radius'],.9)

    def test_resampling_keeps_spatial_dependence(self):
        rng=np.random.default_rng(123);one=rng.normal(size=(3,32,1));data=np.concatenate([one,2*one],axis=2)
        a=bootstrap_profiles(data,8,32,np.random.default_rng(1))
        b=bootstrap_profiles(data,8,32,np.random.default_rng(1))
        np.testing.assert_array_equal(a,b)
        np.testing.assert_allclose(a[:,1],2*a[:,0])

    def test_rejects_invalid_or_uncovered_profiles(self):
        with self.assertRaises(ValueError):measure([float('nan')]*len(GRID))
        with self.assertRaises(ValueError):measure(np.zeros(len(GRID)))
        with self.assertRaises(ValueError):aligned([[{'x_aligned':[-1,1],'density_normalized':[0,1]}]]*3)
        with self.assertRaises(ValueError):estimate(np.zeros((3,31,len(GRID))),32)

    def test_noisy_translation_fixture_is_reproducible(self):
        a,error=synthetic('logistic',1,noise=0,phi=0,count=32)
        b,_=synthetic('logistic',1,noise=0,phi=0,count=32)
        np.testing.assert_array_equal(a,b)
        self.assertLess(error,.001)
        self.assertLess(np.max(abs(measure(a.mean(axis=(0,1)))[0]/truth('logistic')-1)),.002)
