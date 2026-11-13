from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from control_front_shared_noise import ar1,profile,fronts,loading,summarize


class SharedNoiseTest(unittest.TestCase):
    def test_noise_free_shift_recovery(self):
        x=np.linspace(-14,14,561);s=np.array([[0.,.3,-.4]])
        np.testing.assert_allclose(fronts(x,profile(x,s)),s,atol=1e-12)

    def test_common_field_changes_crossing_in_expected_direction(self):
        x=np.linspace(-14,14,561);v=profile(x,np.zeros((1,1)))
        self.assertLess(fronts(x,v+.005*loading(x,'same_sign_common'))[0,0],0)
        self.assertGreater(fronts(x,v+.005*loading(x,'opposite_sign_common'))[0,0],0)
        np.testing.assert_array_equal(loading(x,'independent_nodes'),np.zeros_like(x))

    def test_ar1_recursion_and_reproducibility(self):
        raw=np.random.default_rng(123).normal(size=(2,3))
        actual=ar1(np.random.default_rng(123),(2,3),.7)
        np.testing.assert_array_equal(actual[:,0],raw[:,0])
        np.testing.assert_allclose(actual[:,1],.7*raw[:,0]+np.sqrt(.51)*raw[:,1])
        with self.assertRaises(ValueError):ar1(np.random.default_rng(1),(2,3),1.)

    def test_mcse(self):
        s=summarize([-1.,1.]);self.assertEqual(s['mean'],0);self.assertEqual(s['mcse'],1)
        with self.assertRaises(ValueError):summarize([float('nan'),1])
