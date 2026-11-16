from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from screen_fixed_profile_coverage import estimate,truth
from control_front_shared_noise import profile


class FixedProfileTest(unittest.TestCase):
    def test_known_truth(self):
        x=np.linspace(-14,14,113)
        np.testing.assert_allclose(truth(x,0,32),profile(x,0),atol=1e-14)
        shifted=truth(x,.08,64)
        np.testing.assert_allclose(shifted,truth(x,.08,128),atol=1e-13)
        np.testing.assert_allclose(shifted+shifted[::-1],1,atol=1e-14)
        self.assertGreater(abs(shifted-profile(x,0)).max(),1e-4)

    def test_group_offsets_not_removed_or_counted_as_independent_runs(self):
        raw=np.repeat(np.array([-1.,1.]),3)[:,None,None]*np.ones((6,64,2))
        center,half=estimate(raw,2,3,2.)
        np.testing.assert_array_equal(center,[0,0]);np.testing.assert_allclose(half,[2,2])
        # Repeating correlated windows must not shrink the interval.
        c,h=estimate(np.repeat(raw,2,axis=1),2,3,2.)
        np.testing.assert_array_equal(h,half)

    def test_translation_of_data_preserves_uncertainty(self):
        raw=np.arange(6*4*2).reshape(6,4,2)
        c,h=estimate(raw,3,2,2.);d,k=estimate(raw+100,3,2,2.)
        np.testing.assert_allclose(d,c+100);np.testing.assert_allclose(k,h)

    def test_invalid_groups(self):
        with self.assertRaises(ValueError):estimate(np.zeros((5,64,2)),2,3,2)
        with self.assertRaises(ValueError):estimate(np.full((6,64,2),np.nan),2,3,2)
