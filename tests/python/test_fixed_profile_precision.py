import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from audit_fixed_profile_precision import grouped, band, temporal_diagnostic
from screen_simultaneous_band import critical_value


class FixedProfilePrecisionTests(unittest.TestCase):
    def setUp(self):
        self.sources=[dict(group=f'g{i//3}') for i in range(24)]
        self.raw=np.repeat(np.arange(8.),3)[:,None,None]*np.ones((24,4,2))

    def test_order_and_window_duplication_preserve_group_uncertainty(self):
        expected=grouped(self.raw,self.sources)
        np.testing.assert_array_equal(grouped(self.raw[::-1],self.sources[::-1]),expected)
        np.testing.assert_array_equal(grouped(np.repeat(self.raw,2,axis=1),self.sources),expected)
        mean,half=band(expected,2.)
        np.testing.assert_allclose(mean,3.5)
        np.testing.assert_allclose(half,2*np.std(np.arange(8.),ddof=1)/np.sqrt(8))

    def test_offsets_are_not_removed(self):
        a=grouped(self.raw,self.sources)
        b=grouped(self.raw+7,self.sources)
        np.testing.assert_allclose(b-a,7)
        np.testing.assert_allclose(band(a,2)[1],band(b,2)[1])

    def test_incomplete_groups_and_nonfinite_rejected(self):
        with self.assertRaises(ValueError):grouped(self.raw[:-1],self.sources[:-1])
        self.raw[0,0,0]=np.nan
        with self.assertRaises(ValueError):grouped(self.raw,self.sources)

    def test_joint_family_widens_band(self):
        self.assertGreater(critical_value(560,.025),critical_value(113,.05))

    def test_half_change_retains_translation(self):
        x=np.linspace(-4,4,81)
        first=(1+np.tanh(x))/2
        second=(1+np.tanh(x-.2))/2
        raw=np.tile(np.concatenate([np.tile(first,(32,1)),np.tile(second,(32,1))]),(24,1,1))
        r=temporal_diagnostic(raw,self.sources,x)
        self.assertAlmostEqual(r['first_half_front'],0,places=12)
        self.assertAlmostEqual(r['second_half_front'],.2,places=12)
        np.testing.assert_allclose(r['half_profile_change'],second-first,atol=1e-14)
        with self.assertRaises(ValueError):temporal_diagnostic(raw[:,:32],self.sources,x)
