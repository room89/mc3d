from pathlib import Path
import sys
import unittest
import math
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from screen_simultaneous_band import t7_survival,critical_value,coverage_summary,trial


class SimultaneousBandTest(unittest.TestCase):
    def test_t7_against_integral_and_known_quantile(self):
        self.assertAlmostEqual(t7_survival(0),.5,places=14)
        self.assertAlmostEqual(t7_survival(2.3646242515927844),.025,places=12)
        for x in (1.,3.,7.):
            angle=math.atan2(math.sqrt(7),x);u=np.linspace(0,angle,10001)
            y=np.sin(u)**6
            integral=np.sum((y[1:]+y[:-1])*.5)*(u[1]-u[0])*16/(5*math.pi)
            self.assertAlmostEqual(t7_survival(x),integral,places=8)

    def test_family_error_and_width(self):
        critical=critical_value(113)
        self.assertGreater(critical,critical_value(1))
        self.assertAlmostEqual(226*t7_survival(critical),.05,places=11)
        with self.assertRaises(ValueError):critical_value(0)

    def test_wilson_limits_and_mcse(self):
        r=coverage_summary([True]*2048)
        self.assertEqual(r['coverage'],1);self.assertEqual(r['mcse'],0)
        self.assertLess(r['wilson95'][0],1);self.assertAlmostEqual(r['wilson95'][1],1)
        r=coverage_summary([True,False]);self.assertEqual(r['coverage'],.5)

    def test_trial_seed_reproducible(self):
        b=dict(nodes=9,x_bounds=[-14,14],groups=8,runs_per_group=3,windows=4,shift_sd=.08,temporal_phi=.7,node_noise_sd=.003,common_noise_sd=.003,run_offset_sd=.002,group_offset_sd=.002,cases=['same_sign_common'])
        spec=(0,dict(seed=123),b,critical_value(9))
        self.assertEqual(trial(spec),trial(spec))
        self.assertNotEqual(trial(spec),trial((1,spec[1],b,spec[3])))
