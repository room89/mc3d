import importlib.util
import math
from pathlib import Path
import unittest
spec = importlib.util.spec_from_file_location('argon_vhs', Path(__file__).resolve().parents[2]/'tools/argon_vhs.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)

class ArgonTest(unittest.TestCase):
    def test_reference_and_fit(self):
        self.assertAlmostEqual(a.viscosity(298.15), 22.5666e-6, places=15)
        c,r = a.fit()
        self.assertEqual((c,r), a.fit())
        self.assertTrue(r['accepted'])
        self.assertLessEqual(r['max_relative_error'], .05)
        w=c['vhs_omega'];d=c['vhs_diameter_m']
        mu=15*math.sqrt(a.MASS*a.KB*300/math.pi)/(2*(5-2*w)*(7-2*w)*d*d)
        self.assertAlmostEqual(mu/r['mu_reference_pa_s'],1,places=14)
        with self.assertRaises(ValueError): a.viscosity(20)
