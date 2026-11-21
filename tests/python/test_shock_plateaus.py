import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
    from analyze_shock_plateaus import mean_interval
    from audit_shock_validation import gate
except ImportError:
    np=None

@unittest.skipIf(np is None,'Requires numpy')
class PlateauTest(unittest.TestCase):
    def test_constant_wrong_plateau_is_not_stationary_target(self):
        mean,u,_=mean_interval(np.full((3,32,2),.02),repeats=32)
        self.assertFalse(gate(float(mean[0]),float(u[0]),.01)['passed'])

    def test_between_seed_variation_is_retained(self):
        values=np.zeros((3,32,2));values[0]=-.004;values[2]=.004
        mean,u,_=mean_interval(values,repeats=32)
        self.assertAlmostEqual(float(mean[0]),0.)
        self.assertGreater(u[0],.005)
        self.assertFalse(gate(float(mean[0]),float(u[0]),.01)['passed'])

    def test_invalid_values(self):
        with self.assertRaises(ValueError):mean_interval(np.zeros((2,32,2)),repeats=32)
        with self.assertRaises(ValueError):mean_interval(np.full((3,32,2),np.nan),repeats=32)

if __name__=='__main__':unittest.main()
