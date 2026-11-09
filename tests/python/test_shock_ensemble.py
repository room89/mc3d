import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
    from shock_alignment import GRID
    from shock_ensemble_estimator import estimate,GROUPS
except ImportError:
    np=None

@unittest.skipIf(np is None,'Requires numpy')
class EnsembleTest(unittest.TestCase):
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

    def test_invalid_groups(self):
        with self.assertRaises(ValueError):estimate(self.groups()[:3],repeats=128)
        groups=self.groups();groups[0]['slopes'][1]=float('nan')
        with self.assertRaises(ValueError):estimate(groups,repeats=128)

if __name__=='__main__':unittest.main()
