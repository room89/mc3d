import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from profile_level_bounds import crossing_bounds,width_bounds,monotone_band

class ProfileLevelBoundsTests(unittest.TestCase):
    def test_known_logistic_width_and_translation(self):
        x=np.linspace(-8,8,321);y=(1+np.tanh(x/2))/2
        r=width_bounds(x,y-.005,y+.005);truth=2*np.log(9)
        self.assertLessEqual(r['width'][0],truth);self.assertGreaterEqual(r['width'][1],truth)
        shifted=width_bounds(x+123,y-.005,y+.005)
        np.testing.assert_allclose(shifted['width'],r['width'],atol=1e-13)
        np.testing.assert_allclose(np.array(shifted['crossings'])-123,r['crossings'],atol=1e-13)

    def test_no_subgrid_linear_assumption(self):
        # Same nodal values admit crossings arbitrarily near either endpoint.
        self.assertEqual(crossing_bounds([0,1],[0,1],[0,1],.5),[0,1])
        self.assertEqual(width_bounds([0,1],[0,1],[0,1])['width'],[0,1])

    def test_monotonicity_tightening_and_infeasibility(self):
        _,lo,hi=monotone_band([0,1,2],[0,.4,.2],[.6,.5,1])
        np.testing.assert_allclose(lo,[0,.4,.4]);np.testing.assert_allclose(hi,[.5,.5,1])
        with self.assertRaises(ValueError):monotone_band([0,1],[.8,0],[1,.2])

    def test_insufficient_domain_rejected(self):
        with self.assertRaisesRegex(ValueError,'not bracketed'):
            width_bounds([0,1],[.2,.8],[.3,.9])
        with self.assertRaises(ValueError):width_bounds([1,0],[0,1],[0,1])
        with self.assertRaises(ValueError):width_bounds([0,1],[0,1],[0,1],(.9,.1))

    def test_all_piecewise_linear_profiles_inside_band_are_enclosed(self):
        # Enumerate a finite family, then use its exact envelope as a valid band.
        x=np.array([-3,-1,0,2,4.]);profiles=np.array([
            [0,.02,.2,.8,1],[0,.2,.4,.9,1],[0,.05,.6,.95,1],
            [0,.1,.1,.9,1]])
        r=width_bounds(x,profiles.min(0),profiles.max(0))
        for y in profiles:
            q=[]
            for level in (.1,.9):
                j=np.flatnonzero(y>=level)[0];q.append(float(x[j-1]+(level-y[j-1])/(y[j]-y[j-1])*(x[j]-x[j-1])))
            self.assertLessEqual(r['width'][0],q[1]-q[0]);self.assertGreaterEqual(r['width'][1],q[1]-q[0])
