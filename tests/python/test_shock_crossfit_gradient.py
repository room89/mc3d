import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
    from shock_alignment import GRID, align
    from shock_crossfit_gradient import derivative, estimate_gradient
except ImportError:
    np=None

@unittest.skipIf(np is None,'Requires numpy')
class CrossfitGradientTest(unittest.TestCase):
    def fixture(self):
        y=1/(1+np.exp(-GRID))
        data=np.tile(y,(3,2,1))
        info=[dict(training_seeds=[j for j in range(3) if j!=i],
                   template_x=GRID.tolist(),template_density=y.tolist()) for i in range(3)]
        return data,info

    def test_cubic_derivative(self):
        x,d=derivative(2+3*GRID-2*GRID**2+.2*GRID**3)
        np.testing.assert_allclose(d,3-4*x+.6*x*x,atol=1e-10)

    def test_clean_logistic(self):
        result=estimate_gradient(*self.fixture())
        self.assertLess(abs(result['gradient_thickness']/4-1),.001)

    def test_does_not_maximize_on_heldout(self):
        data,info=self.fixture()
        before=estimate_gradient(data,info)
        data[0]+=np.exp(-((GRID-4)/.15)**2)
        after=estimate_gradient(data,info)
        self.assertEqual(before,after)

    def test_average_slopes_before_inverse(self):
        data,info=self.fixture()
        base=estimate_gradient(data,info)['gradient_thickness']
        data[0]*=2
        self.assertAlmostEqual(estimate_gradient(data,info)['gradient_thickness'],base*.75)

    def test_template_excludes_heldout_in_pipeline(self):
        x=np.linspace(-14,14,281)
        trials=[[dict(x_aligned=x.tolist(),density_normalized=(1/(1+np.exp(-(x-shift)))).tolist())
                 for shift in (-.5,.5)] for _ in range(3)]
        data,info=align(trials)
        before=estimate_gradient(data,info)['folds'][0]['peak_position']
        trials[0][0]['density_normalized']=(1/(1+np.exp(-x/1.1))).tolist()
        data,info=align(trials)
        self.assertEqual(before,estimate_gradient(data,info)['folds'][0]['peak_position'])

    def test_invalid_inputs(self):
        data,info=self.fixture()
        info[0]['training_seeds']=[0,1]
        with self.assertRaises(ValueError):estimate_gradient(data,info)
        data,info=self.fixture();data[0]*=-1
        with self.assertRaises(ValueError):estimate_gradient(data,info)
        with self.assertRaises(ValueError):derivative(np.full(401,np.nan))

if __name__=='__main__':unittest.main()
