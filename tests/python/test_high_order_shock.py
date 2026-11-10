from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_alignment import GRID
from shock_ensemble_estimator import PROFILE_GRID, centered_density
from shock_high_order_estimator import METHOD, derivative5, pooled_point, estimate
from validate_shock_estimator import shape
from calibrate_high_order_shock import cached


class HighOrderShockTest(unittest.TestCase):
    def points(self):
        return [dict(density=shape(PROFILE_GRID/(1+i*.002),'logistic').tolist(),
                     density_cubic=shape(PROFILE_GRID/(1+i*.002),'logistic').tolist()) for i in range(6)]

    def test_polynomial_derivative_and_invalid_data(self):
        for width in (1.5,2.,2.5):
            x,d=derivative5(.5+.1*GRID+1e-6*GRID**5,width)
            np.testing.assert_allclose(d,.1+5e-6*x**4,atol=1e-12,rtol=1e-12)
        with self.assertRaises(ValueError):derivative5(np.full_like(GRID,np.nan),2.)
        with self.assertRaises(ValueError):derivative5(np.ones_like(GRID),.05)

    def test_peak_training_excludes_heldout_subgroup(self):
        points=self.points(); centered=centered_density([p['density_cubic'] for p in points],interpolation='cubic')
        calls=[]
        def position(x,y,width):calls.append(y.copy());return 0.
        with patch('shock_high_order_estimator.peak_position',side_effect=position):
            result=pooled_point(points)
        self.assertEqual(len(calls),24)
        for i in range(18):
            np.testing.assert_allclose(calls[i],np.delete(centered,i%6,axis=0).mean(0))
        self.assertEqual(result['method'],METHOD)

    def test_complete_estimate_and_processing_budget(self):
        y=shape(PROFILE_GRID,'logistic').tolist();g=pooled_point([dict(density=y,density_cubic=y) for _ in range(6)])
        result=estimate([g]*12,repeats=128)
        self.assertEqual(result['status'],'experimental_coverage_unverified')
        self.assertEqual(len(result['density']),len(GRID))
        gradient=result['gradient_thickness']
        self.assertLess(abs(gradient['mean']/4-1),.0001)
        self.assertLessEqual(abs(gradient['mean']-4),gradient['statistical_radius']+gradient['processing_sensitivity'])
        with self.assertRaisesRegex(ValueError,'matching high-order'):
            estimate([dict(g,method='previous_protocol')]*12,repeats=128)

    def test_checkpoint_rejects_changed_protocol_before_running(self):
        import tempfile,json
        with tempfile.TemporaryDirectory() as temp:
            spec=('logistic',0.,0.,.015,0,110000000,64)
            p=Path(temp)/'logistic-0.0-0.0-0.015-0.json'
            p.write_text(json.dumps(dict(protocol={'method':'old'},spec=spec,result={})))
            with self.assertRaisesRegex(ValueError,'protocol mismatch'):
                cached((spec,temp,{'method':METHOD}))


if __name__=='__main__':unittest.main()
