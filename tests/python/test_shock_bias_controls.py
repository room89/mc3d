from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from diagnose_shock_ensemble_bias import quintic_derivative,summarize
from shock_alignment import GRID


class BiasControlTest(unittest.TestCase):
    def test_quintic_control_reproduces_polynomial_derivative(self):
        y=.001*GRID**5+.03*GRID**3+.2*GRID+.5
        x,d=quintic_derivative(y)
        np.testing.assert_allclose(d,.005*x**4+.09*x**2+.2,atol=1e-10)

    def test_paired_error_contrast_removes_common_sampling_noise(self):
        rows=[]
        for i,v in enumerate((-.1,.1)):
            rows.append(dict(repetition=i,relative_errors=dict(oracle_width=v,refined_width=v-.002,
                oracle_cubic=v,refined_cubic=v+.001,refined_quintic_same_peak=v)))
        result=summarize(rows+[dict(repetition=2,error='Unresolved')])
        self.assertEqual(result['failed_count'],1)
        paired=result['summary']['alignment_width']
        self.assertAlmostEqual(paired['mean'],-.002)
        self.assertLess(paired['mean_two_se'],1e-12)


if __name__=='__main__':unittest.main()
