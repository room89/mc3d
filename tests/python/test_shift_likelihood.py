from pathlib import Path
import sys
import math
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from shock_shift_likelihood import marginal_shift_likelihood


class ShiftLikelihoodTest(unittest.TestCase):
    def test_single_shift_matches_gaussian_density(self):
        result=marginal_shift_likelihood([1.,2.],[[0.,0.]],np.eye(2),[.3],[1.])
        self.assertAlmostEqual(result['log_marginal_likelihood'],-math.log(2*math.pi)-2.5)
        self.assertEqual(result['posterior_weights'],[1.])
        self.assertEqual(result['posterior_shift_variance'],0.)

    def test_correlated_covariance_matches_direct_two_point_mixture(self):
        covariance=np.array([[2.,.7],[.7,1.]])
        y=np.array([.4,-.2]); predictions=np.array([[0.,0.],[1.,-.5]])
        prior=np.array([.3,.7]); residual=y-predictions
        q=np.einsum('ij,jk,ik->i',residual,np.linalg.inv(covariance),residual)
        density=np.exp(-.5*q)/(2*np.pi*np.sqrt(np.linalg.det(covariance)))
        expected=prior*density;result=marginal_shift_likelihood(y,predictions,covariance,[-.2,.3],prior)
        self.assertAlmostEqual(result['log_marginal_likelihood'],float(np.log(expected.sum())))
        np.testing.assert_allclose(result['posterior_weights'],expected/expected.sum(),rtol=1e-14)
        reverse=marginal_shift_likelihood(y,predictions[::-1],covariance,[.3,-.2],prior[::-1])
        self.assertAlmostEqual(result['log_marginal_likelihood'],reverse['log_marginal_likelihood'])
        self.assertAlmostEqual(result['posterior_shift_mean'],reverse['posterior_shift_mean'])

    def test_uncertain_shift_retains_mass_and_low_noise_concentrates(self):
        loose=marginal_shift_likelihood([0.],[[0.],[1.]],[[100.]],[-1.,1.],[1.,1.])
        tight=marginal_shift_likelihood([0.],[[0.],[1.]],[[.001]],[-1.,1.],[1.,1.])
        self.assertGreater(loose['posterior_shift_variance'],.99)
        self.assertLess(tight['posterior_shift_variance'],1e-10)
        self.assertTrue(np.isfinite(tight['log_marginal_likelihood']))
        excluded=marginal_shift_likelihood([0.],[[0.],[1.]],[[1.]],[-1.,1.],[0.,1e300])
        self.assertEqual(excluded['posterior_weights'],[0.,1.])

    def test_invalid_covariance_and_prior_are_rejected(self):
        for covariance in ([[0.]], [[-1.]], [[float('nan')]]):
            with self.assertRaises(ValueError):marginal_shift_likelihood([0.],[[0.]],covariance,[0.],[1.])
        for prior in ([0.],[-1.],[float('inf')]):
            with self.assertRaises(ValueError):marginal_shift_likelihood([0.],[[0.]],[[1.]],[0.],prior)
        with self.assertRaises(ValueError):marginal_shift_likelihood([0.,0.],[[0.,0.]],[[1.,.5],[0.,1.]],[0.],[1.])


if __name__ == '__main__':
    unittest.main()
