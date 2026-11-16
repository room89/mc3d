from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_temporal_covariance import estimate


class TemporalCovarianceTest(unittest.TestCase):
    def test_exact_second_moment_and_stationarity_projection(self):
        y=np.array([[[1.,2.],[2.,1.]], [[-1.,-2.],[-2.,-1.]]])
        result=estimate(y,[0,0],np.eye(2))
        np.testing.assert_allclose(result['temporal_covariance'],[[2.5,2],[2,2.5]])
        shifted=y+np.array([0.,1.])
        np.testing.assert_allclose(estimate(shifted,[0,1],np.eye(2))['temporal_covariance'],result['temporal_covariance'])

    def test_known_spatial_whitening(self):
        y=np.array([[[1.,2.],[2.,1.]], [[-1.,-2.],[-2.,-1.]]])
        spatial=np.array([[2.,.3],[.3,1.]])
        colored=y@np.linalg.cholesky(spatial).T
        np.testing.assert_allclose(estimate(colored,[0,0],spatial)['temporal_covariance'],[[2.5,2],[2,2.5]])

    def test_fresh_correlated_noise_recovers_total_covariance(self):
        rng=np.random.default_rng(202611170001)
        temporal=np.array([[.0005,.00042],[.00042,.0005]])
        spatial=np.array([[1.,.3],[.3,2.]])
        noise=(rng.normal(size=(2000,4))@np.linalg.cholesky(np.kron(temporal,spatial)).T).reshape(2000,2,2)
        result=estimate(noise+np.array([0.,1.]),[0,1],spatial)
        relative=np.linalg.norm(np.array(result['temporal_covariance'])-temporal)/np.linalg.norm(temporal)
        self.assertLess(relative,.05)

    def test_invalid_and_degenerate_inputs(self):
        with self.assertRaises(ValueError):estimate(np.zeros((2,2,2)),[0,0],np.eye(2))
        with self.assertRaises(ValueError):estimate(np.ones((2,2,2)),[0,0],np.eye(2))
        with self.assertRaises(ValueError):estimate(np.ones((1,2,2)),[0,0],np.eye(2))
        with self.assertRaises(ValueError):estimate(np.ones((2,2,2)),[0,0],-np.eye(2))


if __name__=='__main__':unittest.main()
