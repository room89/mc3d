from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_hierarchical_likelihood import covariance_factors, marginal_block_likelihood
from shock_shift_likelihood import marginal_shift_likelihood


class HierarchicalLikelihoodTest(unittest.TestCase):
    def test_covariance_has_persistent_and_stationary_parts(self):
        x=np.arange(5)*.1;scale=np.arange(5)+1.
        time,space=covariance_factors(x,scale,4,.1,.2,.8,0.)
        np.testing.assert_allclose(np.diag(time),.05)
        self.assertAlmostEqual(time[0,3],.01+.04*.8**3)
        np.testing.assert_array_equal(space,np.diag(scale**2))
        _,smooth=covariance_factors(x,scale,4,.1,.2,.8,.2)
        np.testing.assert_allclose(np.diag(smooth),scale**2,rtol=1e-14)
        self.assertGreater(smooth[0,1],0)

    def test_joint_factorization_equals_dense_kronecker_likelihood(self):
        time=np.array([[2.,.7],[.7,1.]])
        space=np.array([[1.,.2,.1],[.2,2.,.4],[.1,.4,1.5]])
        y=np.arange(6).reshape(2,3)/5
        predictions=np.array([np.zeros((2,3)),np.ones((2,3))*.5,np.ones((2,3))])
        weights=[.2,.3,.5]
        joint=marginal_block_likelihood(y,predictions,time,space,weights)
        dense=marginal_shift_likelihood(y.ravel(),predictions.reshape(3,6),np.kron(time,space),[-1.,0.,1.],weights)
        self.assertAlmostEqual(joint['log_marginal_likelihood'],dense['log_marginal_likelihood'],places=12)
        np.testing.assert_allclose(joint['posterior_weights'],dense['posterior_weights'],rtol=1e-13)

    def test_independent_time_limit_and_correlated_counterexample(self):
        y=np.array([[1.],[1.]])
        predictions=np.zeros((1,2,1))
        independent=marginal_block_likelihood(y,predictions,np.eye(2),[[1.]],[1.])
        one=marginal_shift_likelihood([1.],[[0.]],[[1.]],[0.],[1.])
        self.assertAlmostEqual(independent['log_marginal_likelihood'],2*one['log_marginal_likelihood'])
        correlated=marginal_block_likelihood(y,predictions,[[1.,.9],[.9,1.]],[[1.]],[1.])
        self.assertGreater(abs(correlated['log_marginal_likelihood']-independent['log_marginal_likelihood']),.1)

    def test_degenerate_covariance_is_not_silently_regularized(self):
        time,space=covariance_factors(np.arange(4),np.ones(4),2,1.,0.,.8,0.)
        with self.assertRaises(ValueError):marginal_block_likelihood(np.zeros((2,4)),np.zeros((1,2,4)),time,space,[1.])
        with self.assertRaises(ValueError):covariance_factors(np.arange(4),np.ones(4),2,.1,.2,1.,0.)


if __name__ == '__main__':unittest.main()
