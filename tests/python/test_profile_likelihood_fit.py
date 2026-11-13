from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from shock_profile_likelihood_fit import KNOTS, design, fit, project, solve_quadratic
from shock_hierarchical_likelihood import marginal_block_likelihood


class ProfileFitTest(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(-10, 10, 41)
        self.profile = 1/(1+np.exp(-KNOTS)); self.profile[[0, -1]] = [0, 1]
        self.T = np.array([[.02, .009], [.009, .02]])
        self.S = np.eye(41)

    def observations(self, shifts):
        return np.array([[design(self.x, KNOTS, s)@self.profile for s in row] for row in shifts])

    def test_projection_and_dense_quadratic(self):
        np.testing.assert_allclose(project([-.2, .8, .4, 1.2]), [0, .6, .6, 1])
        h = np.array([[3., .7], [.7, 2.]])
        exact = np.array([.2, .7]);g = h@exact
        actual, residual, _, ok = solve_quadratic(h, g, [.9, .1])
        self.assertTrue(ok);self.assertLessEqual(residual, 1e-8)
        np.testing.assert_allclose(actual, np.linalg.solve(h, g), atol=2e-8)

    def test_constrained_quadratic_has_weighted_pooled_solution(self):
        # The unconstrained solution (.8, .2) violates monotonicity. The exact
        # active-set solution pools both coordinates at (3*.8 + 1*.2)/4.
        result, _, _, ok = solve_quadratic(np.diag([3., 1.]), [2.4, .2], [0., 1.])
        self.assertTrue(ok)
        np.testing.assert_allclose(result, [.65, .65], atol=2e-8)

    def test_zero_weight_configurations_do_not_identify_missing_nodes(self):
        x = np.linspace(-10, -5, 11)
        y = np.zeros((2, 2, len(x)))
        configs = [[0., 0.], [-10., -10.]]
        with self.assertRaises(ValueError):
            fit(y, x, self.T, np.eye(len(x)), configs, [1., 0.])

    def test_zero_and_known_shifts_recover_profile(self):
        for shifts in (np.zeros((3, 2)), np.array([[-.5, .5], [0., -.5], [.5, 0.]])):
            y = self.observations(shifts)
            result = fit(y, self.x, self.T, self.S, shifts, np.eye(3))
            self.assertEqual(result['status'], 'converged')
            np.testing.assert_allclose(result['profile'], self.profile, atol=1e-6)
            for start in result['starts']:
                self.assertTrue(np.all(np.diff(start['history']) <= 1e-8))
            v = np.array(result['profile'])
            self.assertTrue(np.all(np.diff(v) >= 0));self.assertEqual(v[0], 0);self.assertEqual(v[-1], 1)

    def test_likelihood_matches_existing_joint_kernel(self):
        configs = np.array([[-.05, .05], [.05, -.05]])
        y = self.observations(configs)
        result = fit(y, self.x, self.T, self.S, configs, [.5, .5])
        best = min([s for s in result['starts'] if s['status']=='converged'], key=lambda s:s['objective'])
        predictions = np.array([[design(self.x, KNOTS, s)@best['profile'] for s in row] for row in configs])
        expected = -np.mean([marginal_block_likelihood(row,predictions,self.T,self.S,[.5,.5])['log_marginal_likelihood'] for row in y])
        cov = np.kron(self.T,self.S)
        constant = .5*(len(cov)*np.log(2*np.pi)+np.linalg.slogdet(cov)[1])
        self.assertAlmostEqual(best['objective']+constant, expected, places=9)
        permuted = fit(y[::-1],self.x,self.T,self.S,configs,[.5,.5])
        np.testing.assert_allclose(result['profile'],permuted['profile'],atol=1e-8)

    def test_covariance_support_and_limits(self):
        y = self.observations(np.zeros((2,2)))
        with self.assertRaises(ValueError):fit(y,self.x,-self.T,-self.S,[[0,0]],[1])
        with self.assertRaises(ValueError):fit(y,self.x,np.zeros((2,2)),self.S,[[0,0]],[1])
        with self.assertRaises(ValueError):fit(y[:,:,:4],self.x[:4],self.T,np.eye(4),[[0,0]],[1])
        with self.assertRaises(ValueError):fit(y,self.x,self.T,self.S,[[0,0]],[0])
        result=fit(y,self.x,self.T,self.S,[[0,0]],[1],max_iterations=1)
        self.assertNotEqual(result['status'],'converged')


if __name__ == '__main__':unittest.main()
