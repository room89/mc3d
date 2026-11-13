from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from fit_physical_plateau_mean import design,fit_region


class MeanModelTest(unittest.TestCase):
    def setUp(self):
        self.x=np.linspace(-12,-8,10);self.basis=design(self.x)
        self.beta=np.array([.2,.03,.01,.04]);self.mean=self.basis@self.beta
        self.train=np.broadcast_to(self.mean,(6,64,10)).copy()
        self.held=np.broadcast_to(self.mean+.005,(4,64,10)).copy()
        self.a=[dict(seed=i+1,group=str(i//2)) for i in range(6)]
        self.b=[dict(seed=i+101,group='h'+str(i//2)) for i in range(4)]

    def test_exact_recovery_and_preserved_heldout_offset(self):
        r=fit_region(self.train,self.held,self.x,self.a,self.b)
        np.testing.assert_allclose(r['coefficients'],self.beta,atol=1e-14)
        self.assertAlmostEqual(r['heldout_minus_training_mean'],.005)
        self.assertAlmostEqual(r['summaries']['heldout']['mean']['mean'],.005)
        self.assertAlmostEqual(r['summaries']['heldout']['time_change']['mean'],0.)
        changed=fit_region(self.train,self.held+.03,self.x,self.a,self.b)
        np.testing.assert_array_equal(changed['coefficients'],r['coefficients'])
        np.testing.assert_array_equal(changed['jackknife_coefficient_covariance'],r['jackknife_coefficient_covariance'])

    def test_leave_group_out_does_not_fit_excluded_offset(self):
        self.train[:2]+=.06
        r=fit_region(self.train,self.held,self.x,self.a,self.b)
        # Group 0's own offset cannot enter its excluded-group mean prediction.
        np.testing.assert_allclose(r['leave_group_out_coefficients'][0],self.beta,atol=1e-14)
        np.testing.assert_allclose(r['summaries']['training_leave_group_out']['per_run_mean'][:2],[.06,.06],atol=1e-14)
        self.assertGreater(r['jackknife_coefficient_covariance'][0][0],0.)
        self.assertGreater(r['predicted_mean_se_range'][1],0.)
        uncertainty=r['transfer_with_training_mean_uncertainty']['mean']
        self.assertGreater(uncertainty['combined_descriptive_se'],uncertainty['heldout_group_mcse'])

    def test_overlap_invalid_grid_and_nonfinite_data(self):
        with self.assertRaises(ValueError):fit_region(self.train,self.train,self.x,self.a,self.a)
        with self.assertRaises(ValueError):design(self.x[::-1])
        self.held[0,0,0]=np.nan
        with self.assertRaises(ValueError):fit_region(self.train,self.held,self.x,self.a,self.b)


if __name__=='__main__':unittest.main()
