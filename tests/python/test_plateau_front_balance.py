from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from diagnose_plateau_front_balance import particle_balance, correlation, associations


class FrontBalanceTest(unittest.TestCase):
    def fixture(self):
        bounds=np.column_stack([np.arange(64),np.arange(1,65)])
        text='time;end_time;N_before;N_after;N_in;N_out;Rcount\n'
        text+=''.join(f'{i};{i+1};{100+i};{101+i};2;1;0\n' for i in range(64))
        return text,bounds

    def test_exact_particle_flux_and_population(self):
        text,bounds=self.fixture();r=particle_balance(text,bounds)
        np.testing.assert_array_equal(r['population'],np.arange(100,164))
        self.assertEqual(r['net_particles_by_step_midpoint'],[1.]*64)
        self.assertEqual(r['checked_steps'],64)

    def test_balance_and_coverage_rejected(self):
        text,bounds=self.fixture()
        with self.assertRaises(ValueError):particle_balance(text.replace('100;101;2;1;0','100;102;2;1;0'),bounds)
        with self.assertRaises(ValueError):particle_balance(text,bounds+1)
        with self.assertRaises(ValueError):particle_balance(text.replace('0;1;100','0;0;100'),bounds)

    def test_fractional_step_population_and_discrete_flux(self):
        text,bounds=self.fixture();bounds=np.column_stack([np.arange(64)*.5+.25,np.arange(64)*.5+.75])
        r=particle_balance(text,bounds)
        self.assertAlmostEqual(r['population'][1],100.5)
        self.assertEqual(r['net_particles_by_step_midpoint'][:2],[1.,0.])

    def test_associations_and_constant_series(self):
        t=np.arange(64,dtype=float);y=np.array([t+10,t-10]);z=-2*y
        s=[dict(group='a',seed=1),dict(group='b',seed=2)]
        r=associations(y,z,s)
        self.assertAlmostEqual(r['ensemble_correlation'],-1)
        self.assertAlmostEqual(r['group_curves']['a']['target_half_change'],32)
        self.assertIsNone(correlation(np.ones(64),t))
        with self.assertRaises(ValueError):correlation(t,np.full(64,np.nan))
