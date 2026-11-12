import sys
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_group_pooling import pool_points,crossfit_pool_points
from plan_shock_pooling import assignments

class PoolingTest(unittest.TestCase):
    def test_average_slopes_before_inversion(self):
        r=pool_points([dict(density=[0,1],slopes=[1,2,3]),dict(density=[.2,.8],slopes=[3,4,5])])
        self.assertEqual(r['slopes'],[2,3,4])
        self.assertAlmostEqual(1/r['slopes'][1],1/3)
        self.assertNotAlmostEqual(1/r['slopes'][1],(1/2+1/4)/2)
    def test_reject_mixed_or_invalid_inputs(self):
        a=dict(density=[0,1],slopes=[1,2,3])
        for points in [[],[a,dict(a,slope_cubic=2)],[a,dict(density=[0],slopes=[1,2,3])],[dict(a,density=[float('nan'),1])]]:
            with self.assertRaises(ValueError):pool_points(points)
    def test_identity_for_one_subgroup(self):
        a=dict(density=[0,1],slopes=[1,2,3],density_cubic=[0,1],slope_cubic=2)
        self.assertEqual(pool_points([a]),a)

    def test_physical_seed_plan_is_disjoint_and_preserves_baseline(self):
        old=[[123,456,789]]+[[10000+3*i+j for j in range(3)] for i in range(1,12)]
        groups=assignments(old)
        self.assertEqual([g[0] for g in groups],old)
        self.assertEqual(len({s for g in groups for triplet in g for s in triplet}),216)
        self.assertEqual(groups[0][1],[20000,20001,20002])
        self.assertEqual(groups[-1][-1],[20177,20178,20179])
        with self.assertRaises(ValueError):assignments(old,10000)
        with self.assertRaises(ValueError):assignments(old,2**31-1)
        with self.assertRaises(ValueError):assignments(old[:11])

    def test_pooled_peak_excludes_heldout_profile(self):
        from shock_ensemble_estimator import PROFILE_GRID,centered_density
        from shock_crossfit_gradient import peak_position
        points=[]
        for i in range(6):
            y=1/(1+np.exp(-PROFILE_GRID/(1+.02*i)))
            points.append(dict(density=y.tolist(),density_cubic=y.tolist(),slopes=[.25]*3,slope_cubic=.25))
        seen=[]
        def record(x,y,width):
            seen.append(np.array(y));return peak_position(x,y,width)
        with patch('shock_crossfit_gradient.peak_position',side_effect=record):
            result=crossfit_pool_points(points)
        expected=centered_density([p['density'] for p in points])
        for i in range(6):
            np.testing.assert_allclose(seen[i],np.delete(expected,i,axis=0).mean(0))
        self.assertEqual(result['density'],pool_points(points)['density'])
        self.assertTrue(all(s>0 for s in result['slopes']))
        with self.assertRaises(ValueError):crossfit_pool_points(points[:5])
