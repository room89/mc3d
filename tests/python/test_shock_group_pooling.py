import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_group_pooling import pool_points

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
