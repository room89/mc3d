from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from diagnose_plateau_time import validate_bounds, summarize


class PlateauTimeTest(unittest.TestCase):
    def setUp(self):
        self.bounds = np.column_stack((np.arange(64), np.arange(1, 65)))*.002
        self.sources = [dict(group='a', seed=1), dict(group='b', seed=2)]

    def test_time_conversion_and_offsets(self):
        y = np.array([np.arange(64)*.01+2, np.arange(64)*.01-3])
        r = summarize(y, self.bounds, self.sources)
        self.assertAlmostEqual(r['half_change'], .32)
        self.assertAlmostEqual(r['groups']['a']['mean']-r['groups']['b']['mean'], 5)
        self.assertAlmostEqual(r['block_midpoints_s'][0], .008)
        self.assertLess(r['normalized_time_vs_index_max_error'], 1e-14)
        self.assertEqual(r['positive_group_half_changes'], 2)
        self.assertAlmostEqual(r['leave_one_group_out_half_changes']['a'], .32)

    def test_group_influence(self):
        y = np.zeros((2,64)); y[0,32:] = 1
        r = summarize(y, self.bounds, self.sources)
        self.assertEqual(r['half_change'], .5)
        self.assertEqual(r['leave_one_group_out_half_changes'], {'a': 0., 'b': 1.})
        self.assertEqual(r['groups']['a']['eight_window_means'], [0.]*4+[1.]*4)

    def test_invalid_schedule_and_series(self):
        for a in (self.bounds[::-1], self.bounds[:32], np.full((64,2), np.nan)):
            with self.assertRaises(ValueError): validate_bounds(a)
        b = self.bounds.copy(); b[10] += .001
        with self.assertRaises(ValueError): validate_bounds(b)
        with self.assertRaises(ValueError): summarize(np.zeros((3,64)), self.bounds, self.sources)
        with self.assertRaises(ValueError): summarize(np.full((2,64),np.nan), self.bounds, self.sources)
