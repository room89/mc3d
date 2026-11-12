from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from audit_shock_coverage import exact_interval,tail,case_counts

class CoverageAuditTest(unittest.TestCase):
    def test_exact_endpoints_and_symmetry(self):
        self.assertAlmostEqual(exact_interval(10,10)[0],.025**.1)
        self.assertEqual(exact_interval(10,10)[1],1.)
        self.assertAlmostEqual(exact_interval(0,10)[1],1-.025**.1)
        a=exact_interval(3,15);b=exact_interval(12,15)
        self.assertAlmostEqual(a[0],1-b[1]);self.assertAlmostEqual(a[1],1-b[0])
        wider=exact_interval(3,15,.05/144)
        self.assertLess(wider[0],a[0]);self.assertGreater(wider[1],a[1])

    def test_documented_binomial_tail(self):
        self.assertAlmostEqual(tail(15,3,.1),.18406106910639114)

    def test_rejections_count_as_noncoverage(self):
        c=dict(trial_count=2,runs=[dict(repetition=0,with_processing_covers=[True,True],density_covers=True),dict(repetition=1,error='ambiguous')])
        self.assertEqual(case_counts(c),(2,[1,1,1]))
        c['runs'][1]['repetition']=0
        with self.assertRaises(ValueError):case_counts(c)
