from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from compare_shock_experiment_sources import compare_sources,crossings

class SourceComparisonTest(unittest.TestCase):
    def test_fixed_scale_and_translation_do_not_change_width(self):
        normalized=[dict(x_paper=x,rho_normalized=y) for x,y in [(0,0),(1,.5),(2,1)]]
        dimensional=[dict(x_m=.002*x+.003,rho_normalized=y) for x,y in [(0,0),(1,.5),(2,1)]]
        r=compare_sources(normalized,dimensional,.002)
        self.assertAlmostEqual(r['normalized_source_width_m'],.0032)
        self.assertAlmostEqual(r['dimensional_to_normalized_width_relative_difference'],0)
        self.assertAlmostEqual(r['midpoint_translation_m'],.003)
        self.assertEqual(r['status'],'secondary_source_diagnostic_only')

    def test_missing_support_cannot_be_extrapolated(self):
        with self.assertRaisesRegex(ValueError,'support'):
            crossings([dict(x_m=x,rho_normalized=y) for x,y in [(0,.2),(1,.5),(2,.8)]],'x_m')
