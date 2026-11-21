import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from audit_shock_validation import at,comparison,gate
from analyze_shock_validation import compare, sparta_profile


class AuditTest(unittest.TestCase):
    def test_gate_distinguishes_precision_from_mean_difference(self):
        g=gate(.004,.02,.01)
        self.assertFalse(g['mean_outside_tolerance'])
        self.assertFalse(g['precision_pass'])
        self.assertFalse(g['passed'])
        self.assertTrue(gate(.002,.003,.01)['passed'])
        self.assertTrue(gate(.02,.001,.01)['mean_outside_tolerance'])

    def test_diagnostics_match_original_gate(self):
        a=dict(density=[0.,.5,1.],density_uncertainty=[.001]*3,
               width_10_90=dict(mean=4.6,uncertainty=.01),
               gradient_thickness=dict(mean=4.,uncertainty=.01),processing_uncertainty=.2)
        b=a|{'density':[.001,.52,1.001],'processing_uncertainty':.1}
        audit=comparison(a,b);original=compare(a,b)
        self.assertEqual(audit['density']['means_outside'],1)
        self.assertEqual(audit['density']['passed']==3,original['checks']['density'])
        for key in ('width_10_90','gradient_thickness'):
            self.assertEqual(audit[key]['passed'],original['checks'][key])
        self.assertAlmostEqual(audit['gradient_thickness']['processing_uncertainty'],.3)
        self.assertTrue(audit['width_10_90']['passed'])
        self.assertFalse(audit['gradient_thickness']['passed'])

    def test_interpolation_does_not_extrapolate(self):
        for x,y in [(0,2),(.5,3),(1,4),(2,6)]:
            self.assertEqual(at([0,1,2],[2,4,6],x),y)
        with self.assertRaises(ValueError):at([0,1],[0,1],2)

    def test_sparta_reads_all_windows_for_each_schedule_format(self):
        step = 3
        cells = 'ITEM: CELLS xlo xhi ylo yhi zlo zhi f_avg[1] f_avg[2] f_avg[3] f_avg[4] f_avg[8]\n'
        for i in range(9):
            cells += f'{i / 8} {(i + 1) / 8} 0 1 0 1 {1 + i / 8} 0 0 0 1\n'
        schedules = (
            ({}, 32 * step, 64),
            ({'warmup_windows': 160}, 160 * step, 64),
            ({'burn_steps': 120 * step, 'measurement_steps': 80 * step}, 120 * step, 80),
        )
        for schedule, burn, count in schedules:
            with self.subTest(schedule=schedule), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                ref = {
                    'quick': False,
                    'steps_per_window': step,
                    'source': {'numerics': {'cells_x': 9}, 'statistics_windows_s': [[i, i + 1] for i in range(count)]},
                    'reference': {'left': {'number_density_m3': 1}, 'right': {'number_density_m3': 2}, 'mean_free_path_m': 1},
                    **schedule,
                }
                (folder / 'reference.json').write_text(json.dumps(ref))
                for i in range(1, count + 1):
                    (folder / f'profile.{burn + i * step}').write_text(cells)
                self.assertEqual([window['window'] for window in sparta_profile(folder)], list(range(count)))
