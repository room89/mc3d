import tempfile
import unittest
from pathlib import Path
from benchmark_solver import parse_log, summarize, particle_counts


class BenchmarkTest(unittest.TestCase):
    def test_old_and_new_logs_and_missing_phase(self):
        old = 'text=Time step:1\ntext=dt = 0.03 t = 0\nstopwatch:ClusterStep;total_time:0.02\n'
        old += 'text=Time step:2\ntext=dt = 0.04 t = 0.03\nstopwatch:ClusterStep;total_time:0.04\n'
        result = parse_log(old)
        self.assertEqual(result['steps'], 2)
        self.assertAlmostEqual(result['model_time'], .07)
        self.assertAlmostEqual(result['cluster_total_s'], .06)
        self.assertEqual(result['phase_ms_per_step']['ClusterStep'], 30)
        self.assertNotIn('TimeStep::SyncDt', result['phase_ms_per_step'])
        newer = parse_log(old + 'stopwatch:TimeStep::SyncDt;total_time:1e-3\n')
        self.assertEqual(newer['phase_ms_per_step']['TimeStep::SyncDt'], .5)
        summary = summarize([dict(result, elapsed_s=.1, peak_rss_kib=2048), dict(newer, elapsed_s=.3, peak_rss_kib=4096)])
        self.assertEqual(summary['median_elapsed_s'], .2)
        self.assertEqual(summary['median_peak_rss_kib'], 3072)
        self.assertEqual(summary['range_steps'], [2, 2])
        self.assertIsNone(summary['median_phase_ms_per_step']['TimeStep::SyncDt'])

    def test_bad_logs_are_not_silent_successes(self):
        for log in ('', 'text=Time step:1', 'stopwatch:ClusterStep;total_time:0.1'):
            with self.assertRaises(ValueError):
                parse_log(log)

    def test_counts_use_actual_population(self):
        with tempfile.TemporaryDirectory() as path:
            p = Path(path)
            (p/'initial.csv').write_text('x;y;z;N;ro\n0;0;0;3;1\n1;0;0;7;1\n')
            (p/'final.net').write_text('cluster\n\tcell\n\t\tparticles\t4\n\tcell\n\t\tparticles\t8\n')
            self.assertEqual(particle_counts(p), {'initial_particles': 10, 'final_particles': 12})


if __name__ == '__main__':
    unittest.main()
