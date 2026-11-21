import tempfile
import json
from unittest.mock import patch
import unittest
from pathlib import Path
from benchmark_solver import parse_log, summarize, particle_counts, main


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

    def test_seeds_are_separate_and_warmups_diagnostics_are_excluded(self):
        calls = []

        def fake_run(binary, directory, cfg, rss=False, diagnostic=False):
            calls.append((directory.name, cfg['seed']))
            result = dict(steps=2, model_time=.05, cluster_total_s=.02,
                          elapsed_s=99 if 'warmup' in directory.name or diagnostic else .03,
                          peak_rss_kib=2048, phase_ms_per_step={'ClusterStep': 10})
            if diagnostic:
                result.update(initial_particles=8, final_particles=9)
            return result

        with tempfile.TemporaryDirectory() as output:
            args = ['benchmark_solver.py', '--binary', '/current', '--baseline', '/baseline',
                    '--out', output, '--geometries', 'none', '--threads', '1',
                    '--seeds', '123', '456', '--repeats', '2', '--warmups', '1', '--diagnostics']
            with patch('sys.argv', args), patch('benchmark_solver.run', side_effect=fake_run), patch('builtins.print'):
                main()
            rows = json.loads((Path(output)/'summary.json').read_text())
        self.assertEqual(len(calls), 16)
        self.assertEqual(len(rows), 4)
        self.assertEqual({row['seed'] for row in rows}, {123, 456})
        for row in rows:
            self.assertEqual(len(row['runs']), 2)
            self.assertEqual(row['median_elapsed_s'], .03)
            self.assertEqual(row['diagnostic']['final_particles'], 9)

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
