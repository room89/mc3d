import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2]/'tools/profile_simulation.py'
spec = importlib.util.spec_from_file_location('profile_simulation', SCRIPT)
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


class Parsers(unittest.TestCase):
    def test_ps_cpu_formats(self):
        for text, expected in [('00:02.50', 2.5), ('1:02:03', 3723), ('2-01:02:03', 176523)]:
            self.assertEqual(p.parse_ps('100 '+text), {'rss_bytes': 102400, 'cpu_seconds': expected})
        with self.assertRaises(ValueError):
            p.parse_ps('')

    def test_mac_page_sizes_and_missing_fields(self):
        for size in [4096, 16384]:
            text = f'Mach Virtual Memory Statistics: (page size of {size} bytes)\nPages free: 12.\nSwapins: 3.\nPages occupied by compressor: 5.\n'
            self.assertEqual(p.parse_vm_stat(text), {'free_bytes': 12*size, 'swap_in_bytes': 3*size, 'compressed_bytes': 5*size})
        self.assertEqual(p.parse_swap('vm.swapusage: total = 1.00G used = 0.50G free = 0.50G')['swap_used_bytes'], 2**29)
        with self.assertRaises(ValueError):
            p.parse_swap('not permitted')

    def test_linux_stat_with_unusual_name(self):
        fields = ['0']*22
        fields[0]='S'; fields[7]='7'; fields[9]='2'
        fields[11]='120'; fields[12]='30'; fields[21]='10'
        self.assertEqual(p.parse_linux_stat('123 (a (b) c) '+' '.join(fields), 100, 4096),
                         dict(rss_bytes=40960, cpu_seconds=1.5, minor_faults=7, major_faults=2))

    def test_linux_system(self):
        self.assertEqual(p.parse_meminfo('MemAvailable: 12 kB\nSwapTotal: 20 kB\nSwapFree: 3 kB'),
                         dict(available_bytes=12288, swap_used_bytes=17408))
        self.assertEqual(p.parse_linux_vmstat('pswpin 4\npswpout 2', 4096), dict(swap_in_bytes=16384, swap_out_bytes=8192))
        self.assertEqual(p.parse_psi('some avg10=0.00 total=1000000\nfull avg10=0.00 total=200000'),
                         dict(psi_some_seconds=1., psi_full_seconds=.2))
        self.assertEqual(p.parse_meminfo('MemFree: 1 kB'), dict(free_bytes=1024))

    def test_unavailable_commands_do_not_fail_collection(self):
        with patch.object(p.platform, 'system', return_value='Darwin'), patch.object(p.subprocess, 'run', side_effect=subprocess.TimeoutExpired('ps', 1)):
            sample = p.Sampler(1).sample()
        self.assertIsNone(sample['process']['rss_bytes'])
        self.assertEqual(set(sample['errors']), {'ps', 'vm_stat', 'swapusage'})
        self.assertIn('process.rss_bytes', sample['unavailable'])

    def test_counter_reset_and_missing_values(self):
        a={'process': {'cpu_seconds': 3}, 'system': {'swap_in_bytes': 10}}
        b={'process': {'cpu_seconds': 4}, 'system': {'swap_in_bytes': 2}}
        d,r=p.counter_deltas(a,b)
        self.assertEqual(d['process.cpu_seconds'],1)
        self.assertIsNone(d['system.swap_in_bytes'])
        self.assertIsNone(d['process.major_faults'])
        self.assertEqual(r,['system.swap_in_bytes'])

    def test_log_events(self):
        self.assertEqual(p.parse_event('tskv\ttext=Time step:12',3)['step'],12)
        e=p.parse_event('tskv\ttext=stopwatch:ClusterStep;total_time:1.5',3)
        self.assertEqual(e['estimated_start'],1.5)
        self.assertIsNone(p.parse_event('text=stopwatch:ClusterStep;total_time:1e999',3))
        self.assertIsNone(p.parse_event('unrelated',3))


class Integration(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out=Path(self.temp.name)/'result'

    def command(self, code):
        return [sys.executable,str(SCRIPT),'--out',str(self.out),'--interval','.05','--',sys.executable,'-u','-c',code]

    def test_exit_and_large_partial_log(self):
        code="import sys,time; print('x'*200000); print('tskv\\ttext=Time step:1'); time.sleep(.2); print('tskv\\ttext=stopwatch:ClusterStep;total_time:0.2'); sys.stdout.write('text=Time step:2'); sys.exit(7)"
        done=subprocess.run(self.command(code),capture_output=True,timeout=15)
        self.assertEqual(done.returncode,7,done.stderr)
        result=json.loads((self.out/'summary.json').read_text())
        self.assertEqual(result['completed_steps'],1)
        self.assertEqual(len(result['steps']),1)
        self.assertTrue((self.out/'run.log').read_bytes().endswith(b'text=Time step:2'))
        self.assertGreater((self.out/'run.log').stat().st_size,200000)
        self.assertEqual(json.loads((self.out/'metadata.json').read_text())['child_returncode'],7)

    def test_interrupt_saves_partial_report_and_stops_child(self):
        proc=subprocess.Popen(self.command("import time; print('text=Time step:1'); time.sleep(30)"),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            deadline=time.monotonic()+5
            while time.monotonic()<deadline:
                if (self.out/'events.jsonl').exists() and (self.out/'events.jsonl').stat().st_size:
                    break
                time.sleep(.02)
            proc.send_signal(signal.SIGINT)
            stdout,stderr=proc.communicate(timeout=10)
            self.assertEqual(proc.returncode,130,stderr)
            meta=json.loads((self.out/'metadata.json').read_text())
            self.assertTrue(meta['interrupted'])
            with self.assertRaises(ProcessLookupError):
                os.kill(meta['pid'],0)
            self.assertEqual(json.loads((self.out/'summary.json').read_text())['incomplete_steps'],[1])
        finally:
            if proc.poll() is None:
                proc.kill(); proc.wait()

    def test_start_failure_and_existing_directory(self):
        done=subprocess.run([sys.executable,str(SCRIPT),'--out',str(self.out),'--','/no/such/binary'],capture_output=True,timeout=10)
        self.assertEqual(done.returncode,1)
        self.assertTrue(json.loads((self.out/'metadata.json').read_text())['collector_error'])
        done=subprocess.run(self.command('print(1)'),capture_output=True,timeout=10)
        self.assertEqual(done.returncode,2)

    def test_summary_gap_and_incomplete_step(self):
        self.out.mkdir()
        events=[dict(kind='phase',phase='init',duration_s=1,t=1),dict(kind='step',step=1,t=1.1),
                dict(kind='phase',phase='ClusterStep',duration_s=2,t=3.1),dict(kind='step',step=2,t=3.2)]
        samples=[]
        for t,counter in [(1,10),(4,15)]:
            samples.append(dict(t=t,process=dict.fromkeys(p.PROCESS_KEYS),system={**dict.fromkeys(p.SYSTEM_KEYS),'swap_in_bytes':counter}))
        (self.out/'events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events))
        (self.out/'samples.jsonl').write_text(''.join(json.dumps(s)+'\n' for s in samples))
        result=p.summarize(self.out,1)
        self.assertEqual(result['incomplete_steps'],[2])
        self.assertEqual(result['initialization_events'][0]['phase'],'init')
        metrics=result['slowest_steps'][0]['metrics']
        self.assertTrue(metrics['sampling_gap'])
        self.assertIsNone(metrics['max_observed_rss_bytes'])
        self.assertEqual(metrics['counter_deltas']['system.swap_in_bytes'],5)

    def test_counter_reset_inside_report_window(self):
        self.out.mkdir()
        events = [dict(kind='step', step=1, t=.1),
                  dict(kind='phase', phase='ClusterStep', duration_s=2, t=2.1)]
        samples=[]
        for t, value in [(0,10),(1,1),(3,20)]:
            samples.append(dict(t=t, process=dict.fromkeys(p.PROCESS_KEYS),
                                system={**dict.fromkeys(p.SYSTEM_KEYS), 'swap_in_bytes':value},
                                counter_resets=['system.swap_in_bytes'] if t==1 else []))
        (self.out/'events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events))
        (self.out/'samples.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in samples))
        result=p.summarize(self.out,1)
        self.assertIsNone(result['steps'][0]['metrics']['counter_deltas']['system.swap_in_bytes'])

    def test_signal_exit_code(self):
        done=subprocess.run(self.command('import os,signal; os.kill(os.getpid(),signal.SIGTERM)'),
                            capture_output=True,timeout=10)
        self.assertEqual(done.returncode,128+signal.SIGTERM)

    def test_collector_failure_terminates_child(self):
        with patch.object(p.Sampler,'sample',side_effect=RuntimeError('collector failed')):
            code=p.run([sys.executable,'-c','import time; time.sleep(30)'],self.out,.1)
        self.assertEqual(code,1)
        meta=json.loads((self.out/'metadata.json').read_text())
        self.assertIn('collector failed',meta['collector_error'])
        with self.assertRaises(ProcessLookupError):
            os.kill(meta['pid'],0)


if __name__=='__main__':
    unittest.main()
