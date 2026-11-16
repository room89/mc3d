import json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import configuration
from sparta_shock import export
class WarmupTest(unittest.TestCase):
    def test_only_observation_times_change(self):
        original,ref=configuration('E',123)
        late,_=configuration('E',123,20)
        tau=ref['transit_time_s']
        for a,b in zip(original['statistics_windows_s'],late['statistics_windows_s']):
            for x,y in zip(a,b):self.assertAlmostEqual((y-x)/tau,16)
        self.assertAlmostEqual(late['end_time_s']/tau,28)
        for k in original:
            if k not in ('end_time_s','statistics_windows_s'):self.assertEqual(original[k],late[k])
    def test_sparta_matches_late_schedule(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'case';export(p,'E',123,warmup_transits=20)
            r=json.loads((p/'reference.json').read_text());n=r['steps_per_window']
            self.assertEqual(r['warmup_windows'],160)
            deck=(p/'in.shock').read_text()
            self.assertIn(f'run {160*n}\n',deck);self.assertIn(f'run {64*n}\n',deck)
            self.assertAlmostEqual(224*n*r['timestep'],r['source']['end_time_s'])
    def test_invalid_warmup(self):
        for x in (-1,1.5,float('nan'),True):
            with self.assertRaises(ValueError):configuration('C',123,x)
    def test_resume_rejects_changed_schedule_without_touching_evidence(self):
        from unittest.mock import patch
        from shock_validation import main
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);binary=root/'binary';binary.write_bytes(b'prepare only')
            output=root/'matrix'
            args=['shock_validation','--binary',str(binary),'--output-dir',str(output),'--prepare-only','--variants','E','--warmup-transits','20']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            snapshot=lambda:{str(p.relative_to(output)):p.read_bytes() for p in output.rglob('*') if p.is_file()}
            before=snapshot()
            self.assertEqual(json.loads((output/'validation.json').read_text())['warmup_transits'],20)
            with patch.object(sys,'argv',args[:-1]+['4','--resume']):
                with self.assertRaisesRegex(ValueError,'different warmup'):main()
            self.assertEqual(before,snapshot())
            with patch.object(sys,'argv',args+['--resume']):self.assertEqual(main(),2)
