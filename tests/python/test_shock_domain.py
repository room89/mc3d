import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import configuration,main
from physical_config import prepare
from sparta_shock import export
class DomainTest(unittest.TestCase):
    def test_physical_resolution_weight_and_time_preserved(self):
        a,ra=configuration('E',210001,20);b,rb=configuration('E',210001,20,2)
        self.assertEqual(ra,rb)
        self.assertEqual(b['domain_m'][0],2*a['domain_m'][0])
        self.assertEqual(b['numerics']['cells_x'],2*a['numerics']['cells_x'])
        self.assertEqual(a['statistics_windows_s'],b['statistics_windows_s'])
        self.assertEqual(a['end_time_s'],b['end_time_s'])
        _,ua=prepare(a);_,ub=prepare(b)
        self.assertAlmostEqual(ua['scales']['molecules_per_model_particle'],ub['scales']['molecules_per_model_particle'])
        b['domain_m'][0]/=2;b['numerics']['cells_x']//=2
        self.assertEqual(a,b)
    def test_sparta_step_weight_and_schedule_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a';b=Path(d)/'b'
            export(a,'E',210001,warmup_transits=20);export(b,'E',210001,warmup_transits=20,domain_multiplier=2)
            x=json.loads((a/'reference.json').read_text());y=json.loads((b/'reference.json').read_text())
            for k in ('fnum','timestep','steps_per_window','warmup_windows'):self.assertEqual(x[k],y[k])
            self.assertIn('create_grid 1600 1 1',(b/'in.shock').read_text())
    def test_resume_rejects_domain_change(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);binary=p/'binary';binary.write_bytes(b'prepare only');out=p/'matrix'
            args=['runner','--binary',str(binary),'--output-dir',str(out),'--prepare-only','--variants','E','--domain-multiplier','2']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            before={str(f):f.read_bytes() for f in out.rglob('*') if f.is_file()}
            with patch.object(sys,'argv',args[:-1]+['1','--resume']):
                with self.assertRaisesRegex(ValueError,'different domain'):main()
            self.assertEqual(before,{str(f):f.read_bytes() for f in out.rglob('*') if f.is_file()})
    def test_invalid_multiplier(self):
        for v in (0,3,1.5,True):
            with self.assertRaises(ValueError):configuration('E',123,20,v)
