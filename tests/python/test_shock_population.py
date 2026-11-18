import json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from unittest.mock import patch
from shock_validation import configuration,main
from physical_config import prepare
from sparta_shock import export
class PopulationTest(unittest.TestCase):
    def test_only_population_changes(self):
        a,ra=configuration('E',210001,20,2)
        b,rb=configuration('E',210001,20,2,4)
        self.assertEqual(ra,rb)
        _,ua=prepare(a);_,ub=prepare(b)
        self.assertAlmostEqual(ua['scales']['molecules_per_model_particle']/4,ub['scales']['molecules_per_model_particle'])
        b['numerics']['particles_per_cell']//=4
        self.assertEqual(a,b)
    def test_sparta_weight_and_time(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a';b=Path(d)/'b'
            export(a,'E',210001,warmup_transits=20,domain_multiplier=2)
            export(b,'E',210001,warmup_transits=20,domain_multiplier=2,population_multiplier=4)
            x=json.loads((a/'reference.json').read_text());y=json.loads((b/'reference.json').read_text())
            self.assertEqual(x['fnum']/4,y['fnum'])
            for k in ('timestep','steps_per_window','warmup_windows'):self.assertEqual(x[k],y[k])
    def test_invalid_multiplier(self):
        for v in (0,2,1.5,True):
            with self.assertRaises(ValueError):configuration('E',123,population_multiplier=v)

    def test_resume_rejects_population_change(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);binary=p/'binary';binary.write_bytes(b'prepare only');out=p/'matrix'
            args=['runner','--binary',str(binary),'--output-dir',str(out),'--prepare-only','--variants','E','--population-multiplier','4']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            before={str(f):f.read_bytes() for f in out.rglob('*') if f.is_file()}
            with patch.object(sys,'argv',args[:-1]+['1','--resume']):
                with self.assertRaisesRegex(ValueError,'different population'):main()
            self.assertEqual(before,{str(f):f.read_bytes() for f in out.rglob('*') if f.is_file()})
