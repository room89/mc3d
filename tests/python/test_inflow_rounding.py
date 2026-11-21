import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import configuration,main
from physical_config import prepare

class InflowRoundingTest(unittest.TestCase):
    def test_only_mode_changes(self):
        a,ra=configuration('E',210001,20,2,4)
        b,rb=configuration('E',210001,20,2,4,'stochastic')
        self.assertEqual(ra,rb)
        ca,ua=prepare(a);cb,ub=prepare(b)
        self.assertEqual(ua['scales'],ub['scales'])
        self.assertEqual(cb.pop('hyperfree_inflow_rounding'),'stochastic')
        self.assertEqual(ca,cb)
        b['numerics'].pop('hyperfree_inflow_rounding')
        self.assertEqual(a,b)
    def test_invalid_mode(self):
        for mode in ('unknown',True,None):
            with self.assertRaises(ValueError):configuration('E',123,hyperfree_inflow_rounding=mode)
            data,_=configuration('E',123)
            data['numerics']['hyperfree_inflow_rounding']=mode
            with self.assertRaises(ValueError):prepare(data)
    def test_resume_rejects_mode_change_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binary=root/'binary';binary.write_bytes(b'prepare only');out=root/'case'
            args=['runner','--binary',str(binary),'--output-dir',str(out),'--prepare-only','--variants','E','--hyperfree-inflow-rounding','stochastic']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            cfg=json.loads((out/'E-123'/'solver.json').read_text())
            self.assertEqual(cfg['hyperfree_inflow_rounding'],'stochastic')
            before={str(p):p.read_bytes() for p in out.rglob('*') if p.is_file()}
            with patch.object(sys,'argv',args[:-1]+['accumulated','--resume']):
                with self.assertRaisesRegex(ValueError,'different inflow rounding'):main()
            self.assertEqual(before,{str(p):p.read_bytes() for p in out.rglob('*') if p.is_file()})
