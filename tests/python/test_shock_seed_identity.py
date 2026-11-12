import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import main


class ShockSeedIdentityTest(unittest.TestCase):
    def test_selected_case_is_explicit_and_cannot_expand_on_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);binary=root/'binary';binary.write_bytes(b'preparation-only')
            output=root/'matrix'
            args=['shock_validation','--binary',str(binary),'--sparta-bin',str(binary),
                  '--output-dir',str(output),'--prepare-only','--variants','C']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            report=output/'validation.json';before=report.read_bytes();data=json.loads(before)
            self.assertEqual(data['variants'],['C'])
            self.assertEqual(set(data['runs']),{'C-123','C-456','C-789'})
            self.assertEqual(set(data['sparta_runs']),{'C-123','C-456','C-789'})
            self.assertEqual(data['convergence']['status'],'incomplete')
            with patch.object(sys,'argv',args[:-1]+['D','--resume']):
                with self.assertRaisesRegex(ValueError,'different variants'):main()
            self.assertEqual(before,report.read_bytes())

    def test_custom_seed_preparation_and_resume_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);binary=root/'binary';binary.write_bytes(b'preparation-only')
            output=root/'matrix'
            args=['shock_validation','--binary',str(binary),'--output-dir',str(output),
                  '--prepare-only','--no-time-budget','--seeds','10003','10004','10005']
            with patch.object(sys,'argv',args):self.assertEqual(main(),2)
            report=output/'validation.json';before=report.read_bytes()
            self.assertEqual(json.loads(before)['seeds'],[10003,10004,10005])
            for seed in (10003,10004,10005):
                cfg=json.loads((output/f'C-{seed}'/'solver.json').read_text())
                self.assertEqual(cfg['seed'],seed)
            self.assertFalse((output/'C-123').exists())
            with patch.object(sys,'argv',args[:-3]+['10006','10007','10008','--resume']):
                with self.assertRaisesRegex(ValueError,'different seeds'):main()
            self.assertEqual(before,report.read_bytes())

    def test_duplicate_seeds_rejected_before_creation(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'matrix'
            with patch.object(sys,'argv',['shock_validation','--binary','unused','--output-dir',str(output),'--seeds','1','1','2']):
                with self.assertRaisesRegex(ValueError,'distinct'):main()
            self.assertFalse(output.exists())

if __name__=='__main__':unittest.main()
