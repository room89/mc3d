from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
import analyze_order_subgroup as module
import test_alsmeyer_subgroup as fixtures
from shock_order_sensitivity import METHOD
from shock_validation import digest,save


class OrderReferenceSubgroupTest(unittest.TestCase):
    def test_phases_retain_every_source_window_and_both_orders(self):
        data=fixtures.runs()
        with patch.object(module,'subgroup_point',side_effect=lambda rows:dict(windows=len(rows[0]),windowwise={},meanwise={})) as calculate:
            result=module.subgroup(data,[1,2,3])
        self.assertEqual([len(c.args[0][0]) for c in calculate.call_args_list],[64,32,32,32])
        full,first,last,doubled=[c.args[0] for c in calculate.call_args_list]
        self.assertEqual(first[0],full[0][:32]);self.assertEqual(last[0],full[0][32:])
        self.assertEqual(full[0][0]['density_normalized'][0],-.02)
        np.testing.assert_allclose(doubled[0][0]['density_normalized'],full[0][0]['density_normalized'])
        self.assertEqual(result['source_windows_per_seed'],dict(full=64,first=32,last=32,doubled=64))
        self.assertEqual(result['effective_windows_per_seed'],dict(full=64,first=32,last=32,doubled=32))
        self.assertEqual(result['edges']['window_count'],192)

    def test_complete_verified_matrix_records_new_method_and_source_dependencies(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);proof=fixtures.ReferenceSubgroupTest().fixture(root)
            with patch.object(module,'subgroup',return_value={}) as calculate:
                result=module.analyze(root,proof,root/'result.json')
            self.assertEqual(calculate.call_count,10);self.assertEqual(result['method'],METHOD)
            self.assertEqual(result['status'],'order_subgroup_statistics_only');self.assertEqual(result['physical_acceptance'],'incomplete')
            self.assertEqual(result['sources']['shock_order_sensitivity.py'],digest(Path(module.__file__).with_name('shock_order_sensitivity.py')))
            self.assertIn('analyze_alsmeyer_subgroup.py',result['sources'])

    def test_rehashed_wrong_lineage_fails_before_estimation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);proof=fixtures.ReferenceSubgroupTest().fixture(root);target=root/'mc3d-C.json'
            data=json.loads(target.read_text());data['runs'][0]['source_sha256']='wrong';save(target,data)
            report_path=root/'extraction.json';report=json.loads(report_path.read_text());report['artifacts'][target.name]=digest(target);save(report_path,report)
            with patch.object(module,'subgroup') as calculate,self.assertRaisesRegex(ValueError,'source differs'):
                module.analyze(root,proof,root/'invalid.json')
            calculate.assert_not_called();self.assertFalse((root/'invalid.json').exists())

    def test_data_change_during_estimation_is_not_published(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);proof=fixtures.ReferenceSubgroupTest().fixture(root);target=root/'mc3d-C.json'
            def change(runs,seeds):
                target.write_text(target.read_text()+' ')
                return {}
            with patch.object(module,'subgroup',side_effect=change),self.assertRaisesRegex(ValueError,'inputs changed'):
                module.analyze(root,proof,root/'invalid.json')
            self.assertFalse((root/'invalid.json').exists())

if __name__=='__main__':unittest.main()
