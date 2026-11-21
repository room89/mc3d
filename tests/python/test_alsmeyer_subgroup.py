from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
import analyze_alsmeyer_subgroup as module
from shock_validation import digest,save


def runs():
    x=np.linspace(-20,20,401);y=1/(1+np.exp(-x))
    y[0]=-.02;y[-1]=1.03
    return [dict(seed=s,source_sha256='source',collision_mean_free_path_m=.001,
                 statistics_windows_s=[[i+1,i+2] for i in range(64)],
                 windows=[dict(window=i,front=.01*i,x_aligned=(x-.01*i).tolist(),
                               density_normalized=y.tolist()) for i in range(64)]) for s in (1,2,3)]


class ReferenceSubgroupTest(unittest.TestCase):
    def test_preserves_windows_and_builds_correct_phases(self):
        data=runs()
        with patch.object(module,'group_point',return_value={'point':True}) as point:
            result=module.subgroup(data,[1,2,3])
        self.assertEqual([len(c.args[0][0]) for c in point.call_args_list],[64,32,32,32])
        full,first,last,doubled=[c.args[0] for c in point.call_args_list]
        self.assertEqual(first[0],full[0][:32]);self.assertEqual(last[0],full[0][32:])
        self.assertEqual(full[0][0]['density_normalized'][0],-.02)
        np.testing.assert_allclose(doubled[0][0]['density_normalized'],full[0][0]['density_normalized'])
        self.assertEqual(result['edges']['window_count'],192)
        self.assertEqual(result['collision_mean_free_path_m'],.001)

    def test_incompatible_sampling_and_seed_reuse_rejected(self):
        for change in ('seed','duration','gap','order','scale','grid'):
            data=runs()
            if change=='seed':data[1]['seed']=1
            if change=='duration':data[0]['statistics_windows_s'][0][1]+=.1
            if change=='gap':data[0]['statistics_windows_s'][0]=[.5,1.5]
            if change=='order':data[0]['windows'][1]['window']=0
            if change=='scale':data[1]['collision_mean_free_path_m']*=2
            if change=='grid':
                for w in data[1]['windows']:w['x_aligned']=[v+.1 for v in w['x_aligned']]
            with self.subTest(change=change),self.assertRaises(ValueError):module.prepare_runs(data,[1,2,3])

    def fixture(self,root):
        proof=dict(status='matrix_outputs_verified',tool_sha256=digest(Path(module.__file__).with_name('verify_alsmeyer_outputs.py')),
                   inputs={'seeds':[1,2,3]},cases={})
        artifacts={}
        for key in module.KEYS:
            target=root/(key+'.json');save(target,dict(status='profiles_only',runs=runs()));artifacts[target.name]=digest(target)
            for s in (1,2,3):
                proof['cases'][key.replace('-','/',1)+'-'+str(s)]={'manifest_sha256':'source','execution_sha256':'source'}
        verified=root/'verified.json';save(verified,proof)
        report=dict(status='verified_profiles_extracted',case_count=30,dataset_count=10,
                    verification_sha256=digest(verified),input_verification=proof['inputs'],artifacts=artifacts,
                    sources={n:digest(Path(module.__file__).with_name(n)) for n in module.EXTRACTION_SOURCES})
        save(root/'extraction.json',report);return verified

    def test_verified_matrix_and_tampered_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);verified=self.fixture(root)
            with patch.object(module,'subgroup',return_value={}) as calculate:
                result=module.analyze(root,verified,root/'result.json')
            self.assertEqual(calculate.call_count,10);self.assertEqual(result['physical_acceptance'],'incomplete')
            target=root/'mc3d-C.json';target.write_text(target.read_text()+' ')
            with patch.object(module,'subgroup') as calculate,self.assertRaisesRegex(ValueError,'profile changed'):
                module.analyze(root,verified,root/'invalid.json')
            calculate.assert_not_called();self.assertFalse((root/'invalid.json').exists())

    def test_rehashed_wrong_lineage_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);verified=self.fixture(root);target=root/'mc3d-C.json'
            data=json.loads(target.read_text());data['runs'][0]['source_sha256']='wrong';save(target,data)
            record=root/'extraction.json';report=json.loads(record.read_text());report['artifacts'][target.name]=digest(target);save(record,report)
            with self.assertRaisesRegex(ValueError,'source differs'):module.analyze(root,verified,root/'invalid.json')
            self.assertFalse((root/'invalid.json').exists())


if __name__=='__main__':unittest.main()
