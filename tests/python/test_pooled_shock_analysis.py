import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from reanalyze_pooled_shock import validate_plan,pool_edges,check_coverage


def plan():
    return dict(variant='C',total_groups=12,subgroups_per_group=6,groups=[
        dict(name=f'g{g:02d}',subgroups=[dict(path=f'g{g:02d}-s{s}',seeds=[1+18*g+3*s+k for k in range(3)])
             for s in range(6)]) for g in range(12)])


class PooledAnalysisTest(unittest.TestCase):
    def test_duplicate_seeds_or_folder_aliases_are_rejected(self):
        p=plan();validate_plan(p,Path('/tmp/pooled-analysis-test'))
        bad=copy.deepcopy(p);bad['groups'][1]['subgroups'][0]['seeds']=p['groups'][0]['subgroups'][0]['seeds']
        with self.assertRaisesRegex(ValueError,'distinct'):validate_plan(bad,Path('/tmp'))
        bad=copy.deepcopy(p);bad['groups'][1]['subgroups'][0]['path']='other/../g00-s0'
        with self.assertRaisesRegex(ValueError,'reuse folders'):validate_plan(bad,Path('/tmp'))

    def test_pooling_keeps_outliers_and_worst_boundary(self):
        edges=[dict(full=[.1*i,0],first=[0,0],last=[0,0],boundary_min=20-i,
                    window_outliers=i,window_count=192) for i in range(6)]
        r=pool_edges(edges)
        self.assertAlmostEqual(r['full'][0],.25)
        self.assertEqual(r['window_outliers'],15)
        self.assertEqual(r['window_count'],1152)
        self.assertEqual(r['boundary_min'],15)
        edges[0]['window_count']=96
        with self.assertRaisesRegex(ValueError,'Unequal'):pool_edges(edges)

    def test_earlier_pooling_cannot_qualify_new_peak_method(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'coverage.json'
            p.write_text(json.dumps(dict(protocol=dict(subgroups=6,windows=64,group_count=12,
                bootstrap_repeats=1024,center_profile=True,width_band=True,interpolation_sensitivity=True))))
            with self.assertRaisesRegex(ValueError,'independent pooled peaks'):check_coverage(p,64)
