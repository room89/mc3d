from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from verify_order_screen import expected_protocol,scenario_keys,verify,verify_pair
from calibrate_shock_ensemble import summarize
from shock_validation import save
from validate_shock_estimator import truth


class VerifyOrderScreenTest(unittest.TestCase):
    def row(self,kind,i):
        return dict(repetition=i,estimate=truth(kind).tolist(),statistical_radius=[.01,.01],processing_sensitivity=[.002,.002],relative_error=[0.,0.],
            statistical_covers=[True,True],with_processing_covers=[True,True],density_error=.001,density_radius=.003,
            density_statistical_radius=.002,density_processing_sensitivity=.001,density_covers=True,order_density_difference=.0005,
            order_metric_differences=dict(width_10_90=.001,gradient_thickness=.001))

    def write(self,path,windows=64,start=160000000,complete=False,failure=False):
        protocol=expected_protocol(windows,2,start);folder=path.with_suffix('.parts');folder.mkdir()
        cases=[]
        for case_index,key in enumerate(scenario_keys() if complete else scenario_keys()[:1]):
            rows=[]
            for i in range(2 if complete else 1):
                row=dict(repetition=i,error='Ambiguous minimum') if failure and case_index==0 and i==0 else self.row(key[0],i)
                spec=key+(i,);save(folder/('-'.join(map(str,spec))+'.json'),dict(protocol=protocol,spec=spec,result=row));rows.append(row)
            cases.append(dict(shape=key[0],phi=key[1],spatial=key[2],noise=key[3],motion_amplitude=key[4],**summarize(rows)))
        if complete:save(path,dict(status='screen_with_failures' if failure else 'screen_only',physical_acceptance='incomplete',protocol=protocol,cases=cases))
        return sorted(folder.glob('*.json'))[0]

    def test_partial_is_not_complete_and_final_requires_every_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'report.json';self.write(path)
            result=verify(path,64,2,160000000)
            self.assertEqual(result['status'],'partial_checkpoint_snapshot_verified');self.assertEqual(result['completed_trials'],1);self.assertEqual(result['missing_trials'],143)
            save(path,dict(protocol=result['protocol'],physical_acceptance='incomplete',cases=[]))
            with self.assertRaisesRegex(ValueError,'complete matching'):verify(path,64,2,160000000)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'No screen'):verify(Path(d)/'none.json',64,2,160000000)

    def test_complete_failure_counts_as_noncoverage_and_summary_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'report.json';self.write(path,complete=True,failure=True)
            result=verify(path,64,2,160000000)
            self.assertEqual(result['status'],'complete_screen_evidence_verified');self.assertEqual(result['rejected_trials'],1)
            self.assertEqual(result['completed_scenarios'][0]['density_coverage'],.5)
            data=json.loads(path.read_text());data['cases'][0]['density_coverage']=1.;save(path,data)
            with self.assertRaisesRegex(ValueError,'summary differs'):verify(path,64,2,160000000)

    def test_source_identity_coverage_flags_and_budget_cannot_be_rehashed_away(self):
        for mutation in ('source','flag','budget','identity'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as d:
                path=Path(d)/'report.json';part=self.write(path);data=json.loads(part.read_text())
                if mutation=='source':data['protocol']['sources']['shock_order_sensitivity.py']='0'*64
                elif mutation=='flag':data['result']['density_covers']=False
                elif mutation=='budget':data['result']['density_radius']=.004
                else:data['spec'][-1]=True
                save(part,data)
                with self.assertRaises(ValueError):verify(path,64,2,160000000)

    def test_extra_checkpoint_and_overlapping_seed_ranges_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);full=root/'full.json';half=root/'half.json'
            part=self.write(full);self.write(half,windows=32,start=162000000)
            result=verify_pair(full,half,2,160000000,162000000)
            self.assertEqual(result['physical_acceptance'],'incomplete')
            extra=part.with_name('duplicate.json');extra.write_bytes(part.read_bytes())
            with self.assertRaisesRegex(ValueError,'Unknown, repeated'):verify(full,64,2,160000000)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);full=root/'full.json';half=root/'half.json'
            self.write(full);self.write(half,windows=32,start=160000000)
            with self.assertRaisesRegex(ValueError,'Overlapping'):verify_pair(full,half,2,160000000,160000000)

if __name__=='__main__':unittest.main()
