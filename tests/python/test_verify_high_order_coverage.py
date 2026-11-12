import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from verify_high_order_coverage import validate,verify_pair,SOURCES
from shock_high_order_estimator import METHOD,WIDTHS
from shock_validation import digest
from calibrate_shock_ensemble import summarize
from validate_shock_estimator import truth


def report(windows=64,repetitions=128,start=114000000):
    tools=Path(__file__).resolve().parents[2]/'tools'
    p=dict(method=METHOD,group_count=12,subgroups=6,bootstrap_repeats=1024,windows=windows,seed_start=start,
           noise_levels=[.005,.015,.03],derivative_degree=5,widths=list(WIDTHS),peak_training_degree=3,
           peak_training_width=1.5,primary_interpolation='cubic',sources={n:digest(tools/n) for n in SOURCES})
    cases=[]
    for kind in ('logistic','gompertz'):
        for phi in (0.,.8):
            for spatial in (0.,.2):
                for noise in (.005,.015,.03):
                    rows=[dict(repetition=i,estimate=list(truth(kind)),statistical_radius=[.1,.1],processing_sensitivity=[.01,.01],
                               relative_error=[0.,0.],statistical_covers=[True,True],with_processing_covers=[True,True],
                               density_error=.001,density_radius=.003,density_statistical_radius=.002,
                               density_processing_sensitivity=.001,density_covers=True) for i in range(repetitions)]
                    cases.append(dict(shape=kind,phi=phi,spatial=spatial,noise=noise,**summarize(rows)))
    return dict(status='coverage_study_only',group_count=12,seed_start=start,center_profile=True,protocol=p,cases=cases)


class HighOrderCoverageVerificationTest(unittest.TestCase):
    def test_full_integrity_is_not_qualification(self):
        r=report();result=validate(r,64,128)
        self.assertEqual(result['trial_count'],3072);self.assertNotIn('passed',result)
        with self.assertRaisesRegex(ValueError,'unexpected trials'):validate(report(repetitions=8),64,128)

    def test_duplicate_and_fabricated_coverage_rejected(self):
        r=report();r['cases'][0]['runs'][-1]['repetition']=0
        with self.assertRaisesRegex(ValueError,'unexpected trials'):validate(r,64,128)
        r=report();r['cases'][0]['runs'][0]['with_processing_covers'][0]=False
        with self.assertRaisesRegex(ValueError,'coverage or relative'):validate(r,64,128)
        r=report();r['cases'][0]['density_coverage']=.5
        with self.assertRaisesRegex(ValueError,'summary differs'):validate(r,64,128)

    def test_rejection_counts_as_noncoverage_and_wrong_source_rejected(self):
        r=report();c=r['cases'][0];c['runs'][0]=dict(repetition=0,error='Ambiguous minimum');c.update(summarize(c['runs']))
        r['status']='coverage_study_with_failures';result=validate(r,64,128)
        self.assertEqual(result['rejected_trials'],1);self.assertEqual(result['minimum_density_coverage'],127/128)
        r['protocol']['sources']['shock_high_order_estimator.py']='changed'
        with self.assertRaisesRegex(ValueError,'source mismatch'):validate(r,64,128)

    def test_pair_checks_seed_separation_and_never_accepts_physics(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);paths=[]
            for i,(w,n,s) in enumerate([(64,128,114000000),(32,128,116000000),(64,8,110000000),(32,8,112000000)]):
                p=root/f'{i}.json';p.write_text(json.dumps(report(w,n,s)));paths.append(p)
            result=verify_pair(paths[0],paths[1],paths[2:])
            self.assertEqual(result['physical_acceptance'],'incomplete')
            bad=report(32,128,114000001);paths[1].write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError,'Overlapping'):verify_pair(paths[0],paths[1],paths[2:])


if __name__=='__main__':unittest.main()
