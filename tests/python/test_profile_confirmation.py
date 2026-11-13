from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from confirm_profile_likelihood import prepare, summarize_confirmation, read_rows, check_sources, SEED, SOURCES
from run_profile_likelihood_screen import protocol, generate, digest, save


class ConfirmationTest(unittest.TestCase):
    def make_plan(self):
        p=protocol();p['seed']=SEED;p['critical_value']=2.2414027276049464
        p['planning']={s:dict(fixed_per_stratum=2,target_halfwidth=.1) for s in ('logistic','gompertz')}
        p['specifications']=[dict(shape=s,amplitude=a,noise=n,repetition=r,analytic=False) for s in p['planning'] for a in (.05,.5) for n in ('low','high') for r in range(2)]
        return p

    def test_fixed_planning_formula_and_disjoint_stream_inputs(self):
        old=protocol();scenarios=[]
        for s in ('logistic','gompertz'):
            for a in (.05,.5):
                for n in ('low','high'):
                    scenarios.append(dict(analytic=False,shape=s,amplitude=a,noise=n,trials=16,
                                          methods={'hard':{'mean_ise':1.}},
                                          comparisons={'marginal_minus_hard':{'paired_differences':[-.1,.1]*8}}))
        proof=dict(scenarios=scenarios,plan_sha256='plan',report_sha256='report')
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root);save(folder/'plan.json',old)
            with patch('confirm_profile_likelihood.diagnose',return_value=proof):p=prepare(folder)
        self.assertEqual(p['planning']['gompertz']['fixed_per_stratum'],32)
        self.assertEqual(len(p['specifications']),256)
        spec=p['specifications'][0]
        self.assertFalse(np.array_equal(generate(spec,0,p)[5],generate(spec,0,old)[5]))

    def test_partial_precision_failure_and_advantage_decisions(self):
        p=self.make_plan();rows=[]
        for i,spec in enumerate(p['specifications']):
            rows.append(dict(index=i,specification=spec,results={m:dict(status='converged',integrated_squared_error=1. if m=='marginal' else 2.) for m in p['methods']}))
        self.assertEqual(summarize_confirmation(rows[:-1],p)['decision'],'incomplete')
        self.assertEqual(summarize_confirmation(rows,p)['decision'],'advance_to_noise_model_research_only')
        failure=copy.deepcopy(rows);failure[0]['results']['oracle']['status']='rejected'
        self.assertEqual(summarize_confirmation(failure,p)['decision'],'do_not_advance_fit_failures')
        for row in rows:row['results']['marginal']['integrated_squared_error']=3.
        self.assertEqual(summarize_confirmation(rows,p)['decision'],'no_advantage_for_all_families')
        for row in rows:row['results']['marginal']['integrated_squared_error']=1. if row['index']%2 else 4.
        self.assertEqual(summarize_confirmation(rows,p)['decision'],'inconclusive_precision')
        for v in p['planning'].values():v['target_halfwidth']=10.
        self.assertEqual(summarize_confirmation(rows,p)['decision'],'inconclusive_advantage')

    def test_changed_sources_and_unexpected_checkpoint(self):
        p=self.make_plan();tools=Path(__file__).resolve().parents[2]/'tools'
        p['sources']={n:digest(tools/n) for n in SOURCES}
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root);save(folder/'plan.json',p);check_sources(folder,p)
            bad=copy.deepcopy(p);bad['sources'][SOURCES[0]]='changed'
            with self.assertRaises(ValueError):check_sources(folder,bad)
            save(folder/'trial-9999.json',{})
            with self.assertRaisesRegex(ValueError,'Unexpected checkpoint'):read_rows(folder,p)


if __name__=='__main__':unittest.main()
