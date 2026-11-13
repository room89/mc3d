from pathlib import Path
import copy
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from run_profile_likelihood_screen import protocol, generate, trial, summarize, validate_checkpoint


class ScreenTest(unittest.TestCase):
    def test_protocol_and_independent_streams(self):
        p=protocol();self.assertEqual(len(p['specifications']),208)
        self.assertEqual(sum(not s['analytic'] for s in p['specifications']),192)
        a=generate(p['specifications'][32],32,p);b=generate(p['specifications'][32],32,p)
        np.testing.assert_array_equal(a[5],b[5])
        c=generate(p['specifications'][32],33,p)
        self.assertFalse(np.array_equal(a[5],c[5]))
        self.assertEqual(a[5].shape,(32,2,41))
        self.assertGreater(a[3][0,1],0)

    def test_no_motion_controls_and_checkpoint_tampering(self):
        p=protocol();spec=p['specifications'][0];row=trial(spec,0,p);row['plan_sha256']='test'
        validate_checkpoint(row,0,spec,'test',p)
        self.assertEqual([r['status'] for r in row['results'].values()],['converged']*3)
        np.testing.assert_allclose(row['results']['marginal']['profile'],row['results']['oracle']['profile'])
        bad=copy.deepcopy(row);bad['results']['marginal']['integrated_squared_error']+=1
        with self.assertRaises(ValueError):validate_checkpoint(bad,0,spec,'test',p)
        with self.assertRaises(ValueError):validate_checkpoint(row,0,spec,'other',p)
        self.assertEqual(summarize([row],p)['decision'],'incomplete')

    def test_failures_block_advancement(self):
        p=protocol();rows=[]
        for i,spec in enumerate(p['specifications']):
            rows.append(dict(specification=spec,results={m:dict(status='converged',integrated_squared_error=1 if m=='marginal' else 2) for m in p['methods']}))
        self.assertEqual(summarize(rows,p)['decision'],'advance_to_noise_model_research_only')
        rows[32]['results']['marginal']['status']='iteration_limit'
        report=summarize(rows,p)
        self.assertEqual(report['decision'],'do_not_advance_additional_failures')
        self.assertEqual(report['failures']['marginal'],1)


if __name__=='__main__':unittest.main()
