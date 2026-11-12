from pathlib import Path
import sys
import tempfile
import unittest
import json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from run_profile_likelihood_screen import protocol,trial,save,digest,summarize
from diagnose_profile_likelihood_screen import diagnose,paired


class DiagnosisTest(unittest.TestCase):
    def make_evidence(self,folder):
        p=protocol();p['specifications']=p['specifications'][:2]
        save(folder/'plan.json',p);ph=digest(folder/'plan.json');rows=[]
        for i,s in enumerate(p['specifications']):
            row=trial(s,i,p);row['plan_sha256']=ph;save(folder/f'trial-{i:03d}.json',row);rows.append(row)
        report=summarize(rows,p);report['plan_sha256']=ph
        report['checkpoints']={f'trial-{i:03d}.json':digest(folder/f'trial-{i:03d}.json') for i in range(2)}
        save(folder/'complete.json',report)

    def test_readonly_decomposition_and_start_diagnostics(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root);self.make_evidence(folder)
            before={p.name:p.read_bytes() for p in folder.iterdir()}
            result=diagnose(folder);self.assertEqual(len(result['scenarios']),1)
            s=result['scenarios'][0];self.assertEqual(s['trials'],2)
            self.assertEqual(s['comparisons']['marginal_minus_hard']['mean'],0)
            for method in s['methods'].values():
                self.assertAlmostEqual(method['mean_ise'],method['empirical_mean_error_energy']+method['empirical_scatter_energy'],places=12)
                self.assertEqual(len(method['starts']),2)
            self.assertEqual(before,{p.name:p.read_bytes() for p in folder.iterdir()})
            path=folder/'trial-000.json';path.write_text(path.read_text()+' ')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):diagnose(folder)

    def test_summary_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root);self.make_evidence(folder)
            path=folder/'complete.json';r=json.loads(path.read_text());r['trials']+=1;path.write_text(json.dumps(r))
            with self.assertRaisesRegex(ValueError,'Summary mismatch'):diagnose(folder)

    def test_paired_mcse(self):
        d=paired([1.,2.,3.]);self.assertEqual(d['mean'],2.)
        self.assertAlmostEqual(d['mcse'],1/np.sqrt(3))


if __name__=='__main__':unittest.main()
