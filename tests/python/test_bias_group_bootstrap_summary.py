import json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from summarize_bias_group_bootstrap import summarize,curve_statistics
from shock_validation import digest


class BiasGroupBootstrapSummaryTest(unittest.TestCase):
    def test_complete_evidence_and_rejection_paths(self):
        x=np.arange(64)[:,None]*np.ones((1,401))
        s=curve_statistics(x,np.zeros(401));assert s['empirical_95_percent_supnorm']==60 and s['maximum_supnorm']==63 and s['maximum_absolute_bootstrap_mean_shift']==31.5
        try:curve_statistics(x[:63],np.zeros(401))
        except ValueError:pass
        else:raise AssertionError('accepted incomplete curves')
        with tempfile.TemporaryDirectory() as temp:
         root=Path(temp);prior={'cases':{},'bank_hashes':{'mc3d-training':'training','mc3d-heldout':'heldout'}}
         for name in prior['bank_hashes']:
          for mode in ('windowwise','mean32'):prior['cases'][f'{name}-logistic-{mode}-nonlinear']={'error':[0.]*401}
         pp=root/'prior.json';pp.write_text(json.dumps(prior))
         plan=dict(banks=['mc3d-training','mc3d-heldout'],modes=['windowwise','mean32'],scheduled_replicates=256,replicates_per_case=64,seed_start=190000000,shape='logistic',added_motion=0.,status='exploratory_source_group_bootstrap',bank_hashes=prior['bank_hashes'],prior_evidence_sha256=digest(pp))
         p=root/'plan.json';p.write_text(json.dumps(plan));ph=digest(p);hashes={}
         for ci,(name,mode) in enumerate(( (n,m) for n in plan['banks'] for m in plan['modes'])):
          for i in range(64):
           p=root/f'{name}-{mode}-{i}.json';p.write_text(json.dumps(dict(bank=name,mode=mode,replicate=i,seed=190000000+ci*10000+i,plan_sha256=ph,status='ok',curve=[float(i*(1 if ci<2 else 2))]*401)));hashes[p.name]=digest(p)
         c=root/'complete.json';record=dict(completed_replicates=256,failed_replicates=0,plan_sha256=ph,checkpoint_sha256=hashes);c.write_text(json.dumps(record))
         out=summarize(root,pp);assert out['comparisons']['windowwise']['empirical_95_percent_supnorm']==60
         record['failed_replicates']=1;c.write_text(json.dumps(record))
         try:summarize(root,pp)
         except ValueError as e:assert 'Failed replicates' in str(e)
         else:raise AssertionError('accepted failures')
         record['failed_replicates']=0;c.write_text(json.dumps(record));(root/next(iter(hashes))).write_text('{}')
         try:summarize(root,pp)
         except ValueError as e:assert 'hash mismatch' in str(e)
         else:raise AssertionError('accepted corrupted checkpoint')


if __name__ == '__main__':
    unittest.main()
