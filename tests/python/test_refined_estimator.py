import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
    import numpy as np
    from shock_alignment import shift, shift_many, template
    from experiment_crossfit_gradient import fixture
    from shock_refined_estimator import point, resample_trials
except ImportError:
    np=None

@unittest.skipIf(np is None,'Requires numpy')
class RefinedEstimatorTest(unittest.TestCase):
    def test_batched_shifts_match_scalar(self):
        trials,_=fixture('gompertz',6000000,.8,.2)
        x,y=template(trials,0)
        before=[shift(w,x,y) for w in trials[0]]
        after=shift_many(trials[0],x,y)
        np.testing.assert_allclose([a[0] for a in before],[a[0] for a in after],atol=1e-10,rtol=0)
        np.testing.assert_allclose([a[1]['loss'] for a in before],[a[1]['loss'] for a in after],atol=1e-14,rtol=0)

    def test_different_source_grids_fall_back(self):
        trials,_=fixture('logistic',6000000,0.,0.)
        run=trials[0][:2];x,y=template(trials,0)
        run[0]['x_aligned']=run[0]['x_aligned'][::2]
        run[0]['density_normalized']=run[0]['density_normalized'][::2]
        self.assertEqual(shift_many(run,x,y),[shift(w,x,y) for w in run])

    def test_resampling_keeps_seed_identity_and_blocks(self):
        trials=[[(s,i) for i in range(64)] for s in range(3)]
        sampled=resample_trials(trials,8,np.random.default_rng(1))
        for seed,run in enumerate(sampled):
            self.assertEqual(len(run),64)
            self.assertTrue(all(s==seed for s,i in run))
            for start in range(0,64,8):
                self.assertEqual([i for s,i in run[start:start+8]],[(run[start][1]+j)%64 for j in range(8)])

    def test_noiseless_point(self):
        for kind,exact in [('logistic',4.),('gompertz',np.e),('flat_peak',4.)]:
            trials,_=fixture(kind,6000000,0.,0.,noise=0.)
            r=point(trials,processing=True)
            self.assertLess(abs(r['metrics'][1]/exact-1),.005)
            self.assertEqual(r['seed_metrics'].shape,(3,2))

    def test_invalid_blocks(self):
        with self.assertRaises(ValueError):resample_trials([[1]*64]*3,7,np.random.default_rng(1))

if __name__=='__main__':unittest.main()
