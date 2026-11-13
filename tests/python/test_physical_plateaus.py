from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from diagnose_physical_plateaus import analyze,region_statistics


class PlateauTest(unittest.TestCase):
    def bank(self,offset,group):
        x=np.linspace(-14,14,561);mean=1/(1+np.exp(-x))+offset
        noise=np.zeros((6,64,len(x)))
        return dict(x=x,mean_profile=mean,fluctuations=noise,
                    sources=[dict(group=f'{group}{i//3}',seed=100*group+i+1) for i in range(6)])

    def test_heldout_mean_not_centered_away_and_regions_frozen(self):
        a=self.bank(0.,1);b=self.bank(.002,2)
        result=analyze(a,b)
        for r in result['regions'].values():
            self.assertAlmostEqual(r['heldout_minus_training_mean'],.002,places=12)
            self.assertEqual(r['training']['bounds'],r['heldout']['bounds'])
        changed=analyze(a,self.bank(.005,2))
        self.assertEqual(result['training_center'],changed['training_center'])
        self.assertEqual(result['motion']['training']['region_overlaps'],0)

    def test_time_change_and_run_offsets_retained(self):
        x=np.arange(6.);raw=np.zeros((6,64,6));raw[:,32:,:]=2.
        sources=[dict(group=str(i//3),seed=i) for i in range(6)]
        r=region_statistics(raw,x,np.ones(6,dtype=bool),0.,sources)
        self.assertEqual(r['time_change']['mean'],2.)
        self.assertEqual(r['phases']['full']['second_moment_about_training_mean'],2.)
        raw[:]=3.
        r=region_statistics(raw,x,np.ones(6,dtype=bool),0.,sources)
        self.assertEqual(r['phases']['full']['lag_one_product_about_training_mean'],9.)

    def test_overlapping_sources_rejected(self):
        a=self.bank(0.,1)
        with self.assertRaisesRegex(ValueError,'Overlapping'):analyze(a,a)


if __name__=='__main__':unittest.main()
