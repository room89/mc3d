from pathlib import Path
import sys
import unittest
import copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_alignment import GRID, initial
from shock_ensemble_estimator import centered_density
from shock_order_sensitivity import average_seed_windows, envelope, subgroup_point, pooled_point, estimate
from calibrate_order_sensitivity import moving_fixture
from validate_shock_estimator import shape
from experiment_crossfit_gradient import fixture


class OrderSensitivityTest(unittest.TestCase):
    def record(self, offset=0.):
        return dict(density=shape(GRID-offset,'logistic').tolist(), density_statistical_radius=.002,
                    density_processing_sensitivity=.0003,
                    width_10_90=dict(mean=4.4+offset,statistical_radius=.01,processing_sensitivity=.005),
                    gradient_thickness=dict(mean=4.+offset,statistical_radius=.02,processing_sensitivity=.007))

    def test_envelope_contains_both_density_bands_and_metric_intervals(self):
        a=self.record();b=self.record(.01);b['density_statistical_radius']=.003
        b['gradient_thickness']['statistical_radius']=.03;original=copy.deepcopy(a)
        r=envelope(a,b)
        self.assertEqual(a,original)
        for source in (a,b):
            radius=source['density_statistical_radius']+source['density_processing_sensitivity']
            self.assertLessEqual(float(np.max(abs(np.asarray(source['density'])-r['density'])))+radius,r['density_simultaneous_radius']+1e-14)
            for metric in ('width_10_90','gradient_thickness'):
                needed=abs(source[metric]['mean']-r[metric]['mean'])+source[metric]['statistical_radius']+source[metric]['processing_sensitivity']
                self.assertLessEqual(needed,r[metric]['statistical_radius']+r[metric]['processing_sensitivity']+1e-14)
        self.assertEqual(r['density'],a['density'])
        self.assertEqual(r['status'],'experimental_coverage_unverified')

    def test_envelope_rejects_nonfinite_negative_or_unresolved_bounds(self):
        a=self.record()
        for value in (-.001,float('nan')):
            b=self.record();b['density_statistical_radius']=value
            with self.assertRaises(ValueError):envelope(a,b)
        b=self.record();b['density'][0]=float('inf')
        with self.assertRaises(ValueError):envelope(a,b)
        b=self.record();b['density_statistical_radius']=.11
        with self.assertRaises(ValueError):envelope(a,b)

    def test_mean_preserves_outliers_and_requires_physical_grid_and_complete_seeds(self):
        x=np.linspace(-14,14,281);run=[]
        for i in range(32):
            y=shape(x,'logistic');y[0]=-float(i);y[-1]=float(i+2)
            run.append(dict(x=x.tolist(),density_normalized=y.tolist()))
        reduced=average_seed_windows([run,run,run])
        for seed in reduced:
            self.assertEqual(len(seed),1)
            np.testing.assert_allclose(initial(seed[0])[1],np.mean([w['density_normalized'] for w in run],axis=0))
        with self.assertRaises(ValueError):average_seed_windows([run,run,run[:-1]])
        changed=copy.deepcopy(run);changed[0]['x'][0]+=.01
        with self.assertRaises(ValueError):average_seed_windows([run,run,changed])

    def test_motion_adapter_retains_noise_and_reproduces_original_at_unit_amplitude(self):
        for kind in ('logistic','gompertz'):
            original,_=fixture(kind,150500123,.8,.2,noise=.03)
            unit=moving_fixture(kind,150500123,.8,.2,.03,1.,32)
            still=moving_fixture(kind,150500123,.8,.2,.03,0.,32)
            clean,_=fixture(kind,150500123,.8,.2,noise=0.)
            for s in range(3):
                for i in range(32):
                    x,y,_=initial(original[s][i]);_,cy,_=initial(clean[s][i])
                    np.testing.assert_allclose(initial(unit[s][i])[1],y,atol=2e-14,rtol=1e-12)
                    np.testing.assert_allclose(initial(still[s][i])[1]-shape(x,kind),y-cy,atol=2e-14,rtol=1e-11)
        with self.assertRaises(ValueError):moving_fixture('logistic',1,0.,0.,.03,-1.,32)

    def test_real_subgroup_and_rejection_of_mixed_sampling(self):
        trials=moving_fixture('logistic',150500321,0.,0.,0.,0.,32)
        point=subgroup_point(trials)
        # Independently fitted frames can differ by solver tolerance. Compare
        # the final centered shapes, rather than requiring identical raw shifts.
        a=centered_density(point['windowwise']['density_cubic'],interpolation='cubic')
        b=centered_density(point['meanwise']['density_cubic'],interpolation='cubic')
        np.testing.assert_allclose(a,b,atol=1e-8,rtol=0)
        group=pooled_point([point]*6)
        result=estimate([group]*12,repeats=128,seed=1)
        self.assertLess(result['alignment_order_density_difference'],1e-8)
        bad=copy.deepcopy(group);bad['windows']=64
        with self.assertRaises(ValueError):estimate([group]*11+[bad],repeats=128)
        with self.assertRaises(ValueError):pooled_point([point]*5)


if __name__=='__main__':unittest.main()
