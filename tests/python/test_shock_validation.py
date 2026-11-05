import sys,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import configuration
from physical_config import prepare
from shock_analysis import monotone,metrics,equivalence,interval,pooled
class ShockValidationTest(unittest.TestCase):
    def test_matrix_preserves_weight(self):
        def weight(v):
            d,_=configuration(v,123);c,u=prepare(d);return u['scales']['molecules_per_model_particle']
        self.assertAlmostEqual(weight('A')/weight('B'),1)
        for v in ('D','E','F','G'):self.assertAlmostEqual(weight(v)/weight('C'),1)
    def test_common_time_and_lambda(self):
        a,_=configuration('A',123);f,_=configuration('F',123)
        self.assertEqual(a['statistics_windows_s'],f['statistics_windows_s'])
        ca,_=prepare(a);cf,_=prepare(f)
        self.assertAlmostEqual(ca['kn']*a['domain_m'][0],cf['kn']*f['domain_m'][0])
    def test_profile_translation(self):
        x=[i/100 for i in range(-1000,1001)];y=[1+2/(1+math.exp(-v)) for v in x]
        a=metrics(x,y,1,3,1);b=metrics([v+7 for v in x],y,1,3,1)
        self.assertAlmostEqual(a['front']+7,b['front'])
        self.assertAlmostEqual(a['width_10_90'],2*math.log(9),places=4)
        self.assertAlmostEqual(a['gradient_thickness'],4,places=3)
    def test_uncertainty_is_not_permission(self):
        self.assertFalse(equivalence(.001,.1,.01));self.assertFalse(equivalence(.02,0,.01))
        self.assertTrue(equivalence(.002,.001,.01))
        with self.assertRaises(ValueError):interval([1,2])
    def test_monotone(self):self.assertEqual(monotone([0,2,1,3]),[0,1.5,1.5,3])
    def test_pooled_moments_with_different_origins(self):
        # Two deterministic beams. Local temperatures vanish, pooled T does not.
        rows=[]
        for velocity in (-2.,2.):
            r={'A0':1.,'origin_x':velocity,'origin_y':0.,'origin_z':0.}
            for a in 'xyz':r['A1'+a]=r['A3'+a]=0.
            for ab in ('xx','xy','xz','yy','yz','zz'):r['A2'+ab]=0.
            rows.append(r)
        result=pooled(rows,1,1,1)
        self.assertEqual(result['velocity'],[0.,0.,0.])
        self.assertAlmostEqual(result['temperature'],4/3)
        self.assertEqual(result['pressure'][0][0],8)
        self.assertEqual(result['heat'],[0.,0.,0.])
    def test_constant_noise_is_not_independent_samples(self):
        from shock_analysis import interval
        means=[.9,1.,1.1]
        _,u=interval(means)
        self.assertGreater(u,.2)
    def test_prepare_resume_and_configuration_guard(self):
        import tempfile,subprocess,json
        from pathlib import Path
        script=Path(__file__).resolve().parents[2]/'tools/shock_validation.py'
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'matrix'
            command=[sys.executable,str(script),'--binary',sys.executable,'--output-dir',str(root),'--prepare-only']
            self.assertEqual(subprocess.run(command,capture_output=True).returncode,2)
            report=json.loads((root/'validation.json').read_text())
            self.assertEqual(len(report['runs']),21)
            for key in ('stationarity','convergence','cross_code','experiment'):self.assertEqual(report[key]['status'],'incomplete')
            self.assertEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,2)
            input_file=root/'A-123/input_si.json';data=json.loads(input_file.read_text());data['temperature_k']=301;input_file.write_text(json.dumps(data))
            self.assertNotEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,0)
