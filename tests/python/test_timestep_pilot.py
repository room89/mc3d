import sys,tempfile,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from analyze_timestep_pilot import summarize,completed_reports

class TimestepPilotTests(unittest.TestCase):
    def test_whole_run_offsets_retained(self):
        x=np.linspace(-4,4,81);y=(1+np.tanh(x))/2
        a=np.tile(y,(3,64,1))+np.array([-.01,0,.01])[:,None,None]
        r=summarize(a,x)
        np.testing.assert_allclose(r['between_run_sd'],.01)
        np.testing.assert_allclose(r['mean'],y)
        np.testing.assert_allclose(r['half_changes'],0,atol=1e-14)
    def test_half_changes_preserve_run_identity(self):
        x=np.linspace(-4,4,81);a=np.tile((1+np.tanh(x))/2,(3,64,1))
        a[1,32:]+=.03
        r=summarize(a,x)
        np.testing.assert_allclose(np.array(r['half_changes'])[1],.03)
        with self.assertRaises(ValueError):summarize(a[:,:63],x)
    def test_active_runner_rejected_without_reading_results(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'not terminal'):
                completed_reports(Path(d),{'variants':['C','E']})
    def test_terminal_error_does_not_become_completed(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'C-exit.json').write_text(json.dumps({'returncode':3}))
            with self.assertRaisesRegex(ValueError,'execution error'):
                completed_reports(root,{'variants':['C']})
