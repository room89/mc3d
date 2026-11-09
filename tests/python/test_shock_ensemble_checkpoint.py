import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from calibrate_shock_ensemble import run_cached,summarize


class EnsembleCheckpointTest(unittest.TestCase):
    def test_invalid_estimator_is_saved_and_counts_as_noncoverage(self):
        spec=('gompertz',.8,.2,.03,0,14000000,True)
        with tempfile.TemporaryDirectory() as folder:
            with patch('calibrate_shock_ensemble.run',side_effect=ValueError('Ambiguous minimum')) as compute:
                failed=run_cached((spec,folder,{}))
                self.assertEqual(failed,run_cached((spec,folder,{})))
                compute.assert_called_once()
            valid=dict(statistical_covers=[True,True],with_processing_covers=[True,True],density_covers=True)
            result=summarize([valid,failed])
            self.assertEqual(result['statistical_coverage'],[.5,.5])
            self.assertEqual(result['with_processing_coverage'],[.5,.5])
            self.assertEqual(result['density_coverage'],.5)
            self.assertEqual(result['valid_count'],1)
            self.assertEqual(result['errors'],[failed])
            self.assertEqual(summarize([failed])['statistical_coverage'],[0.,0.])

    def test_resume_reuses_result_and_rejects_changed_protocol(self):
        spec=('gompertz',.8,.2,.03,0,14000000,True)
        protocol={'sources':{'estimator':'original'}}
        with tempfile.TemporaryDirectory() as folder:
            with patch('calibrate_shock_ensemble.run',return_value={'repetition':0}) as compute:
                first=run_cached((spec,folder,protocol))
                self.assertEqual(first,run_cached((spec,folder,protocol)))
                compute.assert_called_once()
                with self.assertRaisesRegex(ValueError,'protocol mismatch'):
                    run_cached((spec,folder,{'sources':{'estimator':'changed'}}))
                changed=spec[:-2]+(15000000,True)
                with self.assertRaisesRegex(ValueError,'protocol mismatch'):
                    run_cached((changed,folder,protocol))
            files=list(Path(folder).iterdir())
            self.assertEqual(len(files),1)
            self.assertEqual(json.loads(files[0].read_text())['result'],first)


if __name__=='__main__':unittest.main()
