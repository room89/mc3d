import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from calibrate_shock_ensemble import run_cached


class EnsembleCheckpointTest(unittest.TestCase):
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
