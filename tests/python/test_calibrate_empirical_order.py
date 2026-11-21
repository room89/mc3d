from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from calibrate_empirical_order import cached,verify_banks
from shock_validation import save,digest


class EmpiricalOrderScreenTest(unittest.TestCase):
    def test_failure_is_retained_and_cached_not_retried(self):
        with tempfile.TemporaryDirectory() as d:
            spec=('mc3d-training','gompertz',.05,0);protocol={'method':'test-only'};task=(spec,d,protocol,'unused')
            with patch('calibrate_empirical_order.trial',side_effect=ValueError('Ambiguous minimum')):
                returned,result=cached(task)
            self.assertEqual(returned,spec);self.assertEqual(result,dict(repetition=0,error='Ambiguous minimum'))
            with patch('calibrate_empirical_order.trial',side_effect=AssertionError('Must not retry rejected checkpoint')):
                self.assertEqual(cached(task),(spec,result))
            with self.assertRaisesRegex(ValueError,'protocol mismatch'):cached((spec,d,{'method':'changed'},'unused'))

    def test_missing_or_replaced_completed_verification_rejected_before_loading(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest=root/'manifest.json';proof=root/'verification.json';save(manifest,dict(banks={}))
            save(proof,dict(status='conditional_empirical_banks_prepared',physical_acceptance='incomplete'))
            with patch('calibrate_empirical_order.load_bank',side_effect=AssertionError('Must fail preflight')):
                with self.assertRaisesRegex(ValueError,'completed physical-bank'):verify_banks(manifest,proof)
            save(proof,dict(status='conditional_empirical_banks_verified',physical_acceptance='incomplete',manifest_sha256=digest(manifest),region_phase_checks=648,loader_sha256='0'*64,banks={}))
            with self.assertRaisesRegex(ValueError,'completed physical-bank'):verify_banks(manifest,proof)

    def test_all_four_verified_bank_identities_required(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest=root/'manifest.json';proof=root/'verification.json';save(manifest,dict(banks={}))
            source=Path(__file__).resolve().parents[2]/'tools/shock_empirical_fluctuations.py'
            save(proof,dict(status='conditional_empirical_banks_verified',physical_acceptance='incomplete',manifest_sha256=digest(manifest),region_phase_checks=648,loader_sha256=digest(source),banks={}))
            with self.assertRaisesRegex(ValueError,'All four'):verify_banks(manifest,proof)

if __name__=='__main__':unittest.main()
