from pathlib import Path
import sys
import unittest
import tempfile
import json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from extract_alsmeyer_profiles import canonical_windows,extract


class ExtractReferenceProfilesTest(unittest.TestCase):
    def rows(self):
        return [dict(front=i*.01,x_aligned=[-1-i*.01,-i*.01,1-i*.01],density_normalized=[-.02,.5,1.03]) for i in range(64)]

    def test_preserves_drift_noise_and_all_windows(self):
        result=canonical_windows(self.rows(),3)
        self.assertEqual(len(result),64);self.assertEqual(result[-1]['window'],63)
        self.assertEqual(result[-1]['front'],.63)
        self.assertEqual(result[0]['density_normalized'],[-.02,.5,1.03])
        self.assertAlmostEqual(result[-1]['x_aligned'][0]+result[-1]['front'],-1.)

    def test_truncation_nonfinite_and_changing_grid_rejected(self):
        with self.assertRaisesRegex(ValueError,'64 complete'):canonical_windows(self.rows()[:-1],3)
        for value in (float('nan'),2.):
            rows=self.rows();rows[-1]['x_aligned'][-1]=value
            with self.assertRaises(ValueError):canonical_windows(rows,3)

    def test_prepared_only_cannot_create_extraction(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);p=root/'verification.json';p.write_text(json.dumps(dict(status='prepared_inputs_verified')))
            with self.assertRaisesRegex(ValueError,'complete output verification'):
                extract(root,root,root/'evidence.json',p,root/'extracted')
            self.assertFalse((root/'extracted').exists())


if __name__=='__main__':unittest.main()
