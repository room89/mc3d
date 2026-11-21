import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from verify_alsmeyer_matrix import verify
from shock_validation import digest


class VerifyReferenceMatrixTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mc3d-reference-verify-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = Path(__file__).resolve().parents[2]
        self.evidence = self.repo/'validation/stage9b/experiment.json'
        self.mc, self.sp = self.root/'mc', self.root/'sp'
        subprocess.run([sys.executable, str(self.repo/'tools/prepare_alsmeyer_matrix.py'),
                        '--reference', str(self.evidence), '--output-dir', str(self.mc)],
                       check=True, capture_output=True)
        subprocess.run([sys.executable, str(self.repo/'tools/export_alsmeyer_sparta.py'),
                        str(self.mc), '--output-dir', str(self.sp)], check=True, capture_output=True)

    def test_valid_inputs_never_imply_execution(self):
        r = verify(self.mc, self.sp, self.evidence)
        self.assertEqual((r['mc3d_cases'], r['sparta_cases']), (21, 9))
        self.assertEqual(r['physical_acceptance'], 'incomplete')
        self.assertFalse(any(self.mc.glob('*/output/manifest.json')))

    def test_rehashed_wrong_normalization_is_still_rejected(self):
        config = self.mc/'C-30001/solver.json'
        c = json.loads(config.read_text()); c['kn'] *= 2; config.write_text(json.dumps(c))
        manifest = self.mc/'preparation.json'; m = json.loads(manifest.read_text())
        next(x for x in m['cases'] if x['folder']=='C-30001')['hashes']['solver.json'] = digest(config)
        manifest.write_text(json.dumps(m))
        export = self.sp/'export.json'; e = json.loads(export.read_text())
        e['preparation_sha256'] = digest(manifest); export.write_text(json.dumps(e))
        with self.assertRaisesRegex(ValueError, 'Normalized configuration'):
            verify(self.mc, self.sp, self.evidence)

    def test_quick_export_cannot_be_treated_as_production(self):
        export = self.sp/'export.json'; e = json.loads(export.read_text()); e['quick'] = True
        export.write_text(json.dumps(e))
        with self.assertRaisesRegex(ValueError, 'quick inputs'):
            verify(self.mc, self.sp, self.evidence)


if __name__ == '__main__': unittest.main()
