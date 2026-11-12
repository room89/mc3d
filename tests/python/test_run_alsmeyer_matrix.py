import json
from pathlib import Path
import shutil
import sys
from threading import Event, Timer
import unittest
import test_alsmeyer_outputs as fixtures
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from run_alsmeyer_matrix import input_hashes, reuse, run_case
from shock_validation import digest


class ReferenceRunnerTest(unittest.TestCase):
    def setUp(self):
        fixture = fixtures.OutputVerificationTest(); fixture.setUp()
        self.addCleanup(fixture.temp.cleanup)
        self.root = fixture.root; fixture.fixture_sparta()
        self.templates = self.root/'templates'; self.templates.mkdir()
        for p in self.root.glob('profile.*'): shutil.move(str(p), self.templates/p.name)
        (self.root/'run.log').unlink()
        self.binary = self.root/'fake_sparta'
        self.inputs = input_hashes('sparta',self.root)

    def executable(self, code):
        self.binary.write_text('#!'+sys.executable+'\n'+code+'\n'); self.binary.chmod(0o755)
        return digest(self.binary)

    def run_fake(self, binary_hash, stop=None):
        return run_case('sparta',self.root,self.binary,binary_hash,self.inputs,'full',stop or Event())

    def test_execute_verify_then_reuse_without_relaunch(self):
        binary_hash = self.executable("from pathlib import Path\nimport shutil\nfor p in Path('templates').glob('profile.*'): shutil.copyfile(p,p.name)\nprint('done')")
        result = self.run_fake(binary_hash)
        self.assertEqual(result['status'],'outputs_verified'); self.assertFalse(result['reused'])
        before = (self.root/'execution.json').read_bytes()
        self.assertTrue(self.run_fake(binary_hash)['reused'])
        self.assertEqual(before,(self.root/'execution.json').read_bytes())
        (self.root/'profile.96').write_text('corrupt')
        with self.assertRaisesRegex(ValueError,'hash mismatch'): self.run_fake(binary_hash)

    def test_nonzero_exit_preserved_and_not_overwritten(self):
        binary_hash = self.executable("raise SystemExit(7)")
        with self.assertRaisesRegex(ValueError,'code 7'): self.run_fake(binary_hash)
        record = json.loads((self.root/'execution.json').read_text())
        self.assertEqual(record['returncode'],7); self.assertEqual(record['status'],'interrupted_or_error')
        before = (self.root/'execution.json').read_bytes()
        with self.assertRaisesRegex(ValueError,'Incomplete attempt'): self.run_fake(binary_hash)
        self.assertEqual(before,(self.root/'execution.json').read_bytes())

    def test_unrecorded_output_and_changed_input_rejected_before_launch(self):
        binary_hash = self.executable("raise AssertionError('must not start')")
        (self.root/'profile.1').write_text('old')
        with self.assertRaisesRegex(ValueError,'Unrecorded output'): self.run_fake(binary_hash)
        (self.root/'profile.1').unlink(); (self.root/'in.shock').write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed before'): self.run_fake(binary_hash)
        self.assertFalse((self.root/'execution.json').exists())

    def test_stop_terminates_child_and_keeps_failed_attempt(self):
        binary_hash = self.executable('import time\ntime.sleep(30)')
        stop = Event(); timer = Timer(.15,stop.set); timer.start(); self.addCleanup(timer.cancel)
        with self.assertRaisesRegex(RuntimeError,'Campaign stopped'): self.run_fake(binary_hash,stop)
        record = json.loads((self.root/'execution.json').read_text())
        self.assertEqual(record['status'],'interrupted_or_error'); self.assertLess(record['returncode'],0)


if __name__ == '__main__': unittest.main()
