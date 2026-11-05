import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "control_simulation", Path(__file__).resolve().parents[2] / "tools/control_simulation.py")
control = importlib.util.module_from_spec(spec)
spec.loader.exec_module(control)


class ControlSimulationTest(unittest.TestCase):
    def test_replace_and_invalid_request_preserves_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "load.control"
            for mode in ("quiet", "pause", "balanced", "full"):
                control.set_mode(path, mode)
                self.assertEqual(path.read_text(), mode + "\n")
            with self.assertRaises(ValueError):
                control.set_mode(path, "unknown")
            self.assertEqual(path.read_text(), "full\n")
            self.assertEqual(list(path.parent.iterdir()), [path])
