"""Fresh output per CTest invocation; preserve failures in CTest output."""
import subprocess,sys,tempfile
from pathlib import Path
with tempfile.TemporaryDirectory(prefix='mc3d-analytical-') as d:
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'tools/check_physics.py'),'--binary',sys.argv[1],'--output-dir',str(Path(d)/'runs'),'--load-mode','full'],check=True)
