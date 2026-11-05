"""Run maintained short examples, including the SI preparation path."""
import json
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
binary=Path(sys.argv[1]).resolve();repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mc3d examples ') as tmp:
    root=Path(tmp);manifests={}
    for config in sorted((repo/'configs/checks').glob('*.json')):
        out=root/config.stem
        r=subprocess.run([str(binary),'--config',str(config),'--output-dir',str(out)],capture_output=True)
        assert r.returncode==0,(config,r.stderr.decode())
        m=json.loads((out/'manifest.json').read_text());assert m['status']=='completed';manifests[config.stem]=m
        assert (out/'initial.dat').is_file() and (out/'final.dat').is_file()
    # Identical HS law, RNG and configuration must produce identical physical data.
    def hashes(m):return {x['path']:x['sha256'] for x in m['outputs'] if x['kind']=='physical'}
    assert hashes(manifests['hs_low_population'])==hashes(manifests['vhs_hs_limit'])
    prep=root/'prepared'
    subprocess.run([sys.executable,str(repo/'tools/physical_config.py'),str(repo/'configs/physical/argon_cube.json'),'--output-dir',str(prep)],check=True,capture_output=True)
    out=root/'physical'
    r=subprocess.run([str(binary),'--config',str(prep/'solver.json'),'--output-dir',str(out)],capture_output=True)
    assert r.returncode==0,r.stderr.decode()
    m=json.loads((out/'manifest.json').read_text());assert m['status']=='completed' and m['units_mode']=='si'
    assert (out/'body_loads.csv').is_file() and (out/'final.dat.surface.csv').is_file()
    with (out/'initial.dat.surface.csv').open() as f:before=list(csv.DictReader(f,delimiter=';'))
    with (out/'final.dat.surface.csv').open() as f:after=list(csv.DictReader(f,delimiter=';'))
    scales=m['units']['scales']
    for a,b in zip(before,after):
        heat=(float(b['absorbed_heat'])-float(a['absorbed_heat']))*scales['energy_j']
        capacity=float(b['heat_capacity'])*scales['heat_capacity_j_k']
        delta=(float(b['temperature'])-float(a['temperature']))*scales['temperature_k']
        assert math.isclose(heat,capacity*delta,rel_tol=1e-9,abs_tol=1e-22),(heat,capacity*delta)
print('Configuration examples passed')
