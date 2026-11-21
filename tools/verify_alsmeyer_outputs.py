"""Verify solver outputs against prepared inputs; no physical acceptance verdict.

These checks establish output completeness and provenance. The caller must also
run verify_alsmeyer_matrix.verify for the complete input matrix. SPARTA requires
an execution record from the launcher because its dump files contain no binary
or input provenance of their own.
"""
import csv
import json
import math
from pathlib import Path
from shock_validation import digest


def read(path):
    return json.loads(path.read_text())


def close(a, b):
    return math.isfinite(float(a)) and math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-14)


def checked_files(root, entries):
    paths = set()
    for entry in entries:
        name = entry['path']; path = Path(name)
        if path.is_absolute() or '..' in path.parts or name in paths or not path.parts:
            raise ValueError('Unsafe or duplicate artifact path')
        target = root/path
        if not target.resolve().is_relative_to(root.resolve()) or digest(target) != entry['sha256']:
            raise ValueError('Artifact hash mismatch: '+name)
        paths.add(name)
    return paths


def mc3d_outputs(folder, binary_sha256):
    folder = Path(folder); out = folder/'output'; manifest = read(out/'manifest.json')
    cfg = read(folder/'solver.json'); units = read(folder/'units.json')
    if manifest['status'] != 'completed' or manifest.get('error') or manifest['binary_sha256'] != binary_sha256:
        raise ValueError('Completed MC3D run with matching binary required')
    if manifest['seed'] != cfg['seed'] or not close(manifest['simulation_time'], cfg['end_time']):
        raise ValueError('MC3D seed or final time mismatch')
    for key, value in cfg.items():
        if key != 'units_file' and manifest['effective_config'].get(key) != value:
            raise ValueError('MC3D effective parameter differs: '+key)
    if manifest['units'] != units:
        raise ValueError('Actual units differ from prepared units')
    checked_files(out, manifest['inputs'])
    for role, filename in [('config', 'solver.json'), ('physical_input', 'input_si.json'), ('units', 'units.json')]:
        entries = [e for e in manifest['inputs'] if e['role'] == role]
        if len(entries) != 1 or entries[0]['sha256'] != digest(folder/filename):
            raise ValueError('Actual MC3D input differs: '+role)
    output_names = checked_files(out, manifest['outputs'])
    bounds = cfg['statistics_windows']; cells = math.prod(cfg['cells_'+a] for a in 'xyz')
    if len(bounds) != 64:
        raise ValueError('64 production windows required')
    for wi, (start, end) in enumerate(bounds):
        filename = f'cell_statistics_{wi}.csv'
        if filename not in output_names:
            raise ValueError('Statistics missing from output manifest')
        with (out/filename).open() as stream:
            rows = list(csv.DictReader(stream, delimiter=';'))
        if len(rows) != cells or {int(r['cell']) for r in rows} != set(range(cells)):
            raise ValueError('Incomplete MC3D cell grid')
        for row in rows:
            if row['complete'] != '1' or row['status'] != 'valid' or int(row['window']) != wi:
                raise ValueError('Incomplete or invalid MC3D window')
            if any(not close(row[k], v) for k, v in [('start', start), ('end', end), ('observed_time', end-start)]):
                raise ValueError('MC3D observation interval differs')
            if any(not math.isfinite(float(v)) for k, v in row.items() if k != 'status'):
                raise ValueError('Nonfinite MC3D moments')
            if any(float(row[k]) <= 0 for k in ('Vgas', 'A0', 'density')) or float(row['Tstat']) < 0:
                raise ValueError('Invalid MC3D population or temperature')
    return dict(status='outputs_verified', windows=64, cells=cells,
                manifest_sha256=digest(out/'manifest.json'), physical_acceptance='incomplete')


def sparta_dump(path, reference, step):
    lines = path.read_text().splitlines()
    fields = ['id', 'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi'] + [f'f_avg[{i}]' for i in range(1, 9)]
    num = reference['source']['numerics']; nx, ny, nz = (num['cells_'+a] for a in 'xyz')
    cells = nx*ny*nz
    if (len(lines) < 9 or lines[:1] != ['ITEM: TIMESTEP'] or int(lines[1]) != step or
            lines[2] != 'ITEM: NUMBER OF CELLS' or int(lines[3]) != cells or
            lines[4] != 'ITEM: BOX BOUNDS oo pp pp' or lines[8].split() != ['ITEM:', 'CELLS', *fields]):
        raise ValueError('SPARTA dump header, time or grid mismatch')
    rows = [list(map(float, line.split())) for line in lines[9:] if line.strip()]
    if len(rows) != cells or any(len(r) != len(fields) or not all(map(math.isfinite, r)) for r in rows):
        raise ValueError('Incomplete or nonfinite SPARTA grid')
    if {r[0] for r in rows} != set(range(1, cells+1)):
        raise ValueError('SPARTA cell identities differ')
    lengths = reference['source']['domain_m']; coords = set()
    for row in rows:
        index = []
        for axis, (length, count) in enumerate(zip(lengths, (nx, ny, nz))):
            lo, hi = row[1+2*axis:3+2*axis]; dx = length/count
            i = round((lo+length/2)/dx)
            if not 0 <= i < count or not close(lo, -length/2+i*dx) or not close(hi, -length/2+(i+1)*dx):
                raise ValueError('SPARTA cell coordinates differ')
            index.append(i)
        coords.add(tuple(index))
        if row[7] <= 0 or row[-1] < 0:
            raise ValueError('Invalid SPARTA density or temperature')
    if len(coords) != cells:
        raise ValueError('Duplicate SPARTA cell coordinates')
    return cells


def sparta_outputs(folder, execution, binary_sha256):
    folder = Path(folder); ref = read(folder/'reference.json')
    if ref['quick'] or execution['returncode'] != 0 or execution['binary_sha256'] != binary_sha256:
        raise ValueError('Successful production SPARTA execution required')
    for name in ('in.shock', 'argon.species', 'argon.vss', 'reference.json'):
        if execution['input_hashes'].get(name) != digest(folder/name):
            raise ValueError('SPARTA execution input differs: '+name)
    step = ref['steps_per_window']
    if ref['burn_steps'] != 32*step or ref['measurement_steps'] != 64*step or ref['dump_period_steps'] != step:
        raise ValueError('SPARTA production sampling differs')
    names = checked_files(folder, execution['outputs'])
    if execution['log_sha256'] != digest(folder/'run.log'):
        raise ValueError('SPARTA execution log differs')
    for wi in range(64):
        filename = f'profile.{(33+wi)*step}'
        if filename not in names:
            raise ValueError('SPARTA measurement missing from execution record')
        cells = sparta_dump(folder/filename, ref, (33+wi)*step)
    return dict(status='outputs_verified', windows=64, cells=cells,
                physical_acceptance='incomplete')


def main():
    import argparse
    from verify_alsmeyer_matrix import verify
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mc3d', type=Path, required=True); p.add_argument('--sparta', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True); p.add_argument('--sparta-bin', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    inputs = verify(a.mc3d, a.sparta, a.evidence)
    binaries = dict(mc3d=digest(a.binary), sparta=digest(a.sparta_bin)); results = {}
    for case in read(a.mc3d/'preparation.json')['cases']:
        name = case['folder']; results['mc3d/'+name] = mc3d_outputs(a.mc3d/name, binaries['mc3d'])
    for case in read(a.sparta/'export.json')['cases']:
        name = case['folder']; folder = a.sparta/name; record = folder/'execution.json'
        result = sparta_outputs(folder, read(record), binaries['sparta'])
        result['execution_sha256'] = digest(record); results['sparta/'+name] = result
    # Revalidate input and executable identities after potentially expensive file reads.
    if inputs != verify(a.mc3d, a.sparta, a.evidence) or binaries != dict(mc3d=digest(a.binary), sparta=digest(a.sparta_bin)):
        raise ValueError('Input or executable changed during verification')
    report = dict(status='matrix_outputs_verified', physical_acceptance='incomplete',
                  inputs=inputs, binaries=binaries, cases=results, tool_sha256=digest(Path(__file__)),
                  limitations=['Output completeness does not establish stationarity, convergence, uncertainty coverage or experimental agreement.',
                               'SPARTA execution provenance is recorded by the launcher, not by the dump format.'])
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x') as f: json.dump(report, f, indent=2, allow_nan=False); f.write('\n')
    print('Verified 21 MC3D and 9 SPARTA outputs; physical acceptance remains incomplete.')


if __name__ == '__main__': main()
