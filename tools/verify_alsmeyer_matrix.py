"""Verify the complete prepared reference matrix; never imply solver completion."""
import argparse
import json
from pathlib import Path
from prepare_alsmeyer_matrix import reference_configuration
from physical_config import prepare
from export_alsmeyer_sparta import render
from shock_validation import VARIANTS, digest


def read(path):
    return json.loads(path.read_text())


def sources(record, names):
    if set(record) != set(names):
        raise ValueError('Incomplete source provenance')
    for name in names:
        if record[name] != digest(Path(__file__).with_name(name)):
            raise ValueError('Source version changed: '+name)


def verify(mc_root, sparta_root, evidence_path):
    mc_root, sparta_root = mc_root.resolve(), sparta_root.resolve()
    evidence = read(evidence_path); prepared = read(mc_root/'preparation.json')
    exported = read(sparta_root/'export.json')
    if prepared['status'] != 'prepared_only' or exported['status'] != 'exported_only' or exported['quick'] is not False:
        raise ValueError('Production preparation/export required; quick inputs are not production')
    if prepared['evidence_sha256'] != digest(evidence_path) or prepared['primary_sha256'] != evidence['primary_source']['sha256']:
        raise ValueError('Primary evidence changed')
    if exported['preparation_sha256'] != digest(mc_root/'preparation.json'):
        raise ValueError('SPARTA refers to a different matrix')
    sources(prepared['sources'], ('prepare_alsmeyer_matrix.py', 'shock_validation.py', 'normal_shock.py', 'physical_config.py', 'argon_vhs.py'))
    sources(exported['sources'], ('export_alsmeyer_sparta.py', 'sparta_preflight.py', 'argon_vhs.py', 'shock_validation.py'))
    seeds = sorted({c['seed'] for c in prepared['cases']})
    identities = {(c['variant'], c['seed']) for c in prepared['cases']}
    if len(seeds) != 3 or len(prepared['cases']) != 21 or identities != {(v, s) for v in VARIANTS for s in seeds}:
        raise ValueError('Complete A-G by three distinct seeds required')
    mc_cases = {}; resolved = set()
    for case in prepared['cases']:
        name = f"{case['variant']}-{case['seed']}"; folder = mc_root/name
        if case['folder'] != name or folder.resolve() in resolved:
            raise ValueError('Case path or independence mismatch')
        resolved.add(folder.resolve())
        for filename in ('input_si.json', 'solver.json', 'units.json', 'reference.json'):
            if digest(folder/filename) != case['hashes'][filename]:
                raise ValueError('MC3D input hash mismatch: '+name+'/'+filename)
        expected_data, expected_ref = reference_configuration(case['variant'], case['seed'], evidence)
        if read(folder/'input_si.json') != expected_data or read(folder/'reference.json') != expected_ref:
            raise ValueError('Reference-domain configuration mismatch: '+name)
        expected_cfg, expected_units = prepare(expected_data)
        cfg, units = read(folder/'solver.json'), read(folder/'units.json')
        if Path(cfg.pop('units_file')).resolve() != (folder/'units.json').resolve() or Path(units.pop('source_file')).resolve() != (folder/'input_si.json').resolve():
            raise ValueError('Input file linkage mismatch: '+name)
        if cfg != expected_cfg or units != expected_units:
            raise ValueError('Normalized configuration or units mismatch: '+name)
        mc_cases[name] = (case, expected_data, expected_ref, units)
    wanted = {f'{v}-{s}' for v in ('C', 'D', 'E') for s in seeds}
    if len(exported['cases']) != 9 or {c['folder'] for c in exported['cases']} != wanted:
        raise ValueError('Complete SPARTA C/D/E by the same seeds required')
    for case in exported['cases']:
        name = case['folder']; folder = sparta_root/name
        mc, data, reference, units = mc_cases[name]
        if case['source_hashes'] != mc['hashes']:
            raise ValueError('SPARTA input lineage mismatch: '+name)
        files, metadata = render(data, reference, units['reference_gas'], quick=False)
        metadata['mc3d_input_hashes'] = mc['hashes']
        for filename in (*files, 'reference.json'):
            if digest(folder/filename) != case['hashes'][filename]:
                raise ValueError('SPARTA input hash mismatch: '+name+'/'+filename)
        if any((folder/f).read_text() != text for f, text in files.items()) or read(folder/'reference.json') != metadata:
            raise ValueError('SPARTA export is not equivalent: '+name)
    return dict(status='prepared_inputs_verified', physical_acceptance='incomplete',
                mc3d_cases=21, sparta_cases=9, seeds=seeds,
                evidence_sha256=digest(evidence_path), preparation_sha256=digest(mc_root/'preparation.json'),
                export_sha256=digest(sparta_root/'export.json'), tool_sha256=digest(Path(__file__)),
                limitations=['Solver execution and output completeness have not been verified.',
                             'Three seeds do not supply the independent-group replication required by the candidate estimator.',
                             'Method qualification, numerical convergence and experimental acceptance remain separate requirements.'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mc3d', type=Path, required=True); p.add_argument('--sparta', type=Path, required=True)
    p.add_argument('--evidence', type=Path, default=Path('validation/stage9b/experiment.json'))
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    result = verify(a.mc3d, a.sparta, a.evidence)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'mc3d_cases', 'sparta_cases', 'physical_acceptance')}))


if __name__ == '__main__': main()
