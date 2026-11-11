"""Prepare a separate Alsmeyer-reference A-G matrix; never execute or accept it."""
import argparse
import copy
import json
from pathlib import Path
from argon_vhs import KB
from physical_config import prepare
from shock_validation import VARIANTS, configuration, digest, save


def reference_configuration(variant, seed, evidence):
    primary = evidence['primary_source']
    source = primary['verified_evidence']['printed_page_499']
    temperature = source['reference_temperature_K']
    pressure = source['reference_pressure_mTorr'] * 1e-3 * 101325 / 760
    length = source['reference_mean_free_path_m']
    if temperature != 300 or pressure <= 0 or length <= 0 or len(primary['sha256']) != 64:
        raise ValueError('Verified 300 K primary reference conditions required')
    # JSON roundtrip removes shared reservoir/initial-state object aliases.
    data, old_ref = json.loads(json.dumps(configuration(variant, seed)))
    density = pressure / (KB * temperature)
    factor = density / data['number_density_m3']
    data['number_density_m3'] = density
    data['domain_m'] = [x / factor for x in data['domain_m']]
    base_length = 40 * length
    data['domain_m'][0] = base_length * (2 if variant == 'F' else 1)
    data['initial_step']['x_m'] /= factor
    for state in list(data['reservoirs'].values()) + [data['initial_step'][key] for key in ('left', 'right')]:
        state['number_density_m3'] *= factor
    tau = base_length / data['speed_m_s']
    data['end_time_s'] = 12 * tau
    data['snapshot_interval_s'] = .2 * tau
    data['statistics_windows_s'] = [[(4+i/8)*tau, (4+(i+1)/8)*tau] for i in range(64)]
    provenance = dict(source_doi=evidence['source_doi'], primary_sha256=primary['sha256'],
                      primary_printed_page=499, experimental_reference_length_m=length,
                      reference_pressure_Pa=pressure,
                      condition_scope='Reference conditions for normalization, not a recovered individual experimental run record.')
    data['provenance'] = provenance
    reference = dict(description=__doc__, mean_free_path_m=old_ref['mean_free_path_m']/factor,
                     transit_time_s=tau, left=copy.deepcopy(data['reservoirs']['x_neg']),
                     right=copy.deepcopy(data['reservoirs']['x_pos']),
                     experimental_reference=provenance,
                     axial_experimental_mean_free_paths=80 if variant == 'F' else 40,
                     sampling_time_policy='All variants use the same base-domain transit time; F changes only domain length and axial cells.')
    return data, reference


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, default=Path('validation/stage9b/experiment.json'))
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--seeds', type=int, nargs=3, default=[30001, 30002, 30003])
    p.add_argument('--variants', nargs='+', choices=list(VARIANTS), default=list(VARIANTS))
    a = p.parse_args()
    if len(set(a.seeds)) != 3 or any(not 0 < s < 2**31 for s in a.seeds):
        p.error('Three distinct positive 31-bit seeds required')
    if len(set(a.variants)) != len(a.variants):
        p.error('Distinct variants required')
    evidence = json.loads(a.reference.read_text())
    cases = [(v, seed, *reference_configuration(v, seed, evidence)) for v in a.variants for seed in a.seeds]
    root = a.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    entries = []
    for variant, seed, data, ref in cases:
        folder = root / f'{variant}-{seed}'; folder.mkdir()
        cfg, units = prepare(data)
        cfg['units_file'] = str(folder / 'units.json')
        units['source_file'] = str(folder / 'input_si.json')
        for name, obj in [('input_si.json', data), ('reference.json', ref), ('solver.json', cfg), ('units.json', units)]:
            save(folder/name, obj)
        entries.append(dict(variant=variant, seed=seed, folder=folder.name, status='prepared',
                            hashes={name: digest(folder/name) for name in ('input_si.json', 'reference.json', 'solver.json', 'units.json')}))
    names = ('prepare_alsmeyer_matrix.py', 'shock_validation.py', 'normal_shock.py', 'physical_config.py', 'argon_vhs.py')
    save(root/'preparation.json', dict(status='prepared_only', physical_acceptance='incomplete',
         primary_sha256=evidence['primary_source']['sha256'], evidence_sha256=digest(a.reference),
         sources={name: digest(Path(__file__).with_name(name)) for name in names}, cases=entries,
         limitations=['No solver was launched.', 'SPARTA export and execution are pending.',
                      'Independent-group replication, convergence and experimental acceptance remain pending.',
                      'These configurations must not be mixed with the engineering C pilot or its pooled campaign.']))
    print(f'Prepared {len(entries)} separate MC3D configurations; no solver launched.')


if __name__ == '__main__':
    main()
