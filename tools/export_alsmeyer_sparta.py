"""Export verified prepared SI cases to SPARTA without changing frozen engineering tools."""
import argparse
import json
import math
from pathlib import Path
from argon_vhs import MASS, KB
from sparta_preflight import rate
from shock_validation import digest, save


def render(data, reference, gas, quick=False):
    if gas['model'] != 'vhs' or gas['mass_kg'] != MASS or gas['reference_temperature_k'] != 300:
        raise ValueError('Pinned argon VHS at 300 K required')
    d, omega = gas['diameter_m'], gas['omega']
    if not math.isfinite(d) or d <= 0 or not math.isfinite(omega) or not .5 <= omega <= 1:
        raise ValueError('Invalid VHS parameters')
    num, lengths = data['numerics'], data['domain_m']
    states = [data['reservoirs'][key] for key in ('x_neg', 'x_pos')]
    if data['geometry'] != {'type': 'none'} or data['initial_step']['x_m'] != 0:
        raise ValueError('Planar centered shock required')
    for i, key in enumerate(('left', 'right')):
        if data['initial_step'][key] != states[i] or reference[key] != states[i]:
            raise ValueError('Initial and reference states differ from reservoirs')
    for axis in ('y', 'z'):
        if any(num[f'boundary_{axis}_{side}'] != 'periodic' for side in ('neg', 'pos')):
            raise ValueError('Periodic transverse boundaries required')
    if any(num[f'boundary_x_{side}'] != 'hyperfree' for side in ('neg', 'pos')):
        raise ValueError('Open reservoir boundaries required')
    tau = reference['transit_time_s']
    expected = [[(4+i/8)*tau, (4+(i+1)/8)*tau] for i in range(64)]
    if data['statistics_windows_s'] != expected or data['end_time_s'] != 12*tau:
        raise ValueError('Require common 4-tau burn and 64 windows through 12 tau')
    weight = data['number_density_m3']*math.prod(lengths) / (
        math.prod(num[f'cells_{axis}'] for axis in 'xyz')*num['particles_per_cell'])
    cell = min(lengths[i]/num[f'cells_{axis}'] for i, axis in enumerate('xyz'))
    parameters = dict(vhs_diameter_m=d, vhs_omega=omega)
    flight_dt = num['cu']*cell/max(abs(s['velocity_m_s'][0])+math.sqrt(2*KB*s['temperature_k']/MASS) for s in states)
    collision_dt = num['collision_courant']/max(s['number_density_m3']*rate(s['temperature_k'], parameters) for s in states)
    perwindow = max(1, math.ceil(tau/8/min(flight_dt, collision_dt)))
    dt = tau/(8*perwindow)
    period = 10 if quick else perwindow
    burn, measure = (10, 20) if quick else (32*perwindow, 64*perwindow)
    lx, ly, lz = lengths
    lines = [f"seed {num['seed']}", 'dimension 3', 'boundary o p p',
             f'create_box {-lx/2} {lx/2} {-ly/2} {ly/2} {-lz/2} {lz/2}',
             f"create_grid {num['cells_x']} {num['cells_y']} {num['cells_z']}",
             'species argon.species Ar', f'global fnum {weight:.17g}']
    for name, state in zip(('up', 'down'), states):
        lines.append(f"mixture {name} Ar nrho {state['number_density_m3']:.17g} temp {state['temperature_k']:.17g} vstream " +
                     ' '.join(map(str, state['velocity_m_s'])))
    for name, lo, hi in [('up', -lx/2, 0), ('down', 0, lx/2)]:
        lines += [f'region {name} block {lo} {hi} {-ly/2} {ly/2} {-lz/2} {lz/2}',
                  f'create_particles {name} n 0 region {name}']
    lines += ['collide vss all argon.vss', 'fix inlet emit/face up xlo', 'fix outlet emit/face down xhi',
              f'timestep {dt:.17g}', 'stats 100', 'stats_style step cpu np ncoll', f'run {burn}',
              'compute moments grid all all nrho u v w usq vsq wsq', 'compute thermal thermal/grid all all temp',
              f'fix avg ave/grid all 1 {period} {period} c_moments[*] c_thermal[*]',
              f'dump profiles grid all {period} profile.* id xlo xhi ylo yhi zlo zhi f_avg[*]',
              'dump_modify profiles format float %.17g', f'run {measure}']
    files = {'argon.species': f'Ar 39.948 {MASS:.17g} 0 0 0 0 0 1 0\n',
             'argon.vss': f'Ar {d:.17g} {omega:.17g} 300 1\n', 'in.shock': '\n'.join(lines)+'\n'}
    metadata = dict(source=data, reference=reference, fnum=weight, timestep=dt,
                    steps_per_window=perwindow, dump_period_steps=period, burn_steps=burn,
                    measurement_steps=measure, quick=quick, flight_dt_limit=flight_dt,
                    collision_dt_limit=collision_dt, physical_acceptance='incomplete')
    return files, metadata


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--variants', nargs='+', choices=['C', 'D', 'E'], default=['C', 'D', 'E'])
    p.add_argument('--seeds', type=int, nargs='+')
    p.add_argument('--quick', action='store_true', help='Separate 30-step syntax smoke inputs, not physical validation')
    a = p.parse_args(); root = a.prepared.resolve(); manifest = root/'preparation.json'
    report = json.loads(manifest.read_text()); manifest_hash = digest(manifest)
    if report['status'] != 'prepared_only' or len(set(a.variants)) != len(a.variants):
        raise ValueError('Distinct variants and a prepared matrix required')
    for name, value in report['sources'].items():
        if digest(Path(__file__).with_name(name)) != value:
            raise ValueError('Preparation source changed: '+name)
    seeds = sorted({c['seed'] for c in report['cases']}) if a.seeds is None else a.seeds
    if len(set(seeds)) != len(seeds) or any(not 0 < s < 2**31 for s in seeds):
        raise ValueError('Distinct positive seeds required')
    cases = []; seen = set()
    for case in report['cases']:
        if case['variant'] not in a.variants or case['seed'] not in seeds: continue
        name = f"{case['variant']}-{case['seed']}"
        if case['folder'] != name or name in seen:
            raise ValueError('Invalid or duplicate case identity')
        seen.add(name); folder = root/name
        for filename in ('input_si.json', 'reference.json', 'units.json', 'solver.json'):
            if digest(folder/filename) != case['hashes'][filename]:
                raise ValueError('Prepared input changed: '+name+'/'+filename)
        data = json.loads((folder/'input_si.json').read_text())
        ref = json.loads((folder/'reference.json').read_text())
        units = json.loads((folder/'units.json').read_text())
        if data['numerics']['seed'] != case['seed'] or ref['experimental_reference']['primary_sha256'] != report['primary_sha256']:
            raise ValueError('Mismatched seed or primary source')
        files, metadata = render(data, ref, units['reference_gas'], a.quick)
        cases.append((name, files, metadata, case['hashes']))
    if len(cases) != len(a.variants)*len(seeds):
        raise ValueError('Selected cases missing from preparation')
    output = a.output_dir.resolve(); output.mkdir(parents=True, exist_ok=False); entries = []
    for name, files, metadata, hashes in cases:
        folder = output/name; folder.mkdir()
        for filename, text in files.items(): (folder/filename).write_text(text)
        metadata['mc3d_input_hashes'] = hashes; save(folder/'reference.json', metadata)
        entries.append(dict(folder=name, source_hashes=hashes,
                            hashes={f: digest(folder/f) for f in [*files, 'reference.json']}))
    if digest(manifest) != manifest_hash: raise ValueError('Preparation changed during export')
    save(output/'export.json', dict(status='exported_only', physical_acceptance='incomplete',
         quick=a.quick, preparation_sha256=manifest_hash,
         sources={n: digest(Path(__file__).with_name(n)) for n in ('export_alsmeyer_sparta.py', 'sparta_preflight.py', 'argon_vhs.py', 'shock_validation.py')}, cases=entries))
    print(f'Exported {len(cases)} SPARTA inputs; no solver launched; quick={a.quick}.')


if __name__ == '__main__': main()
