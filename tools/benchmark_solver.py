#!/usr/bin/env python3
"""Paired solver measurements. Initialization is included only in process wall time.

Parallel trajectories and step counts can differ even with the same seed.
Use --rss on Linux; diagnostic snapshots and warmups never enter measured summaries.
"""
import argparse
import csv
import json
from pathlib import Path
import re
import statistics
import subprocess
import time

NUMBER = r'[0-9.eE+\-]+'


def parse_log(log):
    totals = {}
    for name, duration in re.findall(r'stopwatch:([^;]+);total_time:(' + NUMBER + ')', log):
        totals[name] = totals.get(name, 0.0) + float(duration)
    steps = len(re.findall(r'text=Time step:', log))
    times = re.findall(r'dt = (' + NUMBER + r') t = (' + NUMBER + ')', log)
    if not steps or 'ClusterStep' not in totals or len(times) != steps:
        raise ValueError('Missing steps, ClusterStep or model time in solver log')
    dt, previous_time = map(float, times[-1])
    return dict(steps=steps, model_time=previous_time + dt,
                cluster_total_s=totals['ClusterStep'],
                phase_ms_per_step={k: v * 1000 / steps for k, v in totals.items()})


def summarize(runs):
    phases = set().union(*(r['phase_ms_per_step'] for r in runs))
    result = {'runs': runs, 'median_phase_ms_per_step': {}, 'range_phase_ms_per_step': {}}
    for key in ('elapsed_s', 'cluster_total_s', 'steps', 'model_time', 'peak_rss_kib'):
        values = [r[key] for r in runs if r.get(key) is not None]
        result['median_' + key] = statistics.median(values) if len(values) == len(runs) else None
        result['range_' + key] = [min(values), max(values)] if values else None
    for key in sorted(phases):
        values = [r['phase_ms_per_step'][key] for r in runs if key in r['phase_ms_per_step']]
        result['median_phase_ms_per_step'][key] = statistics.median(values) if len(values) == len(runs) else None
        result['range_phase_ms_per_step'][key] = [min(values), max(values)]
    return result


def particle_counts(directory):
    with (directory / 'initial.csv').open() as stream:
        initial = sum(int(row['N']) for row in csv.DictReader(stream, delimiter=';'))
    final_rows = re.findall(r'^\s*particles\s+(\d+)\s*$', (directory / 'final.net').read_text(), re.M)
    if not final_rows:
        raise ValueError('Missing final particle counts')
    return dict(initial_particles=initial, final_particles=sum(map(int, final_rows)))


def run(binary, directory, cfg, rss=False, diagnostic=False):
    directory.mkdir(parents=True, exist_ok=False)
    config = directory / 'config.json'
    config.write_text(json.dumps(cfg, indent=2))
    command = [str(binary), '--config', str(config)]
    if rss:
        command = ['/usr/bin/time', '-f', '%M', '-o', str(directory / 'rss.txt')] + command
    start = time.perf_counter()
    with (directory / 'run.log').open('w') as log:
        completed = subprocess.run(command, cwd=directory, stdout=log, stderr=log)
    elapsed = time.perf_counter() - start
    if completed.returncode:
        (directory / 'failure.json').write_text(json.dumps({'returncode': completed.returncode}))
        raise RuntimeError(f'Solver failed: {directory}')
    result = parse_log((directory / 'run.log').read_text())
    result['elapsed_s'] = elapsed
    result['peak_rss_kib'] = int((directory / 'rss.txt').read_text().strip()) if rss else None
    if diagnostic:
        result.update(particle_counts(directory))
    (directory / 'result.json').write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--grid', type=int, nargs=3, default=(24, 12, 8))
    parser.add_argument('--particles', type=int, default=60)
    parser.add_argument('--threads', type=int, nargs='+', default=(1, 4))
    parser.add_argument('--geometries', choices=('none', 'cube', 'wedge'), nargs='+', default=('none', 'cube', 'wedge'))
    parser.add_argument('--seeds', type=int, nargs='+', default=(123,))
    parser.add_argument('--warmups', type=int, default=0)
    parser.add_argument('--rss', action='store_true')
    parser.add_argument('--diagnostics', action='store_true')
    args = parser.parse_args()
    if min(args.grid) < 1 or args.particles < 0 or min(args.threads) < 1 or args.repeats < 1 or args.warmups < 0:
        parser.error('Invalid grid, particles, threads, repeats or warmups')
    if len(set(args.seeds)) != len(args.seeds) or any(s < 0 or s > 2**32-1 for s in args.seeds):
        parser.error('Seeds must be unique uint32 values')
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    template = json.loads((Path(__file__).resolve().parents[1] / 'configs/example.json').read_text())
    variants = [('current', args.binary.resolve())]
    if args.baseline:
        variants.insert(0, ('baseline', args.baseline.resolve()))
    results = []
    for geometry in args.geometries:
        for threads in args.threads:
            for seed in args.seeds:
                cfg = dict(template, geometry_type=geometry, thread_pool_size=threads,
                           ncx=args.grid[0], ncy=args.grid[1], ncz=args.grid[2], np=args.particles,
                           seed=seed, end_time=.05, snapshot_interval=1,
                           write_default_snapshot_before_compute=False, write_default_snapshot_after_compute=False,
                           write_speed_before=False, write_speed_after=False, write_times=False,
                           precompute_snapshot='none', final_cell_snapshot='none')
                prefix = f'{geometry}-{threads}-seed{seed}'
                for warmup in range(args.warmups):
                    for name, binary in variants:
                        run(binary, output / f'{prefix}-warmup{warmup}-{name}', cfg)
                runs = {name: [] for name, _ in variants}
                for repeat in range(args.repeats):
                    for name, binary in variants if repeat % 2 == 0 else reversed(variants):
                        runs[name].append(run(binary, output / f'{prefix}-{repeat}-{name}', cfg, args.rss))
                for name, binary in variants:
                    summary = dict(variant=name, geometry=geometry, threads=threads, seed=seed,
                                   grid=args.grid, particles_per_cell=args.particles, **summarize(runs[name]))
                    if args.diagnostics:
                        summary['diagnostic'] = run(binary, output / f'{prefix}-diagnostic-{name}',
                                                   dict(cfg, precompute_snapshot='initial.csv', final_cell_snapshot='final.net'),
                                                   diagnostic=True)
                    results.append(summary)
                    print(json.dumps({k: v for k, v in summary.items() if k != 'runs'}), flush=True)
                    (output / 'summary.json').write_text(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
