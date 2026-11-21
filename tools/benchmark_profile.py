#!/usr/bin/env python3
"""Compare direct/profiled MC3D runs; preserve configs, logs and diagnostics."""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--snapshots', action='store_true')
    a = p.parse_args()
    if a.repeats < 1 or a.threads < 1:
        p.error('repeats and threads must be positive')
    a.out = a.out.resolve()
    a.out.mkdir(parents=True, exist_ok=False)
    config = json.loads(a.config.read_text())
    config.update(seed=123, thread_pool_size=a.threads, end_time=.05,
                  snapshot_interval=1, write_default_snapshot_before_compute=False,
                  write_default_snapshot_after_compute=False, write_speed_before=False,
                  write_speed_after=False, write_times=False, precompute_snapshot='none',
                  final_snapshot='none', final_cell_snapshot='final.net' if a.snapshots else 'none')
    wrapper = Path(__file__).with_name('profile_simulation.py')
    rows = {'direct': [], 'profiled': []}
    for repeat in ([0] if a.snapshots else range(-1, a.repeats)):
        variants = ['direct', 'profiled'] if repeat % 2 == 0 else ['profiled', 'direct']
        for variant in variants:
            directory = a.out / f'{repeat}-{variant}'
            directory.mkdir()
            cfg = directory/'config.json'
            cfg.write_text(json.dumps(config, indent=2)+'\n')
            command = [str(a.binary.resolve()), '--config', str(cfg)]
            if variant == 'profiled':
                command = [sys.executable, str(wrapper), '--out', str(directory/'profile'), '--', *command]
            start = time.monotonic()
            with (directory/'run.log').open('wb') as log:
                subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
            elapsed = time.monotonic()-start
            log_path = directory/('profile/run.log' if variant=='profiled' else 'run.log')
            text = log_path.read_text()
            durations = [float(x) for x in re.findall(r'stopwatch:ClusterStep;total_time:([\d.eE+-]+)', text)]
            if not durations:
                raise RuntimeError('No completed steps: '+str(directory))
            if variant == 'profiled':
                metadata = json.loads((directory/'profile/metadata.json').read_text())
                if metadata['collector_error'] or not metadata['availability']['process']['rss_bytes']:
                    raise RuntimeError('Missing usable process metrics: '+str(directory))
            row = dict(elapsed_s=elapsed, steps=len(durations), mean_step_s=statistics.mean(durations))
            (directory/'result.json').write_text(json.dumps(row, indent=2)+'\n')
            if repeat >= 0:
                rows[variant].append(row)
    if a.snapshots:
        direct = (a.out/'0-direct/final.net').read_bytes()
        profiled = (a.out/'0-profiled/final.net').read_bytes()
        if not direct or direct != profiled:
            raise RuntimeError('Snapshots differ')
    result = {'platform': platform.platform(), 'cpu_count': os.cpu_count(), 'runs': rows,
              'median_elapsed_s': {k:statistics.median(x['elapsed_s'] for x in v) for k,v in rows.items()},
              'median_step_s': {k:statistics.median(x['mean_step_s'] for x in v) for k,v in rows.items()},
              'snapshots_identical': True if a.snapshots else None}
    result['overhead_percent'] = 100*(result['median_elapsed_s']['profiled']/result['median_elapsed_s']['direct']-1)
    (a.out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
