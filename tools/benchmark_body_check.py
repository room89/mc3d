#!/usr/bin/env python3
"""Linux paired body-check benchmark: warmups, phases, RSS and snapshot equivalence."""
import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--config', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--threads', type=int, nargs='+', default=[1, 4, 12])
p.add_argument('--repeats', type=int, default=5)
p.add_argument('--snapshots', action='store_true')
a = p.parse_args()
if a.repeats < 1 or min(a.threads) < 1:
    p.error('repeats and threads must be positive')
out = a.out.resolve()
out.mkdir(parents=True, exist_ok=True)
variants = [('baseline', a.baseline.resolve()), ('current', a.binary.resolve())]
template = json.loads(a.config.read_text())
rows = []
for workers in a.threads:
    runs = {name: [] for name, _ in variants}
    for repeat in range(-1, a.repeats):
        if a.snapshots and repeat != 0:
            continue
        for name, binary in variants if repeat % 2 == 0 else reversed(variants):
            directory = out / f'{workers}-{repeat}-{name}'
            directory.mkdir(exist_ok=False)
            cfg = dict(template, seed=123, thread_pool_size=workers, end_time=.05,
                       snapshot_interval=1, write_default_snapshot_before_compute=False,
                       write_default_snapshot_after_compute=False, write_speed_before=False,
                       write_speed_after=False, write_times=False, precompute_snapshot='none',
                       final_snapshot='none', final_cell_snapshot='final.net' if a.snapshots else 'none')
            config = directory / 'config.json'
            config.write_text(json.dumps(cfg, indent=2))
            start = time.perf_counter()
            with (directory / 'run.log').open('w') as log:
                subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(directory/'rss.txt'),
                                str(binary), '--config', str(config)], cwd=directory,
                               stdout=log, stderr=log, check=True)
            elapsed = time.perf_counter()-start
            text = (directory/'run.log').read_text()
            phases = {}
            for phase, duration in re.findall(r'stopwatch:([^;]+);total_time:([0-9.eE+-]+)', text):
                phases[phase] = phases.get(phase, 0) + float(duration)
            steps = len(re.findall(r'text=Time step:', text))
            if not steps or 'ClusterStep' not in phases:
                raise RuntimeError(f'Incomplete log: {directory}')
            result = dict(steps=steps, elapsed_s=elapsed, peak_rss_kib=int((directory/'rss.txt').read_text()),
                          phase_ms_per_step={k: v*1000/steps for k,v in phases.items()})
            (directory/'result.json').write_text(json.dumps(result, indent=2))
            if repeat >= 0:
                runs[name].append(result)
    if a.snapshots:
        old = (out/f'{workers}-0-baseline/final.net').read_bytes()
        new = (out/f'{workers}-0-current/final.net').read_bytes()
        if not old or old != new:
            raise RuntimeError('Snapshots differ')
        print('Snapshots identical:', len(old), flush=True)
    for name, values in runs.items():
        row = dict(variant=name, threads=workers, runs=values,
                   median_elapsed_s=statistics.median(r['elapsed_s'] for r in values),
                   median_peak_rss_kib=statistics.median(r['peak_rss_kib'] for r in values),
                   median_phase_ms_per_step={k: statistics.median(r['phase_ms_per_step'][k] for r in values)
                                             for k in values[0]['phase_ms_per_step']})
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='runs'}), flush=True)
    (out/'summary.json').write_text(json.dumps(rows, indent=2))
