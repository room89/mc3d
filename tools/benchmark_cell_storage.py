#!/usr/bin/env python3
"""Linux paired empty-grid memory/initialization probe; binaries built from benchmark_cell_storage.cpp."""
import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
rows = []
variants = [('baseline', a.baseline.resolve()), ('current', a.binary.resolve())]
for grid in (40, 80):
    runs = {name: [] for name, _ in variants}
    for repeat in range(5):
        for name, binary in variants if repeat % 2 == 0 else reversed(variants):
            stem = a.out / f'{grid}-{repeat}-{name}'
            rss = stem.with_suffix('.rss')
            with stem.with_suffix('.log').open('w') as log:
                subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(rss), str(binary), str(grid)], stdout=log, stderr=subprocess.STDOUT, check=True)
            text = stem.with_suffix('.log').read_text()
            values = {key: float(re.search(rf'{key}=([0-9.]+)', text)[1]) for key in ('cell_bytes', 'seeded_cell_allocations', 'initialization_ms')}
            values['peak_rss_kib'] = int(rss.read_text().strip())
            runs[name].append(values)
    for name, values in runs.items():
        rows.append(dict(grid=grid, variant=name, runs=values, medians={key: statistics.median(row[key] for row in values) for key in values[0]}))
(a.out / 'summary.json').write_text(json.dumps(rows, indent=2))
print(json.dumps(rows, indent=2))
