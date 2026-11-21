"""Compare statistics output cost and final results using two Release binaries."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import time

from normal_shock import setup
from physical_config import prepare


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    config, _ = prepare(setup()[0])
    # Fixed worker count keeps the exact-output comparison independent of scheduling.
    config.update(end_time=.02, snapshot_interval=0, resolution_diagnostics=False,
                  thread_pool_size=1,
                  statistics_windows=[[i*.02/64, (i+1)*.02/64] for i in range(64)])
    # Use normalized values directly; no unit sidecar is needed for this benchmark.
    config.pop('units_file', None)
    result = {'config': config, 'variants': {}, 'outputs_identical': True}
    expected = None
    for label, binary, interval in [('baseline', args.baseline, 0),
                                    ('cached', args.candidate, 0),
                                    ('cached_sparse', args.candidate, .002)]:
        binary = binary.resolve()
        input_path = root / (label + '.json')
        input_path.write_text(json.dumps(config | {'snapshot_interval': interval}))
        durations = []
        for repeat in range(args.repeats):
            output = root / f'{label}-{repeat}'
            with (root / f'{label}-{repeat}.log').open('w') as log:
                start = time.monotonic()
                subprocess.run([str(binary), '--config', str(input_path),
                                '--output-dir', str(output)], stdout=log,
                               stderr=subprocess.STDOUT, check=True, timeout=300)
                durations.append(time.monotonic()-start)
            paths = [p for pattern in ('*statistics_*.csv', 'conservation*', 'final.dat*')
                     for p in output.glob(pattern)]
            hashes = {p.name: digest(p) for p in paths}
            if expected is None:
                expected = hashes
            if hashes != expected:
                raise AssertionError(f'Final results differ: {label}-{repeat}')
        result['variants'][label] = dict(binary_sha256=digest(binary), seconds=durations,
                                         median_seconds=statistics.median(durations),
                                         compared_files=len(hashes))
    baseline = result['variants']['baseline']['median_seconds']
    for variant in result['variants'].values():
        variant['speedup'] = baseline / variant['median_seconds']
    (root / 'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['variants'], indent=2))


if __name__ == '__main__':
    main()
