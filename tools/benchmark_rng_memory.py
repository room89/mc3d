#!/usr/bin/env python3
"""Compare initialization time and peak RSS of two solver binaries on empty grids.

Each measurement uses a fresh Python helper with exactly one solver child, so
RUSAGE_CHILDREN's peak is not inherited from earlier measurements. RSS includes
cells and runtime overhead, not only RNG. Supported on Linux and macOS.
"""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time


def measure(binary, config, log):
    start = time.perf_counter()
    with Path(log).open("w") as output:
        subprocess.run([binary, "--config", config], cwd=Path(config).parent,
                       stdout=output, stderr=output, check=True)
    elapsed = time.perf_counter() - start
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return dict(elapsed_s=elapsed, peak_rss_bytes=rss if sys.platform == "darwin" else rss * 1024)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--after", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--measure", nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if sys.platform not in ("darwin", "linux"):
        parser.error("RSS measurement supports macOS and Linux")
    if args.measure:
        print(json.dumps(measure(*args.measure)))
        return
    if not all((args.before, args.after, args.out)):
        parser.error("--before, --after and --out are required")
    rows = []
    for side in (32, 64, 80):
        for name, binary in (("before", args.before), ("after", args.after)):
            directory = args.out.resolve() / f"{name}-{side}"
            directory.mkdir(parents=True, exist_ok=True)
            cfg = dict(ncx=side, ncy=side, ncz=side, np=0, seed=123,
                       thread_pool_size=1, geometry_type="none", end_time=0,
                       snapshot_interval=1, write_default_snapshot_before_compute=False,
                       write_default_snapshot_after_compute=False, write_speed_before=False,
                       write_speed_after=False, write_times=False,
                       precompute_snapshot="none", final_cell_snapshot="none")
            config = directory / "config.json"
            config.write_text(json.dumps(cfg, indent=2))
            result = subprocess.check_output(
                [sys.executable, str(Path(__file__).resolve()), "--measure",
                 str(binary.resolve()), str(config), str(directory / "run.log")], text=True)
            row = dict(variant=name, cells=side**3, **json.loads(result))
            rows.append(row)
            print(json.dumps(row), flush=True)
    (args.out / "summary.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
