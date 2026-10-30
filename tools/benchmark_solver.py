#!/usr/bin/env python3
"""Repeat small solver scenarios; keep logs and report phase time per step.

Use the same machine/build configuration for comparisons. Initialization is
included in elapsed time, but not in the ClusterStep phase. Parallel runs may
have different trajectories and step counts even with a fixed seed.
"""
import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, help="Alternate baseline and current binary")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--grid", type=int, nargs=3, default=(24, 12, 8))
    parser.add_argument("--particles", type=int, default=60)
    parser.add_argument("--threads", type=int, nargs="+", default=(1, 4))
    parser.add_argument("--geometries", choices=("none", "cube", "wedge"),
                        nargs="+", default=("none", "cube", "wedge"))
    args = parser.parse_args()
    if min(args.grid) < 1 or args.particles < 0 or min(args.threads) < 1:
        parser.error("grid and threads must be positive; particles must be nonnegative")
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    binary = args.binary.resolve()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).resolve().parents[1] / "configs/example.json"
    variants = [("current", binary)]
    if args.baseline:
        variants.insert(0, ("baseline", args.baseline.resolve()))
    results = []
    for geometry in args.geometries:
        for threads in args.threads:
            variant_runs = {name: [] for name, _ in variants}
            for repeat in range(args.repeats):
                ordered = variants if repeat % 2 == 0 else list(reversed(variants))
                for variant, binary in ordered:
                    directory = output / f"{geometry}-{threads}-{repeat}-{variant}"
                    directory.mkdir(exist_ok=True)
                    cfg = json.loads(source.read_text())
                    cfg.update(geometry_type=geometry, thread_pool_size=threads,
                               ncx=args.grid[0], ncy=args.grid[1], ncz=args.grid[2],
                               np=args.particles, seed=123,
                               end_time=0.05, snapshot_interval=1,
                               write_default_snapshot_before_compute=False,
                               write_default_snapshot_after_compute=False,
                               write_speed_before=False, write_speed_after=False,
                               write_times=False, precompute_snapshot="none",
                               final_cell_snapshot="none")
                    config = directory / "config.json"
                    config.write_text(json.dumps(cfg, indent=2))
                    start = time.perf_counter()
                    with (directory / "run.log").open("w") as log:
                        subprocess.run([str(binary), "--config", str(config)],
                                       cwd=directory, stdout=log, stderr=log, check=True)
                    elapsed = time.perf_counter() - start
                    log = (directory / "run.log").read_text()
                    phases = {}
                    for name, duration in re.findall(
                            r"stopwatch:([^;]+);total_time:([0-9.eE+-]+)", log):
                        phases[name] = phases.get(name, 0.0) + float(duration)
                    steps = len(re.findall(r"text=Time step:", log))
                    if not steps:
                        raise RuntimeError(f"No simulation steps in {directory}")
                    variant_runs[variant].append(dict(elapsed_s=elapsed, steps=steps,
                                     phase_ms_per_step={k: v * 1000 / steps
                                                        for k, v in phases.items()}))
            for variant, runs in variant_runs.items():
                summary = dict(variant=variant, geometry=geometry, threads=threads, runs=runs,
                               median_elapsed_s=statistics.median(r["elapsed_s"] for r in runs),
                               median_phase_ms_per_step={
                                   key: statistics.median(r["phase_ms_per_step"][key] for r in runs)
                                   for key in runs[0]["phase_ms_per_step"]})
                results.append(summary)
                print(json.dumps({k: v for k, v in summary.items() if k != "runs"}), flush=True)
    (output / "summary.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
