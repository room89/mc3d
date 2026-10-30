#!/usr/bin/env python3
"""Render three field maps and optionally a temperature GIF from a solver run."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt

from generate_gif import PLANE_AXES, draw_plot
from snapshot_loader import load_snapshot
from plot_body import render_body


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Run directory containing final.dat and data*.dat")
    parser.add_argument("--out", type=Path, help="Output directory (default: RUN/plots)")
    parser.add_argument("--plane", choices=PLANE_AXES, default="xy")
    parser.add_argument("--slice", type=float, default=0.0,
                        help="Select the cell layer nearest this coordinate (default: 0)")
    parser.add_argument("--gif", action="store_true", help="Also generate a temperature animation")
    parser.add_argument("--vmin", type=float, default=0.0, help="GIF temperature scale minimum")
    parser.add_argument("--vmax", type=float, default=12.0, help="GIF temperature scale maximum")
    args = parser.parse_args()
    if args.vmin >= args.vmax:
        parser.error("--vmin must be smaller than --vmax")
    if args.gif and not list(args.run.glob("data*.dat")):
        parser.error(f"No data*.dat snapshots in {args.run}; omit --gif for PNGs only")
    source = args.run / "final.dat"
    if not source.is_file():
        parser.error(f"Missing {source}; wait for the simulation to finish")
    snapshot = load_snapshot(source)
    axis = PLANE_AXES[args.plane][2]
    coords = snapshot[axis].dropna().unique()
    if not len(coords):
        parser.error(f"No cell coordinates in {source}")
    centre = min(coords, key=lambda value: (abs(value - args.slice), -value))
    layer = snapshot.loc[snapshot[axis] == centre]
    # Keep GIF selection on exactly the same layer, including non-default grids.
    gaps = [abs(value - centre) for value in coords if value != centre]
    thickness = min(gaps) if gaps else 1e-6
    out = args.out or args.run / "plots"
    out.mkdir(parents=True, exist_ok=True)
    print(f"Selected {axis}={centre:g} ({len(layer)} cells)", flush=True)
    for field in ("ro", "T", "vx"):
        fig = draw_plot(layer, args.plane, field, None, True, False, 1, None, "k")
        path = out / f"{field}_{args.plane}.png"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(path, flush=True)
    render_body(args.run, out, args.gif)
    if args.gif:
        subprocess.run([
            sys.executable, str(Path(__file__).with_name("generate_gif.py")),
            "--pattern", str(args.run / "data*.dat"),
            "--plane", args.plane, "--z-value", str(centre),
            "--thickness", str(thickness), "--field", "T", "--show-cbar",
            "--vmin", str(args.vmin), "--vmax", str(args.vmax),
            "--out", str(out / "temperature.gif"),
        ], check=True)


if __name__ == "__main__":
    main()
