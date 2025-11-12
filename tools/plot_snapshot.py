#!/usr/bin/env python3
"""
Utility for visualising MC3D simulation snapshots.

The solver writes CSV-like files (e.g. data.dat, data<t>.dat) with the columns:
    x; y; z; N; ro; T; vx; vy; vz; E

This script filters cells by a chosen plane and renders a scatter plot coloured
by the selected field.  Example:
    python3 tools/plot_snapshot.py --input build/data.dat --plane xy \
        --z-value 0 --thickness 0.05 --field T --show-cbar --out plots/T_xy.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from generate_gif import PLANE_AXES, PLANE_VELOCITY_COMPONENTS, draw_plot, filter_slice
from snapshot_loader import load_snapshot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot MC3D snapshot.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("build/data.dat"),
        help="Path to snapshot file (CSV with ';' separator).",
    )
    parser.add_argument(
        "--plane",
        choices=("xy", "xz", "yz"),
        default="xy",
        help="Plane to display.",
    )
    parser.add_argument(
        "--z-value",
        type=float,
        default=0.0,
        help="Coordinate of the plane centre measured along the excluded axis.",
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=0.05,
        help="Slice thickness around the requested coordinate.",
    )
    parser.add_argument(
        "--field",
        choices=("ro", "T", "vx", "vy", "vz", "E"),
        default="ro",
        help="Field to visualise.",
    )
    parser.add_argument(
        "--show-velocity",
        action="store_true",
        help="Overlay velocity vectors corresponding to the selected plane.",
    )
    parser.add_argument(
        "--quiver-step",
        type=int,
        default=1,
        help="Plot every N-th cell velocity vector (use >1 to declutter).",
    )
    parser.add_argument(
        "--quiver-scale",
        type=float,
        default=None,
        help="Scale factor passed to matplotlib.quiver (smaller values draw longer arrows).",
    )
    parser.add_argument(
        "--quiver-color",
        default="k",
        help="Colour of the velocity arrows.",
    )
    parser.add_argument(
        "--show-cbar",
        action="store_true",
        help="Show colour bar on the plot.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Dots-per-inch for the saved figure (ignored when showing interactively).",
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=None,
        help="Figure size in inches (width height). Defaults to 7 5.5.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional path to save the figure (PNG/PDF/etc). Displays window if omitted.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional custom plot title.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    snapshot = load_snapshot(args.input)
    filtered = filter_slice(snapshot, args.plane, args.z_value, args.thickness)
    figsize = tuple(args.figsize) if args.figsize else None
    fig = draw_plot(
        filtered,
        args.plane,
        args.field,
        args.title,
        args.show_cbar,
        args.show_velocity,
        args.quiver_step,
        args.quiver_scale,
        args.quiver_color,
        figsize=figsize,
    )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()

