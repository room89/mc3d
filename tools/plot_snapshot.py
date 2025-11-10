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
import matplotlib.tri as mtri
import numpy as np
import pandas as pd


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
        "--show-cbar",
        action="store_true",
        help="Show colour bar on the plot.",
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


PLANE_AXES = {
    "xy": ("x", "y", "z"),
    "xz": ("x", "z", "y"),
    "yz": ("y", "z", "x"),
}


def load_snapshot(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Snapshot file not found: {path}")
    return pd.read_csv(path, sep=";", engine="python")


def filter_slice(df: pd.DataFrame, plane: str, centre: float, thickness: float) -> pd.DataFrame:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    half = thickness * 0.5
    mask = np.abs(df[axis_c] - centre) <= half
    filtered = df[mask]
    if filtered.empty:
        raise ValueError(
            f"No cells fall within ±{half:g} of {axis_c}={centre:g}. "
            "Adjust --z-value/--thickness."
        )
    return filtered[[axis_a, axis_b, axis_c, "N", "ro", "T", "vx", "vy", "vz", "E"]]


def draw_plot(df: pd.DataFrame, plane: str, field: str, title: str | None, show_cbar: bool) -> plt.Figure:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    x = df[axis_a].to_numpy()
    y = df[axis_b].to_numpy()
    values = df[field].to_numpy()

    contour = None
    if x.size >= 3:
        try:
            triang = mtri.Triangulation(x, y)
            contour = ax.tricontourf(triang, values, levels=40, cmap="plasma")
            ax.tricontour(triang, values, levels=15, colors="k", linewidths=0.2, alpha=0.3)
        except (ValueError, RuntimeError):
            contour = None

    if contour is None:
        scatter = ax.scatter(x, y, c=values, cmap="plasma", s=40, lw=0.0)
        mappable = scatter
    else:
        mappable = contour

    ax.set_xlabel(axis_a)
    ax.set_ylabel(axis_b)
    ax.set_aspect("equal")
    ax.set_title(title or f"{field} on {plane.upper()} plane ({axis_c}≈{df[axis_c].mean():.3f})")

    if show_cbar:
        fig.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, label=field)

    return fig


def main() -> None:
    args = parse_args()
    snapshot = load_snapshot(args.input)
    filtered = filter_slice(snapshot, args.plane, args.z_value, args.thickness)
    fig = draw_plot(filtered, args.plane, args.field, args.title, args.show_cbar)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=200, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()

