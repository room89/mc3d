#!/usr/bin/env python3
"""
Real-time viewer for MC3D simulation output files.

The solver produces a series of files data<t>.dat containing the same schema as
data.dat.  This script watches a directory for the latest snapshot and updates a
scatter plot every few seconds.
"""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd

DATA_PATTERN = re.compile(r"data([0-9eE\.\+\-]+)\.dat")

PLANE_AXES = {
    "xy": ("x", "y", "z"),
    "xz": ("x", "z", "y"),
    "yz": ("y", "z", "x"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Monitor MC3D simulation results.")
    parser.add_argument(
        "--directory",
        type=Path,
        default=Path("build"),
        help="Directory containing data files (defaults to build/).",
    )
    parser.add_argument(
        "--plane",
        choices=("xy", "xz", "yz"),
        default="xy",
        help="Plane to visualise.",
    )
    parser.add_argument(
        "--coord",
        type=float,
        default=0.0,
        help="Coordinate along the excluded axis (centre of the slice).",
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=0.05,
        help="Thickness of the slice.",
    )
    parser.add_argument(
        "--field",
        choices=("ro", "T", "vx", "vy", "vz", "E"),
        default="ro",
        help="Field to colour by.",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=5.0,
        help="Refresh interval in seconds.",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Render only the latest snapshot and exit.",
    )
    return parser.parse_args()


def list_data_files(directory: Path) -> Iterable[Path]:
    if not directory.exists():
        return []
    files = sorted(
        (p for p in directory.iterdir() if DATA_PATTERN.match(p.name)),
        key=lambda p: p.stat().st_mtime,
    )
    return files


def load_latest(directory: Path) -> pd.DataFrame | None:
    files = list_data_files(directory)
    if not files:
        return None
    try:
        return pd.read_csv(files[-1], sep=";", engine="python")
    except Exception as exc:  # pragma: no cover - best-effort monitoring
        print(f"Failed to read {files[-1]}: {exc}")
        return None


def filter_slice(df: pd.DataFrame, plane: str, centre: float, thickness: float) -> pd.DataFrame:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    half = thickness * 0.5
    mask = np.abs(df[axis_c] - centre) <= half
    return df.loc[mask, [axis_a, axis_b, axis_c, "N", "ro", "T", "vx", "vy", "vz", "E"]]


def update_plot(fig: plt.Figure, ax: plt.Axes, df: pd.DataFrame, plane: str, field: str) -> None:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    ax.clear()
    if df.empty:
        ax.set_title("No cells in slice yet")
        return

    x = df[axis_a].to_numpy()
    y = df[axis_b].to_numpy()
    values = df[field].to_numpy()

    mappable = None
    if x.size >= 3:
        try:
            triang = mtri.Triangulation(x, y)
            contour = ax.tricontourf(triang, values, levels=40, cmap="plasma")
            ax.tricontour(triang, values, levels=12, colors="k", linewidths=0.2, alpha=0.25)
            mappable = contour
        except (ValueError, RuntimeError):
            mappable = None

    if mappable is None:
        scatter = ax.scatter(x, y, c=values, cmap="plasma", s=35, lw=0.0)
        mappable = scatter

    ax.set_xlabel(axis_a)
    ax.set_ylabel(axis_b)
    ax.set_aspect("equal")
    mean_coord = df[axis_c].mean()
    ax.set_title(f"{field} on {plane.upper()} plane ({axis_c}≈{mean_coord:.3f})")

    if hasattr(update_plot, "_colorbar") and update_plot._colorbar is not None:
        update_plot._colorbar.remove()
    update_plot._colorbar = fig.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, label=field)


update_plot._colorbar = None


def main() -> None:
    args = parse_args()
    plt.ion()
    fig, ax = plt.subplots(figsize=(7, 5.5))

    while True:
        df = load_latest(args.directory)
        if df is not None:
            sliced = filter_slice(df, args.plane, args.coord, args.thickness)
            update_plot(fig, ax, sliced, args.plane, args.field)
            fig.canvas.draw_idle()
            fig.canvas.flush_events()
        else:
            ax.set_title("Waiting for MC3D output ...")
            fig.canvas.draw_idle()
            fig.canvas.flush_events()

        if args.once:
            break

        time.sleep(max(0.5, args.interval))

    if not args.once:
        print("Monitoring stopped. Close the window to exit.")
        plt.ioff()
        plt.show()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nStopped by user.")

