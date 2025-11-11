#!/usr/bin/env python3
"""
Generate an animated GIF from a sequence of MC3D snapshot files.

Relies on the plotting utilities from plot_snapshot.py to keep the look
consistent with static figures.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path
from typing import Iterable, List

import numpy as np
from PIL import Image

import matplotlib.pyplot as plt

import plot_snapshot


DATA_ORDER_PATTERN = re.compile(r"([0-9]+(?:\.[0-9]+)?)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate GIF from MC3D snapshots.")
    parser.add_argument(
        "--pattern",
        type=str,
        default="build/data*.dat",
        help="Glob pattern for snapshot files (default: build/data*.dat).",
    )
    parser.add_argument(
        "--plane",
        choices=("xy", "xz", "yz"),
        default="xy",
        help="Plane to visualise.",
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
        help="Field to colour by.",
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
        help="Plot every N-th velocity vector (use >1 to declutter).",
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
        help="Include colour bar on every frame.",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=6.0,
        help="Target frames per second for the GIF animation.",
    )
    parser.add_argument(
        "--frame-step",
        type=int,
        default=1,
        help="Use every N-th file from the sorted list (skip frames to shorten the GIF).",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Limit the number of frames included in the GIF.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=120,
        help="DPI for the rendered frames (controls resolution).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output GIF file path.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Optional lower bound for colour scale (applied to all frames).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Optional upper bound for colour scale (applied to all frames).",
    )
    return parser.parse_args()


def find_snapshot_files(pattern: str) -> List[Path]:
    paths = [Path(p).resolve() for p in glob.glob(pattern)]
    paths = [p for p in paths if p.is_file()]
    if not paths:
        raise FileNotFoundError(f"No files matched pattern: {pattern!r}")

    def sort_key(path: Path) -> tuple[float, str]:
        match = DATA_ORDER_PATTERN.search(path.stem)
        if match:
            try:
                return (float(match.group(1)), path.name)
            except ValueError:
                pass
        return (-1.0, path.name)

    return sorted(paths, key=sort_key)


def dataframe_to_image(fig: plt.Figure, dpi: int) -> Image.Image:
    fig.set_dpi(dpi)
    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    buffer = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8).reshape((height, width, 4))
    # Convert ARGB -> RGBA for Pillow.
    rgba = np.empty_like(buffer)
    rgba[..., 0] = buffer[..., 1]
    rgba[..., 1] = buffer[..., 2]
    rgba[..., 2] = buffer[..., 3]
    rgba[..., 3] = buffer[..., 0]
    image = Image.fromarray(rgba, "RGBA")
    plt.close(fig)
    return image


def apply_color_limits(fig: plt.Figure, vmin: float | None, vmax: float | None) -> None:
    if vmin is None and vmax is None:
        return
    ax = fig.axes[0]
    for collection in ax.collections:
        if hasattr(collection, "set_clim"):
            collection.set_clim(vmin, vmax)
    # Update colorbar if present.
    for mappable in ax.images + ax.collections:
        if hasattr(mappable, "set_clim"):
            mappable.set_clim(vmin, vmax)


def create_gif(args: argparse.Namespace) -> None:
    files = find_snapshot_files(args.pattern)
    if args.frame_step <= 0:
        raise ValueError("--frame-step must be a positive integer.")
    files = files[:: args.frame_step]
    if args.max_frames is not None:
        files = files[: args.max_frames]

    if not files:
        raise ValueError("No files selected for GIF generation after applying filters.")

    frames: List[Image.Image] = []
    for idx, path in enumerate(files, start=1):
        snapshot = plot_snapshot.load_snapshot(path)
        filtered = plot_snapshot.filter_slice(snapshot, args.plane, args.z_value, args.thickness)
        if filtered.empty:
            print(f"[{idx}/{len(files)}] Skipping {path.name}: slice contains no cells.")
            continue

        fig = plot_snapshot.draw_plot(
            filtered,
            args.plane,
            args.field,
            title=f"{path.name}",
            show_cbar=args.show_cbar,
            show_velocity=args.show_velocity,
            quiver_step=args.quiver_step,
            quiver_scale=args.quiver_scale,
            quiver_color=args.quiver_color,
        )

        apply_color_limits(fig, args.vmin, args.vmax)
        frame = dataframe_to_image(fig, dpi=args.dpi)
        frames.append(frame)
        print(f"[{idx}/{len(files)}] Added frame from {path.name}")

    if not frames:
        raise RuntimeError("No frames were generated. Check slice parameters or input pattern.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    duration_ms = 1000.0 / max(args.fps, 0.1)
    frames[0].save(
        args.out,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        disposal=2,
    )
    print(f"GIF written to {args.out} ({len(frames)} frames, {args.fps} fps)")


def main() -> None:
    args = parse_args()
    create_gif(args)


if __name__ == "__main__":
    main()

