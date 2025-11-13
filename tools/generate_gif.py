#!/usr/bin/env python3
"""
Generate an animated GIF from a sequence of MC3D snapshot files.

Relies on the plotting utilities from plot_snapshot.py to keep the look
consistent with static figures.
"""

from __future__ import annotations

import argparse
import gc
import glob
import re
from pathlib import Path
from typing import List, Iterator

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from PIL import Image

from snapshot_loader import load_snapshot

# Try to import imageio for streaming GIF writing (more memory efficient)
try:
    import imageio
    HAS_IMAGEIO = True
except ImportError:
    HAS_IMAGEIO = False

PLANE_AXES = {
    "xy": ("x", "y", "z"),
    "xz": ("x", "z", "y"),
    "yz": ("y", "z", "x"),
}

PLANE_VELOCITY_COMPONENTS = {
    "xy": ("vx", "vy"),
    "xz": ("vx", "vz"),
    "yz": ("vy", "vz"),
}

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
        "--figsize",
        type=float,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=None,
        help="Figure size in inches for each frame (width height). Defaults to 7 5.5.",
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

def filter_slice(
    df: pd.DataFrame, plane: str, centre: float, thickness: float, *, allow_empty: bool = False
) -> pd.DataFrame:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    half = thickness * 0.5
    mask = np.abs(df[axis_c] - centre) <= half
    columns = [axis_a, axis_b, axis_c, "N", "ro", "T", "vx", "vy", "vz", "E"]
    filtered = df.loc[mask, columns]
    if filtered.empty and not allow_empty:
        raise ValueError(
            f"No cells fall within ±{half:g} of {axis_c}={centre:g}. Adjust --z-value/--thickness."
        )
    return filtered


def draw_plot(
    df: pd.DataFrame,
    plane: str,
    field: str,
    title: str | None,
    show_cbar: bool,
    show_velocity: bool,
    quiver_step: int,
    quiver_scale: float | None,
    quiver_color: str,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    axis_a, axis_b, axis_c = PLANE_AXES[plane]
    fig, ax = plt.subplots(figsize=figsize or (7, 5.5))
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

    if show_velocity:
        if quiver_step <= 0:
            raise ValueError("--quiver-step must be a positive integer.")
        vel_a, vel_b = PLANE_VELOCITY_COMPONENTS[plane]
        u = df[vel_a].to_numpy()
        v = df[vel_b].to_numpy()

        step = slice(None, None, quiver_step)
        quiver_kwargs = {
            "angles": "xy",
            "scale_units": "xy",
            "color": quiver_color,
            "pivot": "mid",
            "linewidths": 0.4,
        }
        if quiver_scale is not None:
            quiver_kwargs["scale"] = quiver_scale

        ax.quiver(x[step], y[step], u[step], v[step], **quiver_kwargs)

    return fig


def dataframe_to_image(fig: plt.Figure, dpi: int, use_rgb: bool = True) -> Image.Image:
    """
    Convert matplotlib figure to PIL Image.

    Args:
        fig: Matplotlib figure to convert
        dpi: Resolution for rendering
        use_rgb: If True, convert to RGB (saves ~25% memory). If False, use RGBA.

    Returns:
        PIL Image object
    """
    fig.set_dpi(dpi)
    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    buffer = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8).reshape((height, width, 4))

    if use_rgb:
        # Convert ARGB -> RGB for Pillow (more memory efficient, no alpha channel needed for GIF)
        rgb = np.empty((height, width, 3), dtype=np.uint8)
        rgb[..., 0] = buffer[..., 1]  # R
        rgb[..., 1] = buffer[..., 2]  # G
        rgb[..., 2] = buffer[..., 3]  # B
        image = Image.fromarray(rgb, "RGB")
        # Explicitly delete intermediate arrays to free memory immediately
        del rgb
    else:
        # Convert ARGB -> RGBA for Pillow.
        rgba = np.empty_like(buffer)
        rgba[..., 0] = buffer[..., 1]
        rgba[..., 1] = buffer[..., 2]
        rgba[..., 2] = buffer[..., 3]
        rgba[..., 3] = buffer[..., 0]
        image = Image.fromarray(rgba, "RGBA")
        del rgba

    # Free memory from buffer and close figure
    del buffer
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


def generate_frames(
    files: List[Path],
    args: argparse.Namespace,
    figsize: tuple[float, float] | None,
) -> Iterator[Image.Image]:
    """
    Generator that yields frames one at a time to avoid storing all in memory.
    """
    use_rgb = True

    for idx, path in enumerate(files, start=1):
        snapshot = load_snapshot(path)
        filtered = filter_slice(
            snapshot, args.plane, args.z_value, args.thickness, allow_empty=True
        )

        # Free snapshot data immediately after filtering
        del snapshot

        if filtered.empty:
            print(f"[{idx}/{len(files)}] Skipping {path.name}: slice contains no cells.")
            continue

        fig = draw_plot(
            filtered,
            args.plane,
            args.field,
            title=f"{path.name}",
            show_cbar=args.show_cbar,
            show_velocity=args.show_velocity,
            quiver_step=args.quiver_step,
            quiver_scale=args.quiver_scale,
            quiver_color=args.quiver_color,
            figsize=figsize,
        )

        apply_color_limits(fig, args.vmin, args.vmax)
        frame = dataframe_to_image(fig, dpi=args.dpi, use_rgb=use_rgb)

        # Free filtered data after creating the plot
        del filtered

        print(f"[{idx}/{len(files)}] Generated frame from {path.name}")

        # Periodic garbage collection for large frame counts
        if idx % 10 == 0:
            gc.collect()

        yield frame


def create_gif_streaming(args: argparse.Namespace, frames: Iterator[Image.Image]) -> None:
    """
    Create GIF using imageio with streaming (memory efficient).
    """
    args.out.parent.mkdir(parents=True, exist_ok=True)
    duration = 1.0 / max(args.fps, 0.1)

    # Convert PIL Images to numpy arrays on-the-fly and write directly
    print("Writing GIF with streaming (memory efficient)...")
    frame_count = 0

    with imageio.get_writer(
        args.out,
        mode="I",
        duration=duration,
        loop=0,
    ) as writer:
        for frame in frames:
            # Convert PIL Image to numpy array
            frame_array = np.array(frame)
            writer.append_data(frame_array)
            frame_count += 1
            # Free frame immediately after writing
            del frame
            if frame_count % 10 == 0:
                gc.collect()

    print(f"GIF written to {args.out} ({frame_count} frames, {args.fps} fps)")


def create_gif_legacy(args: argparse.Namespace, frames: List[Image.Image]) -> None:
    """
    Create GIF using PIL (requires all frames in memory).
    """
    if not frames:
        raise RuntimeError("No frames were generated. Check slice parameters or input pattern.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    duration_ms = 1000.0 / max(args.fps, 0.1)

    num_frames = len(frames)
    print(f"Writing GIF with {num_frames} frames (all frames in memory)...")
    frames[0].save(
        args.out,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        disposal=2,
    )

    # Clear frames from memory after writing
    del frames
    gc.collect()

    print(f"GIF written to {args.out} ({num_frames} frames, {args.fps} fps)")


def create_gif(args: argparse.Namespace) -> None:
    files = find_snapshot_files(args.pattern)
    if args.frame_step <= 0:
        raise ValueError("--frame-step must be a positive integer.")
    files = files[:: args.frame_step]
    if args.max_frames is not None:
        files = files[: args.max_frames]

    if not files:
        raise ValueError("No files selected for GIF generation after applying filters.")

    figsize = tuple(args.figsize) if args.figsize else None

    # Use streaming approach if imageio is available (much more memory efficient)
    if HAS_IMAGEIO:
        frames_generator = generate_frames(files, args, figsize)
        create_gif_streaming(args, frames_generator)
    else:
        # Fallback to legacy method (requires all frames in memory)
        print("Warning: imageio not available. Using legacy method (all frames in memory).")
        print("Install imageio for memory-efficient streaming: pip install imageio")
        frames: List[Image.Image] = []
        for frame in generate_frames(files, args, figsize):
            frames.append(frame)
        create_gif_legacy(args, frames)


def main() -> None:
    args = parse_args()
    create_gif(args)


if __name__ == "__main__":
    main()

