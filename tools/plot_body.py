#!/usr/bin/env python3
"""Plot polygon temperatures and integral body loads from an MC3D run."""
from __future__ import annotations

import argparse
import io
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
from PIL import Image

VERTICES = [f"{axis}{i}" for i in (1, 2, 3) for axis in "xyz"]


def read_numeric(path: Path, columns: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=";")
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if frame.empty or not np.isfinite(frame[columns].to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: empty or nonfinite data")
    return frame


def surface_frames(run: Path) -> list[pd.DataFrame]:
    # Prefer the final snapshot when a timed snapshot has the same timestamp.
    paths = sorted(run.glob("*.dat.surface.csv"), key=lambda p: (p.name == "final.dat.surface.csv", p.name))
    frames = {}
    for path in paths:
        frame = read_numeric(path, ["time", "polygon", *VERTICES,
                                    "temperature", "heat_capacity", "absorbed_heat"])
        if frame.time.nunique() != 1 or frame.polygon.duplicated().any():
            raise ValueError(f"{path}: inconsistent time or duplicate polygon IDs")
        if (frame.temperature <= 0).any() or (frame.heat_capacity <= 0).any():
            raise ValueError(f"{path}: nonpositive temperature or heat capacity")
        frames[float(frame.time.iloc[0])] = frame.sort_values("polygon").reset_index(drop=True)
    result = [frames[t] for t in sorted(frames)]
    if result:
        reference = result[0][["polygon", *VERTICES, "heat_capacity"]].to_numpy()
        for frame in result[1:]:
            values = frame[["polygon", *VERTICES, "heat_capacity"]].to_numpy()
            if values.shape != reference.shape or not np.array_equal(values, reference):
                raise ValueError("Surface geometry or heat capacities change between snapshots")
    return result


def surface_figure(frame: pd.DataFrame, norm: Normalize):
    vertices = frame[VERTICES].to_numpy().reshape(-1, 3, 3)
    fig = plt.figure(figsize=(9, 6), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    mesh = Poly3DCollection(vertices, cmap="inferno", norm=norm,
                            edgecolors="#333333", linewidths=0.35)
    mesh.set_array(frame.temperature.to_numpy())
    ax.add_collection3d(mesh)
    points = vertices.reshape(-1, 3)
    lo, hi = points.min(axis=0), points.max(axis=0)
    span = hi - lo
    padding = max(float(span.max()) * 0.04, 1e-9)
    for setter, low, high in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), lo, hi):
        setter(low-padding, high+padding)
    ax.set_box_aspect(np.maximum(span, padding))
    ax.view_init(elev=25, azim=-125)
    ax.set(xlabel="x", ylabel="y", zlabel="z",
           title=f"Surface temperature · t = {frame.time.iloc[0]:.5g}")
    fig.colorbar(mesh, ax=ax, shrink=0.7, pad=0.08, label="Wall temperature (solver units)")
    return fig


def save(fig, path: Path) -> None:
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print(path, flush=True)


def plot_surfaces(frames: list[pd.DataFrame], out: Path, animate: bool) -> None:
    low = min(f.temperature.min() for f in frames)
    high = max(f.temperature.max() for f in frames)
    if high == low:
        low, high = low - max(abs(low)*0.01, 0.01), high + max(abs(high)*0.01, 0.01)
    norm = Normalize(low, high)
    save(surface_figure(frames[-1], norm), out / "surface_temperature.png")
    times = [f.time.iloc[0] for f in frames]
    means = [np.average(f.temperature, weights=f.heat_capacity) for f in frames]
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout="constrained")
    axes[0].fill_between(times, [f.temperature.min() for f in frames],
                         [f.temperature.max() for f in frames], color="tab:orange", alpha=0.2,
                         label="Polygon min–max")
    axes[0].plot(times, means, ".-", color="tab:orange", label="Heat-capacity-weighted mean")
    axes[0].set(ylabel="Wall temperature", title="Surface heating (solver units)")
    axes[0].legend()
    axes[1].plot(times, [f.absorbed_heat.sum() for f in frames], ".-", color="tab:red")
    axes[1].set(xlabel="Time", ylabel="Net absorbed heat")
    for ax in axes:
        ax.grid(alpha=0.25)
    save(fig, out / "surface_history.png")
    if animate:
        images = []
        try:
            for frame in frames:
                fig = surface_figure(frame, norm)
                buffer = io.BytesIO()
                fig.savefig(buffer, format="png", dpi=100)
                plt.close(fig)
                buffer.seek(0)
                with Image.open(buffer) as source:
                    images.append(source.convert("RGB"))
                buffer.close()
            # Playback time is proportional to simulation time; hold final frame.
            gaps = np.diff(times)
            durations = [max(20, round(4000*gap/(times[-1]-times[0]))) for gap in gaps]
            durations.append(800)
            path = out / "surface_temperature.gif"
            images[0].save(path, save_all=True, append_images=images[1:],
                           duration=durations, loop=0)
            print(path, flush=True)
        finally:
            for img in images:
                img.close()


def plot_loads(path: Path, out: Path) -> None:
    components = ["Fx", "Fy", "Fz"]
    frame = read_numeric(path, ["time", "dt", *components, *[f"mean_{c}" for c in components]])
    if (frame.dt <= 0).any() or (np.diff(frame.time) <= 0).any():
        raise ValueError(f"{path}: time must increase and dt must be positive")
    groups = [(components, "Body force (solver units)", "body_forces.png")]
    coefficients = ["Cx", "Cy", "Cz", "Cd", "Cl"]
    if any(c in frame.columns for c in coefficients):
        frame = read_numeric(path, [*frame.columns, *coefficients, *[f"mean_{c}" for c in coefficients],
                                    "reference_area", "dynamic_pressure", "flow_angle"])
        groups.append((coefficients, "Force coefficients", "body_coefficients.png"))
    for columns, title, filename in groups:
        fig, axes = plt.subplots(len(columns), 1, figsize=(10, 2*len(columns)+1),
                                  sharex=True, layout="constrained")
        for ax, column in zip(axes, columns):
            ax.plot(frame.time, frame[column], color="tab:blue", alpha=0.55,
                    linewidth=0.9, label="Step average")
            ax.plot(frame.time, frame[f"mean_{column}"], color="tab:orange",
                    linewidth=1.8, label="Time average from start")
            ax.axhline(0, color="gray", linewidth=0.5)
            ax.set_ylabel(column)
            ax.grid(alpha=0.25)
        axes[0].set_title(title)
        axes[0].legend(loc="best", ncol=2)
        axes[-1].set_xlabel("Time")
        if filename == "body_coefficients.png":
            ref = frame.iloc[-1]
            fig.suptitle(f"Aref = {ref.reference_area:g}, q∞ = {ref.dynamic_pressure:g}, "
                         f"α = {ref.flow_angle:g} rad · Cx/Cy/Cz: global axes; Cd/Cl: flow axes", fontsize=10)
        save(fig, out / filename)


def render_body(run: Path, out: Path, animate: bool = False) -> bool:
    frames = surface_frames(run)
    loads = run / "body_loads.csv"
    if not frames and not loads.exists():
        return False
    out.mkdir(parents=True, exist_ok=True)
    if frames:
        plot_surfaces(frames, out, animate)
    if loads.exists():
        plot_loads(loads, out)
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--out", type=Path, help="Default: RUN/plots")
    parser.add_argument("--gif", action="store_true", help="Animate surface temperature with a fixed color scale")
    args = parser.parse_args()
    try:
        if not render_body(args.run, args.out or args.run / "plots", args.gif):
            parser.error("No *.dat.surface.csv or body_loads.csv found")
    except (ValueError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
