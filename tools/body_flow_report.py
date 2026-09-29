#!/usr/bin/env python3
"""Render a cube-flow diagnostic and its no-body control (CSV snapshots).

Requires numpy, pandas, matplotlib and pillow. Cell averages are plotted without
spatial interpolation. This is a qualitative diagnostic, not solver validation.
"""

import argparse
import io
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image


def snapshots(directory):
    timed = []
    for path in directory.glob("data*.dat"):
        match = re.fullmatch(r"data([0-9.eE+-]+)\.dat", path.name)
        if match:
            timed.append((float(match[1]), path))
    if not timed or not (directory / "initial.dat").exists():
        raise ValueError(f"Missing initial or timed snapshots in {directory}")
    return [(0.0, directory / "initial.dat"), *sorted(timed)]


def cube_bounds(cfg):
    x = cfg["geometry_cube_x"]
    return np.array([x, -cfg["geometry_cube_width"] / 2,
                     -cfg["geometry_cube_height"] / 2]), np.array([
        x + cfg["geometry_cube_length"], cfg["geometry_cube_width"] / 2,
        cfg["geometry_cube_height"] / 2])


def load_case(directory, cfg, has_body):
    series, slices = [], {"xy": [], "xz": []}
    lo, hi = cube_bounds(cfg)
    reference_coords = None
    dx = cfg["lx"] / cfg["ncx"]
    for t, path in snapshots(directory):
        frame = pd.read_csv(path, sep=";").sort_values(["x", "y", "z"])
        if not np.isfinite(frame.to_numpy()).all():
            raise ValueError(f"Non-finite data in {path}")
        if (frame.N < 0).any() or (frame["T"] < -1e-10).any():
            raise ValueError(f"Negative count/temperature in {path}")
        coords = frame[["x", "y", "z"]].to_numpy()
        if reference_coords is None:
            if len(frame) != cfg["ncx"] * cfg["ncy"] * cfg["ncz"]:
                raise ValueError(f"Grid size does not match config in {path}")
            for axis in "xyz":
                count, length = cfg["nc" + axis], cfg["l" + axis]
                expected = -length / 2 + (np.arange(count) + 0.5) * length / count
                actual = np.sort(frame[axis].unique())
                if len(actual) != count or not np.allclose(actual, expected):
                    raise ValueError(f"Grid coordinates do not match config in {path}")
            reference_coords = coords
        elif not np.array_equal(coords, reference_coords):
            raise ValueError(f"Grid changed in {path}")
        n = frame.N.sum()
        if n <= 0:
            raise ValueError(f"Empty gas domain in {path}")
        inside = ((coords > lo + 1e-9) & (coords < hi - 1e-9)).all(axis=1)
        series.append(dict(time=t, count=float(n),
            ux=float((frame.N * frame.vx).sum() / n),
            temperature=float((frame.N * frame["T"]).sum() / n),
            energy=float((frame.N * frame.E).sum()),
            interior_count=float(frame.loc[inside, "N"].sum()) if has_body else 0))
        for plane, normal in [("xy", "z"), ("xz", "y")]:
            # Two central layers, when grid spacing is isotropic.
            slab = frame.loc[frame[normal].abs() < dx + 1e-9].copy()
            for field in ["vx", "vy", "vz", "T"]:
                slab[field] *= slab.N
            grouped = slab.groupby(list(plane), sort=True)
            summed = grouped[["N", "vx", "vy", "vz", "T"]].sum()
            summed["samples"] = grouped.size()
            slices[plane].append(summed)
    log = (directory / "run.log").read_text() if (directory / "run.log").exists() else ""
    if "Computation is over." not in log:
        raise ValueError(f"Simulation completion not confirmed in {directory}/run.log")
    return pd.DataFrame(series), slices, log.count("Particles inside body at time")


def average(frames, cfg):
    total = frames[0].copy()
    for frame in frames[1:]:
        total = total.add(frame)
    result = total.copy()
    result["density"] = total.N / total.samples / cfg["np"]
    for field in ["vx", "vy", "vz", "T"]:
        result[field] = total[field] / total.N.replace(0, np.nan)
    result["ux"] = result.vx / (cfg["s"] * np.sqrt(2 * cfg["temperature"]))
    result["temperature"] = result["T"] / cfg["temperature"]
    result["nT"] = result.density * result.temperature
    return result


def grid(frame, plane, field):
    table = frame[field].unstack(plane[0])
    return table.columns.to_numpy(), table.index.to_numpy(), table.to_numpy()


def map_panel(ax, frame, plane, field, cfg, body, limits, title, arrows=False):
    x, y, values = grid(frame, plane, field)
    cmap = {"density": "viridis", "temperature": "inferno", "ux": "coolwarm",
            "nT": "magma"}[field]
    image = ax.pcolormesh(x, y, values, shading="nearest", cmap=cmap,
                         vmin=limits[field][0], vmax=limits[field][1], rasterized=True)
    lo, hi = cube_bounds(cfg)
    side = 1 if plane == "xy" else 2
    if arrows:
        _, _, u = grid(frame, plane, "vx")
        _, _, v = grid(frame, plane, "v" + plane[1])
        xx, yy = np.meshgrid(x, y)
        if body:
            mask = (xx > lo[0]) & (xx < hi[0]) & (yy > lo[side]) & (yy < hi[side])
            u, v = np.where(mask, np.nan, u), np.where(mask, np.nan, v)
        ax.quiver(xx[::4, ::5], yy[::4, ::5], u[::4, ::5], v[::4, ::5],
                  color="#172338", scale=35, width=0.0025)
    if body:
        ax.add_patch(Rectangle((lo[0], lo[side]), hi[0] - lo[0], hi[side] - lo[side],
                              facecolor="#283342", edgecolor="white", linewidth=1.2, zorder=5))
    ax.set(xlabel="x", ylabel=plane[1], title=title, aspect="equal")
    ax.set_xlim(-cfg["lx"] / 2, cfg["lx"] / 2)
    extent = cfg["ly"] if plane == "xy" else cfg["lz"]
    ax.set_ylim(-extent / 2, extent / 2)
    return image


def render(body_dir, control_dir, output, cfg):
    if cfg["geometry_type"] != "cube" or cfg["s"] <= 0 or cfg["temperature"] <= 0:
        raise ValueError("This report requires a cube and a positive X inflow")
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titleweight": "bold", "figure.facecolor": "#f7f9fc"})
    body, body_slices, warnings = load_case(body_dir, cfg, True)
    control, control_slices, _ = load_case(control_dir, cfg, False)
    for name, history in [("body", body), ("control", control)]:
        history.to_csv(output / f"{name}_history.csv", index=False)
        if abs(history.time.iloc[-1] - cfg["end_time"]) > 1e-6:
            raise ValueError(f"Incomplete {name} run")
    start = 2 * cfg["end_time"] / 3
    body_late = np.flatnonzero(body.time.to_numpy() >= start)
    control_late = np.flatnonzero(control.time.to_numpy() >= start)
    means = {
        "body": {p: average([body_slices[p][i] for i in body_late], cfg) for p in body_slices},
        "control": {p: average([control_slices[p][i] for i in control_late], cfg) for p in control_slices}}
    limits = {"density": (0, 3), "temperature": (0.5, 3), "ux": (-0.25, 1.25), "nT": (0, 5)}
    titles = {"density": "Плотность n / n₀", "temperature": "Температура T / T₀",
              "ux": "Продольная скорость uₓ / U₀", "nT": "Величина nT / (n₀T₀)"}
    period = f"усреднение t = {body.time.iloc[body_late[0]]:.2f}–{body.time.iloc[-1]:.2f}"
    for plane in ["xy", "xz"]:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
        fig.suptitle(f"ОБТЕКАНИЕ КУБА · срез {plane.upper()}\n{period}; поток слева направо",
                     fontsize=17, fontweight="bold")
        for ax, field in zip(axes.flat, titles):
            im = map_panel(ax, means["body"][plane], plane, field, cfg, True,
                           limits, titles[field], arrows=field == "ux")
            fig.colorbar(im, ax=ax, shrink=0.78, extend="both")
        fig.savefig(output / f"fields_{plane}.png", dpi=160)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    fig.suptitle(f"КОНТРОЛЬ БЕЗ ТЕЛА И КУБ · XY\n{period}; одинаковые шкалы",
                 fontsize=17, fontweight="bold")
    for row, field in enumerate(["density", "ux"]):
        for col, case in enumerate(["control", "body"]):
            ax = axes[row, col]
            im = map_panel(ax, means[case]["xy"], "xy", field, cfg, case == "body", limits,
                           ("Без тела · " if case == "control" else "Куб · ") + titles[field])
            fig.colorbar(im, ax=ax, shrink=0.78, extend="both")
    fig.savefig(output / "control_comparison.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout="constrained")
    u0 = cfg["s"] * np.sqrt(2 * cfg["temperature"])
    for case, name, color in [(body, "Куб", "#cf633d"), (control, "Без тела", "#2879a2")]:
        for ax, values in zip(axes.flat, [case["count"] / case["count"].iloc[0],
                case.ux / u0, case.temperature / cfg["temperature"], case.energy / case["count"]]):
            ax.plot(case.time, values, label=name, color=color, linewidth=2)
            ax.axvspan(start, cfg["end_time"], color="#cbd5e1", alpha=0.2)
            ax.grid(alpha=0.2)
            ax.set_xlabel("Время t")
    for ax, title in zip(axes.flat, ["Число частиц N / Nнач", "Средняя скорость ⟨uₓ⟩ / U₀",
                                   "Средняя температура ⟨T⟩ / T₀", "Кинетическая энергия на частицу"]):
        ax.set_title(title)
    axes[0, 0].legend()
    fig.suptitle("ДИАГНОСТИКА · открытая система: N, импульс и энергия могут меняться",
                 fontsize=15, fontweight="bold")
    fig.savefig(output / "history.png", dpi=160)
    plt.close(fig)

    for plane, field, filename in [("xy", "density", "density_xy.gif"),
                                    ("xz", "ux", "velocity_xz.gif")]:
        frames = []
        for i, t in enumerate(body.time):
            first = max(0, i - 2)
            frame = average(body_slices[plane][first:i + 1], cfg)
            fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
            fig.suptitle("ОБТЕКАНИЕ КУБА · " + titles[field], fontsize=16, fontweight="bold")
            im = map_panel(ax, frame, plane, field, cfg, True, limits,
                f"{plane.upper()} · t = {t:.3f} · среднее по {i-first+1} снимкам", arrows=field == "ux")
            fig.colorbar(im, ax=ax, shrink=0.85, extend="both")
            buffer = io.BytesIO()
            fig.savefig(buffer, format="png", dpi=110)
            plt.close(fig)
            buffer.seek(0)
            with Image.open(buffer) as image:
                frames.append(image.convert("RGB").copy())
            buffer.close()
        frames[0].save(output / filename, save_all=True, append_images=frames[1:],
                       duration=150, loop=0)
        for frame in frames:
            frame.close()

    summary = {"config": cfg, "late_average_start": start,
        "body_interior_warning_count": warnings,
        "max_particles_in_interior_cells": float(body.interior_count.max()),
        "finite_nonnegative_checks": "passed", "cases": {}}
    for case, name in [(body, "body"), (control, "control")]:
        late = case.loc[case.time >= start]
        summary["cases"][name] = {
            "initial_count": float(case["count"].iloc[0]),
            "final_count": float(case["count"].iloc[-1]),
            "late_count_min": float(late["count"].min()),
            "late_count_max": float(late["count"].max()),
            "late_mean_ux_over_u0": float(late.ux.mean() / u0),
            "late_mean_temperature_over_t0": float(late.temperature.mean() / cfg["temperature"])}
    mean = means["body"]["xy"].reset_index()
    lo, hi = cube_bounds(cfg)
    length = hi[0] - lo[0]
    for name, a, b in [("upstream", lo[0] - 0.75 * length, lo[0]),
                       ("wake", hi[0], hi[0] + 1.5 * length)]:
        region = mean.loc[mean.x.between(a, b) & (mean.y.abs() < hi[1])]
        summary[name] = {"density_over_n0": float(region.density.mean()),
                         "ux_over_u0": float(np.average(region.ux, weights=region.N))}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--body", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    render(args.body, args.control, args.out, json.loads(args.config.read_text()))
