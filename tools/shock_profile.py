#!/usr/bin/env python3
"""Compare central upstream profiles on two cube-flow grids.

The Mx=1 crossing is a flow indicator, not a precise kinetic shock surface.
Requires numpy, pandas, matplotlib. Input directories contain config.json,
initial.dat, timed data*.dat files and run.log from configs/hypersonic_cube.json.
"""
import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def crossing(x, y, value=1):
    for i in range(len(x) - 1):
        if y[i] > value >= y[i + 1]:
            return float(x[i] + (value - y[i]) * (x[i+1] - x[i]) / (y[i+1] - y[i]))
    return None


def profiles(directory):
    cfg = json.loads((directory / "config.json").read_text())
    files = []
    for path in directory.glob("data*.dat"):
        match = re.fullmatch(r"data([0-9.eE+-]+)\.dat", path.name)
        if match:
            files.append((float(match[1]), path))
    files.sort()
    if (directory / "final.dat").exists() and (not files or files[-1][0] < cfg["end_time"] - 1e-6):
        files.append((cfg["end_time"], directory / "final.dat"))
    if not files or abs(files[-1][0] - cfg["end_time"]) > 1e-6:
        raise ValueError(f"Incomplete snapshots in {directory}")
    log = (directory / "run.log").read_text()
    if "Computation is over." not in log or "Particles inside body at time" in log:
        raise ValueError(f"Incomplete run or body penetration in {directory}")
    history, late, front_history = [], [], []
    wall = cfg["geometry_cube_x"]
    u0 = cfg["s"] * np.sqrt(2 * cfg["temperature"])
    width = cfg["geometry_cube_width"]
    for t, path in [(0, directory / "initial.dat"), *files]:
        frame = pd.read_csv(path, sep=";")
        if not np.isfinite(frame.to_numpy()).all() or (frame.N < 0).any() or (frame["T"] < 0).any():
            raise ValueError(f"Invalid fields in {path}")
        n = frame.N.sum()
        mean_v = np.array([(frame.N*frame[v]).sum()/n for v in ("vx", "vy", "vz")])
        global_temperature = (2*(frame.N*frame.E).sum()/n - mean_v @ mean_v)/3
        history.append(dict(t=t, N=int(n), ux=float((frame.N*frame.vx).sum()/n),
            T=float((frame.N*frame["T"]).sum()/n), global_T=float(global_temperature)))
        if t < 2 * cfg["end_time"] / 3:
            continue
        slab = frame.loc[(frame.y.abs() < width / 4) & (frame.z.abs() < width / 4)].copy()
        slab["uxN"] = slab.vx * slab.N
        slab["TN"] = slab["T"] * slab.N
        groups = slab.groupby("x", sort=True)
        p = groups[["N", "uxN", "TN"]].sum()
        p["samples"] = groups.size()
        late.append(p)
        upstream = p.loc[p.index < wall]
        mach = upstream.uxN / upstream.N / np.sqrt((5/3) * upstream.TN / upstream.N)
        location = crossing(upstream.index.to_numpy(), mach.to_numpy())
        if location is not None:
            front_history.append(dict(t=t, distance=wall-location))
    pooled = late[0].copy()
    for p in late[1:]:
        pooled = pooled.add(p)
    pooled["density"] = pooled.N / pooled.samples / cfg["np"]
    pooled["ux"] = pooled.uxN / pooled.N.replace(0, np.nan)
    pooled["temperature"] = pooled.TN / pooled.N.replace(0, np.nan)
    pooled["mach"] = pooled.ux / np.sqrt((5/3) * pooled.temperature)
    up = pooled.loc[pooled.index < wall]
    location = crossing(up.index.to_numpy(), up.mach.to_numpy())
    # Maximum density gradient is reported separately: a wall Knudsen layer
    # can dominate it, so it must not be silently labelled the shock position.
    x = up.index.to_numpy()
    gradient = np.diff(up.density.to_numpy()) / np.diff(x)
    i = int(np.argmax(gradient))
    metrics = dict(dx=cfg["lx"]/cfg["ncx"], grid=[cfg[k] for k in ("ncx","ncy","ncz")],
        mach_one_distance=None if location is None else wall-location,
        density_gradient_distance=float(wall-(x[i]+x[i+1])/2),
        initial_count=history[0]["N"], final_count=history[-1]["N"],
        final_ux_over_u0=history[-1]["ux"]/u0, final_temperature=history[-1]["T"],
        final_global_temperature=history[-1]["global_T"])
    if cfg["geometry_type"] == "none":
        metrics.pop("density_gradient_distance")
        metrics.pop("mach_one_distance")
    if front_history:
        metrics["instantaneous_mach_one_distance_p10_p90"] = np.quantile(
            [p["distance"] for p in front_history], [0.1,0.9]).tolist()
    return cfg, pooled, pd.DataFrame(history), pd.DataFrame(front_history), metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results = {name: profiles(getattr(args,name)) for name in ("coarse","fine","control")}
    cfg = results["fine"][0]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"figure.facecolor":"#f7f9fc"})
    fig, axes = plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for name, label, color in [("coarse","80×40×20","#2b7b9b"),("fine","160×40×20","#d56a40")]:
        config,p,history,front,metrics = results[name]
        p.to_csv(args.out/f"{name}_profile.csv")
        history.to_csv(args.out/f"{name}_history.csv",index=False)
        front.to_csv(args.out/f"{name}_front_history.csv",index=False)
        wall=config["geometry_cube_x"]
        p=p.loc[(p.index < wall) & (p.index > wall-0.4)]
        d=wall-p.index.to_numpy()
        for ax,field in zip(axes.flat[:3],["density","ux","mach"]):
            values=p[field].to_numpy()
            if field=="ux": values=values/(config["s"]*np.sqrt(2*config["temperature"]))
            ax.plot(d,values,".-",label=label,color=color)
            ax.set_xlim(0,0.4)
            ax.set_xlabel("Расстояние перед лобовой гранью")
            ax.grid(alpha=0.2)
        if not front.empty:
            axes[1,1].plot(front.t,front.distance,".-",color=color,label=label)
    for ax,title in zip(axes.flat[:3],["Плотность n/n₀","Скорость uₓ/U₀","Индикатор Mₓ = uₓ / √(5T/3)"]):
        ax.set_title(title)
    axes[1,0].axhline(1,color="black",linestyle="--",linewidth=1)
    axes[1,1].set(title="Расстояние до Mₓ=1 во времени",xlabel="Время t",ylabel="Расстояние от грани")
    axes[1,1].grid(alpha=0.2)
    axes[0,0].legend()
    fig.suptitle("ФРОНТ ПЕРЕД КУБОМ · сравнение разрешения вдоль потока\n"
                 "Mₓ=1 — диагностический индикатор, не точная ударная поверхность",fontsize=15)
    fig.savefig(args.out/"shock_profiles.png",dpi=160)
    plt.close(fig)
    fig, axes=plt.subplots(1,3,figsize=(13,4),layout="constrained")
    for name,label in [("control","Без тела"),("coarse","Куб · 80"),("fine","Куб · 160")]:
        config,_,h,_,_=results[name]
        u0=config["s"]*np.sqrt(2*config["temperature"])
        for ax,v in zip(axes,[h.N/h.N.iloc[0],h.ux/u0,h.global_T/config["temperature"]]):
            ax.plot(h.t,v,label=label)
            ax.grid(alpha=0.2); ax.set_xlabel("t")
    for ax,title in zip(axes,["N/Nнач","⟨uₓ⟩/U₀","Глобальная дисперсия / T₀"]): ax.set_title(title)
    axes[0].legend()
    fig.suptitle("Все шесть границ: hyperfree · контроль однородного газа")
    fig.savefig(args.out/"reservoir_check.png",dpi=160)
    plt.close(fig)
    summary={name:r[4] for name,r in results.items()}
    (args.out/"shock_summary.json").write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
