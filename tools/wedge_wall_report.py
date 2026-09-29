#!/usr/bin/env python3
"""Compare wedge grids and cube wall temperatures; requires numpy/pandas/matplotlib/pillow.

Each input directory contains its effective config.json, initial.dat, data*.dat,
final.dat and run.log. Density is N/np per full cell; cut cells are excluded
from front detection. The oblique shock relation is an inviscid 2D reference,
not an acceptance criterion for this finite-span rarefied flow.
"""
import argparse
import io
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from PIL import Image


def reference(cfg):
    gamma = 5/3
    mach = cfg["s"]*np.sqrt(2/gamma)
    # CreateWedge uses height=length*sin(alpha), not length*tan(alpha).
    theta = np.arctan(np.sin(cfg["geometry_wedge_alpha"]))
    beta = np.linspace(np.arcsin(1/mach)+1e-7, np.pi/2-1e-7, 10000)
    rhs = 2/np.tan(beta)*(mach**2*np.sin(beta)**2-1)/(mach**2*(gamma+np.cos(2*beta))+2)
    indices = np.flatnonzero(rhs >= np.tan(theta))
    if not len(indices):
        raise ValueError("No attached weak-shock reference for these parameters")
    i = indices[0]
    b = float(np.interp(np.tan(theta), rhs[i-1:i+1], beta[i-1:i+1]))
    mn2 = mach**2*np.sin(b)**2
    density = (gamma+1)*mn2/((gamma-1)*mn2+2)
    pressure = (2*gamma*mn2-(gamma-1))/(gamma+1)
    return dict(mach=mach, theta_degrees=float(np.degrees(theta)),
                beta_degrees=float(np.degrees(b)), density_ratio=density,
                temperature_ratio=pressure/density)


def load(directory):
    cfg = json.loads((directory/"config.json").read_text())
    log = (directory/"run.log").read_text()
    if "Computation is over." not in log or "Particles inside body at time" in log:
        raise ValueError(f"Incomplete run or penetration warning: {directory}")
    files = []
    for f in directory.glob("data*.dat"):
        match = re.fullmatch(r"data([0-9.eE+-]+)\.dat",f.name)
        if match: files.append((float(match[1]),f))
    files.sort()
    if not files: raise ValueError(f"No timed snapshots: {directory}")
    if files[-1][0] < cfg["end_time"]-1e-6:
        files.append((cfg["end_time"],directory/"final.dat"))
    frames, times = [], []
    for t,path in [(0,directory/"initial.dat"),*files]:
        f = pd.read_csv(path,sep=";")
        if not np.isfinite(f.to_numpy()).all() or (f.N<0).any() or (f["T"]<0).any():
            raise ValueError(f"Invalid fields: {path}")
        if len(f) != cfg["ncx"]*cfg["ncy"]*cfg["ncz"]:
            raise ValueError(f"Incorrect grid: {path}")
        f = f.loc[f.z.abs() < cfg["lz"]/cfg["ncz"]].copy()
        f["TN"] = f.N*f["T"]
        f["uN"] = f.N*f.vx
        g = f.groupby(["x","y"])
        p = g[["N","TN","uN"]].sum()
        p["samples"] = g.size()
        frames.append(p); times.append(t)
    late = [p for p,t in zip(frames,times) if t>=2*cfg["end_time"]/3]
    return cfg, times, frames, fields(sum_frames(late),cfg)


def sum_frames(frames):
    total = frames[0].copy()
    for f in frames[1:]: total = total.add(f)
    return total


def fields(p,cfg):
    p=p.copy()
    p["density"] = p.N/p.samples/cfg["np"]
    p["temperature"] = p.TN/p.N.replace(0,np.nan)/cfg["temperature"]
    p["ux"] = p.uN/p.N.replace(0,np.nan)/(cfg["s"]*np.sqrt(2*cfg["temperature"]))
    p["pressure"] = p.density*p.temperature
    return p


def transect(p,cfg,distance):
    target=cfg["geometry_wedge_x"]+distance
    xs=p.index.get_level_values("x").unique().to_numpy()
    # Average a fixed physical strip, reducing noise at the same location on both grids.
    chosen=xs[np.abs(xs-target)<0.025-1e-9]
    strip=p.loc[p.index.get_level_values("x").isin(chosen)]
    raw=strip.groupby("y")[["N","TN","uN","samples"]].sum()
    q=fields(raw,cfg)
    wall=distance*np.sin(cfg["geometry_wedge_alpha"])
    return q.loc[q.index>wall],wall


def front_distance(q,wall,cfg):
    # Exclude cut and adjacent cells. Threshold defines a diagnostic contour,
    # not the shock maximum nor its exact kinetic position.
    q=q.loc[q.index>wall+1.5*cfg["ly"]/cfg["ncy"]]
    y=q.index.to_numpy(); n=q.density.to_numpy()
    crossings=[]
    for i in range(len(y)-1):
        if n[i]>=1.5>n[i+1]:
            crossings.append(y[i]+(1.5-n[i])*(y[i+1]-y[i])/(n[i+1]-n[i]))
    return None if not crossings else float(max(crossings)-wall)


def panel(ax,p,cfg,field,title):
    v=p[field].unstack("x")
    limits={"density":(0,6),"temperature":(0,6),"ux":(0,1.2)}
    im=ax.pcolormesh(v.columns,v.index,v.to_numpy(),shading="nearest",
                     vmin=limits[field][0],vmax=limits[field][1],cmap="viridis")
    x=cfg["geometry_wedge_x"]; length=cfg["geometry_wedge_length"]
    h=length*np.sin(cfg["geometry_wedge_alpha"])
    ax.add_patch(Polygon([[x,0],[x+length,h],[x+length,-h]],facecolor="#263445",edgecolor="white"))
    ref=reference(cfg); b=np.radians(ref["beta_degrees"])
    for sign in [-1,1]: ax.plot([x,x+length],[0,sign*length*np.tan(b)],"--",color="white",lw=1)
    ax.set(xlim=(x-.15,x+length+.15),ylim=(-.35,.35),xlabel="x",ylabel="y",title=title)
    return im


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ["coarse","fine","cold","warm","out"]: parser.add_argument("--"+key,type=Path,required=True)
    a=parser.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    runs={k:load(getattr(a,k)) for k in ["coarse","fine","cold","warm"]}
    common=["lx","ly","lz","np","s","temperature","kn","cu","alpha","end_time"]
    common += [f"boundary_{axis}_{side}" for axis in "xyz" for side in ["neg","pos"]]
    for left,right,geometry in [("coarse","fine","wedge"),("cold","warm","cube")]:
        l,r=runs[left][0],runs[right][0]
        keys=common+[k for k in l if k.startswith(f"geometry_{geometry}_")]
        if geometry=="cube": keys += ["ncx","ncy","ncz"]
        else: keys += ["wall_temperature","ncz"]
        if l["geometry_type"]!=geometry or r["geometry_type"]!=geometry or l["alpha"]!=0:
            raise ValueError(f"Expected an aligned {geometry} pair")
        if any(l.get(k,1 if k=="wall_temperature" else None)!=r.get(k,1 if k=="wall_temperature" else None) for k in keys):
            raise ValueError(f"Non-comparable configurations: {left}, {right}")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"figure.facecolor":"#f7f9fc"})
    cfg,times,frames,p=runs["fine"]
    fig,axes=plt.subplots(1,3,figsize=(15,5),layout="constrained")
    for ax,field,title in zip(axes,["density","temperature","ux"],["Плотность N/np","Температура T/T₀","Скорость uₓ/U₀"]):
        im=panel(ax,p,cfg,field,title);fig.colorbar(im,ax=ax,shrink=.7,extend="max")
    fig.suptitle("КЛИН · среднее за последнюю треть расчёта\nПунктир: невязкий двумерный ориентир; плотность в пересечённых телом ячейках занижена")
    fig.savefig(a.out/"wedge_fields.png",dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout="constrained")
    summary={"reference":reference(cfg),"wedge":{},"cube":{}}
    for name,color in [("coarse","#2879a2"),("fine","#d4683c")]:
        c,_,_,mean=runs[name]; summary["wedge"][name]={}
        label="×".join(str(c[k]) for k in ["ncx","ncy","ncz"])
        for j,d in enumerate([.2,.3,.4]):
            q,wall=transect(mean,c,d)
            q.to_csv(a.out/f"{name}_transect_{d}.csv")
            summary["wedge"][name][str(d)]={"wall_y":wall,"density_1_5_gap":front_distance(q,wall,c)}
            for i,field in enumerate(["density","temperature"]):
                axes[i,j].plot(q.index-wall,q[field],".-",label=label,color=color)
                title="Плотность N/np" if field=="density" else "Температура T/T₀"
                axes[i,j].set(xlim=(0,.22),xlabel="y − yстенки",title=f"x − xвершины = {d}; {title}")
                axes[i,j].grid(alpha=.2)
    axes[0,0].legend();fig.suptitle("КЛИН · поперечные профили на двух сетках; координата 0 — поверхность")
    fig.savefig(a.out/"wedge_profiles.png",dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout="constrained")
    for name in ["cold","warm"]:
        c,_,_,mean=runs[name]
        strip=mean.loc[np.abs(mean.index.get_level_values("y"))<.05]
        q=fields(strip.groupby("x")[["N","TN","uN","samples"]].sum(),c)
        q=q.loc[(q.index<c["geometry_cube_x"]) & (q.index>c["geometry_cube_x"]-.3)]
        q.to_csv(a.out/f"cube_{name}_profile.csv")
        summary["cube"][name]={"wall_temperature":c.get("wall_temperature",1),
            "nearest_cell_density":float(q.density.iloc[-1]),
            "nearest_cell_temperature":float(q.temperature.iloc[-1]),
            "nearest_cell_nT":float(q.pressure.iloc[-1])}
        for ax,field in zip(axes,["density","temperature","pressure"]):
            ax.plot(c["geometry_cube_x"]-q.index,q[field],".-",label=f'Tw={c.get("wall_temperature",1)}')
            title={"density":"Плотность n/n₀","temperature":"Температура T/T₀","pressure":"Величина nT/(n₀T₀)"}[field]
            ax.set(xlabel="Расстояние перед кубом",title=title);ax.grid(alpha=.2)
    axes[0].legend();fig.suptitle("КУБ · влияние температуры диффузной стенки; одинаковая сетка 80×40×20")
    fig.savefig(a.out/"cube_wall_temperature.png",dpi=160);plt.close(fig)
    for field in ["density","temperature"]:
        images=[]
        for i,t in enumerate(times):
            f=fields(sum_frames(frames[max(0,i-2):i+1]),cfg)
            fig,ax=plt.subplots(figsize=(10,6),layout="constrained")
            title="плотность N/np" if field=="density" else "температура T/T₀"
            im=panel(ax,f,cfg,field,f"Клин · {title} · t={t:.3f} · среднее до 3 снимков")
            fig.colorbar(im,ax=ax,shrink=.8,extend="max")
            buf=io.BytesIO();fig.savefig(buf,format="png",dpi=100);plt.close(fig);buf.seek(0)
            with Image.open(buf) as image: images.append(image.convert("RGB").copy())
            buf.close()
        images[0].save(a.out/f"wedge_{field}.gif",save_all=True,append_images=images[1:],duration=160,loop=0)
        for image in images:image.close()
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=="__main__": main()
