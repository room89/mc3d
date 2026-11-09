"""Experimental mean-profile estimator; does not replace stage-9b acceptance.

Fixed before reanalysis: aligned grid [-10,10] at 0.05 lambda; local cubic
full width 1 lambda, sensitivity at 0.5 lambda; circular blocks 8/16 windows;
hierarchical seed/block percentile envelope plus a conservative seed-t envelope.
Only whole profiles are resampled, retaining spatial dependence.
Requires numpy. No independent-time-sample or calibrated-coverage claim.
"""
import math
import numpy as np
from shock_analysis import monotone

GRID=np.linspace(-10,10,401)


def measure(profiles, grid=GRID, fit_width=1.):
    y=np.atleast_2d(np.asarray(profiles,dtype=float))
    x=np.asarray(grid,dtype=float);dx=x[1]-x[0]
    if y.shape[1]!=len(x) or not np.isfinite(y).all() or dx<=0 or not np.allclose(np.diff(x),dx):
        raise ValueError('Finite profiles on a uniform increasing grid required')
    radius=int(round(fit_width/(2*dx)))
    if radius<2 or 2*radius+1>len(x):raise ValueError('Unresolved fit width')
    t=np.arange(-radius,radius+1)*dx/fit_width
    design=np.column_stack([t**k for k in range(4)])
    weights=np.linalg.pinv(design)[1]/fit_width
    result=[]
    for row in y:
        m=monotone(row)
        if not m[0]<.1<.9<m[-1]:raise ValueError('Missing plateau crossings')
        width=np.interp(.9,m,x)-np.interp(.1,m,x)
        slope=np.convolve(row,weights[::-1],mode='valid').max()
        if slope<=0:raise ValueError('Unresolved gradient')
        result.append([width,1/slope])
    return np.asarray(result)


def aligned(trials,grid=GRID):
    if len(trials)<3 or len({len(r) for r in trials})!=1:raise ValueError('Equal-length runs, at least three seeds required')
    result=[]
    for run in trials:
        rows=[]
        for window in run:
            x=np.asarray(window['x_aligned']);y=np.asarray(window['density_normalized'])
            if len(x)!=len(y) or not np.isfinite(x).all() or not np.all(np.diff(x)>0) or grid[0]<x[0] or grid[-1]>x[-1]:
                raise ValueError('Comparison grid not covered by valid coordinates')
            rows.append(np.interp(grid,x,y))
        result.append(rows)
    return np.asarray(result)


def bootstrap_profiles(data,block,repeats,rng):
    seeds,count,points=data.shape
    if count%block:raise ValueError('Block must divide the window count')
    weights=np.zeros((repeats,seeds*count))
    for row in weights:
        for seed in rng.integers(0,seeds,size=seeds):
            starts=rng.integers(0,count,size=count//block)
            indices=(starts[:,None]+np.arange(block))%count
            row+=np.bincount(seed*count+indices.ravel(),minlength=seeds*count)
    weights/=seeds*count
    return weights@data.reshape(seeds*count,points)


def radius(samples,point):
    lo,hi=np.quantile(samples,[.025,.975],axis=0)
    return np.maximum(np.abs(lo-point),np.abs(hi-point))


def estimate(data,repeats=512,seed=20261107):
    data=np.asarray(data,dtype=float)
    if data.ndim!=3 or data.shape[0]!=3 or data.shape[1]<32 or data.shape[1]%16 or data.shape[2]!=len(GRID) or not np.isfinite(data).all():
        raise ValueError('Expected three finite runs with a multiple of 16 windows')
    if repeats<32:raise ValueError('Need at least 32 bootstrap replicates')
    point_profile=data.mean(axis=(0,1));point=measure(point_profile)[0]
    alternative=measure(point_profile,fit_width=.5)[0]
    seed_profiles=data.mean(axis=1);seed_metrics=measure(seed_profiles)
    # Offset accounts for noncommutation of averaging and nonlinear measurement.
    seed_u=4.303*seed_metrics.std(axis=0,ddof=1)/math.sqrt(3)+abs(seed_metrics.mean(axis=0)-point)
    density_seed_u=4.303*seed_profiles.std(axis=0,ddof=1)/math.sqrt(3)
    rng=np.random.default_rng(seed);blocks={};total=seed_u.copy();density_total=density_seed_u.copy()
    for block in (8,16):
        sampled=bootstrap_profiles(data,block,repeats,rng);m=measure(sampled)
        u=radius(m,point);density_u=radius(sampled,point_profile)
        simultaneous=float(np.quantile(np.max(abs(sampled-point_profile),axis=1),.95))
        total=np.maximum(total,u);density_total=np.maximum(density_total,density_u)
        blocks[str(block)]=dict(metric_radius=u.tolist(),density_max_pointwise_radius=float(density_u.max()),density_simultaneous_radius=simultaneous)
    keys=('width_10_90','gradient_thickness')
    result={k:dict(mean=float(point[i]),statistical_radius=float(total[i]),seed_t_radius=float(seed_u[i]),processing_sensitivity=float(abs(point[i]-alternative[i])) if i else 0.) for i,k in enumerate(keys)}
    result.update(blocks=blocks,seed_metrics=seed_metrics.tolist(),density=point_profile.tolist(),density_pointwise_radius=density_total.tolist(),
                  density_simultaneous_radius=max(float(density_seed_u.max()),*(b['density_simultaneous_radius'] for b in blocks.values())))
    return result
