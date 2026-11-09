"""Experimental inference from 12 independent three-seed simulation groups."""
import math
import numpy as np
from shock_alignment import GRID,align
from shock_analysis import monotone,crossing
from shock_crossfit_gradient import estimate_gradient
from shock_mean_profile import measure,radius

GROUPS=12
PROFILE_GRID=np.linspace(-11,11,441)
# NIST/SEMATECH Student t table, df=11, two-sided 95%, rounded upward.
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
T95=2.201


def group_point(trials,center_profile=False):
    data,diagnostics=align(trials,grid=PROFILE_GRID if center_profile else GRID,refine_passes=1)
    core=data[:,:,20:-20] if center_profile else data
    return dict(density=data.mean((0,1)).tolist(),
                slopes=[1/estimate_gradient(core,diagnostics,width=w)['gradient_thickness'] for w in (1.,1.5,2.)])


def centered_density(profiles):
    result=[]
    for y in np.atleast_2d(profiles):
        center=crossing(PROFILE_GRID,monotone(y),.5)
        points=GRID+center
        if points[0]<PROFILE_GRID[0] or points[-1]>PROFILE_GRID[-1]:
            raise ValueError('Centered density outside measured support')
        result.append(np.interp(points,PROFILE_GRID,y))
    return np.asarray(result)


def estimate(groups,repeats=2048,seed=20261109,center_profile=False):
    if len(groups)!=GROUPS or repeats<128:
        raise ValueError('Exactly 12 independent groups and at least 128 resamples required')
    profiles=np.asarray([g['density'] for g in groups],float)
    slopes=np.asarray([g['slopes'] for g in groups],float)
    grid=PROFILE_GRID if center_profile else GRID
    if profiles.shape!=(GROUPS,len(grid)) or slopes.shape!=(GROUPS,3) or not np.isfinite(profiles).all() or not np.isfinite(slopes).all() or np.any(slopes<=0):
        raise ValueError('Invalid group profiles or slopes')
    profile=profiles.mean(0);slope=slopes.mean(0)
    point=np.array([measure(profile,grid=grid)[0,0],1/slope[1]])
    group_width=measure(profiles,grid=grid)[:,0]
    width_u=T95*group_width.std(ddof=1)/math.sqrt(GROUPS)+abs(group_width.mean()-point[0])
    slope_u=T95*slopes[:,1].std(ddof=1)/math.sqrt(GROUPS)
    if slope[1]<=slope_u:raise ValueError('Unresolved mean slope')
    gradient_u=max(abs(1/(slope[1]-slope_u)-point[1]),abs(1/(slope[1]+slope_u)-point[1]))
    seed_u=np.array([width_u,gradient_u])
    weights=np.random.default_rng(seed).multinomial(GROUPS,[1/GROUPS]*GROUPS,size=repeats)/GROUPS
    sampled_profiles=weights@profiles;sampled_slopes=weights@slopes
    sampled=np.column_stack([measure(sampled_profiles,grid=grid)[:,0],1/sampled_slopes[:,1]])
    total=np.maximum(seed_u,radius(sampled,point))
    if center_profile:
        profiles=centered_density(profiles)
        sampled_profiles=centered_density(sampled_profiles)
        profile=centered_density(profile)[0]
    density_t=T95*profiles.std(0,ddof=1)/math.sqrt(GROUPS)
    density_boot=float(np.quantile(np.max(abs(sampled_profiles-profile),axis=1),.95))
    density_u=max(float(density_t.max()),density_boot)
    proc=max(abs(1/slope[0]-point[1]),abs(1/slope[2]-point[1]))
    result={k:dict(mean=float(point[i]),statistical_radius=float(total[i]),group_t_radius=float(seed_u[i]),
                   processing_sensitivity=0. if i==0 else float(proc)) for i,k in enumerate(('width_10_90','gradient_thickness'))}
    result.update(status='experimental_coverage_unverified',group_count=GROUPS,center_profile=center_profile,
        density=profile.tolist(),density_simultaneous_radius=density_u,
        density_bootstrap_radius=density_boot,density_max_t_radius=float(density_t.max()),
        group_widths=group_width.tolist(),group_slopes=slopes.tolist())
    return result
