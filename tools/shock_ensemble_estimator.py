"""Experimental inference from 12 independent three-seed simulation groups."""
import math
import numpy as np
from shock_alignment import GRID,align,initial
from shock_analysis import monotone,crossing
from shock_crossfit_gradient import estimate_gradient
from shock_mean_profile import measure,radius

GROUPS=12
PROFILE_GRID=np.linspace(-11,11,441)
# NIST/SEMATECH Student t table, df=11, two-sided 95%, rounded upward.
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
T95=2.201


def cubic_profile(x,y,points):
    x=np.asarray(x);y=np.asarray(y);points=np.asarray(points)
    if len(x)<4 or len(y)!=len(x) or not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(points).all() or np.any(np.diff(x)<=0):
        raise ValueError('Finite increasing source grid required')
    if np.any(points<x[0]) or np.any(points>x[-1]):raise ValueError('Cubic interpolation outside measured support')
    indices=np.clip(np.searchsorted(x,points)-1,1,len(x)-3)[:,None]+np.arange(-1,3)
    nodes=x[indices];result=np.zeros(len(points))
    for i in range(4):
        weight=np.ones(len(points))
        for j in range(4):
            if j!=i:weight*=(points-nodes[:,j])/(nodes[:,i]-nodes[:,j])
        result+=weight*y[indices[:,i]]
    return result


def group_point(trials,center_profile=False,interpolation_sensitivity=False):
    data,diagnostics=align(trials,grid=PROFILE_GRID if center_profile else GRID,refine_passes=1)
    core=data[:,:,20:-20] if center_profile else data
    result=dict(density=data.mean((0,1)).tolist(),
                slopes=[1/estimate_gradient(core,diagnostics,width=w)['gradient_thickness'] for w in (1.,1.5,2.)])
    if interpolation_sensitivity:
        grid=PROFILE_GRID if center_profile else GRID;alternative=[]
        for windows,info in zip(trials,diagnostics):
            rows=[]
            for window,fit in zip(windows,info['windows']):
                x,y,_=initial(window)
                rows.append(cubic_profile(x,y,grid+fit['front']))
            alternative.append(rows)
        alternative=np.asarray(alternative);core=alternative[:,:,20:-20] if center_profile else alternative
        result['density_cubic']=alternative.mean((0,1)).tolist()
        result['slope_cubic']=1/estimate_gradient(core,diagnostics,width=1.5)['gradient_thickness']
    return result


def centered_density(profiles,interpolation='linear'):
    if interpolation not in ('linear','cubic'):raise ValueError('Unknown centering interpolation')
    result=[]
    for y in np.atleast_2d(profiles):
        monotonic=monotone(y)
        center=crossing(PROFILE_GRID,monotonic,.5)
        if interpolation=='cubic':
            i=int(np.searchsorted(monotonic,.5,side='right'))-1
            lo,hi=PROFILE_GRID[i],PROFILE_GRID[i+1]
            for _ in range(32):
                mid=(lo+hi)/2
                if cubic_profile(PROFILE_GRID,monotonic,np.array([mid]))[0]<.5:lo=mid
                else:hi=mid
            center=(lo+hi)/2
        points=GRID+center
        if points[0]<PROFILE_GRID[0] or points[-1]>PROFILE_GRID[-1]:
            raise ValueError('Centered density outside measured support')
        result.append(np.interp(points,PROFILE_GRID,y) if interpolation=='linear' else cubic_profile(PROFILE_GRID,y,points))
    return np.asarray(result)


def width_density_interval(profile,uncertainty,grid=GRID):
    """Invert a simultaneous density band for the 10--90% width.

    If the band contains the monotone true profile, these crossing bounds
    contain its width. This does not establish coverage of the input band.
    """
    if not np.isfinite(uncertainty) or not 0<=uncertainty<.1:
        raise ValueError('Unresolved density band for width crossings')
    y=monotone(profile)
    lo=crossing(grid,y,.9-uncertainty)-crossing(grid,y,.1+uncertainty)
    hi=crossing(grid,y,.9+uncertainty)-crossing(grid,y,.1-uncertainty)
    return [lo,hi]


def estimate(groups,repeats=2048,seed=20261109,center_profile=False,width_band=False,interpolation_sensitivity=False):
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
    density_statistical=density_u;density_processing=0.;gradient_interpolation=0.;width_processing=0.
    if interpolation_sensitivity:
        alternatives=np.asarray([g['density_cubic'] for g in groups],float)
        alternative_slopes=np.asarray([g['slope_cubic'] for g in groups],float)
        if alternatives.shape!=(GROUPS,len(grid)) or alternative_slopes.shape!=(GROUPS,) or not np.isfinite(alternatives).all() or not np.isfinite(alternative_slopes).all() or np.any(alternative_slopes<=0):
            raise ValueError('Invalid interpolation alternatives')
        alternative=alternatives.mean(0)
        width_processing=abs(measure(alternative,grid=grid)[0,0]-point[0])
        if center_profile:alternative=centered_density(alternative,interpolation='cubic')[0]
        density_processing=float(np.max(abs(alternative-profile)))
        gradient_interpolation=abs(1/alternative_slopes.mean()-point[1])
        density_u+=density_processing
    width_interval=None
    if width_band:
        comparison_grid=GRID if center_profile else grid
        statistical_interval=width_density_interval(profile,density_statistical,comparison_grid)
        total[0]=max(total[0],*(abs(v-point[0]) for v in statistical_interval))
        width_interval=width_density_interval(profile,density_u,comparison_grid)
        width_processing=max(width_processing,0.,max(abs(v-point[0]) for v in width_interval)-total[0])
    proc=max(abs(1/slope[0]-point[1]),abs(1/slope[2]-point[1]))
    result={k:dict(mean=float(point[i]),statistical_radius=float(total[i]),group_t_radius=float(seed_u[i]),
                   processing_sensitivity=float(width_processing) if i==0 else float(proc+gradient_interpolation)) for i,k in enumerate(('width_10_90','gradient_thickness'))}
    result.update(status='experimental_coverage_unverified',group_count=GROUPS,center_profile=center_profile,
        density=profile.tolist(),density_simultaneous_radius=density_u,
        density_bootstrap_radius=density_boot,density_max_t_radius=float(density_t.max()),
        group_widths=group_width.tolist(),group_slopes=slopes.tolist())
    if width_band:result['width_density_interval']=width_interval
    if interpolation_sensitivity:result.update(density_statistical_radius=density_statistical,
        density_processing_sensitivity=density_processing,gradient_interpolation_sensitivity=float(gradient_interpolation))
    return result
