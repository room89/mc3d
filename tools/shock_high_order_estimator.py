"""Experimental cubic-profile/quintic-gradient candidate; coverage unqualified.

The frozen earlier estimator is unchanged. Twelve independent groups each pool
six three-seed subgroups. Peak selection always excludes the evaluated subgroup.
"""
import math
import numpy as np
from shock_alignment import GRID
from shock_crossfit_gradient import peak_position
from shock_ensemble_estimator import GROUPS, T95, PROFILE_GRID, centered_density, width_density_interval
from shock_mean_profile import measure, radius

METHOD = 'cubic_quintic_crossfit_v1'
WIDTHS = (1.5, 2., 2.5)


def derivative5(profile, width):
    y = np.asarray(profile,float); dx = GRID[1]-GRID[0]
    if y.shape != GRID.shape or not np.isfinite(y).all() or not math.isfinite(width) or width<=0:
        raise ValueError('Finite profile and positive width required')
    half = int(round(width/(2*dx)))
    if half<3 or 2*half+1>len(GRID): raise ValueError('Unresolved quintic derivative')
    t = np.arange(-half,half+1)*dx/width
    weights = np.linalg.pinv(np.column_stack([t**k for k in range(6)]))[1]/width
    return GRID[half:-half], np.convolve(y,weights[::-1],mode='valid')


def pooled_point(points):
    if len(points)!=6: raise ValueError('Six independent subgroups required')
    primary = np.asarray([p['density_cubic'] for p in points],float)
    alternative = np.asarray([p['density'] for p in points],float)
    if primary.shape!=(6,len(PROFILE_GRID)) or alternative.shape!=primary.shape or not np.isfinite(primary).all() or not np.isfinite(alternative).all():
        raise ValueError('Finite extended subgroup profiles required')
    centered = centered_density(primary,interpolation='cubic')
    linear = centered_density(alternative)
    def selected(profiles,width):
        values=[]
        for heldout in range(6):
            training=np.delete(profiles,heldout,axis=0).mean(0)
            position=peak_position(GRID,training,width=1.5)
            x,d=derivative5(profiles[heldout],width)
            if not x[0]<=position<=x[-1]: raise ValueError('Peak outside evaluation support')
            value=float(np.interp(position,x,d))
            if value<=0: raise ValueError('Nonpositive heldout slope')
            values.append(value)
        return float(np.mean(values))
    return dict(method=METHOD,density=primary.mean(0).tolist(),density_linear=alternative.mean(0).tolist(),
                slopes=[selected(centered,w) for w in WIDTHS],slope_linear=selected(linear,2.))


def estimate(groups,repeats=1024,seed=20261110):
    if len(groups)!=GROUPS or repeats<128 or any(g.get('method')!=METHOD for g in groups):
        raise ValueError('Twelve groups with matching high-order protocol required')
    profiles=np.asarray([g['density'] for g in groups],float)
    alternatives=np.asarray([g['density_linear'] for g in groups],float)
    slopes=np.asarray([g['slopes'] for g in groups],float)
    alternative_slopes=np.asarray([g['slope_linear'] for g in groups],float)
    if (profiles.shape!=(GROUPS,len(PROFILE_GRID)) or alternatives.shape!=profiles.shape or
            slopes.shape!=(GROUPS,3) or alternative_slopes.shape!=(GROUPS,) or
            any(not np.isfinite(a).all() for a in (profiles,alternatives,slopes,alternative_slopes)) or
            np.any(slopes<=0) or np.any(alternative_slopes<=0)):
        raise ValueError('Invalid high-order group statistics')
    profile=profiles.mean(0); slope=slopes.mean(0)
    point=np.array([measure(profile,grid=PROFILE_GRID)[0,0],1/slope[1]])
    widths=measure(profiles,grid=PROFILE_GRID)[:,0]
    width_u=T95*widths.std(ddof=1)/math.sqrt(GROUPS)+abs(widths.mean()-point[0])
    slope_u=T95*slopes[:,1].std(ddof=1)/math.sqrt(GROUPS)
    if slope[1]<=slope_u: raise ValueError('Unresolved mean slope')
    gradient_u=max(abs(1/(slope[1]-slope_u)-point[1]),abs(1/(slope[1]+slope_u)-point[1]))
    group_u=np.array([width_u,gradient_u])
    weights=np.random.default_rng(seed).multinomial(GROUPS,[1/GROUPS]*GROUPS,size=repeats)/GROUPS
    sampled_profiles=weights@profiles; sampled_slopes=weights@slopes
    sampled=np.column_stack([measure(sampled_profiles,grid=PROFILE_GRID)[:,0],1/sampled_slopes[:,1]])
    total=np.maximum(group_u,radius(sampled,point))
    centered=centered_density(profiles,interpolation='cubic')
    sampled_centered=centered_density(sampled_profiles,interpolation='cubic')
    profile=centered_density(profile,interpolation='cubic')[0]
    density_t=T95*centered.std(0,ddof=1)/math.sqrt(GROUPS)
    density_boot=float(np.quantile(np.max(abs(sampled_centered-profile),axis=1),.95))
    density_statistical=max(float(density_t.max()),density_boot)
    alternative=alternatives.mean(0)
    width_processing=abs(measure(alternative,grid=PROFILE_GRID)[0,0]-point[0])
    alternate_centered=centered_density(alternative)[0]
    density_processing=float(np.max(abs(alternate_centered-profile)))
    density_u=density_statistical+density_processing
    gradient_interpolation=abs(1/alternative_slopes.mean()-point[1])
    statistical_interval=width_density_interval(profile,density_statistical,GRID)
    total[0]=max(total[0],*(abs(v-point[0]) for v in statistical_interval))
    interval=width_density_interval(profile,density_u,GRID)
    width_processing=max(width_processing,0.,max(abs(v-point[0]) for v in interval)-total[0])
    processing=max(abs(1/slope[0]-point[1]),abs(1/slope[2]-point[1]))
    result={k:dict(mean=float(point[i]),statistical_radius=float(total[i]),group_t_radius=float(group_u[i]),
                   processing_sensitivity=float(width_processing if i==0 else processing+gradient_interpolation))
            for i,k in enumerate(('width_10_90','gradient_thickness'))}
    result.update(status='experimental_coverage_unverified',method=METHOD,group_count=GROUPS,center_profile=True,
                  density=profile.tolist(),density_simultaneous_radius=density_u,
                  density_statistical_radius=density_statistical,density_processing_sensitivity=density_processing,
                  density_bootstrap_radius=density_boot,density_max_t_radius=float(density_t.max()),
                  gradient_interpolation_sensitivity=float(gradient_interpolation),
                  gradient_width_sensitivity=float(processing),width_density_interval=interval,
                  group_widths=widths.tolist(),group_slopes=slopes.tolist())
    return result
