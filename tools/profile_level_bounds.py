"""Outer crossing/width bounds conditional on a simultaneous nodal profile band.

Assume a continuous nondecreasing target. No interpolation, smoothing or alignment
is imposed between nodes. These operations cannot establish coverage of the input.
"""
import numpy as np


def monotone_band(x, lower, upper):
    x,lo,hi=(np.asarray(v,float) for v in (x,lower,upper))
    if x.ndim!=1 or len(x)<2 or lo.shape!=x.shape or hi.shape!=x.shape:
        raise ValueError('Matching one-dimensional bands required')
    if not np.isfinite([x,lo,hi]).all() or np.any(np.diff(x)<=0) or np.any(lo>hi):
        raise ValueError('Finite ordered nodes and valid band required')
    lo=np.maximum.accumulate(lo)
    hi=np.minimum.accumulate(hi[::-1])[::-1]
    if np.any(lo>hi):raise ValueError('No nondecreasing profile satisfies this band')
    return x,lo,hi


def crossing_bounds(x,lower,upper,level):
    """Closed outer bounds for inf{x: f(x)>=level}, allowing arbitrary subgrid shape."""
    if not np.isfinite(level):raise ValueError('Finite level required')
    x,lo,hi=monotone_band(x,lower,upper)
    below=np.flatnonzero(hi<level)
    above=np.flatnonzero(lo>=level)
    if not len(below) or not len(above):
        raise ValueError('Crossing not bracketed for every profile; extend domain or reduce uncertainty')
    return [float(x[below[-1]]),float(x[above[0]])]


def width_bounds(x,lower,upper,levels=(.1,.9)):
    if len(levels)!=2 or not np.isfinite(levels).all() or levels[0]>=levels[1]:
        raise ValueError('Two increasing finite levels required')
    a=crossing_bounds(x,lower,upper,levels[0])
    b=crossing_bounds(x,lower,upper,levels[1])
    return dict(levels=list(levels),crossings=[a,b],
        width=[max(0.,b[0]-a[1]),b[1]-a[0]],
        interpretation='Outer set conditional on coverage and monotonicity; not a new confidence level')
