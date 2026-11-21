"""Experimental support-bounded global translation fit; coverage unqualified.

For linear interpolation, Huber loss is convex inside each interpolation cell.
Minimize every cell, then retain the existing separated-minimum ambiguity rule.
This module is deliberately not wired into the frozen production estimator.
"""
import numpy as np
from shock_alignment import PreparedWindow

METHOD = 'support_bounded_cellwise_huber_v1'


def _vector(values, name):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 4 or not np.isfinite(values).all():
        raise ValueError(f'Finite one-dimensional {name} with four nodes required')
    return values


def global_shift(window, tx, ty, output_grid):
    if isinstance(window, PreparedWindow):
        x, y = window.x, window.y
    else:
        x = np.asarray(window['x'] if 'x' in window else window['x_aligned'], dtype=float)
        if 'x' not in window:
            x = x + window.get('front', 0.)
        y = window['density_normalized']
    x, y = _vector(x, 'source grid'), _vector(y, 'source density')
    tx, ty = _vector(tx, 'template grid'), _vector(ty, 'template density')
    output_grid = _vector(output_grid, 'output grid')
    if len(x) != len(y) or len(tx) != len(ty):
        raise ValueError('Grid and density dimensions differ')
    if any(np.any(np.diff(a) <= 0) for a in (x, tx, output_grid)):
        raise ValueError('Increasing source, template and output grids required')
    mask = (ty >= .1) & (ty <= .9)
    u, target = tx[mask], ty[mask]
    if len(u) < 4:
        raise ValueError('Unresolved template')
    lower = max(x[0]-u[0], x[0]-output_grid[0])
    upper = min(x[-1]-u[-1], x[-1]-output_grid[-1])
    if upper-lower <= 2e-4:
        raise ValueError('Insufficient translation support')
    # Collapse roundoff copies of the same cell boundary on regular grids.
    knots = np.unique(np.round((x[:, None]-u).reshape(-1), 12))
    edges = np.r_[lower, knots[(knots > lower) & (knots < upper)], upper]
    left, right = edges[:-1].copy(), edges[1:].copy()
    mid = (left+right)/2
    points = u[None, :]+mid[:, None]
    indices = np.clip(np.searchsorted(x, points, side='right')-1, 0, len(x)-2)
    slopes = (y[indices+1]-y[indices])/(x[indices+1]-x[indices])
    residuals = y[indices]+(points-x[indices])*slopes-target
    for _ in range(40):
        probe = (left+right)/2
        residual = residuals+(probe-mid)[:, None]*slopes
        derivative = np.mean(np.clip(residual, -.02, .02)*slopes, axis=1)
        left = np.where(derivative < 0, probe, left)
        right = np.where(derivative < 0, right, probe)
    candidates = (left+right)/2
    residual = np.abs(residuals+(candidates-mid)[:, None]*slopes)
    costs = np.mean(np.where(residual <= .02, .5*residual**2, .02*(residual-.01)), axis=1)
    best = int(costs.argmin()); front = float(candidates[best]); loss = float(costs[best])
    if min(front-lower, upper-front) <= 1e-4:
        raise ValueError('Boundary minimum on available support')
    # Include neighboring equal-cost cells when a minimum lies on a knot.
    local = (costs <= np.r_[np.inf, costs[:-1]]+1e-14) & (costs <= np.r_[costs[1:], np.inf]+1e-14)
    separated = local & (np.abs(candidates-front) >= .25)
    competing = float(costs[separated].min()) if np.any(separated) else None
    if competing is not None and competing-loss <= max(1e-8, .01*loss):
        raise ValueError('Ambiguous minimum')
    check = np.abs(np.interp(u+front, x, y)-target)
    actual = float(np.mean(np.where(check <= .02, .5*check**2, .02*(check-.01))))
    if abs(actual-loss) > 1e-12:
        raise ValueError('Cellwise objective verification failed')
    return front, dict(method=METHOD, status='experimental_coverage_unverified',
                       front=front, loss=actual, support_bounds=[float(lower), float(upper)],
                       interpolation_cells=len(candidates), competing_minimum_loss=competing)
