"""Unqualified interval envelope for two orders of averaging and alignment.

The windowwise high-order estimator remains the reported point estimate. The
second order averages each seed on its physical grid before fitting translation.
The envelope contains both estimators' intervals; this algebraic property does
not establish coverage when both estimators are biased. Independent calibration
and physical precision checks remain mandatory.
"""
import copy
import math
import numpy as np
import shock_high_order_estimator as high
from shock_alignment import initial, GRID
from shock_ensemble_estimator import group_point, width_density_interval

METHOD = 'alignment_order_envelope_v1'


def average_seed_windows(trials):
    if len(trials) != 3 or len({len(run) for run in trials}) != 1 or len(trials[0]) not in (32, 64):
        raise ValueError('Three equal seeds with 32 or 64 windows required')
    result = []
    for run in trials:
        raw = [initial(w) for w in run]
        x = raw[0][0]
        if any(len(v[0]) != len(x) or not np.allclose(v[0], x, rtol=1e-14, atol=1e-12) for v in raw):
            raise ValueError('Fixed physical grid required before averaging')
        result.append([dict(x=x.tolist(), density_normalized=np.mean([v[1] for v in raw], axis=0).tolist())])
    return result


def subgroup_point(trials):
    averaged = average_seed_windows(trials)
    options = dict(center_profile=True, interpolation_sensitivity=True)
    return dict(windows=len(trials[0]), windowwise=group_point(trials, **options),
                meanwise=group_point(averaged, **options))


def pooled_point(points):
    if len(points) != 6 or len({p['windows'] for p in points}) != 1 or points[0]['windows'] not in (32, 64):
        raise ValueError('Six matching 32/64-window subgroups required')
    return dict(method=METHOD, windows=points[0]['windows'],
                windowwise=high.pooled_point([p['windowwise'] for p in points]),
                meanwise=high.pooled_point([p['meanwise'] for p in points]))


def envelope(primary, alternative):
    """Contain both input intervals without assuming either has true coverage."""
    for result in (primary, alternative):
        density = np.asarray(result['density'], float)
        if density.shape != GRID.shape or not np.isfinite(density).all():
            raise ValueError('Finite centered density required')
        fields = [result['density_statistical_radius'], result['density_processing_sensitivity']]
        for name in ('width_10_90', 'gradient_thickness'):
            metric = result[name]
            if not math.isfinite(metric['mean']) or metric['mean'] <= 0:
                raise ValueError('Positive finite metric required')
            fields.extend([metric['statistical_radius'], metric['processing_sensitivity']])
        if any(not math.isfinite(v) or v < 0 for v in fields):
            raise ValueError('Finite nonnegative radii required')
    result = dict(density=list(primary['density']), group_count=12, center_profile=True)
    differences = {}
    for name in ('width_10_90', 'gradient_thickness'):
        a, b = primary[name], alternative[name]
        difference = abs(b['mean'] - a['mean'])
        result[name] = dict(mean=a['mean'])
        result[name]['statistical_radius'] = max(a['statistical_radius'], b['statistical_radius'])
        result[name]['processing_sensitivity'] = max(a['processing_sensitivity'], difference + b['processing_sensitivity'])
        differences[name] = difference
    difference = float(np.max(abs(np.asarray(alternative['density']) - primary['density'])))
    statistical = max(primary['density_statistical_radius'], alternative['density_statistical_radius'])
    processing = max(primary['density_processing_sensitivity'], difference + alternative['density_processing_sensitivity'])
    total = statistical + processing
    interval = width_density_interval(result['density'], total, GRID)
    width = result['width_10_90']
    width['processing_sensitivity'] = max(width['processing_sensitivity'],
        max(abs(v - width['mean']) for v in interval) - width['statistical_radius'])
    result.update(method=METHOD, status='experimental_coverage_unverified',
        density_statistical_radius=statistical, density_processing_sensitivity=processing,
        density_simultaneous_radius=total, width_density_interval=interval,
        alignment_order_density_difference=difference, alignment_order_metric_differences=differences,
        primary_estimate=copy.deepcopy(primary), alternative_estimate=copy.deepcopy(alternative),
        limitation='Interval containment is algebraic, not a coverage or physical acceptance claim.')
    return result


def estimate(groups, repeats=1024, seed=20261111):
    if len(groups) != 12 or any(g.get('method') != METHOD for g in groups):
        raise ValueError('Twelve matching order-envelope groups required')
    if len({g['windows'] for g in groups}) != 1 or groups[0]['windows'] not in (32, 64):
        raise ValueError('Matching window counts required')
    primary = high.estimate([g['windowwise'] for g in groups], repeats=repeats, seed=seed)
    alternative = high.estimate([g['meanwise'] for g in groups], repeats=repeats, seed=seed)
    result = envelope(primary, alternative)
    result['windows'] = groups[0]['windows']
    return result
