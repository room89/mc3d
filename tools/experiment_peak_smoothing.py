"""Paired factorial diagnostic of peak selection with exact front alignment."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from decompose_gradient_bias import decompose
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_analysis import monotone
from shock_crossfit_gradient import derivative, peak_position
from validate_shock_estimator import truth

METHODS = ('raw', 'smooth', 'monotone', 'smooth_monotone')


def training_profile(profile, method):
    y = np.asarray(profile, float)
    if y.shape != GRID.shape or not np.isfinite(y).all() or method not in METHODS:
        raise ValueError('Finite training profile and known method required')
    if method in ('smooth', 'smooth_monotone'):
        offsets = np.arange(-15, 16) * .05
        kernel = np.exp(-.5 * (offsets / .25)**2)
        kernel /= kernel.sum()
        y = np.convolve(np.pad(y, 15, mode='edge'), kernel, 'valid')
    if method in ('monotone', 'smooth_monotone'):
        y = np.asarray(monotone(y))
    return y


def evaluate(data, kind, exact):
    data = np.asarray(data, float)
    if data.ndim != 3 or data.shape[0] != 3 or data.shape[1] < 1 or data.shape[2] != len(GRID) or not np.isfinite(data).all():
        raise ValueError('Three finite aligned seed arrays required')
    true_peak = math.log(math.log(2)) if kind == 'gompertz' else 0.
    slopes = {m: [] for m in (*METHODS, 'known_peak')}
    positions = {m: [] for m in METHODS}
    for seed in range(3):
        training = data[[j for j in range(3) if j != seed]].mean((0, 1))
        x, heldout = derivative(data[seed].mean(0))
        slopes['known_peak'].append(float(np.interp(true_peak, x, heldout)))
        for method in METHODS:
            position = peak_position(GRID, training_profile(training, method))
            positions[method].append(position - true_peak)
            slopes[method].append(float(np.interp(position, x, heldout)))
    result = {m: dict(normalized_slope=float(np.mean(v)) * exact,
                     peak_errors=positions.get(m, [0.] * 3)) for m, v in slopes.items()}
    if any(v['normalized_slope'] <= 0 for v in result.values()):
        raise ValueError('Nonpositive slope')
    return result


def summarize(rows):
    result = {}
    for method in (*METHODS, 'known_peak'):
        slopes = [r[method]['normalized_slope'] for r in rows]
        errors = np.asarray([r[method]['peak_errors'] for r in rows])
        result[method] = dict(**decompose(slopes), peak_mean_error=float(errors.mean()),
                             peak_rmse=float(np.sqrt(np.mean(errors**2))))
    contrasts = {}
    for a, b in [('smooth', 'raw'), ('monotone', 'raw'),
                 ('smooth_monotone', 'smooth'), ('smooth_monotone', 'raw'),
                 ('smooth_monotone', 'known_peak')]:
        delta = np.asarray([1/r[a]['normalized_slope'] - 1/r[b]['normalized_slope'] for r in rows])
        contrasts[a+' minus '+b] = dict(mean=float(delta.mean()),
            mc_se=float(delta.std(ddof=1)/math.sqrt(len(delta))))
    return dict(methods=result, paired_thickness_contrasts=contrasts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--repeats', type=int, default=256)
    a = p.parse_args()
    if a.repeats < 2:
        p.error('At least two repetitions required')
    if a.output.exists():
        raise FileExistsError(a.output)
    result = dict(status='diagnostic_only', seed_start=3600000,
                  repeats=a.repeats, noise=.03, numpy_version=np.__version__,
                  smoothing_sigma=.25, derivative_width=1., cases=[], noiseless=[])
    for kind in ('logistic', 'gompertz', 'flat_peak'):
        exact = 4. if kind == 'flat_peak' else truth(kind)[1]
        clean = [evaluate(fixture(kind, 3600000+i, 0., 0., noise=0.)[1], kind, exact) for i in range(8)]
        result['noiseless'].append(dict(shape=kind, **summarize(clean)))
        for phi in (0., .8):
            for spatial in (0., .2):
                rows, errors = [], []
                for i in range(a.repeats):
                    try:
                        data = fixture(kind, 3600000+i, phi, spatial)[1]
                        rows.append(dict(repetition=i, **evaluate(data, kind, exact)))
                    except ValueError as e:
                        errors.append(dict(repetition=i, error=str(e)))
                c = dict(shape=kind, phi=phi, spatial=spatial, errors=errors,
                         summary=summarize(rows) if len(rows) >= 2 else {}, runs=rows)
                result['cases'].append(c)
                print(kind, phi, spatial, len(rows), len(errors), flush=True)
    if any(c['errors'] for c in result['cases']):
        result['status'] = 'diagnostic_with_errors'
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
