"""Isolate training-front alignment and template recentering on Gompertz."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from experiment_peak_smoothing import training_profile
from shock_alignment import GRID, template
from shock_analysis import crossing
from shock_crossfit_gradient import derivative, peak_position
from decompose_gradient_bias import decompose

METHODS = ('exact_training', 'exact_training_recentered', 'crossing_training_recentered')


def evaluate(trials, oracle):
    slopes = {m: [] for m in METHODS}
    positions = {m: [] for m in METHODS}
    true_peak = math.log(math.log(2))
    for seed in range(3):
        y = training_profile(oracle[[j for j in range(3) if j != seed]].mean((0, 1)), 'smooth_monotone')
        center = crossing(GRID, y, .5)
        tx, ty = template(trials, seed)
        templates = ((GRID, y), (GRID-center, y), (tx, ty))
        x, d = derivative(oracle[seed].mean(0))
        for method, (px, py) in zip(METHODS, templates):
            peak = peak_position(px, py)
            if not x[0] <= peak <= x[-1]:
                raise ValueError('Peak outside measurement coverage')
            slopes[method].append(float(np.interp(peak, x, d)))
            positions[method].append(peak-true_peak)
    result = {m: dict(normalized_slope=float(np.mean(slopes[m]))*math.e,
                     peak_errors=positions[m]) for m in METHODS}
    if any(v['normalized_slope'] <= 0 for v in result.values()):
        raise ValueError('Nonpositive slope')
    return result


def summarize(rows):
    methods = {}
    for m in METHODS:
        errors = np.asarray([r[m]['peak_errors'] for r in rows])
        methods[m] = dict(**decompose([r[m]['normalized_slope'] for r in rows]),
                          peak_mean_error=float(errors.mean()),
                          peak_rmse=float(np.sqrt(np.mean(errors**2))))
    contrasts = {}
    for a, b in [(METHODS[1], METHODS[0]), (METHODS[2], METHODS[1])]:
        delta = np.asarray([1/r[a]['normalized_slope']-1/r[b]['normalized_slope'] for r in rows])
        contrasts[a+' minus '+b] = dict(mean=float(delta.mean()),
            mc_se=float(delta.std(ddof=1)/math.sqrt(len(delta))))
    return dict(methods=methods, paired_thickness_contrasts=contrasts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=256)
    a = p.parse_args()
    if a.repeats < 2:
        p.error('At least two repetitions required')
    if a.output.exists():
        raise FileExistsError(a.output)
    result = dict(status='diagnostic_only', shape='gompertz', seed_start=3600000,
                  repeats=a.repeats, noise=.03, numpy_version=np.__version__, cases=[])
    clean = [evaluate(*fixture('gompertz', 3600000+i, 0., 0., noise=0.)) for i in range(8)]
    result['noiseless'] = summarize(clean)
    for phi in (0., .8):
        for spatial in (0., .2):
            rows, errors = [], []
            for i in range(a.repeats):
                try:
                    rows.append(dict(repetition=i, **evaluate(*fixture('gompertz', 3600000+i, phi, spatial))))
                except ValueError as e:
                    errors.append(dict(repetition=i, error=str(e)))
            c = dict(phi=phi, spatial=spatial, errors=errors,
                     summary=summarize(rows) if len(rows) >= 2 else {}, runs=rows)
            result['cases'].append(c)
            print(phi, spatial, len(rows), len(errors), flush=True)
    if any(c['errors'] for c in result['cases']):
        result['status'] = 'diagnostic_with_errors'
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
