"""Separate slope error and reciprocal nonlinearity; diagnostic, not acceptance."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID, align
from shock_crossfit_gradient import derivative, estimate_gradient, peak_position
from validate_shock_estimator import shape, truth


def decompose(slopes):
    """Exact empirical identity for slopes normalized by the analytic slope."""
    s = np.asarray(slopes, float)
    if s.ndim != 1 or len(s) < 2 or not np.isfinite(s).all() or np.any(s <= 0):
        raise ValueError('At least two finite positive slopes required')
    mean = float(s.mean())
    reciprocal = 1 / s
    slope_component = 1 / mean - 1
    return dict(mean_relative_slope=mean - 1,
                slope_sd=float(s.std(ddof=1)),
                thickness_bias=float(reciprocal.mean() - 1),
                thickness_mc_se=float(reciprocal.std(ddof=1) / math.sqrt(len(s))),
                inverse_mean_slope_bias=slope_component,
                reciprocal_jensen_gap=float(reciprocal.mean() - 1 / mean))


def slopes_for(trials, oracle, kind, exact):
    position = math.log(math.log(2)) if kind == 'gompertz' else 0.
    x, d = derivative(oracle.mean((0, 1)))
    known = float(np.interp(position, x, d))
    selected = []
    for seed in range(3):
        training = oracle[[j for j in range(3) if j != seed]].mean((0, 1))
        peak = peak_position(GRID, training)
        x, d = derivative(oracle[seed].mean(0))
        selected.append(float(np.interp(peak, x, d)))
    data, diagnostics = align(trials)
    candidate = estimate_gradient(data, diagnostics)
    return dict(oracle_known_peak=known * exact,
                oracle_selected_peak=float(np.mean(selected)) * exact,
                template_selected_peak=exact / candidate['gradient_thickness'])


def analytic_control(kind, exact):
    # Express the analytic profile in coordinates centered at its 50% crossing.
    center = -math.log(math.log(2)) if kind == 'gompertz' else 0.
    if kind == 'flat_peak':
        z = np.clip(GRID / 2.5, -1, 1)
        y = .5 + 5 / 8 * (z - z**5 / 5)
    else:
        y = shape(GRID + center, kind)
    x, d = derivative(y)
    normalized = float(np.interp(-center, x, d)) * exact
    return dict(normalized_slope=normalized, thickness_bias=1 / normalized - 1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=128)
    a = p.parse_args()
    if a.repeats < 2:
        p.error('Need at least two repetitions')
    if a.output.exists():
        raise FileExistsError(a.output)
    result = dict(status='diagnostic_only', seed_start=2400000,
                  repeats=a.repeats, numpy_version=np.__version__,
                  noise=.03, noiseless=[], cases=[])
    for kind in ('logistic', 'gompertz', 'flat_peak'):
        exact = 4. if kind == 'flat_peak' else truth(kind)[1]
        clean = []
        for i in range(8):
            trials, oracle = fixture(kind, 2400000+i, 0., 0., noise=0.)
            clean.append(slopes_for(trials, oracle, kind, exact))
        result['noiseless'].append(dict(shape=kind, direct=analytic_control(kind, exact),
            translated={k: decompose([r[k] for r in clean]) for k in clean[0]}))
        for phi in (0., .8):
            for spatial in (0., .2):
                rows, errors = [], []
                for i in range(a.repeats):
                    try:
                        trials, oracle = fixture(kind, 2400000+i, phi, spatial)
                        rows.append(dict(repetition=i, **slopes_for(trials, oracle, kind, exact)))
                    except ValueError as e:
                        errors.append(dict(repetition=i, error=str(e)))
                summary = {k: decompose([r[k] for r in rows]) for k in
                           ('oracle_known_peak', 'oracle_selected_peak', 'template_selected_peak')} if len(rows) >= 2 else {}
                case = dict(shape=kind, phi=phi, spatial=spatial,
                            summary=summary, errors=errors, runs=rows)
                result['cases'].append(case)
                print(json.dumps({k: v for k, v in case.items() if k != 'runs'}), flush=True)
    if any(c['errors'] for c in result['cases']):
        result['status'] = 'diagnostic_with_errors'
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
