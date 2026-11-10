"""Pointwise interval diagnostic against digitized primary data; never acceptance."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from shock_alignment import GRID


def interval(value):
    if len(value) != 2 or not all(math.isfinite(x) for x in value) or value[0] > value[1]:
        raise ValueError('Finite ordered interval required')
    return value


def compare(primary, summary, model_length):
    if not math.isfinite(model_length) or model_length <= 0:
        raise ValueError('Positive model reference length required')
    density = np.asarray(summary['density'], float)
    if summary.get('center_profile') is not True or density.shape != GRID.shape or not np.isfinite(density).all():
        raise ValueError('Centered estimator output on the 401-point comparison grid required')
    radius = summary['density_simultaneous_radius']
    scale = primary['reference_length_m'] / model_length
    length_error = primary['stated_relative_length_accuracy']
    density_error = primary['stated_normalized_density_accuracy']
    if not all(math.isfinite(x) for x in (radius, scale, length_error, density_error)) or radius < 0 or scale <= 0 or not 0 <= length_error < 1 or density_error < 0:
        raise ValueError('Invalid accuracy or scale')
    results = []
    for point in primary['points']:
        x = interval(point['x_paper_interval'])
        y = interval(point['rho_interval'])
        coordinates = [a * scale * b for a, b in itertools.product(x, (1-length_error, 1+length_error))]
        lo, hi = min(coordinates), max(coordinates)
        if lo < GRID[0] or hi > GRID[-1]:
            raise ValueError('Experimental point exceeds model support; extrapolation forbidden')
        # Include interior vertices: a noisy sampled profile need not be monotone.
        positions = np.r_[lo, GRID[(GRID > lo) & (GRID < hi)], hi]
        values = np.interp(positions, GRID, density)
        model = [float(values.min()-radius), float(values.max()+radius)]
        experiment = [y[0]-density_error, y[1]+density_error]
        gap = max(experiment[0]-model[1], model[0]-experiment[1], 0.)
        results.append(dict(pixel_row=point['pixel_row'], model_x_interval=[lo, hi],
                            model_density_interval=model, experimental_density_interval=experiment,
                            nonoverlap_gap=float(gap)))
    if not results:
        raise ValueError('Experimental points required')
    return dict(status='pointwise_diagnostic_only', physical_acceptance='incomplete',
                model_reference_length_m=model_length, reference_length_ratio=scale,
                point_count=len(results), nonoverlap_count=sum(p['nonoverlap_gap'] > 0 for p in results),
                maximum_nonoverlap_gap=max(p['nonoverlap_gap'] for p in results), points=results,
                limitations=['Pointwise overlap is necessary, not sufficient, for a common global scale and origin.',
                             'Raster samples are correlated; no confidence level or independent-trial count is inferred.',
                             'Unmeasured gaps and tails, scan distortion and manual trace reproducibility remain unresolved.',
                             'Numerical convergence and method qualification must be assessed separately.',
                             'Experimental and digitization intervals are propagated without fitting spatial scale.'])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--primary', type=Path, required=True)
    p.add_argument('--analysis', type=Path, required=True)
    p.add_argument('--variant', default='C')
    p.add_argument('--model-reference-length-m', type=float, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    source = json.loads(a.analysis.read_text())
    # The summary has no explicit coordinate array; verify its producer convention.
    names = ('shock_ensemble_estimator.py', 'shock_alignment.py')
    for name in names:
        if source['protocol']['sources'].get(name) != digest(Path(__file__).with_name(name)):
            raise ValueError('Estimator source mismatch: '+name)
    result = compare(json.loads(a.primary.read_text()), source['summaries'][a.variant], a.model_reference_length_m)
    result.update(variant=a.variant, source_hashes={'primary': digest(a.primary), 'analysis': digest(a.analysis)},
                  tool_sha256=digest(Path(__file__)))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'point_count', 'nonoverlap_count', 'maximum_nonoverlap_gap')}))


if __name__ == '__main__':
    main()
