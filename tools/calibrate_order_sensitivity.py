"""Independent screen for alignment-order intervals, including small front motion."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import initial, GRID
from shock_order_sensitivity import METHOD, subgroup_point, pooled_point, estimate
from calibrate_shock_ensemble import summarize
from shock_validation import digest, save
from validate_shock_estimator import truth, shape
from verify_high_order_coverage import SOURCES as BASE_SOURCES

SOURCES = tuple(BASE_SOURCES) + ('shock_order_sensitivity.py', 'calibrate_order_sensitivity.py', 'verify_high_order_coverage.py')


def moving_fixture(kind, seed, phi, spatial, noise, amplitude, windows):
    if kind not in ('logistic', 'gompertz') or amplitude < 0 or not math.isfinite(amplitude):
        raise ValueError('Supported shape and nonnegative motion amplitude required')
    if windows not in (32, 64):
        raise ValueError('32 or 64 windows required')
    raw, _ = fixture(kind, seed, phi, spatial, noise=noise)
    clean, _ = fixture(kind, seed, phi, spatial, noise=0.)
    result = []
    for run, clean_run in zip(raw, clean):
        rows = []
        for window, clean_window in zip(run[:windows], clean_run[:windows]):
            x, y, _ = initial(window)
            cx, cy, _ = initial(clean_window)
            if not np.allclose(x, cx, rtol=1e-14, atol=1e-12):
                raise ValueError('Fixture grid mismatch')
            i = int(abs(cy - .5).argmin())
            inverse = math.log(cy[i] / (1 - cy[i])) if kind == 'logistic' else -math.log(-math.log(cy[i]))
            shift = x[i] - inverse
            np.testing.assert_allclose(shape(x - shift, kind), cy, rtol=1e-12, atol=1e-14)
            values = shape(x - amplitude * shift, kind) + y - cy
            rows.append(dict(x=x.tolist(), density_normalized=values.tolist()))
        result.append(rows)
    return result


def trial(spec, protocol):
    kind, phi, spatial, noise, amplitude, repetition = spec
    start = protocol['seed_start']; windows = protocol['windows']; groups = []
    for g in range(12):
        points = [subgroup_point(moving_fixture(kind, start + repetition * 100 + g * 6 + j,
                    phi, spatial, noise, amplitude, windows)) for j in range(6)]
        groups.append(pooled_point(points))
    result = estimate(groups, repeats=1024, seed=start + 1000000 + repetition)
    exact = truth(kind); center = -math.log(math.log(2)) if kind == 'gompertz' else 0.
    values = np.array([result[k]['mean'] for k in ('width_10_90', 'gradient_thickness')])
    u = np.array([result[k]['statistical_radius'] for k in ('width_10_90', 'gradient_thickness')])
    proc = np.array([result[k]['processing_sensitivity'] for k in ('width_10_90', 'gradient_thickness')])
    error = float(np.max(abs(np.asarray(result['density']) - shape(GRID + center, kind))))
    return dict(repetition=repetition, estimate=values.tolist(), statistical_radius=u.tolist(),
        processing_sensitivity=proc.tolist(), relative_error=(values/exact-1).tolist(),
        statistical_covers=(abs(values-exact)<=u).tolist(), with_processing_covers=(abs(values-exact)<=u+proc).tolist(),
        density_error=error, density_radius=result['density_simultaneous_radius'],
        density_covers=bool(error <= result['density_simultaneous_radius']),
        density_statistical_radius=result['density_statistical_radius'], density_processing_sensitivity=result['density_processing_sensitivity'],
        order_density_difference=result['alignment_order_density_difference'], order_metric_differences=result['alignment_order_metric_differences'])


def cached(task):
    spec, folder, protocol = task
    path = Path(folder) / ('-'.join(map(str, spec)) + '.json')
    if path.exists():
        saved = json.loads(path.read_text())
        if saved['protocol'] != protocol or saved['spec'] != list(spec):
            raise ValueError('Checkpoint protocol mismatch')
        return spec, saved['result']
    try:
        result = trial(spec, protocol)
    except ValueError as error:
        result = dict(repetition=spec[-1], error=str(error))
    save(path, dict(protocol=protocol, spec=spec, result=result))
    return spec, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed-start', type=int, required=True)
    parser.add_argument('--windows', type=int, choices=(32,64), default=64)
    parser.add_argument('--replicates', type=int, default=8)
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.replicates < 2 or args.jobs < 1 or args.seed_start <= 0 or args.replicates*100 >= 1000000:
        parser.error('Invalid seed or replication parameters')
    if args.output.exists():
        raise FileExistsError(args.output)
    folder = args.output.with_suffix('.parts'); folder.mkdir(parents=True, exist_ok=args.resume)
    protocol = dict(method=METHOD, status='screen_only', group_count=12, subgroups=6, bootstrap_repeats=1024,
        windows=args.windows, seed_start=args.seed_start, replicates=args.replicates,
        motion_amplitudes=[0., .05, 1.], sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES},
        limitation='Screen does not qualify coverage or model physical heteroscedastic noise; fresh complete validation remains required.')
    keys = [(k,t,s,n,a) for k in ('logistic','gompertz') for t in (0.,.8) for s in (0.,.2) for n in (.005,.015,.03) for a in protocol['motion_amplitudes']]
    tasks = [(key+(i,), str(folder), protocol) for key in keys for i in range(args.replicates)]
    rows = {key:[] for key in keys}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(cached, task) for task in tasks]
        for number, future in enumerate(as_completed(futures), 1):
            spec, result = future.result(); rows[tuple(spec[:-1])].append(result)
            print('completed', number, 'of', len(tasks), 'errors', sum('error' in r for values in rows.values() for r in values), flush=True)
    if any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()):
        raise ValueError('Frozen source changed')
    cases = [dict(shape=k[0],phi=k[1],spatial=k[2],noise=k[3],motion_amplitude=k[4],
        **summarize(sorted(rows[k],key=lambda r:r['repetition']))) for k in keys]
    save(args.output, dict(status='screen_with_failures' if any(c['errors'] for c in cases) else 'screen_only',
        physical_acceptance='incomplete', protocol=protocol, cases=cases))

if __name__ == '__main__':
    main()
