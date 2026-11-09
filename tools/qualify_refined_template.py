"""Fixed one-pass refinement screen on fresh paired synthetic samples."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID, align
from shock_crossfit_gradient import derivative, estimate_gradient, peak_position
from shock_mean_profile import measure
from validate_shock_estimator import truth


def run_case(spec):
    kind, phi, spatial, repeats, seed_start, fit_width, noise, scope = spec
    exact = 4. if kind == 'flat_peak' else truth(kind)[1]
    if kind=='flat_peak':
        lo,hi=0.,1.
        for _ in range(60):
            z=(lo+hi)/2
            if .5+5/8*(z-z**5/5)<.9:lo=z
            else:hi=z
        exact_width=5*(lo+hi)/2
    else:exact_width=truth(kind)[0]
    true_peak = math.log(math.log(2)) if kind == 'gompertz' else 0.
    rows, errors = [], []
    for i in range(repeats):
        try:
            trials, oracle = fixture(kind, seed_start+i, phi, spatial,noise=noise)
            row = dict(repetition=i)
            for passes, name in [(0, 'baseline'), (1, 'refined')]:
                data, diagnostics = align(trials, refine_passes=passes)
                width = fit_width if passes else 1.
                estimate = estimate_gradient(data, diagnostics, width=width)
                peaks, controlled_slopes = [], []
                for seed, info in enumerate(diagnostics):
                    pos = peak_position(info['template_x'], info['template_density'],width=width)
                    peaks.append(pos-true_peak)
                    x, d = derivative(oracle[seed].mean(0),width=width)
                    controlled_slopes.append(float(np.interp(pos, x, d)))
                if min(controlled_slopes) <= 0:
                    raise ValueError('Nonpositive controlled slope')
                metrics = measure(data.mean((0,1)),fit_width=width)[0]
                row[name] = dict(gradient_bias=estimate['gradient_thickness']/exact-1,
                    max_gradient_bias=float(metrics[1]/exact-1),
                    controlled_gradient_bias=float(1/np.mean(controlled_slopes)/exact-1),
                    peak_errors=peaks, width=float(metrics[0]),width_bias=float(metrics[0]/exact_width-1))
            rows.append(row)
        except ValueError as e:
            errors.append(dict(repetition=i,error=str(e)))
    summaries = {}
    for name in ('baseline', 'refined'):
        summary = {}
        for field in ('gradient_bias','max_gradient_bias','controlled_gradient_bias','width_bias'):
            v = np.array([r[name][field] for r in rows])
            if len(v) >= 2:
                summary[field] = dict(mean=float(v.mean()), sd=float(v.std(ddof=1)),
                    mc_se=float(v.std(ddof=1)/math.sqrt(len(v))),
                    bias_plus_two_se=float(abs(v.mean())+2*v.std(ddof=1)/math.sqrt(len(v))))
        peaks = np.array([r[name]['peak_errors'] for r in rows])
        summary['peak_rmse'] = float(np.sqrt(np.mean(peaks**2))) if len(rows) else None
        summary['peak_mean_error'] = float(peaks.mean()) if len(rows) else None
        summaries[name] = summary
    passed = False
    if not errors and len(rows) == repeats:
        b = summaries['baseline']['gradient_bias'];q = summaries['refined']['gradient_bias']
        passed = q['bias_plus_two_se'] <= .015 and q['sd'] <= 1.25*b['sd']
        if scope=='metrics':passed = passed and summaries['refined']['width_bias']['bias_plus_two_se']<=.005
    contrasts = {}
    if len(rows) >= 2:
        for field in ('gradient_bias', 'controlled_gradient_bias'):
            v = np.array([r['refined'][field]-r['baseline'][field] for r in rows])
            contrasts[field] = dict(mean=float(v.mean()),mc_se=float(v.std(ddof=1)/math.sqrt(len(v))))
    return dict(shape=kind,phi=phi,spatial=spatial,errors=errors,
                passed=bool(passed),summary=summaries,paired_changes=contrasts,runs=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeats',type=int,default=64)
    p.add_argument('--jobs',type=int,default=3)
    p.add_argument('--fit-width',type=float,choices=[1.,1.5],default=1.)
    p.add_argument('--seed-start',type=int,default=4800000)
    p.add_argument('--noise',type=float,choices=[.015,.03],default=.03)
    p.add_argument('--scope',choices=['gradient','metrics'],default='gradient')
    a=p.parse_args()
    if a.repeats<2 or a.jobs<1:p.error('Need at least two repetitions and one worker')
    if a.output.exists():raise FileExistsError(a.output)
    specs=[(kind,phi,spatial,a.repeats,a.seed_start,a.fit_width,a.noise,a.scope) for kind in ('logistic','gompertz','flat_peak') for phi in (0.,.8) for spatial in (0.,.2)]
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for case in pool.map(run_case,specs):
            cases.append(case)
            print(case['shape'],case['phi'],case['spatial'],case['passed'],case['summary']['refined'],flush=True)
    result=dict(status='screen_passed' if all(c['passed'] for c in cases) else 'screen_not_qualified',
                seed_start=a.seed_start,repeats=a.repeats,noise=a.noise,scope=a.scope,refine_passes=1,fit_width=a.fit_width,baseline_fit_width=1.,numpy_version=np.__version__,
                criteria=dict(gradient_bias_plus_two_se_max=.015,sd_ratio_max=1.25,width_bias_plus_two_se_max=.005 if a.scope=='metrics' else None),cases=cases)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
