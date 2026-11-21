"""Verify complete high-order calibration evidence; never certify coverage."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from calibrate_shock_ensemble import summarize
from shock_high_order_estimator import METHOD, WIDTHS
from shock_validation import digest, save
from validate_shock_estimator import truth

SOURCES=('calibrate_high_order_shock.py','shock_high_order_estimator.py','calibrate_shock_ensemble.py',
         'shock_ensemble_estimator.py','shock_alignment.py','shock_crossfit_gradient.py','shock_mean_profile.py',
         'shock_analysis.py','experiment_crossfit_gradient.py','validate_shock_estimator.py','shock_validation.py')


def validate(report,windows,repetitions):
    p=report['protocol'];start=p['seed_start']
    expected=dict(method=METHOD,group_count=12,subgroups=6,bootstrap_repeats=1024,windows=windows,
                  seed_start=start,noise_levels=[.005,.015,.03],derivative_degree=5,widths=list(WIDTHS),
                  peak_training_degree=3,peak_training_width=1.5,primary_interpolation='cubic',
                  sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES})
    if p!=expected or type(start)!=int or start<0 or windows not in (32,64):
        raise ValueError('High-order protocol or source mismatch')
    if report['group_count']!=12 or report['seed_start']!=start or report['center_profile'] is not True:
        raise ValueError('Report metadata mismatch')
    cases=report['cases'];wanted={(k,t,s,n) for k in ('logistic','gompertz') for t in (0.,.8) for s in (0.,.2) for n in (.005,.015,.03)}
    if len(cases)!=24 or {(c['shape'],c['phi'],c['spatial'],c['noise']) for c in cases}!=wanted:
        raise ValueError('Complete 24-case matrix required')
    failures=0
    for c in cases:
        rows=c['runs'];ids=[r['repetition'] for r in rows]
        if (c['trial_count']!=repetitions or len(rows)!=repetitions or any(type(i)!=int for i in ids) or
                set(ids)!=set(range(repetitions))):
            raise ValueError('Missing, repeated or unexpected trials')
        for r in rows:
            if 'error' in r:
                if set(r)!={'repetition','error'} or not isinstance(r['error'],str) or not r['error']:
                    raise ValueError('Invalid rejected-trial record')
                continue
            arrays=[np.asarray(r[k],float) for k in ('estimate','statistical_radius','processing_sensitivity','relative_error')]
            if any(a.shape!=(2,) or not np.isfinite(a).all() for a in arrays) or any(np.any(a<0) for a in arrays[1:3]):
                raise ValueError('Invalid interval or estimate')
            value,u,processing,relative=arrays;exact=truth(c['shape'])
            stat=(abs(value-exact)<=u).tolist();total=(abs(value-exact)<=u+processing).tolist()
            if (r['statistical_covers']!=stat or r['with_processing_covers']!=total or
                    any(type(v)!=bool for k in ('statistical_covers','with_processing_covers') for v in r[k]) or
                    not np.allclose(relative,value/exact-1,rtol=1e-12,atol=1e-14)):
                raise ValueError('Metric coverage or relative error inconsistent')
            numbers=[r[k] for k in ('density_error','density_radius','density_statistical_radius','density_processing_sensitivity')]
            if any(not math.isfinite(v) or v<0 for v in numbers):raise ValueError('Invalid density error or radius')
            error,radius,statistical,processing=numbers
            if (type(r['density_covers'])!=bool or r['density_covers']!=(error<=radius) or
                    not math.isclose(radius,statistical+processing,rel_tol=1e-12,abs_tol=1e-14)):
                raise ValueError('Density coverage or budget inconsistent')
        summary=summarize(rows)
        if any(c.get(k)!=v for k,v in summary.items() if k!='runs'):
            raise ValueError('Stored summary differs from trial outcomes')
        failures+=len(summary['errors'])
    expected_status='coverage_study_with_failures' if failures else 'coverage_study_only'
    if report['status']!=expected_status:raise ValueError('Study status hides rejected trials')
    return dict(windows=windows,trial_count=24*repetitions,rejected_trials=failures,
                minimum_metric_coverage=[min(c['with_processing_coverage'][i] for c in cases) for i in (0,1)],
                minimum_density_coverage=min(c['density_coverage'] for c in cases),
                fixture_seed_range=[start,start+(repetitions-1)*100+71],
                bootstrap_seed_range=[start+1000000,start+1000000+repetitions-1])


def verify_pair(full,half,screens):
    if len(screens)!=2:raise ValueError('Both screening reports required')
    paths=[Path(full),Path(half),*map(Path,screens)];hashes=[digest(p) for p in paths]
    summaries=[validate(json.loads(p.read_text()),w,n) for p,w,n in zip(paths,(64,32,64,32),(128,128,8,8))]
    ranges=[(i,s[k]) for i,s in enumerate(summaries) for k in ('fixture_seed_range','bootstrap_seed_range')]
    for index,(i,a) in enumerate(ranges):
        for j,b in ranges[index+1:]:
            if max(a[0],b[0])<=min(a[1],b[1]):raise ValueError('Overlapping calibration or screening seed ranges')
    if hashes!=[digest(p) for p in paths]:raise ValueError('Coverage report changed during verification')
    return dict(status='complete_matching_evidence',method=METHOD,physical_acceptance='incomplete',
                full=summaries[0],half=summaries[1],screening=summaries[2:],report_sha256=hashes,
                tool_sha256=digest(Path(__file__)),
                limitations=['Completion, protocol identity and internal consistency are checked; nominal coverage is not certified.',
                             'Seed separation does not replace scientific assessment or exclude undocumented model selection.',
                             'Synthetic coverage does not establish numerical convergence or experimental agreement.'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coverage',type=Path,required=True)
    p.add_argument('--half-coverage',type=Path,required=True)
    p.add_argument('--screen',type=Path,nargs=2,required=True,metavar=('FULL','HALF'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    result=verify_pair(a.coverage,a.half_coverage,a.screen)
    a.output.parent.mkdir(parents=True,exist_ok=True);save(a.output,result)
    print('Complete matching calibration evidence; no qualification verdict.')


if __name__=='__main__':main()
