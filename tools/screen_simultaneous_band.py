"""Independent fixed-grid simultaneous-band screen with frozen Bonferroni t law."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import numpy as np
from screen_fixed_profile_coverage import PROTOCOL as BASE, truth, estimate
from control_front_shared_noise import ar1, profile, loading
from shock_validation import digest

PROTOCOL=BASE.with_name('SIMULTANEOUS_BAND_PROTOCOL.json')


def t7_survival(x):
    # t=sqrt(7)*tan(theta); integral becomes cos(theta)^6.
    # Complementary angle avoids subtraction from one in the upper tail.
    if not math.isfinite(x) or x<0: raise ValueError('Finite nonnegative t required')
    angle=math.atan2(math.sqrt(7),x);s=math.sin(angle);c=math.cos(angle)
    integral=angle
    for n in (2,4,6): integral=(n-1)*integral/n-s**(n-1)*c/n
    return 16*integral/(5*math.pi)


def critical_value(nodes,alpha=.05):
    if nodes<1 or not 0<alpha<1: raise ValueError('Valid family size and alpha required')
    target=alpha/(2*nodes);lo=0.;hi=1.
    while t7_survival(hi)>target: hi*=2
    for _ in range(70):
        mid=(lo+hi)/2
        if t7_survival(mid)>target:lo=mid
        else:hi=mid
    return hi


def coverage_summary(values):
    v=np.asarray(values,float);n=len(v);p=float(v.mean());z=1.959963984540054
    radius=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    center=(p+z*z/(2*n))/(1+z*z/n)
    return dict(coverage=p,mcse=float(v.std(ddof=1)/np.sqrt(n)),wilson95=[center-radius,center+radius])


def trial(spec):
    i,p,b,critical=spec
    rng=np.random.default_rng(np.random.SeedSequence([p['seed'],i]))
    x=np.linspace(*b['x_bounds'],b['nodes']);target=truth(x,b['shift_sd'],64)
    shape=(b['groups']*b['runs_per_group'],b['windows'])
    shifts=b['shift_sd']*ar1(rng,shape,b['temporal_phi'])
    node=b['node_noise_sd']*ar1(rng,shape+(len(x),),b['temporal_phi'])
    common=b['common_noise_sd']*ar1(rng,shape,b['temporal_phi'])
    offsets=b['run_offset_sd']*rng.normal(size=(shape[0],1,1))
    offsets+=np.repeat(b['group_offset_sd']*rng.normal(size=(b['groups'],1,1)),b['runs_per_group'],axis=0)
    base=profile(x,shifts)+node+offsets;result={}
    for case in b['cases']:
        mean,half=estimate(base+common[...,None]*loading(x,case),b['groups'],b['runs_per_group'],critical)
        errors=mean-target
        result[case]=dict(covered=bool(np.all(abs(errors)<=half)),mean_half_width=float(half.mean()),
            max_abs_error=float(abs(errors).max()),selected_errors=errors[[len(x)//2,int(np.argmin(abs(x-10)))]].tolist())
    return result


def screen(workers=1):
    p=json.loads(PROTOCOL.read_text());b=json.loads(BASE.read_text());hashes={PROTOCOL.name:digest(PROTOCOL),BASE.name:digest(BASE)}
    if p['nodes']!=b['nodes'] or p['groups']!=b['groups'] or b['groups']!=8:raise ValueError('Protocol mismatch')
    if workers<1:raise ValueError('Positive worker count required')
    critical=critical_value(p['nodes'],p['family_alpha'])
    specs=[(i,p,b,critical) for i in range(p['trials'])]
    if workers==1:rows=list(map(trial,specs))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(trial,specs,chunksize=8))
    summaries={}
    for case in b['cases']:
        s=coverage_summary([r[case]['covered'] for r in rows]);criterion=p['screen_criterion']
        s.update(mean_half_width=float(np.mean([r[case]['mean_half_width'] for r in rows])),
            criterion_met=s['wilson95'][0]>=criterion['minimum_wilson_lower_bound'] and s['mcse']<=criterion['maximum_coverage_mcse'])
        summaries[case]=s
    if any(digest(PROTOCOL.with_name(n))!=h for n,h in hashes.items()):raise ValueError('Protocol changed')
    return dict(status='synthetic_fixed_grid_band_screen',physical_acceptance='incomplete',protocol=p,base_protocol=b,
        protocol_sha256=hashes,critical_value=critical,width_ratio_to_pointwise=critical/b['critical_value'],numpy_version=np.__version__,
        summaries=summaries,trials=rows,all_screen_criteria_met=all(s['criterion_met'] for s in summaries.values()),
        sources={n:digest(Path(__file__).with_name(n)) for n in ('screen_simultaneous_band.py','screen_fixed_profile_coverage.py','control_front_shared_noise.py')},
        limitations=['Only the fixed grid is covered; no continuous-profile or gradient-thickness guarantee.',
            'Bonferroni controls family error if individual t intervals are valid; this screen probes illustrative nonlinear group means only.',
            'Conservative width is reported; no physical noise calibration or physical campaign qualification.'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--workers',type=int,default=1);a=parser.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    result=screen(a.workers)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(critical=result['critical_value'],summaries=result['summaries'])))
