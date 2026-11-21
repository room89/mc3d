"""Fresh synthetic screening/qualification for the separate high-order candidate."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_ensemble_estimator import group_point
from shock_high_order_estimator import METHOD, WIDTHS, pooled_point, estimate
from calibrate_shock_ensemble import summarize
from shock_validation import digest, save
from validate_shock_estimator import truth, shape


def trial(spec):
    kind,phi,spatial,noise,i,seed_start,windows=spec
    groups=[]
    for g in range(12):
        points=[group_point([r[:windows] for r in fixture(kind,seed_start+i*100+g*6+j,phi,spatial,noise=noise)[0]],
                            center_profile=True,interpolation_sensitivity=True) for j in range(6)]
        groups.append(pooled_point(points))
    r=estimate(groups,repeats=1024,seed=seed_start+1000000+i)
    exact=truth(kind);center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    values=np.array([r[k]['mean'] for k in ('width_10_90','gradient_thickness')])
    u=np.array([r[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')])
    proc=np.array([r[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    error=float(np.max(abs(np.array(r['density'])-shape(GRID+center,kind))))
    return dict(repetition=i,estimate=values.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),
                relative_error=(values/exact-1).tolist(),statistical_covers=(abs(values-exact)<=u).tolist(),
                with_processing_covers=(abs(values-exact)<=u+proc).tolist(),density_error=error,
                density_radius=r['density_simultaneous_radius'],density_covers=bool(error<=r['density_simultaneous_radius']),
                density_statistical_radius=r['density_statistical_radius'],density_processing_sensitivity=r['density_processing_sensitivity'])


def cached(task):
    spec,folder,protocol=task;kind,phi,spatial,noise,i,*_=spec
    path=Path(folder)/f'{kind}-{phi}-{spatial}-{noise}-{i}.json'
    if path.exists():
        saved=json.loads(path.read_text())
        if saved['protocol']!=protocol or saved['spec']!=list(spec): raise ValueError('Checkpoint protocol mismatch')
        return saved['result']
    try: result=trial(spec)
    except ValueError as error: result=dict(repetition=i,error=str(error))
    save(path,dict(protocol=protocol,spec=spec,result=result));return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--replicates',type=int,default=8);p.add_argument('--jobs',type=int,default=5)
    p.add_argument('--windows',type=int,choices=(32,64),default=64);p.add_argument('--seed-start',type=int,required=True)
    p.add_argument('--noise-levels',nargs='+',type=float,default=[.005,.015,.03]);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.replicates<2 or a.jobs<1: p.error('Invalid replication counts')
    if not a.noise_levels or len(set(a.noise_levels))!=len(a.noise_levels) or any(not math.isfinite(n) or n<=0 for n in a.noise_levels):
        p.error('Distinct finite positive noise levels required')
    if a.output.exists():raise FileExistsError(a.output)
    folder=a.output.with_suffix('.parts');folder.mkdir(parents=True,exist_ok=a.resume)
    names=('calibrate_high_order_shock.py','shock_high_order_estimator.py','calibrate_shock_ensemble.py',
           'shock_ensemble_estimator.py','shock_alignment.py','shock_crossfit_gradient.py','shock_mean_profile.py',
           'shock_analysis.py','experiment_crossfit_gradient.py','validate_shock_estimator.py','shock_validation.py')
    protocol=dict(method=METHOD,group_count=12,subgroups=6,bootstrap_repeats=1024,windows=a.windows,
                  seed_start=a.seed_start,noise_levels=a.noise_levels,derivative_degree=5,widths=list(WIDTHS),
                  peak_training_degree=3,peak_training_width=1.5,primary_interpolation='cubic',
                  sources={n:digest(Path(__file__).with_name(n)) for n in names})
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for kind in ('logistic','gompertz'):
            for phi in (0.,.8):
                for spatial in (0.,.2):
                    for noise in a.noise_levels:
                        tasks=[((kind,phi,spatial,noise,i,a.seed_start,a.windows),str(folder),protocol) for i in range(a.replicates)]
                        rows=list(pool.map(cached,tasks));case=dict(shape=kind,phi=phi,spatial=spatial,noise=noise,**summarize(rows));cases.append(case)
                        print(kind,phi,spatial,noise,case['valid_count'],case['with_processing_coverage'],case['density_coverage'],flush=True)
    if any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()):raise ValueError('Source changed')
    save(a.output,dict(status='coverage_study_with_failures' if any(c['errors'] for c in cases) else 'coverage_study_only',
                       group_count=12,seed_start=a.seed_start,center_profile=True,protocol=protocol,cases=cases))


if __name__=='__main__':main()
