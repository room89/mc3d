"""Pilot of independent-group inference; no production simulations."""
import argparse,hashlib,json,math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_ensemble_estimator import group_point,estimate,GROUPS
from shock_group_pooling import pool_points
from validate_shock_estimator import truth,shape


def run(spec):
    kind,phi,spatial,noise,i,seed_start,center_profile=spec[:7]
    width_band=spec[7] if len(spec)>7 else False
    interpolation_sensitivity=spec[8] if len(spec)>8 else False
    windows=spec[9] if len(spec)>9 else 64
    subgroups=spec[10] if len(spec)>10 else 1
    groups=[]
    for g in range(GROUPS):
        points=[group_point([r[:windows] for r in fixture(kind,seed_start+i*100+g*subgroups+j,phi,spatial,noise=noise)[0]],
            center_profile=center_profile,interpolation_sensitivity=interpolation_sensitivity) for j in range(subgroups)]
        groups.append(pool_points(points))
    r=estimate(groups,repeats=1024,seed=seed_start+1000000+i,center_profile=center_profile,width_band=width_band,
        interpolation_sensitivity=interpolation_sensitivity)
    exact=truth(kind);center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    values=np.array([r[k]['mean'] for k in ('width_10_90','gradient_thickness')])
    u=np.array([r[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')])
    proc=np.array([r[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    error=float(np.max(abs(np.array(r['density'])-shape(GRID+center,kind))))
    return dict(repetition=i,estimate=values.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),
                relative_error=(values/exact-1).tolist(),statistical_covers=(abs(values-exact)<=u).tolist(),
                with_processing_covers=(abs(values-exact)<=u+proc).tolist(),density_error=error,
                density_radius=r['density_simultaneous_radius'],density_covers=bool(error<=r['density_simultaneous_radius']),
                density_statistical_radius=r.get('density_statistical_radius',r['density_simultaneous_radius']),
                density_processing_sensitivity=r.get('density_processing_sensitivity',0.))


def run_cached(task):
    spec,folder,protocol=task
    kind,phi,spatial,noise,i,_,_=spec[:7]
    checkpoint=Path(folder)/f'{kind}-{phi}-{spatial}-{noise}-{i}.json'
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        if saved['protocol']!=protocol or saved['spec']!=list(spec):
            raise ValueError('Checkpoint protocol mismatch')
        return saved['result']
    try:
        row=run(spec)
    except ValueError as error:
        # An unresolved estimator is a failed trial, never an omitted observation.
        row=dict(repetition=i,error=str(error))
    temporary=checkpoint.with_suffix('.tmp')
    temporary.write_text(json.dumps(dict(protocol=protocol,spec=spec,result=row),allow_nan=False)+'\n')
    temporary.replace(checkpoint)
    return row


def summarize(rows):
    if not rows:raise ValueError('At least one trial required')
    return dict(trial_count=len(rows),valid_count=sum('error' not in r for r in rows),
        errors=[r for r in rows if 'error' in r],
        statistical_coverage=np.mean([r.get('statistical_covers',[False,False]) for r in rows],axis=0).tolist(),
        with_processing_coverage=np.mean([r.get('with_processing_covers',[False,False]) for r in rows],axis=0).tolist(),
        density_coverage=float(np.mean([r.get('density_covers',False) for r in rows])),runs=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--replicates',type=int,default=32);p.add_argument('--jobs',type=int,default=3)
    p.add_argument('--seed-start',type=int,default=12000000)
    p.add_argument('--center-profile',action='store_true')
    p.add_argument('--matrix',action='store_true')
    p.add_argument('--focused',action='store_true',help='Four low-noise/high-spatial-noise diagnostic cases')
    p.add_argument('--width-band',action='store_true',help='Also bound width by inversion of the simultaneous density band')
    p.add_argument('--noise-levels',nargs='+',type=float)
    p.add_argument('--interpolation-sensitivity',action='store_true')
    p.add_argument('--windows',type=int,choices=[32,64],default=64)
    p.add_argument('--subgroups',type=int,choices=[1,6],default=1,help='Independent three-seed subgroups averaged within each of 12 groups')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if a.noise_levels is None:a.noise_levels=[.015,.03] if a.matrix or a.focused else [.03]
    if a.replicates<2 or a.jobs<1:p.error('Invalid replication counts')
    if a.matrix and a.focused:p.error('Choose either --matrix or --focused')
    if any(not math.isfinite(n) or n<=0 for n in a.noise_levels) or len(set(a.noise_levels))!=len(a.noise_levels):
        p.error('Noise levels must be distinct, finite and positive')
    if a.output.exists():raise FileExistsError(a.output)
    folder=a.output.with_suffix('.parts')
    if folder.exists() and not a.resume:raise FileExistsError(folder)
    folder.mkdir(parents=True,exist_ok=True)
    protocol=dict(group_count=GROUPS,bootstrap_repeats=1024,seed_start=a.seed_start,
        center_profile=a.center_profile,width_band=a.width_band,interpolation_sensitivity=a.interpolation_sensitivity,
        noise_levels=a.noise_levels,windows=a.windows,subgroups=a.subgroups,sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('calibrate_shock_ensemble.py','shock_group_pooling.py','shock_ensemble_estimator.py','shock_alignment.py',
                         'shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py',
                         'experiment_crossfit_gradient.py','validate_shock_estimator.py')})
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        specs=[(k,p,s,n) for k in ('logistic','gompertz') for p in (0.,.8) for s in (0.,.2) for n in a.noise_levels] if a.matrix else [(k,.8,.2,n) for k in ('logistic','gompertz') for n in a.noise_levels]
        if a.focused:specs=[(k,0.,s,n) for k in ('logistic','gompertz') for s,n in ((.2,max(a.noise_levels)),(0.,min(a.noise_levels)))]
        for kind,phi,spatial,noise in specs:
            rows=[]
            tasks=[((kind,phi,spatial,noise,i,a.seed_start,a.center_profile,a.width_band,a.interpolation_sensitivity,a.windows,a.subgroups),str(folder),protocol) for i in range(a.replicates)]
            for row in pool.map(run_cached,tasks):
                rows.append(row);print(kind,phi,spatial,noise,len(rows),'/',a.replicates,flush=True)
            cases.append(dict(shape=kind,noise=noise,phi=phi,spatial=spatial,**summarize(rows)))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    status='coverage_study_with_failures' if any(c['errors'] for c in cases) else 'coverage_study_only'
    a.output.write_text(json.dumps(dict(status=status,group_count=GROUPS,seed_start=a.seed_start,
        center_profile=a.center_profile,protocol=protocol,cases=cases),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
