"""Pilot of independent-group inference; no production simulations."""
import argparse,hashlib,json,math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_ensemble_estimator import group_point,estimate,GROUPS
from validate_shock_estimator import truth,shape


def run(spec):
    kind,phi,spatial,noise,i,seed_start,center_profile=spec
    groups=[group_point(fixture(kind,seed_start+i*100+g,phi,spatial,noise=noise)[0],center_profile=center_profile) for g in range(GROUPS)]
    r=estimate(groups,repeats=1024,seed=seed_start+1000000+i,center_profile=center_profile)
    exact=truth(kind);center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    values=np.array([r[k]['mean'] for k in ('width_10_90','gradient_thickness')])
    u=np.array([r[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')])
    proc=np.array([r[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    error=float(np.max(abs(np.array(r['density'])-shape(GRID+center,kind))))
    return dict(repetition=i,estimate=values.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),
                relative_error=(values/exact-1).tolist(),statistical_covers=(abs(values-exact)<=u).tolist(),
                with_processing_covers=(abs(values-exact)<=u+proc).tolist(),density_error=error,
                density_radius=r['density_simultaneous_radius'],density_covers=bool(error<=r['density_simultaneous_radius']))


def run_cached(task):
    spec,folder,protocol=task
    kind,phi,spatial,noise,i,_,_=spec
    checkpoint=Path(folder)/f'{kind}-{phi}-{spatial}-{noise}-{i}.json'
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        if saved['protocol']!=protocol or saved['spec']!=list(spec):
            raise ValueError('Checkpoint protocol mismatch')
        return saved['result']
    row=run(spec)
    temporary=checkpoint.with_suffix('.tmp')
    temporary.write_text(json.dumps(dict(protocol=protocol,spec=spec,result=row),allow_nan=False)+'\n')
    temporary.replace(checkpoint)
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--replicates',type=int,default=32);p.add_argument('--jobs',type=int,default=3)
    p.add_argument('--seed-start',type=int,default=12000000)
    p.add_argument('--center-profile',action='store_true')
    p.add_argument('--matrix',action='store_true')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if a.replicates<2 or a.jobs<1:p.error('Invalid replication counts')
    if a.output.exists():raise FileExistsError(a.output)
    folder=a.output.with_suffix('.parts')
    if folder.exists() and not a.resume:raise FileExistsError(folder)
    folder.mkdir(parents=True,exist_ok=True)
    protocol=dict(group_count=GROUPS,bootstrap_repeats=1024,seed_start=a.seed_start,
        center_profile=a.center_profile,sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('calibrate_shock_ensemble.py','shock_ensemble_estimator.py','shock_alignment.py',
                         'shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py',
                         'experiment_crossfit_gradient.py','validate_shock_estimator.py')})
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        specs=[(k,p,s,n) for k in ('logistic','gompertz') for p in (0.,.8) for s in (0.,.2) for n in (.015,.03)] if a.matrix else [(k,.8,.2,.03) for k in ('logistic','gompertz')]
        for kind,phi,spatial,noise in specs:
            rows=[]
            tasks=[((kind,phi,spatial,noise,i,a.seed_start,a.center_profile),str(folder),protocol) for i in range(a.replicates)]
            for row in pool.map(run_cached,tasks):
                rows.append(row);print(kind,phi,spatial,noise,len(rows),'/',a.replicates,flush=True)
            cases.append(dict(shape=kind,noise=noise,phi=phi,spatial=spatial,
                statistical_coverage=np.mean([r['statistical_covers'] for r in rows],axis=0).tolist(),
                with_processing_coverage=np.mean([r['with_processing_covers'] for r in rows],axis=0).tolist(),
                density_coverage=float(np.mean([r['density_covers'] for r in rows])),runs=rows))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(status='coverage_study_only',group_count=GROUPS,seed_start=a.seed_start,
        center_profile=a.center_profile,protocol=protocol,cases=cases),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
