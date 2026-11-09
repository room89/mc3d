"""Pilot of independent-group inference; no production simulations."""
import argparse,json,math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_ensemble_estimator import group_point,estimate,GROUPS
from validate_shock_estimator import truth,shape


def run(spec):
    kind,i=spec
    groups=[group_point(fixture(kind,12000000+i*100+g,.8,.2)[0]) for g in range(GROUPS)]
    r=estimate(groups,repeats=1024,seed=13000000+i)
    exact=truth(kind);center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    values=np.array([r[k]['mean'] for k in ('width_10_90','gradient_thickness')])
    u=np.array([r[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')])
    proc=np.array([r[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    error=float(np.max(abs(np.array(r['density'])-shape(GRID+center,kind))))
    return dict(repetition=i,estimate=values.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),
                relative_error=(values/exact-1).tolist(),statistical_covers=(abs(values-exact)<=u).tolist(),
                with_processing_covers=(abs(values-exact)<=u+proc).tolist(),density_error=error,
                density_radius=r['density_simultaneous_radius'],density_covers=bool(error<=r['density_simultaneous_radius']))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--replicates',type=int,default=32);p.add_argument('--jobs',type=int,default=3);a=p.parse_args()
    if a.replicates<2 or a.jobs<1:p.error('Invalid replication counts')
    if a.output.exists():raise FileExistsError(a.output)
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for kind in ('logistic','gompertz'):
            rows=[]
            for row in pool.map(run,[(kind,i) for i in range(a.replicates)]):
                rows.append(row);print(kind,len(rows),'/',a.replicates,flush=True)
            cases.append(dict(shape=kind,noise=.03,phi=.8,spatial=.2,
                statistical_coverage=np.mean([r['statistical_covers'] for r in rows],axis=0).tolist(),
                with_processing_coverage=np.mean([r['with_processing_covers'] for r in rows],axis=0).tolist(),
                density_coverage=float(np.mean([r['density_covers'] for r in rows])),runs=rows))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(dict(status='pilot_only',group_count=GROUPS,seed_start=12000000,cases=cases),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
