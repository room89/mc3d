"""Coverage diagnostic for full-procedure intervals, using independent synthetic runs."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import hashlib
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID
from shock_refined_estimator import estimate
from validate_shock_estimator import shape,truth


def run(args):
    kind,phi,spatial,i,repeats,folder,protocol=args
    checkpoint=Path(folder)/f'{kind}-{phi}-{spatial}-{i}.json'
    if checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        if saved['protocol']!=protocol:raise ValueError('Checkpoint protocol mismatch')
        return saved['result']
    trials,_=fixture(kind,8400000+i,phi,spatial)
    r=estimate(trials,repeats=repeats,seed=9400000+i*1000,jobs=1)
    exact=truth(kind)
    center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    density=shape(GRID+center,kind)
    metric=np.array([r[k]['mean'] for k in ('width_10_90','gradient_thickness')])
    u=np.array([r[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')])
    proc=np.array([r[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    blocks=np.max([b['metric_radius'] for b in r['blocks'].values()],axis=0)
    density_error=float(np.max(abs(np.asarray(r['density'])-density)))
    row=dict(repetition=i,estimate=metric.tolist(),statistical_radius=u.tolist(),
                processing_sensitivity=proc.tolist(),relative_error=(metric/exact-1).tolist(),
                statistical_covers=(abs(metric-exact)<=u).tolist(),
                with_processing_covers=(abs(metric-exact)<=u+proc).tolist(),
                block_only_covers=(abs(metric-exact)<=blocks).tolist(),
                density_error=density_error,density_radius=r['density_simultaneous_radius'],
                density_covers=bool(density_error<=r['density_simultaneous_radius']))
    temporary=checkpoint.with_suffix('.tmp')
    temporary.write_text(json.dumps(dict(protocol=protocol,result=row),allow_nan=False)+'\n')
    temporary.replace(checkpoint)
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--replicates',type=int,default=32)
    p.add_argument('--bootstrap-repeats',type=int,default=128)
    p.add_argument('--jobs',type=int,default=3)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--matrix',action='store_true',help='All logistic/Gompertz correlation scenarios')
    a=p.parse_args()
    if a.replicates<2 or a.bootstrap_repeats<32 or a.jobs<1:p.error('Invalid replication counts')
    if a.output.exists():raise FileExistsError(a.output)
    specs=[('gompertz',.8,.2)] if not a.matrix else [(k,p,s) for k in ('logistic','gompertz') for p in (0.,.8) for s in (0.,.2)]
    folder=a.output.with_suffix('.parts')
    if folder.exists() and not a.resume:raise FileExistsError(folder)
    folder.mkdir(parents=True,exist_ok=True)
    protocol=dict(bootstrap_repeats=a.bootstrap_repeats,seed_start=8400000,
        sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
                 ('shock_alignment.py','shock_crossfit_gradient.py','shock_refined_estimator.py','shock_mean_profile.py','experiment_crossfit_gradient.py')})
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for kind,phi,spatial in specs:
            rows=[]
            for row in pool.map(run,[(kind,phi,spatial,i,a.bootstrap_repeats,str(folder),protocol) for i in range(a.replicates)]):
                rows.append(row)
                print(kind,phi,spatial,'replicate',row['repetition']+1,'/',a.replicates,flush=True)
            cases.append(dict(shape=kind,phi=phi,spatial=spatial,
                statistical_coverage=np.mean([r['statistical_covers'] for r in rows],axis=0).tolist(),
                with_processing_coverage=np.mean([r['with_processing_covers'] for r in rows],axis=0).tolist(),
                block_only_coverage=np.mean([r['block_only_covers'] for r in rows],axis=0).tolist(),
                density_coverage=float(np.mean([r['density_covers'] for r in rows])),runs=rows))
    result=dict(status='coverage_diagnostic_only',noise=.03,seed_start=8400000,
                replicates=a.replicates,bootstrap_repeats=a.bootstrap_repeats,protocol=protocol,cases=cases,
                limitation='Finite coverage counts; not proof of exact 95% coverage or full physical qualification')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
