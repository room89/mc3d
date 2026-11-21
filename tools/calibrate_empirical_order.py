"""Conditional empirical stress screen for the frozen alignment-order candidate.

The historical banks include front motion and run offsets. This screen does not
create independent physical runs or qualify the physical noise law.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import math
from pathlib import Path
import numpy as np
from calibrate_order_sensitivity import SOURCES as ORDER_SOURCES
from calibrate_shock_ensemble import summarize
from shock_empirical_fluctuations import load_bank,fixture
from shock_order_sensitivity import METHOD,subgroup_point,pooled_point,estimate
from shock_alignment import GRID
from shock_validation import digest,save
from validate_shock_estimator import truth,shape

SOURCES=tuple(ORDER_SOURCES)+('shock_empirical_fluctuations.py','calibrate_empirical_order.py')
BANK_NAMES=('mc3d-training','mc3d-heldout','sparta-training','sparta-heldout')
_BANK_CACHE={}


def verify_banks(manifest_path,verification_path):
    manifest_path=Path(manifest_path);verification_path=Path(verification_path)
    manifest_hash=digest(manifest_path);proof_hash=digest(verification_path)
    manifest=json.loads(manifest_path.read_text());proof=json.loads(verification_path.read_text())
    if (proof.get('status')!='conditional_empirical_banks_verified' or proof.get('physical_acceptance')!='incomplete' or
        proof.get('manifest_sha256')!=manifest_hash or proof.get('region_phase_checks')!=648 or
        proof.get('loader_sha256')!=digest(Path(__file__).with_name('shock_empirical_fluctuations.py'))):
        raise ValueError('Matching completed physical-bank verification required')
    if set(proof['banks'])!=set(BANK_NAMES) or set(manifest['banks'])!=set(BANK_NAMES):
        raise ValueError('All four verified banks required')
    banks={}
    for name in BANK_NAMES:
        bank=load_bank(manifest_path,name);entry=manifest['banks'][name];checked=proof['banks'][name]
        if (checked['sha256']!=entry['sha256'] or checked['shape']!=entry['shape'] or
            checked['sources']!=entry['sources'] or bank['bank_sha256']!=checked['sha256']):
            raise ValueError('Bank identity differs from verified evidence')
        banks[name]=dict(sha256=bank['bank_sha256'],shape=list(bank['fluctuations'].shape),sources=bank['sources'])
    if digest(manifest_path)!=manifest_hash or digest(verification_path)!=proof_hash:
        raise ValueError('Bank evidence changed during preflight')
    return dict(manifest_sha256=manifest_hash,verification_sha256=proof_hash,banks=banks)


def trial(spec,protocol,manifest_path):
    bank_name,kind,amplitude,i=spec;key=(str(manifest_path),bank_name,protocol['bank_evidence']['manifest_sha256'])
    if key not in _BANK_CACHE:
        bank=load_bank(manifest_path,bank_name)
        if (bank['manifest_sha256']!=key[2] or bank['bank_sha256']!=protocol['bank_evidence']['banks'][bank_name]['sha256']):
            raise ValueError('Bank changed before trial')
        _BANK_CACHE[key]=bank
    bank=_BANK_CACHE[key];groups=[];start=protocol['seed_start'];windows=protocol['windows']
    for g in range(12):
        points=[]
        for j in range(6):
            runs,_=fixture(kind,start+i*100+g*6+j,bank,windows,amplitude)
            points.append(subgroup_point(runs))
        groups.append(pooled_point(points))
    result=estimate(groups,repeats=1024,seed=start+1000000+i)
    exact=truth(kind);center=-math.log(math.log(2)) if kind=='gompertz' else 0.
    values=np.array([result[k]['mean'] for k in ('width_10_90','gradient_thickness')]);u=np.array([result[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')]);proc=np.array([result[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
    error=float(np.max(abs(np.asarray(result['density'])-shape(GRID+center,kind))))
    return dict(repetition=i,estimate=values.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),relative_error=(values/exact-1).tolist(),statistical_covers=(abs(values-exact)<=u).tolist(),with_processing_covers=(abs(values-exact)<=u+proc).tolist(),density_error=error,density_radius=result['density_simultaneous_radius'],density_covers=bool(error<=result['density_simultaneous_radius']),density_statistical_radius=result['density_statistical_radius'],density_processing_sensitivity=result['density_processing_sensitivity'],order_density_difference=result['alignment_order_density_difference'],order_metric_differences=result['alignment_order_metric_differences'])


def cached(task):
    spec,folder,protocol,manifest_path=task;path=Path(folder)/('-'.join(map(str,spec))+'.json')
    if path.exists():
        saved=json.loads(path.read_text())
        if saved['protocol']!=protocol or saved['spec']!=list(spec):raise ValueError('Checkpoint protocol mismatch')
        return spec,saved['result']
    try:result=trial(spec,protocol,manifest_path)
    except ValueError as error:result=dict(repetition=spec[-1],error=str(error))
    save(path,dict(protocol=protocol,spec=spec,result=result));return spec,result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--verified',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--seed-start',type=int,required=True)
    parser.add_argument('--windows',type=int,choices=(32,64),default=64);parser.add_argument('--replicates',type=int,default=8)
    parser.add_argument('--jobs',type=int,default=1);parser.add_argument('--prepare-only',action='store_true');parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.replicates<2 or args.replicates*100>=1000000 or args.jobs<1 or args.seed_start<=0:parser.error('Invalid seed or replication counts')
    if args.output.exists():raise FileExistsError(args.output)
    evidence=verify_banks(args.manifest,args.verified)
    protocol=dict(method=METHOD,status='conditional_empirical_stress',group_count=12,subgroups=6,bootstrap_repeats=1024,windows=args.windows,seed_start=args.seed_start,replicates=args.replicates,motion_amplitudes=[0.,.05,1.],bank_evidence=evidence,sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES},limitations=['Conditional screen on finite historical banks, not new independent physical realizations.','Historical fluctuations include front motion and run offsets; separate bank centering does not test physical mean agreement.','Cyclic windows preserve native order with one wrap boundary; no per-run mean subtraction.','Complete independent parametric calibration and actual physical precision/convergence remain required.'])
    folder=args.output.with_suffix('.parts');plan_path=args.output.with_suffix('.plan.json')
    if plan_path.exists():
        if not args.resume or json.loads(plan_path.read_text())['protocol']!=protocol:raise ValueError('Prepared empirical protocol changed or --resume missing')
    folder.mkdir(parents=True,exist_ok=args.resume)
    keys=[(b,k,a) for b in BANK_NAMES for k in ('logistic','gompertz') for a in protocol['motion_amplitudes']]
    plan=dict(status='empirical_screen_prepared',physical_acceptance='incomplete',protocol=protocol,scheduled_trials=len(keys)*args.replicates)
    if not plan_path.exists():save(plan_path,plan)
    if args.prepare_only:
        print('Prepared',plan['scheduled_trials'],'conditional empirical trials; none executed.');return
    rows={key:[] for key in keys};tasks=[(key+(i,),str(folder),protocol,str(args.manifest.resolve())) for key in keys for i in range(args.replicates)]
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(cached,t) for t in tasks]
        for number,future in enumerate(as_completed(futures),1):
            spec,result=future.result();rows[tuple(spec[:-1])].append(result)
            print('completed',number,'of',len(tasks),'errors',sum('error' in r for values in rows.values() for r in values),flush=True)
    if verify_banks(args.manifest,args.verified)!=evidence or any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()):raise ValueError('Empirical screen input changed')
    cases=[dict(bank=k[0],shape=k[1],motion_amplitude=k[2],**summarize(sorted(rows[k],key=lambda r:r['repetition']))) for k in keys]
    save(args.output,dict(status='empirical_stress_with_failures' if any(c['errors'] for c in cases) else 'empirical_stress_only',physical_acceptance='incomplete',protocol=protocol,cases=cases))

if __name__=='__main__':main()
