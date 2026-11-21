"""Experimental full-procedure block resampling; coverage must be calibrated."""
from concurrent.futures import ProcessPoolExecutor
import math
import numpy as np
from shock_alignment import align, prepare_trials
from shock_crossfit_gradient import estimate_gradient
from shock_mean_profile import measure, radius


def point(trials, processing=False):
    data, diagnostics = align(trials, refine_passes=1)
    profile = data.mean((0, 1))
    gradient = estimate_gradient(data, diagnostics, width=1.5)
    metrics = np.array([measure(profile)[0, 0], gradient['gradient_thickness']])
    per_seed = np.array([[measure(run.mean(0))[0, 0], 1/fold['slope']]
                         for run, fold in zip(data, gradient['folds'])])
    result = dict(metrics=metrics, density=profile, seed_metrics=per_seed,
                  seed_density=data.mean(1))
    if processing:
        result['processing_sensitivity'] = np.array([0., max(abs(
            estimate_gradient(data, diagnostics, width=w)['gradient_thickness']-metrics[1])
            for w in (1., 2.))])
    return result


def resample_trials(trials, block, rng):
    # Keep original seed identities: a held-out seed must never reappear in training.
    if len(trials)!=3 or len({len(run) for run in trials})!=1:
        raise ValueError('Three equal-length seeds required')
    count=len(trials[0])
    if block<1 or count<block or count%block:
        raise ValueError('Block must divide the window count')
    result=[]
    for run in trials:
        starts=rng.integers(0,count,size=count//block)
        indices=((starts[:,None]+np.arange(block))%count).ravel()
        result.append([run[i] for i in indices])
    return result


def replicate(args):
    trials,block,seed=args
    estimate=point(resample_trials(trials,block,np.random.default_rng(seed)))
    return estimate['metrics'],estimate['density']


def estimate(trials,repeats=128,seed=20261109,jobs=1):
    if repeats<32 or jobs<1:raise ValueError('At least 32 repetitions and one worker required')
    trials=prepare_trials(trials)
    center=point(trials,processing=True)
    seed_u=4.303*center['seed_metrics'].std(0,ddof=1)/math.sqrt(3)
    seed_u+=abs(center['seed_metrics'].mean(0)-center['metrics'])
    density_seed_u=4.303*center['seed_density'].std(0,ddof=1)/math.sqrt(3)
    total=seed_u.copy();density_total=float(density_seed_u.max());density_pointwise=density_seed_u.copy();blocks={}
    # Templates, training shifts, held-out shifts and derivative peaks are rebuilt
    # from original windows in every replicate. No interpolation is compounded.
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for block in (8,16):
            tasks=((trials,block,seed+block*100000+i) for i in range(repeats))
            samples=list(pool.map(replicate,tasks)) if jobs>1 else [replicate(t) for t in tasks]
            metrics=np.asarray([s[0] for s in samples]);density=np.asarray([s[1] for s in samples])
            u=radius(metrics,center['metrics'])
            du=float(np.quantile(np.max(abs(density-center['density']),axis=1),.95))
            total=np.maximum(total,u);density_total=max(density_total,du)
            density_pointwise=np.maximum(density_pointwise,radius(density,center['density']))
            blocks[str(block)]=dict(metric_radius=u.tolist(),density_simultaneous_radius=du,
                                    bootstrap_metrics=metrics.tolist())
    result={k:dict(mean=float(center['metrics'][i]),statistical_radius=float(total[i]),
                   seed_t_radius=float(seed_u[i]),processing_sensitivity=float(center['processing_sensitivity'][i]))
            for i,k in enumerate(('width_10_90','gradient_thickness'))}
    result.update(status='experimental_coverage_unverified',
        resampling='Stratified circular blocks; original seed identities preserved; full estimator rebuilt',
        seed_envelope='Conservative seed-t heuristic; cross-fold dependence requires coverage calibration',
        blocks=blocks,seed_metrics=center['seed_metrics'].tolist(),density=center['density'].tolist(),
        density_simultaneous_radius=density_total,density_pointwise_radius=density_pointwise.tolist())
    return result


def main():
    import argparse,hashlib,json
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeats',type=int,default=128)
    p.add_argument('--jobs',type=int,default=3)
    p.add_argument('--seed',type=int,default=20261109)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    raw=a.profiles.read_bytes()
    result=estimate(json.loads(raw),repeats=a.repeats,jobs=a.jobs,seed=a.seed)
    result.update(source_name=a.profiles.name,source_sha256=hashlib.sha256(raw).hexdigest(),
                  repeats=a.repeats,rng_seed=a.seed,numpy_version=np.__version__)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
