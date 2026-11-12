"""Frozen, sequential synthetic screen; never qualifies physical acceptance."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import time
import numpy as np
from shock_profile_likelihood_fit import METHOD, KNOTS, design, fit

SOURCE_NAMES = ('shock_profile_likelihood_fit.py', 'run_profile_likelihood_screen.py')
BASE_SEED = 202611120000


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    # Exclusivity prevents accidentally replacing evidence from a previous run.
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def protocol():
    specs = [dict(shape=s, amplitude=a, noise=n, repetition=r, analytic=False)
             for s, a, n, r in itertools.product(('logistic','gompertz'),(0.,.05,.5),('low','high'),range(16))]
    specs += [dict(shape=s, amplitude=.05, noise=n, repetition=r, analytic=True)
              for s, n, r in itertools.product(('logistic','gompertz'),('low','high'),range(4))]
    return dict(method=METHOD, physical_acceptance='incomplete', seed=BASE_SEED,
                independent_runs=32, windows=2, grid=[-10,10,41], knots=KNOTS.tolist(),
                temporal_phi=.8, spatial_rho=.5, run_sigma={'low':.01,'high':.02},
                window_sigma={'low':.02,'high':.04}, max_iterations=200,
                methods=['marginal','hard','oracle'], specifications=specs,
                sources={n:digest(Path(__file__).with_name(n)) for n in SOURCE_NAMES},
                numpy_version=np.__version__,
                decision='On matched models, require lower mean ISE for both shapes at nonzero motion, no additional failures; uncertainty crossing zero is inconclusive.',
                limitations=['Known covariance and shift law; synthetic blocks only.',
                             'No coverage qualification, fitted noise or physical acceptance.',
                             'Analytic controls are outside the piecewise-linear model and excluded from the matched-model decision.'])


def curve(x, kind):
    return 1/(1+np.exp(-x)) if kind=='logistic' else np.exp(-np.exp(-x))


def generate(spec, index, p):
    x=np.linspace(*p['grid']); times=p['windows'];runs=p['independent_runs']
    a=spec['amplitude']; values=[0.] if a==0 else [-a,0.,a]
    configs=np.array(list(itertools.product(values,repeat=times)))
    streams=np.random.SeedSequence([p['seed'],index]).spawn(2)
    shifts_rng,noise_rng=[np.random.default_rng(s) for s in streams]
    labels=shifts_rng.integers(len(configs),size=runs)
    lag=np.abs(np.arange(times)[:,None]-np.arange(times))
    temporal=p['run_sigma'][spec['noise']]**2+p['window_sigma'][spec['noise']]**2*p['temporal_phi']**lag
    scale=1+2*(x-x[0])/(x[-1]-x[0])
    spatial=scale[:,None]*scale[None,:]*p['spatial_rho']**np.abs(np.arange(len(x))[:,None]-np.arange(len(x)))
    truth=curve(KNOTS,spec['shape']);truth[[0,-1]]=[0,1]
    predictions=np.array([[curve(x-s,spec['shape']) if spec['analytic'] else design(x,KNOTS,s)@truth for s in row] for row in configs])
    factor=np.linalg.cholesky(np.kron(temporal,spatial))
    y=predictions[labels]+(noise_rng.normal(size=(runs,times*len(x)))@factor.T).reshape(runs,times,len(x))
    return x,configs,labels,temporal,spatial,y,truth


def trial(spec, index, p):
    x,configs,labels,temporal,spatial,y,truth=generate(spec,index,p)
    grid=np.linspace(-8,8,321)
    exact=curve(grid,spec['shape']) if spec['analytic'] else np.interp(grid,KNOTS,truth)
    results={};started=time.monotonic()
    for method in p['methods']:
        weights=np.eye(len(configs))[labels] if method=='oracle' else np.ones(len(configs))
        try:
            result=fit(y,x,temporal,spatial,configs,weights,hard=method=='hard',max_iterations=p['max_iterations'])
            error=np.interp(grid,KNOTS,result['profile'])-exact
            result.update(signed_error=error.tolist(),maximum_absolute_error=float(abs(error).max()),
                          integrated_squared_error=float(np.trapezoid(error**2,grid)))
        except ValueError as error:
            result=dict(status='rejected',error=str(error),physical_acceptance='incomplete')
        results[method]=result
    return dict(index=index,specification=spec,results=results,evaluation_grid=grid.tolist(),
                truth=exact.tolist(),true_configuration_indices=labels.tolist(),
                elapsed_seconds=time.monotonic()-started,physical_acceptance='incomplete')


def validate_checkpoint(row, index, spec, plan_hash, p):
    if (row.get('plan_sha256') != plan_hash or row.get('index') != index or
            row.get('specification') != spec or set(row.get('results', {})) != set(p['methods'])):
        raise ValueError('Checkpoint identity mismatch')
    grid = np.linspace(-8, 8, 321)
    truth = curve(KNOTS, spec['shape']); truth[[0, -1]] = [0, 1]
    expected = curve(grid, spec['shape']) if spec['analytic'] else np.interp(grid, KNOTS, truth)
    if not np.array_equal(row['evaluation_grid'], grid) or not np.array_equal(row['truth'], expected):
        raise ValueError('Checkpoint evaluation grid or truth changed')
    for result in row['results'].values():
        if result['status'] == 'rejected':
            if not result.get('error'): raise ValueError('Missing rejection reason')
            continue
        if result['status'] not in ('converged', 'ambiguous', 'iteration_limit', 'inner_iteration_limit', 'objective_increased'):
            raise ValueError('Unknown fit status')
        profile = np.asarray(result['profile'], float)
        if profile.shape != KNOTS.shape or not np.isfinite(profile).all() or np.any(np.diff(profile) < 0) or profile[0] != 0 or profile[-1] != 1:
            raise ValueError('Invalid checkpoint profile')
        error = np.interp(grid, KNOTS, profile)-expected
        if (not np.array_equal(result['signed_error'], error) or
                result['maximum_absolute_error'] != float(abs(error).max()) or
                result['integrated_squared_error'] != float(np.trapezoid(error**2, grid))):
            raise ValueError('Checkpoint metric mismatch')
        if len(result['starts']) != 3: raise ValueError('Missing start diagnostics')
    return row


def summarize(rows,p):
    matched=[r for r in rows if not r['specification']['analytic']]
    failures={m:sum(r['results'][m]['status']!='converged' for r in matched) for m in p['methods']}
    families={}
    for shape in ('logistic','gompertz'):
        selected=[r for r in matched if r['specification']['shape']==shape and r['specification']['amplitude']>0]
        pairs=[r for r in selected if all(r['results'][m]['status']=='converged' for m in ('marginal','hard'))]
        groups={}
        for r in pairs:
            s=r['specification'];key=(s['amplitude'],s['noise'])
            groups.setdefault(key,[]).append(r['results']['marginal']['integrated_squared_error']-r['results']['hard']['integrated_squared_error'])
        complete=len(pairs)==len(selected) and len(selected)==64
        mean=float(np.mean([v for values in groups.values() for v in values])) if pairs else None
        se=float(np.sqrt(sum(np.var(v,ddof=1)/len(v) for v in groups.values()))/4) if complete else None
        families[shape]=dict(planned=64,available=len(selected),paired_successes=len(pairs),
                             mean_ise_difference=mean,stratified_mcse=se,
                             descriptive_normal_interval=[mean-1.96*se,mean+1.96*se] if se is not None else None)
    complete=len(rows)==len(p['specifications'])
    decision='incomplete'
    if complete:
        if failures['marginal']>failures['hard']:decision='do_not_advance_additional_failures'
        elif any(v['paired_successes']!=64 for v in families.values()):decision='inconclusive_failed_pairs'
        elif any(v['mean_ise_difference']>=0 for v in families.values()):decision='no_demonstrated_improvement'
        elif any(v['descriptive_normal_interval'][1]>=0 for v in families.values()):decision='inconclusive_monte_carlo_uncertainty'
        else:decision='advance_to_noise_model_research_only'
    return dict(status='complete' if complete else 'partial',physical_acceptance='incomplete',
                trials=len(rows),matched_trials=len(matched),failures=failures,families=families,decision=decision,
                analytic_control_trials=sum(r['specification']['analytic'] for r in rows),
                limitations=p['limitations']+['Paired summaries require both fits to converge; failures are counted separately.',
                                             'Intervals describe Monte Carlo uncertainty of this synthetic comparison, not profile confidence bands.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--max-new-trials',type=int)
    args=parser.parse_args()
    if args.max_new_trials is not None and args.max_new_trials<1:parser.error('Positive trial limit required')
    p=protocol();folder=args.directory;plan=folder/'plan.json'
    if folder.exists():
        if not args.resume or not plan.exists() or json.loads(plan.read_text())!=p:
            raise ValueError('Existing directory requires identical plan and --resume')
    else:
        folder.mkdir(parents=True);save(plan,p)
    plan_hash=digest(plan)
    if args.prepare_only:print('Prepared',len(p['specifications']),'trials',flush=True);return
    # One process only. Failure to lower priority must not silently run full priority.
    if hasattr(os,'nice'):os.nice(19)
    rows=[];new=0
    for i,spec in enumerate(p['specifications']):
        path=folder/f'trial-{i:03d}.json'
        if path.exists():
            row=json.loads(path.read_text())
            validate_checkpoint(row,i,spec,plan_hash,p)
        else:
            if args.max_new_trials is not None and new>=args.max_new_trials:break
            if protocol()!=p or digest(plan)!=plan_hash:raise ValueError('Sources or plan changed')
            row=trial(spec,i,p);row['plan_sha256']=plan_hash
            validate_checkpoint(row,i,spec,plan_hash,p);save(path,row);new+=1
        rows.append(row)
        print('completed',len(rows),'/',len(p['specifications']),{m:r['status'] for m,r in row['results'].items()},flush=True)
    if protocol()!=p or digest(plan)!=plan_hash:raise ValueError('Sources or plan changed during screen')
    report=summarize(rows,p);report['plan_sha256']=plan_hash
    report['checkpoints']={f'trial-{r["index"]:03d}.json':digest(folder/f'trial-{r["index"]:03d}.json') for r in rows}
    target=folder/('complete.json' if report['status']=='complete' else f'partial-{len(rows):03d}.json')
    if not target.exists():save(target,report)
    print(json.dumps(report['families']),flush=True)


if __name__=='__main__':main()
