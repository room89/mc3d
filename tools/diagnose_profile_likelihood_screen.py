"""Read-only post-hoc diagnosis of a complete frozen synthetic screen.

Scenario intervals are descriptive, unadjusted for multiple comparisons. No
new trials, fitted corrections, threshold changes or qualification decisions.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from run_profile_likelihood_screen import digest, validate_checkpoint, summarize


def paired(values):
    values = np.asarray(values, float)
    mean = float(values.mean())
    se = float(values.std(ddof=1)/np.sqrt(len(values)))
    return dict(paired_differences=values.tolist(), mean=mean, mcse=se, descriptive_interval=[mean-1.96*se, mean+1.96*se])


def diagnose(folder):
    folder = Path(folder); plan_path = folder/'plan.json'; report_path = folder/'complete.json'
    plan_hash, report_hash = digest(plan_path), digest(report_path)
    p, report = json.loads(plan_path.read_text()), json.loads(report_path.read_text())
    if report.get('status') != 'complete' or report.get('plan_sha256') != plan_hash:
        raise ValueError('Matching completed report required')
    hashes = report['checkpoints']; rows = []
    if set(hashes) != {f'trial-{i:03d}.json' for i in range(len(p['specifications']))}:
        raise ValueError('Incomplete checkpoint inventory')
    for i, spec in enumerate(p['specifications']):
        path = folder/f'trial-{i:03d}.json'
        if digest(path) != hashes[path.name]: raise ValueError('Checkpoint hash mismatch')
        row = json.loads(path.read_text()); validate_checkpoint(row,i,spec,plan_hash,p); rows.append(row)
    expected = summarize(rows,p)
    if any(report[k] != v for k,v in expected.items()): raise ValueError('Summary mismatch')
    scenarios = []
    for analytic, shape, amplitude, noise in itertools.product((False,True),('logistic','gompertz'),(0.,.05,.5),('low','high')):
        selected = [r for r in rows if (r['specification']['analytic'],r['specification']['shape'],r['specification']['amplitude'],r['specification']['noise']) == (analytic,shape,amplitude,noise)]
        if not selected: continue
        if any(z['status'] != 'converged' for r in selected for z in r['results'].values()):
            raise ValueError('This diagnostic requires all methods converged; failures cannot be silently omitted')
        grid = np.array(selected[0]['evaluation_grid']); methods = {}
        for method in p['methods']:
            results = [r['results'][method] for r in selected]
            errors = np.array([z['signed_error'] for z in results]); mean_error = errors.mean(axis=0)
            bias_energy = float(np.trapezoid(mean_error**2,grid))
            scatter_energy = float(np.trapezoid(np.mean((errors-mean_error)**2,axis=0),grid))
            ise = float(np.mean([z['integrated_squared_error'] for z in results]))
            if not np.isclose(ise,bias_energy+scatter_energy,rtol=1e-12,atol=1e-15):
                raise ValueError('Empirical error decomposition mismatch')
            starts=[]
            for r,z in zip(selected,results):
                converged=[s for s in z['starts'] if s['status']=='converged']
                profiles=np.array([s['profile'] for s in converged]); objectives=np.array([s['objective'] for s in converged])
                spread=float(np.max(np.ptp(profiles,axis=0))) if len(profiles) else None
                gap=float(np.ptp(objectives)) if len(objectives) else None
                starts.append(dict(index=r['index'],converged_starts=len(converged),
                                   objective_gap=gap,max_profile_spread=spread))
            methods[method]=dict(mean_ise=ise,empirical_mean_error_energy=bias_energy,
                                empirical_scatter_energy=scatter_energy,
                                mean_signed_error=mean_error.tolist(),
                                mean_error_peak_x=float(grid[np.argmax(abs(mean_error))]),
                                starts=starts,max_start_profile_spread=max(s['max_profile_spread'] for s in starts),
                                max_start_objective_gap=max(s['objective_gap'] for s in starts))
        comparisons={}
        for a,b in [('marginal','hard'),('marginal','oracle'),('hard','oracle')]:
            comparisons[a+'_minus_'+b]=paired([r['results'][a]['integrated_squared_error']-r['results'][b]['integrated_squared_error'] for r in selected])
        scenarios.append(dict(analytic=analytic,shape=shape,amplitude=amplitude,noise=noise,
                              trials=len(selected),trial_indices=[r['index'] for r in selected],methods=methods,comparisons=comparisons))
    if digest(plan_path)!=plan_hash or digest(report_path)!=report_hash or any(digest(folder/name)!=h for name,h in hashes.items()):
        raise ValueError('Evidence changed during diagnosis')
    return dict(status='posthoc_diagnostic_only',physical_acceptance='incomplete',
                plan_sha256=plan_hash,report_sha256=report_hash,script_sha256=digest(Path(__file__)),
                original_decision=report['decision'],scenarios=scenarios,
                limitations=['No new independent evidence or changed screen decision.',
                             'Scenario intervals are descriptive and not multiplicity adjusted.',
                             'Empirical mean-error energy includes finite-sample variance: it is not an unbiased estimate of squared bias.',
                             'Agreement of three starts does not prove global optimality.',
                             'Oracle differences measure loss with unknown shifts in this specified synthetic model, not a physical error source.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=diagnose(args.directory)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    for s in result['scenarios']:
        print(s['analytic'],s['shape'],s['amplitude'],s['noise'],
              {m:round(v['mean_ise'],8) for m,v in s['methods'].items()},
              s['comparisons']['marginal_minus_hard'])


if __name__=='__main__':main()
