"""Independent, fixed-size confirmation of the frozen profile estimator.

The pilot sizes this experiment only. Confirmation observations never enter the
pilot summaries, estimator settings or sample-size calculation. No early success
stopping or automatic sample-size extension is permitted.
"""
import argparse
import itertools
import json
import math
import os
from pathlib import Path
from statistics import NormalDist
import numpy as np
from diagnose_profile_likelihood_screen import diagnose
from run_profile_likelihood_screen import protocol, trial, validate_checkpoint, digest, save

SEED = 202611140000
SOURCES = ('shock_profile_likelihood_fit.py', 'run_profile_likelihood_screen.py',
           'diagnose_profile_likelihood_screen.py', 'confirm_profile_likelihood.py')


def prepare(pilot):
    proof = diagnose(pilot)
    old = json.loads((Path(pilot)/'plan.json').read_text())
    if old['seed'] == SEED: raise ValueError('Confirmation must use a different seed family')
    p = protocol()
    if any(old['sources'][name] != p['sources'][name] for name in p['sources']):
        raise ValueError('Estimator or generator differs from the pilot')
    z = NormalDist().inv_cdf(1-.05/(2*2))  # two families, two-sided Bonferroni
    sizes = {}; planning = {}
    for shape in ('logistic','gompertz'):
        cells = [s for s in proof['scenarios'] if not s['analytic'] and s['shape']==shape and s['amplitude']>0]
        if len(cells)!=4 or any(s['trials']!=16 for s in cells):
            raise ValueError('Four complete 16-trial pilot strata required per shape')
        reference = float(np.mean([s['methods']['hard']['mean_ise'] for s in cells]))
        target = .025*reference  # precision planning, not a physical tolerance
        variance_sum = sum(float(np.var(s['comparisons']['marginal_minus_hard']['paired_differences'],ddof=1)) for s in cells)
        required = math.ceil(z*z*variance_sum/(16*target*target))
        count = max(32,16*math.ceil(required/16))
        sizes[shape] = count
        planning[shape] = dict(pilot_hard_mean_ise=reference,target_halfwidth=target,
                               pilot_stratum_variance_sum=variance_sum,required_per_stratum=required,
                               fixed_per_stratum=count)
    p.update(seed=SEED,stage='independent_confirmation',pilot_plan_sha256=proof['plan_sha256'],
             pilot_report_sha256=proof['report_sha256'],planning=planning,critical_value=z,
             sources={name:digest(Path(__file__).with_name(name)) for name in SOURCES},
             specifications=[dict(shape=s,amplitude=a,noise=n,repetition=r,analytic=False)
                             for s in sizes for a,n in itertools.product((.05,.5),('low','high')) for r in range(sizes[s])],
             decision='Only after all fixed trials: complete converged pairs, no failures, both precision targets met and both adjusted interval upper endpoints below zero.',
             limitations=['Known synthetic covariance/shift law; not density-band coverage or physical qualification.',
                          'Normal Monte Carlo intervals with two-family Bonferroni adjustment are asymptotic, not exact.',
                          'Pilot precision planning is approximate; missed precision yields inconclusive, never automatic extension.',
                          'No new estimator fitting choices, exclusions, optional stopping or pilot/confirmation pooling.'])
    return p


def summarize_confirmation(rows,p):
    failures={m:sum(r['results'][m]['status']!='converged' for r in rows) for m in p['methods']}
    families={}
    for shape,planning in p['planning'].items():
        n=planning['fixed_per_stratum'];selected=[r for r in rows if r['specification']['shape']==shape]
        groups=[];metrics={m:[] for m in p['methods']};successes=0
        for amplitude,noise in itertools.product((.05,.5),('low','high')):
            cell=[r for r in selected if r['specification']['amplitude']==amplitude and r['specification']['noise']==noise]
            pairs=[r for r in cell if all(r['results'][m]['status']=='converged' for m in p['methods'])]
            successes+=len(pairs)
            differences=[r['results']['marginal']['integrated_squared_error']-r['results']['hard']['integrated_squared_error'] for r in pairs]
            groups.append(differences)
            for m in metrics:metrics[m].extend(r['results'][m]['integrated_squared_error'] for r in pairs)
        complete=len(selected)==4*n and all(len(v)==n for v in groups)
        mean=se=interval=None;precise=False
        if complete:
            mean=float(np.mean([np.mean(v) for v in groups]))
            se=float(np.sqrt(sum(np.var(v,ddof=1)/n for v in groups))/4)
            radius=p['critical_value']*se;interval=[mean-radius,mean+radius]
            precise=bool(radius<=planning['target_halfwidth'])
        families[shape]=dict(planned=4*n,available=len(selected),complete_pairs=successes,
                             mean_ise={m:float(np.mean(v)) if v else None for m,v in metrics.items()},
                             mean_difference=mean,mcse=se,adjusted_interval=interval,
                             precision_target_met=precise,target_halfwidth=planning['target_halfwidth'])
    complete=len(rows)==len(p['specifications']);decision='incomplete'
    if complete:
        if any(failures.values()):decision='do_not_advance_fit_failures'
        elif not all(v['precision_target_met'] for v in families.values()):decision='inconclusive_precision'
        elif any(v['adjusted_interval'][0]>=0 for v in families.values()):decision='no_advantage_for_all_families'
        elif any(v['adjusted_interval'][1]>=0 for v in families.values()):decision='inconclusive_advantage'
        else:decision='advance_to_noise_model_research_only'
    return dict(status='complete' if complete else 'partial',physical_acceptance='incomplete',
                trials=len(rows),planned=len(p['specifications']),failures=failures,
                families=families,decision=decision,limitations=p['limitations'])


def read_rows(folder,p):
    ph=digest(folder/'plan.json');rows=[];hashes={}
    allowed={f'trial-{i:04d}.json' for i in range(len(p['specifications']))}
    if any(path.name not in allowed for path in folder.glob('trial-*.json')):
        raise ValueError('Unexpected checkpoint')
    for i,spec in enumerate(p['specifications']):
        path=folder/f'trial-{i:04d}.json'
        if not path.exists():continue
        h=digest(path);row=json.loads(path.read_text());validate_checkpoint(row,i,spec,ph,p)
        rows.append(row);hashes[path.name]=h
    if digest(folder/'plan.json')!=ph or any(digest(folder/name)!=h for name,h in hashes.items()):
        raise ValueError('Evidence changed during read-back')
    return rows,hashes


def check_sources(folder,p):
    if json.loads((folder/'plan.json').read_text())!=p or any(digest(Path(__file__).with_name(n))!=h for n,h in p['sources'].items()):
        raise ValueError('Frozen plan or sources changed')
    if p['numpy_version']!=np.__version__:raise ValueError('NumPy version changed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--pilot',type=Path)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--audit-only',action='store_true')
    parser.add_argument('--max-new-trials',type=int)
    args=parser.parse_args();folder=args.directory;plan=folder/'plan.json'
    if args.max_new_trials is not None and args.max_new_trials<1:parser.error('Positive limit required')
    if plan.exists():
        if not(args.resume or args.audit_only):raise ValueError('Existing plan requires --resume or --audit-only')
        p=json.loads(plan.read_text())
    else:
        if args.pilot is None or args.resume or args.audit_only:parser.error('New experiment requires --pilot')
        p=prepare(args.pilot);folder.mkdir(parents=True,exist_ok=False);save(plan,p)
    check_sources(folder,p);ph=digest(plan)
    if args.prepare_only:print('Frozen',len(p['specifications']),'trials',p['planning'],flush=True);return
    rows,hashes=read_rows(folder,p)
    if args.audit_only:
        report=json.loads((folder/'complete.json').read_text());expected=summarize_confirmation(rows,p)
        if report.get('plan_sha256')!=ph or report.get('checkpoints')!=hashes or any(report[k]!=v for k,v in expected.items()):
            raise ValueError('Completed report mismatch')
        if expected['status']!='complete':raise ValueError('Incomplete evidence')
        check_sources(folder,p)
        print('Verified',len(rows),'trials:',expected['decision'],flush=True);return
    os.nice(19);completed={r['index'] for r in rows};new=0
    for i,spec in enumerate(p['specifications']):
        if i in completed:continue
        if args.max_new_trials is not None and new>=args.max_new_trials:break
        check_sources(folder,p)
        row=trial(spec,i,p);row['plan_sha256']=ph;validate_checkpoint(row,i,spec,ph,p)
        save(folder/f'trial-{i:04d}.json',row);rows.append(row);new+=1
        print('completed',len(rows),'/',len(p['specifications']),{m:z['status'] for m,z in row['results'].items()},flush=True)
    check_sources(folder,p);rows,hashes=read_rows(folder,p);report=summarize_confirmation(rows,p)
    report.update(plan_sha256=ph,checkpoints=hashes)
    target=folder/('complete.json' if report['status']=='complete' else f'partial-{len(rows):04d}.json')
    if target.exists():
        if json.loads(target.read_text())!=report:raise ValueError('Existing report differs')
    else:save(target,report)
    print(json.dumps(report['families']),report['decision'],flush=True)


if __name__=='__main__':main()
