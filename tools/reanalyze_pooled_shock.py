"""Analyze verified six-subgroup C data with independent pooled peak selection.

Never declares complete physical acceptance: C alone does not establish grid,
time-step, population and domain convergence or experimental agreement.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
from reanalyze_shock_campaign import process_group,plateau_checks
from reanalyze_shock_mean import compare
from shock_group_pooling import crossfit_pool_points
from shock_ensemble_estimator import estimate
from shock_validation import digest,save


def validate_plan(plan,root):
    if plan.get('variant')!='C' or plan.get('total_groups')!=12 or plan.get('subgroups_per_group')!=6:
        raise ValueError('Require twelve groups of six C subgroups')
    groups=plan['groups']
    if len(groups)!=12 or [g['name'] for g in groups]!=[f'g{i:02d}' for i in range(12)]:
        raise ValueError('Invalid group identities')
    seeds=[];folders=[]
    for g in groups:
        if len(g['subgroups'])!=6:raise ValueError('Six subgroups required')
        for item in g['subgroups']:
            if len(item['seeds'])!=3:raise ValueError('Three seeds per subgroup required')
            seeds.extend(item['seeds']);folders.append((root/item['path']).resolve())
    if len(set(seeds))!=216 or any(type(s)!=int or not 0<s<2**31 for s in seeds):
        raise ValueError('Subgroups must have 216 distinct positive seeds')
    if len(set(folders))!=72:raise ValueError('Subgroups reuse folders')
    return groups


def check_coverage(path,windows):
    report=json.loads(path.read_text());p=report['protocol']
    if p.get('subgroups')!=6 or p.get('windows')!=windows or p.get('group_count')!=12 or p.get('bootstrap_repeats')!=1024:
        raise ValueError('Wrong pooled coverage dimensions')
    if not all(p.get(k) is True for k in ('crossfit_pooled_peak','center_profile','width_band','interpolation_sensitivity')):
        raise ValueError('Coverage does not use independent pooled peaks')
    names=('calibrate_shock_ensemble.py','shock_group_pooling.py','shock_ensemble_estimator.py','shock_alignment.py',
           'shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py','experiment_crossfit_gradient.py','validate_shock_estimator.py')
    for name in names:
        if p['sources'].get(name)!=digest(Path(__file__).with_name(name)):
            raise ValueError('Coverage source changed: '+name)
    required={(k,t,s,n) for k in ('logistic','gompertz') for t in (0.,.8) for s in (0.,.2) for n in (.005,.015,.03)}
    cases=report['cases']
    if len(cases)!=24 or {(c['shape'],c['phi'],c['spatial'],c['noise']) for c in cases}!=required:
        raise ValueError('Full 24-case coverage matrix required')
    for c in cases:
        if c['trial_count']<128 or len(c['runs'])!=c['trial_count'] or len({r['repetition'] for r in c['runs']})!=c['trial_count']:
            raise ValueError('Incomplete or repeated coverage trials')
    return dict(sha256=digest(path),status=report['status'],windows=windows,
        minimum_total_metric_coverage=np.min([c['with_processing_coverage'] for c in cases],axis=0).tolist(),
        minimum_density_coverage=min(c['density_coverage'] for c in cases),
        rejected_trials=sum(len(c['errors']) for c in cases),
        limitation='Matching and completion checks are not automatic qualification of nominal coverage')


def pool_edges(edges):
    if len(edges)!=6:raise ValueError('Six edge summaries required')
    if len({e['window_count'] for e in edges})!=1:raise ValueError('Unequal subgroup sampling')
    result={p:np.mean([e[p] for e in edges],axis=0).tolist() for p in ('full','first','last')}
    result.update(boundary_min=min(e['boundary_min'] for e in edges),
        window_outliers=sum(e['window_outliers'] for e in edges),window_count=sum(e['window_count'] for e in edges))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--coverage',type=Path,required=True);p.add_argument('--half-coverage',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.jobs<1:p.error('Positive worker count required')
    if a.output.exists():raise FileExistsError(a.output)
    plan_hash=digest(a.plan);plan=json.loads(a.plan.read_text());root=a.plan.resolve().parent
    groups=validate_plan(plan,root)
    baseline=root/plan['baseline_campaign']
    if digest(baseline)!=plan['baseline_campaign_sha256']:raise ValueError('Baseline campaign changed')
    # Fail before expensive analysis if any planned subgroup is not completed.
    for g in groups:
        for item in g['subgroups']:
            file=root/item['path']/'validation.json'
            if not file.exists():raise ValueError('Planned subgroup has no execution report: '+g['name'])
            if 'validation_sha256' in item and digest(file)!=item['validation_sha256']:
                raise ValueError('Existing subgroup report changed')
            from shock_independent_campaign import execution_complete
            report=json.loads(file.read_text())
            if report.get('seeds',[123,456,789])!=item['seeds'] or not execution_complete(report,['C'],item['seeds']):
                raise ValueError('Incomplete or mismatched subgroup: '+g['name'])
    coverage=[check_coverage(a.coverage,64),check_coverage(a.half_coverage,32)]
    names=('reanalyze_pooled_shock.py','reanalyze_shock_campaign.py','reanalyze_shock_mean.py',
           'shock_group_pooling.py','shock_ensemble_estimator.py','shock_alignment.py','shock_crossfit_gradient.py',
           'shock_mean_profile.py','shock_analysis.py','analyze_shock_validation.py','audit_shock_validation.py',
           'shock_independent_campaign.py','shock_validation.py','normal_shock.py','physical_config.py','argon_vhs.py')
    protocol=dict(plan_sha256=plan_hash,coverage=coverage,bootstrap_repeats=1024,resampling_seed=20261109,
        sources={n:digest(Path(__file__).with_name(n)) for n in names})
    parts=a.output.with_suffix('.parts')
    if parts.exists() and not a.resume:raise FileExistsError(parts)
    parts.mkdir(parents=True,exist_ok=True)
    tasks=[(g['name']+f'-s{i}',str((root/item['path']).resolve()),item['seeds'],['C'],
            plan['binary_sha256'],plan['sparta_sha256'],str(parts),protocol)
           for g in groups for i,item in enumerate(g['subgroups'])]
    results={};hashes={}
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for name,result,inputs in pool.map(process_group,tasks):
            results[name]=result;hashes[name]=inputs;print(name,'analyzed',flush=True)
    summaries={};halves={};windows={};plateaus={}
    for key in ('C','sparta-C'):
        phases={}
        for phase in ('full','first','last','doubled'):
            pooled=[crossfit_pool_points([results[g['name']+f'-s{i}'][key][phase] for i in range(6)]) for g in groups]
            phases[phase]=estimate(pooled,repeats=1024,seed=20261109,center_profile=True,width_band=True,interpolation_sensitivity=True)
        summaries[key]=phases['full'];halves[key]=compare(phases['first'],phases['last'])
        windows[key]=compare(phases['full'],phases['doubled'])
        plateaus[key]=plateau_checks([pool_edges([results[g['name']+f'-s{i}'][key]['edges'] for i in range(6)]) for g in groups])
    if digest(a.plan)!=plan_hash or digest(baseline)!=plan['baseline_campaign_sha256']:
        raise ValueError('Plan or baseline changed during analysis')
    if any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()):
        raise ValueError('Analysis source changed')
    save(a.output,dict(status='candidate_pilot_diagnostics_only',physical_acceptance='incomplete',protocol=protocol,
        variants=['C'],input_hashes=hashes,summaries=summaries,comparisons={'C vs sparta-C':compare(summaries['C'],summaries['sparta-C'])},
        half_comparisons=halves,window_length_comparisons=windows,plateaus=plateaus))


if __name__=='__main__':main()
