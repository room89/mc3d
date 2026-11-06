"""Explain existing shock gates without modifying data, tolerances or verdicts."""
import argparse
from bisect import bisect_left
import hashlib
import json
from pathlib import Path
import statistics as st

from shock_analysis import equivalence, interval, metrics


def gate(difference, uncertainty, tolerance):
    return dict(difference=difference, uncertainty=uncertainty, tolerance=tolerance,
                precision_limit=tolerance/2,
                mean_outside_tolerance=abs(difference)>tolerance,
                precision_pass=uncertainty<=tolerance/2,
                coverage_pass=abs(difference)+uncertainty<=tolerance,
                passed=equivalence(difference, uncertainty, tolerance))


def comparison(a, b):
    density=[gate(x-y,u+v,.01) for x,y,u,v in zip(
        a['density'],b['density'],a['density_uncertainty'],b['density_uncertainty'])]
    result={'density':dict(points=len(density),passed=sum(g['passed'] for g in density),
        means_outside=sum(g['mean_outside_tolerance'] for g in density),
        precision_fail=sum(not g['precision_pass'] for g in density),
        max_absolute_difference=max(abs(g['difference']) for g in density),
        max_uncertainty=max(g['uncertainty'] for g in density))}
    for key in ('width_10_90','gradient_thickness'):
        a0,b0=a[key],b[key]
        statistical=a0['uncertainty']+b0['uncertainty']
        processing=a['processing_uncertainty']+b['processing_uncertainty'] if key=='gradient_thickness' else 0
        result[key]=gate(a0['mean']-b0['mean'],statistical+processing,.03*b0['mean'])
        result[key].update(statistical_uncertainty=statistical,processing_uncertainty=processing)
    return result


def at(x,y,p):
    if p<x[0] or p>x[-1]:raise ValueError('Point outside profile')
    i=max(1,bisect_left(x,p))
    return y[i-1]+(p-x[i-1])/(x[i]-x[i-1])*(y[i]-y[i-1])


def uncertainty_sources(trials, key):
    return dict(between_seeds=interval([st.mean(w[key] for w in run) for run in trials])[1],
        blocks={str(block):[interval([st.mean(w[key] for w in run[i:i+block])
            for i in range(0,len(run),block)])[1] for run in trials] for block in (8,16)})


def trial_audit(trials, summary):
    if len(trials)!=3 or any(len(run)!=64 for run in trials):
        raise ValueError('Expected three seeds with 64 production windows each')
    result={'seeds':[],'uncertainty_sources':{},'aligned_mean_profile_diagnostic':{}}
    for seed,run in zip((123,456,789),trials):
        item=dict(seed=seed,plateau_failed_windows=sum(w['plateau_error']>.01 for w in run),
            plateau_max=max(w['plateau_error'] for w in run),
            boundary_failed_windows=sum(w['boundary_distance']<10 for w in run),
            boundary_min=min(w['boundary_distance'] for w in run),halves={})
        # Signed edge averages distinguish transient/window noise from a persistent offset.
        edges=[]
        for w in run:
            n=max(1,len(w['density_normalized'])//10);y=w['density_normalized']
            edges.append((st.mean(y[:n]),st.mean(y[-n:])-1))
        item['mean_plateau_offsets']=[st.mean(p[i] for p in edges) for i in (0,1)]
        for key in ('width_10_90','gradient_thickness'):
            first=[st.mean(w[key] for w in run[i:i+8]) for i in range(0,32,8)]
            last=[st.mean(w[key] for w in run[i:i+8]) for i in range(32,64,8)]
            a,ua=interval(first);b,ub=interval(last)
            item['halves'][key]=gate(a-b,ua+ub,.03*summary[key]['mean'])
        result['seeds'].append(item)
    for key in ('width_10_90','gradient_thickness'):
        result['uncertainty_sources'][key]=uncertainty_sources(trials,key)
    # Exploratory estimand only: align windows, average density, then measure.
    # Do not substitute it into the existing acceptance results.
    grid=[-10+i*.05 for i in range(401)];profiles=[];per_seed=[]
    for run in trials:
        y=[st.mean(at(w['x_aligned'],w['density_normalized'],x) for w in run) for x in grid]
        profiles.append(y);m=metrics(grid,y,0,1,1)
        per_seed.append({k:m[k] for k in ('width_10_90','gradient_thickness','processing_uncertainty')})
    y=[st.mean(p[i] for p in profiles) for i in range(len(grid))]
    m=metrics(grid,y,0,1,1)
    result['aligned_mean_profile_diagnostic']=dict(per_seed=per_seed,
        all_seeds={k:m[k] for k in ('width_10_90','gradient_thickness','processing_uncertainty')},
        warning='Exploratory estimator change, not an acceptance result')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('series',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();source=a.series/'validation.json';raw=source.read_bytes();report=json.loads(raw)
    result={'source_report_sha256':hashlib.sha256(raw).hexdigest(),
            'interpretation':'Diagnostics only; original gates and report are unchanged',
            'binary_sha256':report['binary_sha256'],'sparta_sha256':report['sparta_sha256'],
            'profile_sha256':{},'mc3d':{},'comparisons':{},'cross_code':{}}
    for name,summary in report['summaries'].items():
        data=(a.series/f'profiles-{name}.json').read_bytes()
        result['profile_sha256'][name]=hashlib.sha256(data).hexdigest()
        trials=json.loads(data)
        result['mc3d'][name]=trial_audit(trials,summary)
    for pair in report['convergence']['comparisons']:
        left,right=pair.split('-');result['comparisons'][pair]=comparison(report['summaries'][left],report['summaries'][right])
    s=report['sparta_summaries']
    for label,left,right in [('grid',s['D'],s['C']),('time',s['E'],s['C']),('mc3d',report['summaries']['C'],s['C'])]:
        result['cross_code'][label]=comparison(left,right)
    for section, original in [('comparisons', report['convergence']['comparisons']),
                              ('cross_code', report['cross_code']['comparisons'])]:
        for label, detail in result[section].items():
            checks={key:(value['passed']==value['points'] if key=='density' else value['passed'])
                    for key,value in detail.items()}
            if checks!=original[label]['checks']:
                raise ValueError('Audit disagrees with saved gates: '+label)
    if a.output.exists():raise FileExistsError(a.output)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if source.read_bytes()!=raw:raise RuntimeError('Source report changed during audit')
    print(a.output)


if __name__=='__main__':main()
