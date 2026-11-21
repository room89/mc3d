"""Verify complete or partial order-envelope screen evidence, never qualify it."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from calibrate_order_sensitivity import SOURCES
from calibrate_shock_ensemble import summarize
from shock_order_sensitivity import METHOD
from shock_validation import digest, save
from validate_shock_estimator import truth

LIMITATION = 'Screen does not qualify coverage or model physical heteroscedastic noise; fresh complete validation remains required.'


def expected_protocol(windows, repetitions, start):
    if windows not in (32,64) or type(repetitions)!=int or not 2<=repetitions<10000 or type(start)!=int or start<=0:
        raise ValueError('Invalid expected screen parameters')
    return dict(method=METHOD,status='screen_only',group_count=12,subgroups=6,bootstrap_repeats=1024,
        windows=windows,seed_start=start,replicates=repetitions,motion_amplitudes=[0.,.05,1.],
        sources={name:digest(Path(__file__).with_name(name)) for name in SOURCES},limitation=LIMITATION)


def scenario_keys():
    return [(k,t,s,n,a) for k in ('logistic','gompertz') for t in (0.,.8) for s in (0.,.2)
            for n in (.005,.015,.03) for a in (0.,.05,1.)]


def validate_result(row, spec):
    if type(row.get('repetition'))!=int or row['repetition']!=spec[-1]:
        raise ValueError('Trial repetition mismatch')
    if 'error' in row:
        if set(row)!={'repetition','error'} or not isinstance(row['error'],str) or not row['error']:
            raise ValueError('Invalid retained failure')
        return
    names=('estimate','statistical_radius','processing_sensitivity','relative_error')
    values,u,proc,relative=[np.asarray(row[name],float) for name in names]
    if any(a.shape!=(2,) or not np.isfinite(a).all() for a in (values,u,proc,relative)) or np.any(values<=0) or np.any(u<0) or np.any(proc<0):
        raise ValueError('Invalid metric estimate or radius')
    exact=truth(spec[0])
    for key,expected in [('statistical_covers',(abs(values-exact)<=u).tolist()),
                         ('with_processing_covers',(abs(values-exact)<=u+proc).tolist())]:
        actual=row[key]
        if len(actual)!=2 or any(type(v)!=bool for v in actual) or actual!=expected:
            raise ValueError('Metric coverage flag inconsistent')
    if not np.allclose(relative,values/exact-1,rtol=1e-12,atol=1e-14):
        raise ValueError('Relative metric error inconsistent')
    error,radius,statistical,processing=[row[k] for k in ('density_error','density_radius','density_statistical_radius','density_processing_sensitivity')]
    if any(isinstance(v,bool) or not math.isfinite(v) or v<0 for v in (error,radius,statistical,processing)):
        raise ValueError('Invalid density error or radius')
    if not math.isclose(radius,statistical+processing,rel_tol=1e-12,abs_tol=1e-14):
        raise ValueError('Density uncertainty budget inconsistent')
    if type(row['density_covers'])!=bool or row['density_covers']!=(error<=radius):
        raise ValueError('Density coverage flag inconsistent')
    difference=row['order_density_difference'];metric_differences=row['order_metric_differences']
    if not math.isfinite(difference) or not 0<=difference<=processing+1e-14:
        raise ValueError('Invalid order density difference')
    if set(metric_differences)!={'width_10_90','gradient_thickness'}:
        raise ValueError('Missing order metric difference')
    for i,name in enumerate(('width_10_90','gradient_thickness')):
        value=metric_differences[name]
        if not math.isfinite(value) or not 0<=value<=proc[i]+1e-14:
            raise ValueError('Order metric difference exceeds processing allowance')


def verify(report_path, windows, repetitions, start):
    report_path=Path(report_path);folder=report_path.with_suffix('.parts')
    p=expected_protocol(windows,repetitions,start);keys=scenario_keys()
    wanted={key+(i,) for key in keys for i in range(repetitions)}
    report_hash=digest(report_path) if report_path.exists() else None
    paths=sorted(folder.glob('*.json'));hashes={};rows={};seen=set()
    if not folder.is_dir() or not paths:
        raise ValueError('No screen checkpoint evidence available')
    for path in paths:
        before=digest(path);saved=json.loads(path.read_text());spec=saved['spec']
        if (not isinstance(spec,list) or len(spec)!=6 or type(spec[0])!=str or
            any(type(v)!=float for v in spec[1:5]) or type(spec[5])!=int):
            raise ValueError('Invalid scenario identity types')
        spec=tuple(spec)
        if spec not in wanted or spec in seen or path.name!='-'.join(map(str,spec))+'.json':
            raise ValueError('Unknown, repeated or misnamed checkpoint')
        if saved['protocol']!=p:
            raise ValueError('Screen protocol or frozen source mismatch')
        validate_result(saved['result'],spec)
        seen.add(spec);rows[spec]=saved['result'];hashes[path.name]=before
    if report_hash is not None:
        report=json.loads(report_path.read_text())
        if seen!=wanted or report.get('protocol')!=p or report.get('physical_acceptance')!='incomplete':
            raise ValueError('Final report lacks the complete matching checkpoint matrix')
        cases=report['cases'];reported_keys=[(c['shape'],c['phi'],c['spatial'],c['noise'],c['motion_amplitude']) for c in cases]
        if len(cases)!=72 or len(set(reported_keys))!=72 or set(reported_keys)!=set(keys):
            raise ValueError('Missing or repeated final scenarios')
        for key,case in zip(reported_keys,cases):
            expected=summarize([rows[key+(i,)] for i in range(repetitions)])
            if any(case.get(k)!=v for k,v in expected.items()):
                raise ValueError('Final summary differs from checkpoints')
        status='screen_with_failures' if any('error' in r for r in rows.values()) else 'screen_only'
        if report['status']!=status:
            raise ValueError('Final status hides rejected trials')
    if any(digest(folder/name)!=h for name,h in hashes.items()):
        raise ValueError('Checkpoint changed during verification')
    if report_hash is not None and digest(report_path)!=report_hash:
        raise ValueError('Report changed during verification')
    if expected_protocol(windows,repetitions,start)!=p:
        raise ValueError('Frozen source changed during verification')
    rejected=sum('error' in r for r in rows.values())
    completed_cases=[]
    for key in keys:
        if all(key+(i,) in rows for i in range(repetitions)):
            summary=summarize([rows[key+(i,)] for i in range(repetitions)])
            completed_cases.append(dict(shape=key[0],phi=key[1],spatial=key[2],noise=key[3],motion_amplitude=key[4],
                trial_count=repetitions,rejected_trials=len(summary['errors']),metric_coverage=summary['with_processing_coverage'],density_coverage=summary['density_coverage']))
    return dict(status='complete_screen_evidence_verified' if report_hash is not None else 'partial_checkpoint_snapshot_verified',
        physical_acceptance='incomplete',protocol=p,report_sha256=report_hash,checkpoint_sha256=hashes,
        scheduled_trials=len(wanted),completed_trials=len(rows),rejected_trials=rejected,missing_trials=len(wanted-seen),
        completed_scenarios=completed_cases,fixture_seed_range=[start,start+(repetitions-1)*100+71],
        bootstrap_seed_range=[start+1000000,start+1000000+repetitions-1],
        limitations=['Snapshot excludes files committed after directory enumeration; a live process may continue.',
                    'Coverage flags are recomputed from stored metric estimates and density errors; expensive trajectories are not rerun.',
                    'Complete screen evidence is not full coverage qualification or physical acceptance.'])


def verify_pair(full,half,repetitions,full_start,half_start):
    a=verify(full,64,repetitions,full_start);b=verify(half,32,repetitions,half_start)
    ranges=[r[k] for r in (a,b) for k in ('fixture_seed_range','bootstrap_seed_range')]
    for i,left in enumerate(ranges):
        for right in ranges[i+1:]:
            if max(left[0],right[0])<=min(left[1],right[1]):
                raise ValueError('Overlapping fixture or bootstrap seed ranges')
    if a['protocol']['sources']!=b['protocol']['sources']:
        raise ValueError('Pair uses different source snapshots')
    complete=all(r['status']=='complete_screen_evidence_verified' for r in (a,b))
    return dict(status='complete_screen_pair_verified' if complete else 'partial_screen_pair_snapshot_verified',physical_acceptance='incomplete',full=a,half=b,
                verifier_sha256=digest(Path(__file__)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',type=Path,required=True);parser.add_argument('--half',type=Path,required=True)
    parser.add_argument('--full-seed-start',type=int,required=True);parser.add_argument('--half-seed-start',type=int,required=True)
    parser.add_argument('--replicates',type=int,default=8);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=verify_pair(args.full,args.half,args.replicates,args.full_seed_start,args.half_seed_start)
    save(args.output,result)
    print('Verified screen evidence snapshot; no qualification verdict.')

if __name__=='__main__':main()
