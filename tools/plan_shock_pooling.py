"""Prepare a sixfold C campaign plan without starting or reusing new solver runs."""
import argparse
import json
import os
from pathlib import Path
from shock_validation import digest,save
from shock_independent_campaign import verify_baseline


def assignments(existing,seed_start=20000):
    if len(existing)!=12 or any(len(s)!=3 for s in existing):
        raise ValueError('Twelve three-seed baseline groups required')
    groups=[]
    for i,seeds in enumerate(existing):
        groups.append([list(seeds)]+[[seed_start+15*i+3*j+k for k in range(3)] for j in range(5)])
    flat=[s for group in groups for triplet in group for s in triplet]
    if len(flat)!=216 or len(set(flat))!=216 or any(type(s)!=int or not 0<s<2**31 for s in flat):
        raise ValueError('Require 216 distinct positive 31-bit seeds')
    return groups


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('baseline',type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--seed-start',type=int,default=20000);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    baseline=a.baseline.resolve();source=baseline/'campaign.json';source_hash=digest(source);campaign=json.loads(source.read_text())
    if campaign['status']!='execution_complete' or campaign['protocol']['variants']!=['C']:
        raise ValueError('Completed C-only baseline campaign required')
    cp=campaign['protocol'];old=campaign['groups']
    if set(old)!={f'g{i:02d}' for i in range(12)}:raise ValueError('Invalid baseline group identities')
    seeds=assignments([old[f'g{i:02d}']['seeds'] for i in range(12)],a.seed_start)
    if cp['seeds']!=[s[0] for s in seeds]:raise ValueError('Baseline seeds differ from protocol')
    if len({(baseline/g['path']).resolve() for g in old.values()})!=12:
        raise ValueError('Baseline groups reuse the same folder')
    if digest(baseline/old['g00']['path']/'validation.json')!=cp['baseline_report_sha256']:
        raise ValueError('Original baseline report changed')
    root=a.output.resolve().parent;groups=[]
    for i,subgroups in enumerate(seeds):
        name=f'g{i:02d}';folder=(baseline/old[name]['path']).resolve()
        report=json.loads((folder/'validation.json').read_text())
        if report.get('seeds',[123,456,789])!=subgroups[0]:raise ValueError('Baseline run seeds differ')
        verify_baseline(folder,report,['C'],subgroups[0],cp['binary_sha256'],cp['sparta_sha256'])
        items=[dict(path=os.path.relpath(folder,root),seeds=subgroups[0],status='existing_verified',
                    validation_sha256=digest(folder/'validation.json'))]
        for j,triplet in enumerate(subgroups[1:],1):
            items.append(dict(path=f'{name}-s{j}',seeds=triplet,status='planned'))
        groups.append(dict(name=name,subgroups=items))
    if digest(source)!=source_hash:raise ValueError('Baseline campaign changed')
    root.mkdir(parents=True,exist_ok=True)
    save(a.output,dict(status='plan_only',physical_acceptance='incomplete',variant='C',
        baseline_campaign=os.path.relpath(source,root),baseline_campaign_sha256=source_hash,
        binary=os.path.relpath((baseline/cp['binary']).resolve(),root),binary_sha256=cp['binary_sha256'],
        sparta_binary=os.path.relpath((baseline/cp['sparta_binary']).resolve(),root),sparta_sha256=cp['sparta_sha256'],
        required_coverage=dict(subgroups=6,windows=[64,32],cases_per_matrix=24,trials_per_case=128),
        total_groups=12,subgroups_per_group=6,unique_seed_count=216,new_solver_runs=360,groups=groups))
    print('Prepared plan: 216 unique seeds, 360 new solver runs; none started')


if __name__=='__main__':main()
