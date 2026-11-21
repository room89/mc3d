"""Plan independent reference-domain replication; do not launch or accept runs."""
import argparse
import os
from pathlib import Path
from shock_high_order_estimator import METHOD
from shock_validation import digest, save
from verify_alsmeyer_matrix import verify


def assignments(existing, seed_start=40000):
    if len(existing)!=3:raise ValueError('One existing three-seed subgroup required')
    triplets=[list(existing)]+[[seed_start+3*i+j for j in range(3)] for i in range(71)]
    flat=[s for triplet in triplets for s in triplet]
    if len(set(flat))!=216 or any(type(s)!=int or not 0<s<2**31 for s in flat):
        raise ValueError('Require 216 distinct positive 31-bit seeds')
    return [triplets[6*i:6*i+6] for i in range(12)]


def plan(mc,sparta,evidence,binary,sparta_binary,output,seed_start=40000):
    mc,sparta,evidence,binary,sparta_binary,output=map(lambda p:Path(p).resolve(),(mc,sparta,evidence,binary,sparta_binary,output))
    if output.exists():raise FileExistsError(output)
    proof=verify(mc,sparta,evidence);seeds=assignments(proof['seeds'],seed_start);root=output.parent
    groups=[]
    for i,subgroups in enumerate(seeds):
        items=[]
        for j,triplet in enumerate(subgroups):
            existing=i==0 and j==0;name=f'g{i:02d}-s{j}'
            items.append(dict(name=name,seeds=triplet,
                              mc3d=os.path.relpath(mc,root) if existing else name+'/mc3d',
                              sparta=os.path.relpath(sparta,root) if existing else name+'/sparta',
                              status='existing_inputs_verified' if existing else 'inputs_not_prepared'))
        groups.append(dict(name=f'g{i:02d}',subgroups=items))
    names=('plan_alsmeyer_replication.py','verify_alsmeyer_matrix.py','prepare_alsmeyer_matrix.py',
           'export_alsmeyer_sparta.py','shock_high_order_estimator.py','physical_config.py',
           'shock_validation.py','normal_shock.py','argon_vhs.py','sparta_preflight.py')
    result=dict(status='plan_only',physical_acceptance='incomplete',method=METHOD,
                reference_evidence=os.path.relpath(evidence,root),reference_evidence_sha256=digest(evidence),
                existing_input_verification=proof,
                binaries={kind:dict(path=os.path.relpath(path,root),sha256=digest(path)) for kind,path in
                          [('mc3d',binary),('sparta',sparta_binary)]},
                sources={n:digest(Path(__file__).with_name(n)) for n in names},
                total_groups=12,subgroups_per_group=6,seeds_per_subgroup=3,unique_seed_count=216,
                mc3d_variants=list('ABCDEFG'),sparta_variants=list('CDE'),
                planned_solver_runs=dict(mc3d=1512,sparta=648,total=2160),
                additional_input_matrices=71,groups=groups,
                required_before_execution=['Complete independent calibration and scientific assessment of the selected estimator.',
                                           'Verify pinned binaries and every prepared reference-domain input matrix.'],
                execution_order=['Run the existing prepared 30-case reference matrix as a pilot.',
                                 'Inspect outputs, runtime, memory, stationarity and observed noise before expanding to the remaining groups.'],
                limitations=['This plan creates no solver inputs or output data beyond the plan itself.',
                             'Existing inputs are reused once; existing solver completion is not inferred.',
                             'Engineering-domain data cannot be substituted for reference-domain groups.',
                             'Seed identity is shared between variants for paired comparisons, but never repeated between independent subgroups.',
                             'The planned replication does not guarantee that precision or convergence criteria will pass.'])
    if proof!=verify(mc,sparta,evidence):raise ValueError('Prepared reference inputs changed')
    root.mkdir(parents=True,exist_ok=True);save(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('mc3d','sparta','evidence','binary','sparta-bin','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seed-start',type=int,default=40000);a=p.parse_args()
    r=plan(a.mc3d,a.sparta,a.evidence,a.binary,a.sparta_bin,a.output,a.seed_start)
    print(r['unique_seed_count'],'unique seeds;',r['planned_solver_runs']['total'],'planned runs; no solver launched.')


if __name__=='__main__':main()
