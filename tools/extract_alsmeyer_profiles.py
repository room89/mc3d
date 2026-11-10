"""Extract verified reference-matrix profiles; never imply physical acceptance."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from analyze_shock_validation import profile
from reanalyze_shock_mean import sparta_density
from shock_validation import digest, save
from verify_alsmeyer_matrix import verify
from verify_alsmeyer_outputs import mc3d_outputs, sparta_outputs

SOURCES=('extract_alsmeyer_profiles.py','verify_alsmeyer_matrix.py','verify_alsmeyer_outputs.py',
         'analyze_shock_validation.py','reanalyze_shock_mean.py','shock_analysis.py','shock_validation.py')


def read(path): return json.loads(path.read_text())


def canonical_windows(rows, cells_x):
    if len(rows)!=64: raise ValueError('Exactly 64 complete windows required')
    result=[]; grid=None
    for wi,row in enumerate(rows):
        x=np.asarray(row['x_aligned'],float); y=np.asarray(row['density_normalized'],float);front=float(row['front'])
        if x.shape!=(cells_x,) or y.shape!=x.shape or not np.isfinite(x).all() or not np.isfinite(y).all() or not math.isfinite(front) or np.any(np.diff(x)<=0):
            raise ValueError('Finite profile on the complete axial grid required')
        raw=x+front
        if grid is None: grid=raw
        elif not np.allclose(raw,grid,rtol=1e-12,atol=1e-11): raise ValueError('Physical grid changes between windows')
        # Preserve fluctuations outside [0,1], every window and the physical front.
        result.append(dict(window=wi,front=front,x_aligned=x.tolist(),density_normalized=y.tolist()))
    return result


def extract(mc_root,sparta_root,evidence,verified,output):
    mc_root=Path(mc_root);sparta_root=Path(sparta_root);output=Path(output)
    if output.exists():raise FileExistsError(output)
    source_hashes={name:digest(Path(__file__).with_name(name)) for name in SOURCES}
    verification=read(verified);verified_hash=digest(verified)
    if verification['status']!='matrix_outputs_verified' or verification['tool_sha256']!=digest(Path(__file__).with_name('verify_alsmeyer_outputs.py')):
        raise ValueError('Current complete output verification required')
    proof=verify(mc_root,sparta_root,evidence)
    if proof!=verification['inputs']:raise ValueError('Input verification identity differs')
    cases={};datasets={};snapshots={}
    def remember(path):snapshots[path]=digest(path)
    remember(verified);remember(evidence);remember(mc_root/'preparation.json');remember(sparta_root/'export.json')
    for kind,root,manifest in [('mc3d',mc_root,'preparation.json'),('sparta',sparta_root,'export.json')]:
        for case in read(root/manifest)['cases']:
            name=case['folder'];folder=root/name;identity=kind+'/'+name
            ref=read(folder/'reference.json');remember(folder/'reference.json')
            if kind=='mc3d':
                actual=mc3d_outputs(folder,verification['binaries'][kind])
                run_manifest=read(folder/'output/manifest.json');remember(folder/'output/manifest.json')
                for entry in run_manifest['inputs']+run_manifest['outputs']:remember(folder/'output'/entry['path'])
                for filename in ('solver.json','input_si.json','units.json'):remember(folder/filename)
                rows=profile(folder);cfg=read(folder/'solver.json');nx=cfg['cells_x'];seed=cfg['seed']
                reference=ref;intervals=read(folder/'input_si.json')['statistics_windows_s']
                data_source=actual['manifest_sha256']
            else:
                execution=read(folder/'execution.json');remember(folder/'execution.json');remember(folder/'run.log')
                actual=sparta_outputs(folder,execution,verification['binaries'][kind])
                actual['execution_sha256']=digest(folder/'execution.json')
                for filename in execution['input_hashes']:remember(folder/filename)
                for entry in execution['outputs']:remember(folder/entry['path'])
                rows=sparta_density(folder,{});nx=ref['source']['numerics']['cells_x'];seed=ref['source']['numerics']['seed']
                reference=ref['reference'];intervals=ref['source']['statistics_windows_s'];data_source=actual['execution_sha256']
            if actual!=verification['cases'].get(identity):raise ValueError('Output verification differs: '+identity)
            variant=name.split('-')[0];key=kind+'-'+variant
            windows=canonical_windows(rows,nx)
            item=dict(seed=seed,windows=windows,source_sha256=data_source,
                      collision_mean_free_path_m=reference['mean_free_path_m'],statistics_windows_s=intervals)
            datasets.setdefault(key,[]).append(item);cases[identity]=actual
    if set(cases)!=set(verification['cases']):raise ValueError('Output verification contains a different case matrix')
    for key,items in datasets.items():
        if sorted(i['seed'] for i in items)!=proof['seeds']:raise ValueError('Missing or duplicate independent seed: '+key)
        items.sort(key=lambda i:i['seed'])
    if any(digest(path)!=value for path,value in snapshots.items()):raise ValueError('Input changed during extraction')
    if proof!=verify(mc_root,sparta_root,evidence):raise ValueError('Prepared matrix changed during extraction')
    if any(digest(Path(__file__).with_name(n))!=h for n,h in source_hashes.items()):raise ValueError('Extraction source changed')
    output.mkdir(parents=True,exist_ok=False);artifacts={}
    for key,items in datasets.items():
        target=output/(key+'.json')
        save(target,dict(status='profiles_only',physical_acceptance='incomplete',
                         coordinates='x and front in collision mean free paths; density=(n-n1)/(n2-n1)',
                         runs=items))
        artifacts[target.name]=digest(target)
    report=dict(status='verified_profiles_extracted',physical_acceptance='incomplete',
                verification_sha256=verified_hash,input_verification=proof,sources=source_hashes,
                case_count=len(cases),dataset_count=len(datasets),artifacts=artifacts,
                limitations=['Three seeds per variant are one subgroup, not twelve independent groups with six subgroups each.',
                             'No uncertainty, stationarity, numerical convergence or experimental agreement verdict is calculated.',
                             'MC3D extraction also enforces the existing conservation checks; SPARTA has separate execution provenance.'])
    save(output/'extraction.json',report);return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('mc3d','sparta','evidence','verified','output-dir'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=extract(a.mc3d,a.sparta,a.evidence,a.verified,a.output_dir)
    print(result['status'],result['case_count'],'cases; physical acceptance remains incomplete')


if __name__=='__main__':main()
