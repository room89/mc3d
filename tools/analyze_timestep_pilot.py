"""Descriptive whole-run analysis of the frozen C/E pilot; never acceptance."""
import argparse
import json
from pathlib import Path
import numpy as np
from shock_validation import digest, verify_run
from analyze_shock_validation import profile
from reanalyze_shock_mean import sparta_density
from shock_analysis import crossing, monotone


def summarize(raw, x):
    a=np.asarray(raw,float);x=np.asarray(x,float)
    if a.shape!=(3,64,len(x)) or not np.isfinite(a).all() or np.any(np.diff(x)<=0):
        raise ValueError('Three complete finite 64-window runs on increasing nodes required')
    means=a.mean(axis=1)
    return dict(run_means=means.tolist(),mean=means.mean(axis=0).tolist(),
        between_run_sd=means.std(axis=0,ddof=1).tolist(),
        half_changes=(a[:,32:].mean(axis=1)-a[:,:32].mean(axis=1)).tolist(),
        run_fronts=[float(crossing(x,monotone(y),.5)) for y in means])


def completed_reports(root, protocol):
    reports={}
    for v in protocol['variants']:
        terminal=root/(v+'-exit.json')
        if not terminal.exists():raise ValueError('Runner is not terminal: '+v)
        if json.loads(terminal.read_text())['returncode'] not in (0,1,2):
            raise ValueError('Runner execution error: '+v)
        r=json.loads((root/v/'validation.json').read_text())
        if r['seeds']!=protocol['seeds'] or r['variants']!=[v]:raise ValueError('Identity mismatch')
        if r['binary_sha256']!=protocol['binary_sha256'] or r['sparta_sha256']!=protocol['sparta_sha256']:
            raise ValueError('Binary mismatch')
        keys={f'{v}-{s}' for s in protocol['seeds']}
        for field in ('runs','sparta_runs'):
            if set(r[field])!=keys or any(r[field][k]['status']!='passed' for k in keys):
                raise ValueError('Incomplete execution')
        reports[v]=r
    return reports


def analyze(root):
    root=Path(root).resolve();p=json.loads((root/'protocol.json').read_text())
    launch=json.loads((root/'launch.json').read_text())
    if digest(root/'protocol.json')!=launch['protocol_sha256']:raise ValueError('Protocol changed')
    if p['variants']!=['C','E'] or p['seeds']!=[210001,210002,210003]:raise ValueError('Frozen pilot required')
    reports=completed_reports(root,p)  # Never inspect partial physical results.
    hashes={}
    for f in root.rglob('*'):
        if f.is_file() and f.suffix not in ('.pyc',) and '__pycache__' not in f.parts:
            hashes[str(f.relative_to(root))]=digest(f)
    for n,h in launch['sources'].items():
        if digest(root/'runtime-tools'/n)!=h:raise ValueError('Runtime source changed')
    results={};grid=None
    for v,r in reports.items():
        for code in ('mc3d','sparta'):
            runs=[]
            for seed in p['seeds']:
                key=f'{v}-{seed}';folder=root/v/(key if code=='mc3d' else 'sparta-'+key)
                if code=='mc3d':
                    if not verify_run(folder/'output',p['binary_sha256']):raise ValueError('Incomplete manifest')
                    windows=profile(folder)
                else:
                    for name,h in r['sparta_runs'][key]['hashes'].items():
                        if digest(folder/name)!=h:raise ValueError('SPARTA output changed')
                    windows=sparta_density(folder,{})
                raw=[]
                for w in windows:
                    x=np.asarray(w['x_aligned'])+w['front']
                    mask=(x>=-14)&(x<=14);x=x[mask]
                    if grid is None:grid=x
                    if x.shape!=grid.shape or not np.allclose(x,grid,rtol=0,atol=1e-12):raise ValueError('Grid mismatch')
                    raw.append(np.asarray(w['density_normalized'])[mask])
                runs.append(raw)
            results[code+'-'+v]=summarize(runs,grid)
    comparisons={}
    for a,b in [('mc3d-C','sparta-C'),('mc3d-E','sparta-E'),('mc3d-C','mc3d-E'),('sparta-C','sparta-E')]:
        d=np.array(results[a]['mean'])-results[b]['mean']
        comparisons[a+' vs '+b]=dict(difference=d.tolist(),max_abs_difference=float(abs(d).max()))
    if any(digest(root/n)!=h for n,h in hashes.items()):raise ValueError('Evidence changed during analysis')
    return dict(status='descriptive_precision_pilot',physical_acceptance='incomplete',
        x=grid.tolist(),seeds=p['seeds'],summaries=results,comparisons=comparisons,input_sha256=hashes,
        analysis_sources={n:digest(Path(__file__).with_name(n)) for n in ('analyze_timestep_pilot.py','analyze_shock_validation.py','reanalyze_shock_mean.py','shock_analysis.py')},
        limitations=['Three independent runs per case, not independent windows or paired solver paths.',
            'Between-run SD is descriptive; no qualified confidence band or sample-size guarantee.',
            'Engineering-domain pilot, not experimental-condition confirmation or thickness qualification.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    if a.output.resolve().is_relative_to(a.root.resolve()):raise ValueError('Write analysis outside immutable campaign')
    r=analyze(a.root)
    with a.output.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
