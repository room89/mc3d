"""Read-back and deterministic input reconstruction for the plateau screen."""
import argparse
import json
from pathlib import Path
import numpy as np
from run_profile_likelihood_screen import generate,digest
from shock_profile_likelihood_fit import KNOTS
from shock_temporal_covariance import estimate


def audit(folder):
    folder=Path(folder);plan_path=folder/'plan.json';report_path=folder/'complete.json'
    p=json.loads(plan_path.read_text());r=json.loads(report_path.read_text());ph=digest(plan_path);rh=digest(report_path)
    if r['plan_sha256']!=ph or p['seed']==p['calibration_seed']:raise ValueError('Plan or independent seed identity mismatch')
    if any(digest(Path(__file__).with_name(n))!=h for n,h in p['sources'].items()):raise ValueError('Frozen source mismatch')
    expected={f'trial-{i:03d}.json' for i in range(len(p['specifications']))}
    if set(r['checkpoints'])!=expected or {f.name for f in folder.glob('trial-*.json')}!=expected:raise ValueError('Trial inventory mismatch')
    rows=[];differences=[]
    for i,spec in enumerate(p['specifications']):
        path=folder/f'trial-{i:03d}.json'
        if digest(path)!=r['checkpoints'][path.name]:raise ValueError('Checkpoint hash mismatch')
        row=json.loads(path.read_text())
        if row['index']!=i or row['specification']!=spec or row['plan_sha256']!=ph:raise ValueError('Trial identity mismatch')
        x,_,_,t,s,_,truth=generate(spec,i,p)
        calibration=dict(p,seed=p['calibration_seed']);cx,_,_,_,cs,cy,_=generate(spec,i,calibration)
        mask=np.abs(cx)>=8+spec['amplitude'];covariance=estimate(cy[:,:,mask],(cx[mask]>0).astype(float),cs[np.ix_(mask,mask)])
        if covariance!=row['covariance'] or t.tolist()!=row['true_temporal_covariance']:raise ValueError('Reconstructed covariance mismatch')
        error=float(np.linalg.norm(np.array(covariance['temporal_covariance'])-t)/np.linalg.norm(t))
        if error!=row['relative_covariance_error']:raise ValueError('Covariance error mismatch')
        grid=np.linspace(-8,8,321)
        for result in row['results'].values():
            profile=np.asarray(result['profile']);e=np.interp(grid,KNOTS,profile)-np.interp(grid,KNOTS,truth)
            if not np.isfinite(profile).all() or np.any(np.diff(profile)<0) or profile[0]!=0 or profile[-1]!=1:raise ValueError('Invalid fitted profile')
            if float(np.trapezoid(e**2,grid))!=result['ise']:raise ValueError('ISE mismatch')
        rows.append(row);differences.append(row['results']['estimated']['ise']-row['results']['known']['ise'])
    failures={m:sum(z['results'][m]['status']!='converged' for z in rows) for m in ('known','estimated')}
    means={m:float(np.mean([z['results'][m]['ise'] for z in rows])) for m in failures}
    if r['trials']!=len(rows) or r['failures']!=failures or r['mean_ise']!=means or r['mean_relative_covariance_error']!=float(np.mean([z['relative_covariance_error'] for z in rows])):raise ValueError('Aggregate mismatch')
    if digest(report_path)!=rh or digest(plan_path)!=ph or any(digest(folder/n)!=h for n,h in r['checkpoints'].items()):raise ValueError('Evidence changed')
    return dict(status='verified_diagnostic_only',physical_acceptance='incomplete',trials=len(rows),
                mean_paired_ise_difference=float(np.mean(differences)),paired_mcse=float(np.std(differences,ddof=1)/np.sqrt(len(rows))),
                report_sha256=rh,plan_sha256=ph,failures=failures)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,required=True)
    print(json.dumps(audit(parser.parse_args().directory),indent=2))
