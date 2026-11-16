"""Independent synthetic plateau calibration; diagnostic only, no acceptance."""
import argparse
import itertools
import json
import os
from pathlib import Path
import numpy as np
from run_profile_likelihood_screen import protocol,generate,digest,save
from shock_profile_likelihood_fit import fit,KNOTS
from shock_temporal_covariance import estimate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();folder=args.directory
    p=protocol();p.update(seed=202611150000,stage='independent_plateau_calibration',
        specifications=[dict(shape=s,amplitude=a,noise=n,repetition=r,analytic=False) for s,a,n,r in itertools.product(('logistic','gompertz'),(.05,.5),('low','high'),range(16))],
        sources={n:digest(Path(__file__).with_name(n)) for n in ('shock_temporal_covariance.py','screen_temporal_covariance.py','run_profile_likelihood_screen.py','shock_profile_likelihood_fit.py')})
    p['calibration_seed']=202611160000
    if folder.exists():
        if not args.resume or json.loads((folder/'plan.json').read_text())!=p:raise ValueError('Matching plan and resume required')
    else:folder.mkdir(parents=True);save(folder/'plan.json',p)
    ph=digest(folder/'plan.json');os.nice(19);rows=[]
    for i,spec in enumerate(p['specifications']):
        if any(digest(Path(__file__).with_name(n))!=h for n,h in p['sources'].items()):raise ValueError('Source changed')
        path=folder/f'trial-{i:03d}.json'
        if path.exists():
            row=json.loads(path.read_text())
            if row['plan_sha256']!=ph or row['index']!=i or row['specification']!=spec:raise ValueError('Checkpoint mismatch')
        else:
            x,configs,_,true_t,spatial,y,truth=generate(spec,i,p)
            calibration=dict(p,seed=p['calibration_seed'])
            cx,_,_,_,cs,cy,_=generate(spec,i,calibration)
            # Exact plateaus for this synthetic model, including every allowed shift.
            mask=np.abs(cx)>=8+spec['amplitude'];means=(cx[mask]>0).astype(float)
            covariance=estimate(cy[:,:,mask],means,cs[np.ix_(mask,mask)])
            learned=np.array(covariance['temporal_covariance']);results={}
            grid=np.linspace(-8,8,321);exact=np.interp(grid,KNOTS,truth)
            for name,t in [('known',true_t),('estimated',learned)]:
                z=fit(y,x,t,spatial,configs,np.ones(len(configs)))
                error=np.interp(grid,KNOTS,z['profile'])-exact
                z['ise']=float(np.trapezoid(error**2,grid));results[name]=z
            row=dict(index=i,specification=spec,plan_sha256=ph,covariance=covariance,
                     true_temporal_covariance=true_t.tolist(),results=results,
                     relative_covariance_error=float(np.linalg.norm(learned-true_t)/np.linalg.norm(true_t)))
            save(path,row)
        rows.append(row);print('completed',len(rows),'/',len(p['specifications']),flush=True)
    failures={m:sum(r['results'][m]['status']!='converged' for r in rows) for m in ('known','estimated')}
    report=dict(status='diagnostic_complete',physical_acceptance='incomplete',trials=len(rows),failures=failures,
        mean_relative_covariance_error=float(np.mean([r['relative_covariance_error'] for r in rows])),
        mean_ise={m:float(np.mean([r['results'][m]['ise'] for r in rows])) for m in failures},
        plan_sha256=ph,checkpoints={f'trial-{i:03d}.json':digest(folder/f'trial-{i:03d}.json') for i in range(len(rows))},
        limitation='Known spatial covariance and exact synthetic plateaus; no physical noise calibration or interval coverage qualification. Mean errors include unconverged outputs if any; failures must be inspected.')
    target=folder/'complete.json'
    if target.exists():
        if json.loads(target.read_text())!=report:raise ValueError('Report mismatch')
    else:save(target,report)
    print(json.dumps({k:v for k,v in report.items() if k!='checkpoints'}),flush=True)


if __name__=='__main__':main()
