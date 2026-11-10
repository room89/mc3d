"""Paired oracle controls for ensemble bias; never a physical acceptance gate."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from experiment_crossfit_gradient import fixture
from shock_alignment import GRID,align
from shock_crossfit_gradient import derivative,estimate_gradient,peak_position
from shock_mean_profile import measure
from validate_shock_estimator import shape,truth


def quintic_derivative(y,width=1.5):
    dx=GRID[1]-GRID[0];radius=int(round(width/(2*dx)))
    t=np.arange(-radius,radius+1)*dx/width
    weights=np.linalg.pinv(np.column_stack([t**k for k in range(6)]))[1]/width
    return GRID[radius:-radius],np.convolve(y,weights[::-1],mode='valid')


def fixed_slope(y,position,degree):
    x,d=derivative(y,width=1.5) if degree==3 else quintic_derivative(y)
    return float(np.interp(position,x,d))


def run(spec):
    kind,phi,spatial,noise,repetition,seed_start=spec[:6]
    subgroups=spec[6] if len(spec)>6 else 1
    window_count=spec[7] if len(spec)>7 else 64
    position=math.log(math.log(2)) if kind=='gompertz' else 0.
    profiles={'oracle':[],'refined':[]};slopes={name:[] for name in
        ('oracle_cubic','oracle_quintic','oracle_selected_peak','refined_cubic','refined_quintic_same_peak')}
    for group in range(12*subgroups):
        trials,oracle=fixture(kind,seed_start+repetition*100+group,phi,spatial,noise=noise)
        trials=[r[:window_count] for r in trials];oracle=oracle[:,:window_count]
        data,diagnostics=align(trials,refine_passes=1)
        profiles['oracle'].append(oracle.mean((0,1)))
        profiles['refined'].append(data.mean((0,1)))
        for degree,name in ((3,'oracle_cubic'),(5,'oracle_quintic')):
            slopes[name].append(fixed_slope(oracle.mean((0,1)),position,degree))
        selected=[]
        for heldout in range(3):
            training=oracle[[s for s in range(3) if s!=heldout]].mean((0,1))
            peak=peak_position(GRID,training,width=1.5)
            selected.append(fixed_slope(oracle[heldout].mean(0),peak,3))
        slopes['oracle_selected_peak'].append(float(np.mean(selected)))
        slopes['refined_cubic'].append(1/estimate_gradient(data,diagnostics,width=1.5)['gradient_thickness'])
        heldout=[]
        for windows,info in zip(data,diagnostics):
            peak=peak_position(info['template_x'],info['template_density'],width=1.5)
            heldout.append(fixed_slope(windows.mean(0),peak,5))
        slopes['refined_quintic_same_peak'].append(float(np.mean(heldout)))
    exact=truth(kind)
    errors={name+'_width':float(measure(np.mean(rows,axis=0))[0,0]/exact[0]-1)
            for name,rows in profiles.items()}
    errors.update({name:float(1/np.mean(values)/exact[1]-1) for name,values in slopes.items()})
    return dict(repetition=repetition,relative_errors=errors)


def cached(task):
    spec,folder,protocol=task
    kind,phi,spatial,noise,i,_=spec[:6]
    path=Path(folder)/f'{kind}-{phi}-{spatial}-{noise}-{i}.json'
    if path.exists():
        saved=json.loads(path.read_text())
        if saved['protocol']!=protocol or saved['spec']!=list(spec):raise ValueError('Checkpoint mismatch')
        return saved['result']
    try:result=run(spec)
    except ValueError as error:result=dict(repetition=i,error=str(error))
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(dict(protocol=protocol,spec=spec,result=result),allow_nan=False)+'\n')
    tmp.replace(path)
    return result


def summarize(rows):
    valid=[r['relative_errors'] for r in rows if 'error' not in r]
    measures={}
    if len(valid)>=2:
        for name in valid[0]:
            v=np.asarray([r[name] for r in valid])
            measures[name]=dict(mean=float(v.mean()),sd=float(v.std(ddof=1)),
                mean_two_se=float(2*v.std(ddof=1)/math.sqrt(len(v))))
        for name,left,right in [('alignment_width','refined_width','oracle_width'),
                                ('alignment_and_peak_gradient','refined_cubic','oracle_cubic'),
                                ('derivative_order_gradient','refined_quintic_same_peak','refined_cubic')]:
            v=np.asarray([r[left]-r[right] for r in valid])
            measures[name]=dict(mean=float(v.mean()),sd=float(v.std(ddof=1)),
                mean_two_se=float(2*v.std(ddof=1)/math.sqrt(len(v))))
        if all('oracle_selected_peak' in r for r in valid):
            for name,left,right in [('peak_selection_only','oracle_selected_peak','oracle_cubic'),
                                    ('alignment_given_peak_selection','refined_cubic','oracle_selected_peak')]:
                v=np.asarray([r[left]-r[right] for r in valid])
                measures[name]=dict(mean=float(v.mean()),sd=float(v.std(ddof=1)),
                    mean_two_se=float(2*v.std(ddof=1)/math.sqrt(len(v))))
    return dict(valid_count=len(valid),failed_count=len(rows)-len(valid),summary=measures,runs=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--replicates',default=64,type=int)
    p.add_argument('--jobs',default=4,type=int);p.add_argument('--seed-start',default=18000000,type=int)
    p.add_argument('--subgroups',type=int,choices=[1,6],default=1)
    p.add_argument('--windows',type=int,choices=[32,64],default=64)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.replicates<2 or a.jobs<1:p.error('At least two repetitions and one worker required')
    if a.output.exists():raise FileExistsError(a.output)
    folder=a.output.with_suffix('.parts')
    if folder.exists() and not a.resume:raise FileExistsError(folder)
    folder.mkdir(parents=True,exist_ok=True)
    names=('diagnose_shock_ensemble_bias.py','experiment_crossfit_gradient.py','shock_alignment.py',
           'shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py','validate_shock_estimator.py')
    protocol=dict(seed_start=a.seed_start,group_count=12,subgroups=a.subgroups,windows=a.windows,sources={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names})
    cases=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for kind,phi,spatial,noise in [('logistic',0.,.2,.03),('gompertz',0.,.2,.03),
                                      ('logistic',0.,0.,.015),('gompertz',0.,0.,.015)]:
            rows=[]
            for row in pool.map(cached,[((kind,phi,spatial,noise,i,a.seed_start,a.subgroups,a.windows),str(folder),protocol) for i in range(a.replicates)]):
                rows.append(row);print(kind,phi,spatial,noise,len(rows),'/',a.replicates,flush=True)
            cases.append(dict(shape=kind,phi=phi,spatial=spatial,noise=noise,**summarize(rows)))
    result=dict(status='diagnostic_only',protocol=protocol,cases=cases,
        limitation='Oracle controls use known synthetic front positions; not an estimator for physical data. Bias summaries are conditional on successful evaluation; failures remain explicit.')
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
