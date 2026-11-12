"""Separate mean-plateau and boundary diagnostics; original stationarity gate unchanged."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from audit_shock_validation import gate
from reanalyze_shock_mean import sparta_density


def mean_interval(values,repeats=2048):
    values=np.asarray(values,float)
    if values.ndim!=3 or values.shape[0]!=3 or values.shape[1]%16 or values.shape[1]<16 or values.shape[2]!=2 or not np.isfinite(values).all():
        raise ValueError('Three finite paired plateau series required')
    mean=values.mean((0,1));u=4.303*values.mean(1).std(0,ddof=1)/math.sqrt(3)
    rng=np.random.default_rng(20261109);count=values.shape[1];details={}
    for block in (8,16):
        sampled=[]
        for _ in range(repeats):
            sums=[]
            for run in values:
                starts=rng.integers(0,count,count//block)
                ids=((starts[:,None]+np.arange(block))%count).ravel()
                sums.append(run[ids].mean(0))
            sampled.append(np.mean(sums,axis=0))
        lo,hi=np.quantile(sampled,[.025,.975],axis=0)
        radius=np.maximum(abs(lo-mean),abs(hi-mean));u=np.maximum(u,radius)
        details[str(block)]=radius.tolist()
    return mean,u,details


def analyze(trials):
    errors=[];boundaries=[]
    for run in trials:
        rows=[]
        for w in run:
            y=np.asarray(w['density_normalized']);x=np.asarray(w['x_aligned'])
            n=max(1,len(y)//10)
            rows.append([float(y[:n].mean()),float(y[-n:].mean()-1)])
            boundaries.append(float(min(-x[0],x[-1])))
        errors.append(rows)
    errors=np.asarray(errors)
    mean,u,blocks=mean_interval(errors)
    first,fu,_=mean_interval(errors[:,:errors.shape[1]//2])
    last,lu,_=mean_interval(errors[:,errors.shape[1]//2:])
    return dict(status='diagnostic_only',criterion='Mean plateau and half comparison use original 1% tolerance and half-tolerance precision requirement; not the legacy every-window gate',
                mean_plateaus={k:gate(float(mean[i]),float(u[i]),.01) for i,k in enumerate(('upstream','downstream'))},
                half_changes={k:gate(float(last[i]-first[i]),float(fu[i]+lu[i]),.01) for i,k in enumerate(('upstream','downstream'))},
                block_radii=blocks,minimum_boundary_distance=min(boundaries),boundary_pass=min(boundaries)>=10,
                windows_outside_one_percent=int(np.sum(np.max(abs(errors),axis=2)>.01)),
                total_windows=int(errors.shape[0]*errors.shape[1]))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('series',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    report=json.loads((a.series/'validation.json').read_text())
    seeds=tuple(report.get('seeds',(123,456,789)))
    if len(seeds)!=3 or len(set(seeds))!=3:raise ValueError('Invalid seed identities')
    result=dict(status='diagnostic_only',input_hashes={},variants={})
    for name in list('ABCDEFG')+['sparta-C','sparta-D','sparta-E']:
        if name.startswith('sparta'):
            trials=[sparta_density(a.series/f'{name}-{s}',result['input_hashes']) for s in seeds]
        else:
            path=a.series/f'profiles-{name}.json';raw=path.read_bytes();result['input_hashes'][path.name]=hashlib.sha256(raw).hexdigest();trials=json.loads(raw)
        result['variants'][name]=analyze(trials);print(name,result['variants'][name],flush=True)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
