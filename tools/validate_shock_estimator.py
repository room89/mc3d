"""Synthetic calibration of the experimental estimator, with fixed random seeds."""
import argparse,json,math
from pathlib import Path
import numpy as np
from shock_analysis import monotone,crossing
from shock_mean_profile import GRID,measure,estimate


def shape(x,kind,scale=1):
    z=np.asarray(x)/scale
    if kind=='logistic':return 1/(1+np.exp(-z))
    return np.exp(-np.exp(-z))


def truth(kind,scale=1):
    if kind=='logistic':return np.array([2*math.log(9)*scale,4*scale])
    return np.array([math.log(math.log(.1)/math.log(.9))*scale,math.e*scale])


def synthetic(kind,seed,noise=.015,phi=0.,count=64,dx=.1,oracle=False):
    rng=np.random.default_rng(seed);x=np.arange(-14,14+dx/2,dx);trials=[];errors=[]
    for _ in range(3):
        rows=[];previous=rng.normal(0,noise,len(x))
        for _ in range(count):
            previous=phi*previous+math.sqrt(1-phi*phi)*rng.normal(0,noise,len(x))
            shift=rng.uniform(-1,1);y=shape(x-shift,kind)+previous
            front=crossing(x,monotone(y),.5)
            exact=shift if kind=='logistic' else shift-math.log(math.log(2))
            errors.append(front-exact)
            rows.append(np.interp(GRID,x-(exact if oracle else front),y))
        trials.append(rows)
    return np.asarray(trials),float(np.sqrt(np.mean(np.array(errors)**2)))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--replicates',type=int,default=16);a=p.parse_args()
    if a.replicates<2:p.error('At least two replicates')
    result={'numpy_version':np.__version__,'noiseless_relative_tolerance':.01,'status':'experimental_calibration_only','bootstrap_repeats':128,'replicates_per_case':a.replicates,'noiseless':[],'noisy':[]}
    for kind in ('logistic','gompertz'):
        for scale in (.25,1.,2.):
            for dx in (.05,.1):
                x=np.arange(-10,10+dx/2,dx);m=measure(shape(x,kind,scale),x)[0];exact=truth(kind,scale)
                result['noiseless'].append(dict(shape=kind,scale=scale,dx=dx,estimate=m.tolist(),truth=exact.tolist(),relative_error=(m/exact-1).tolist(),resolved=bool(np.max(abs(m/exact-1))<.01)))
    for kind in ('logistic','gompertz'):
        for phi in (0.,.8):
            rows=[];exact=truth(kind)
            for i in range(a.replicates):
                data,alignment_rms=synthetic(kind,10000+i,phi=phi)
                fit=estimate(data,repeats=128,seed=20000+i)
                m=np.array([fit[k]['mean'] for k in ('width_10_90','gradient_thickness')]);u=np.array([fit[k]['statistical_radius'] for k in ('width_10_90','gradient_thickness')]);proc=np.array([fit[k]['processing_sensitivity'] for k in ('width_10_90','gradient_thickness')])
                oracle_data,_=synthetic(kind,10000+i,phi=phi,oracle=True)
                oracle_m=measure(oracle_data.mean(axis=(0,1)))[0]
                rows.append(dict(oracle_alignment_estimate=oracle_m.tolist(),estimate=m.tolist(),statistical_radius=u.tolist(),processing_sensitivity=proc.tolist(),alignment_rms=alignment_rms,statistical_covers=(abs(m-exact)<=u).tolist(),with_processing_covers=(abs(m-exact)<=u+proc).tolist()))
            result['noisy'].append(dict(shape=kind,phi=phi,noise=.015,truth=exact.tolist(),mean_relative_bias=(np.mean([v['estimate'] for v in rows],axis=0)/exact-1).tolist(),oracle_alignment_relative_bias=(np.mean([v['oracle_alignment_estimate'] for v in rows],axis=0)/exact-1).tolist(),statistical_coverage=np.mean([v['statistical_covers'] for v in rows],axis=0).tolist(),with_processing_coverage=np.mean([v['with_processing_covers'] for v in rows],axis=0).tolist(),runs=rows))
            print(kind,phi,result['noisy'][-1]['mean_relative_bias'],flush=True)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    if a.output.exists():raise FileExistsError(a.output)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
