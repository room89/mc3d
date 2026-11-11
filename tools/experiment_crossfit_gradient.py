"""Paired pilot of independent gradient selection; not production acceptance."""
import argparse,json,math
from pathlib import Path
import numpy as np
from shock_alignment import GRID,align
from shock_analysis import monotone,crossing
from shock_mean_profile import aligned,measure
from shock_crossfit_gradient import derivative,estimate_gradient
from validate_shock_estimator import shape,truth


def fixture(kind,seed,phi,spatial,noise=.03):
    rng=np.random.default_rng(seed);x=np.linspace(-14,14,281);trials=[];exact=[]
    def sample():
        v=rng.normal(size=len(x))
        if spatial:
            z=np.arange(-6,7)*.1;k=np.exp(-.5*(z/spatial)**2);k/=np.sqrt(sum(k*k))
            v=np.convolve(np.pad(v,6,mode='reflect'),k,'valid')
        return noise*v
    for _ in range(3):
        run=[];previous=sample()
        for _ in range(64):
            previous=phi*previous+math.sqrt(1-phi*phi)*sample();shift=rng.uniform(-1,1)
            if kind=='flat_peak':
                z=np.clip((x-shift)/2.5,-1,1);clean=.5+5/8*(z-z**5/5)
            else:clean=shape(x-shift,kind)
            y=clean+previous;c=crossing(x,monotone(y),.5)
            exact.append(shift-math.log(math.log(2)) if kind=='gompertz' else shift)
            run.append(dict(x_aligned=(x-c).tolist(),front=c,density_normalized=y.tolist()))
        trials.append(run)
    oracle=np.asarray([[np.interp(GRID+exact[s*64+i],x,w['density_normalized']) for i,w in enumerate(run)] for s,run in enumerate(trials)])
    return trials,oracle


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,default=32);a=p.parse_args()
    if a.repeats<2:p.error('Need at least two repetitions')
    if a.output.exists():raise FileExistsError(a.output)
    cases=[]
    for kind in ('logistic','gompertz','flat_peak'):
        exact=4. if kind=='flat_peak' else truth(kind)[1]
        for phi in (0.,.8):
            for spatial in (0.,.2):
                rows=[];errors=[]
                for i in range(a.repeats):
                    try:
                        trials,oracle=fixture(kind,1200000+i,phi,spatial);data,diagnostics=align(trials)
                        independent=estimate_gradient(data,diagnostics)
                        ox,od=derivative(oracle.mean((0,1)))
                        true_position=math.log(math.log(2)) if kind=='gompertz' else 0.
                        rows.append(dict(crossfit=independent['gradient_thickness']/exact-1,
                            template_max=measure(data.mean((0,1)))[0,1]/exact-1,
                            crossing50_max=measure(aligned(trials).mean((0,1)))[0,1]/exact-1,
                            oracle_max=measure(oracle.mean((0,1)))[0,1]/exact-1,
                            oracle_known_peak=1/np.interp(true_position,ox,od)/exact-1,
                            folds=independent['folds']))
                    except ValueError as e:errors.append(dict(repetition=i,error=str(e)))
                summary={}
                for method in ('crossfit','template_max','crossing50_max','oracle_max','oracle_known_peak'):
                    v=np.array([r[method] for r in rows]);sd=float(v.std(ddof=1)) if len(v)>1 else None
                    summary[method]=dict(mean_bias=float(v.mean()) if len(v) else None,sd=sd,bias_plus_two_se=float(abs(v.mean())+2*sd/math.sqrt(len(v))) if sd is not None else None)
                q=summary['crossfit'];old=summary['template_max']
                passed=not errors and q['bias_plus_two_se']<=.015 and q['sd']<=1.25*old['sd']
                c=dict(shape=kind,phi=phi,spatial=spatial,noise=.03,repeats=a.repeats,passed=bool(passed),summary=summary,errors=errors,runs=rows);cases.append(c)
                print(kind,phi,spatial,passed,q,flush=True)
    result=dict(status='pilot_passed' if all(c['passed'] for c in cases) else 'pilot_not_qualified',seed_start=1200000,criteria=dict(bias_plus_two_se_max=.015,sd_ratio_max=1.25),cases=cases)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
