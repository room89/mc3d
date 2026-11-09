"""Fixed synthetic qualification of cross-fitted alignment; no physical runs."""
import argparse,json
from pathlib import Path
import numpy as np
from shock_alignment import align
from shock_mean_profile import GRID,aligned,measure
from shock_analysis import crossing,monotone
from validate_shock_estimator import shape,truth

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();cases=[]
 for kind in ('logistic','gompertz'):
  for noise in (.015,.03):
   for phi in (0.,.8):
    for spatial in (0.,.2):
     runs=[];errors=[];count=128 if noise==.015 and spatial==0 else 64
     for repeat in range(count):
      rng=np.random.default_rng(900000+repeat);x=np.linspace(-14,14,281);trials=[];exact=[];old_front=[]
      def sample():
       v=rng.normal(size=len(x))
       if spatial:
        z=np.arange(-6,7)*.1;k=np.exp(-.5*(z/spatial)**2);k/=np.sqrt(sum(k*k));v=np.convolve(np.pad(v,6,mode='reflect'),k,'valid')
       return v*noise
      for seed in range(3):
       windows=[];previous=sample()
       for wi in range(64):
        previous=phi*previous+np.sqrt(1-phi*phi)*sample();d=rng.uniform(-1,1);y=shape(x-d,kind)+previous;c=crossing(x,monotone(y),.5)
        exact.append(d if kind=='logistic' else d-np.log(np.log(2)));old_front.append(c)
        windows.append(dict(x_aligned=(x-c).tolist(),front=c,density_normalized=y.tolist()))
       trials.append(windows)
      try:
       data,diagnostics=align(trials);new_front=np.array([w['front'] for s in diagnostics for w in s['windows']]);exact=np.array(exact)
       new_error=new_front-exact;old_error=np.array(old_front)-exact
       oracle=np.array([[np.interp(GRID+exact[s*64+i],x,w['density_normalized']) for i,w in enumerate(run)] for s,run in enumerate(trials)])
       runs.append(dict(new=(measure(data.mean((0,1)))[0]/truth(kind)-1).tolist(),old=(measure(aligned(trials).mean((0,1)))[0]/truth(kind)-1).tolist(),oracle=(measure(oracle.mean((0,1)))[0]/truth(kind)-1).tolist(),new_rmse=float(np.std(new_error)),old_rmse=float(np.std(old_error))))
      except ValueError as e:errors.append(dict(repeat=repeat,error=str(e)))
     values=np.array([v['new'] for v in runs]);bounds=abs(values.mean(0))+2*values.std(0,ddof=1)/np.sqrt(len(values));ratio=np.mean([v['new_rmse'] for v in runs])/np.mean([v['old_rmse'] for v in runs]);limit=.75 if noise==.015 and spatial==0 else 1.
     row=dict(shape=kind,noise=noise,phi=phi,spatial=spatial,repeats=count,errors=errors,bias_plus_two_se=bounds.tolist(),rmse_ratio=float(ratio),qualified=bool(not errors and np.all(bounds<=[.005,.015]) and ratio<=limit),runs=runs);cases.append(row);print({k:v for k,v in row.items() if k not in ('runs','errors')},flush=True)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(dict(status='qualified' if all(c['qualified'] for c in cases) else 'alignment_not_qualified',cases=cases),indent=2)+'\n')
if __name__=='__main__':main()
