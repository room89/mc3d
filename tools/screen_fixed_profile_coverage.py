"""Frozen fixed-coordinate mean and group-level interval coverage screen."""
import argparse
import json
from pathlib import Path
import numpy as np
from control_front_shared_noise import ar1, loading, profile
from shock_validation import digest

PROTOCOL=Path(__file__).resolve().parents[1]/'validation/stage9b/FIXED_PROFILE_PROTOCOL.json'


def truth(x, sd, order):
    nodes, weights=np.polynomial.hermite.hermgauss(order)
    return weights@profile(x,np.sqrt(2)*sd*nodes)/np.sqrt(np.pi)


def estimate(raw, groups, runs_per_group, critical):
    a=np.asarray(raw,float)
    if a.ndim!=3 or a.shape[0]!=groups*runs_per_group or groups<2 or not np.isfinite(a).all():
        raise ValueError('Finite complete independent groups required')
    means=a.mean(axis=1).reshape(groups,runs_per_group,a.shape[-1]).mean(axis=1)
    center=means.mean(axis=0); se=means.std(axis=0,ddof=1)/np.sqrt(groups)
    return center,critical*se


def screen():
    p=json.loads(PROTOCOL.read_text());ph=digest(PROTOCOL)
    x=np.linspace(*p['x_bounds'],p['nodes']); target=truth(x,p['shift_sd'],p['truth_quadrature_orders'][0])
    discrepancy=float(np.max(abs(target-truth(x,p['shift_sd'],p['truth_quadrature_orders'][1]))))
    if discrepancy>p['truth_max_difference']: raise ValueError('Truth quadrature not converged')
    shape=(p['groups']*p['runs_per_group'],p['windows']); errors={c:[] for c in p['cases']};covered={c:[] for c in p['cases']}; widths={c:[] for c in p['cases']}
    for i in range(p['trials']):
        rng=np.random.default_rng(np.random.SeedSequence([p['seed'],i]))
        shifts=p['shift_sd']*ar1(rng,shape,p['temporal_phi'])
        node=p['node_noise_sd']*ar1(rng,shape+(len(x),),p['temporal_phi'])
        common=p['common_noise_sd']*ar1(rng,shape,p['temporal_phi'])
        offsets=p['run_offset_sd']*rng.normal(size=(shape[0],1,1))
        offsets=offsets+np.repeat(p['group_offset_sd']*rng.normal(size=(p['groups'],1,1)),p['runs_per_group'],axis=0)
        base=profile(x,shifts)+node+offsets
        for c in p['cases']:
            raw=base+common[...,None]*loading(x,c)
            mean,half=estimate(raw,p['groups'],p['runs_per_group'],p['critical_value'])
            error=mean-target
            errors[c].append(error);covered[c].append(abs(error)<=half);widths[c].append(half)
    summary={};details={}
    for c in p['cases']:
        e=np.asarray(errors[c]);h=np.asarray(covered[c]);w=np.asarray(widths[c]);fraction=h.mean(axis=1)
        def rate(v):
            return dict(value=float(v.mean()),mcse=float(v.std(ddof=1)/np.sqrt(len(v))))
        summary[c]=dict(bias=e.mean(axis=0).tolist(),bias_mcse=(e.std(axis=0,ddof=1)/np.sqrt(p['trials'])).tolist(),
            pointwise_coverage=h.mean(axis=0).tolist(),pointwise_coverage_mcse=(h.std(axis=0,ddof=1)/np.sqrt(p['trials'])).tolist(),
            mean_half_width=w.mean(axis=0).tolist(),mean_node_coverage=rate(fraction),
            whole_grid_coverage_of_pointwise_intervals=rate(h.all(axis=1)),
            selected_nodes={str(x[j]):dict(coverage=rate(h[:,j]),bias=rate(e[:,j])) for j in (len(x)//2,int(np.argmin(abs(x-10))))})
        details[c]=dict(errors=e.tolist(),covered=h.tolist(),half_widths=w.tolist())
    if digest(PROTOCOL)!=ph: raise ValueError('Protocol changed')
    return dict(status='exploratory_fixed_profile_coverage',physical_acceptance='incomplete',protocol=p,protocol_sha256=ph,
        numpy_version=np.__version__,x=x.tolist(),truth=target.tolist(),truth_quadrature_difference=discrepancy,
        summaries=summary,trials=details,sources={n:digest(Path(__file__).with_name(n)) for n in ('screen_fixed_profile_coverage.py','control_front_shared_noise.py')},
        limitations=['Illustrative Gaussian laws only; no physical calibration or new physical data.',
            'Eight independent groups, each with three runs sharing a persistent group offset; whole run windows stay together.',
            'Intervals are pointwise nominal 95 percent; whole-grid coverage is reported but is not targeted at 95 percent.',
            'Coverage precision uses independent trials, not spatial nodes or windows as replicates.',
            'The target is a shift-averaged laboratory-frame profile, not an intrinsic shock shape or its gradient thickness.',
            'No training is needed for the frozen arithmetic mean; no model selection or tuning on trial outcomes.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=screen()
    with a.output.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
