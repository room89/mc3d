"""Fixed synthetic shared-field control; not a physical noise calibration."""
import argparse
import json
from pathlib import Path
import numpy as np
from shock_analysis import monotone, crossing
from shock_validation import digest
from diagnose_plateau_front_balance import correlation

PROTOCOL = dict(seed=20261114, trials=32, runs=3, windows=64, nodes=560,
    x_bounds=[-13.975,13.975], shift_sd=.08, node_noise_sd=.003,
    common_noise_sd=.003, temporal_phi=.7,
    cases=['independent_nodes','same_sign_common','opposite_sign_common'])


def ar1(rng, shape, phi):
    if len(shape)<2 or not 0<=phi<1: raise ValueError('Stationary run/window array required')
    a=rng.normal(size=shape)
    for i in range(1,shape[1]): a[:,i]=phi*a[:,i-1]+np.sqrt(1-phi*phi)*a[:,i]
    return a


def profile(x, shifts):
    return .5*(1+np.tanh((np.asarray(x)-np.asarray(shifts)[...,None])/2))


def fronts(x, values):
    return np.array([[crossing(x,monotone(w),.5) for w in run] for run in values])


def loading(x, case):
    if case=='independent_nodes': return np.zeros_like(x)
    if case=='same_sign_common': return np.ones_like(x)
    if case=='opposite_sign_common': return np.tanh((x-4)/2)
    raise ValueError('Unknown fixed case')


def trial(seed):
    p=PROTOCOL; rng=np.random.default_rng(seed); shape=(p['runs'],p['windows'])
    x=np.linspace(*p['x_bounds'],p['nodes']); mask=(x>=8)&(x<=12)
    shift=p['shift_sd']*ar1(rng,shape,p['temporal_phi']); truth=profile(x,shift)
    # Matched independent measurement replicate at the SAME true shifts.
    node=[p['node_noise_sd']*ar1(rng,shape+(len(x),),p['temporal_phi']) for _ in range(2)]
    common=[p['common_noise_sd']*ar1(rng,shape,p['temporal_phi']) for _ in range(2)]
    result={}
    for case in p['cases']:
        observed=[truth+node[i]+common[i][...,None]*loading(x,case) for i in range(2)]
        estimated=[fronts(x,v) for v in observed]
        # Oracle residual isolates measurement noise, including finite plateau tail.
        noise=(observed[0]-truth)[:,:,mask].mean(axis=2).mean(axis=0)
        plateau=(observed[0][:,:,mask].mean(axis=2)-profile(x,0)[mask].mean()).mean(axis=0)
        true=shift.mean(axis=0); same=estimated[0].mean(axis=0); independent=estimated[1].mean(axis=0)
        result[case]=dict(oracle_noise_same_front=correlation(noise,same),
            oracle_noise_independent_front=correlation(noise,independent),
            oracle_noise_true_front=correlation(noise,true),
            oracle_noise_front_error=correlation(noise,same-true),
            fixed_profile_residual_same_front=correlation(plateau,same),
            front_rmse=float(np.sqrt(np.mean((estimated[0]-shift)**2))))
    return result


def summarize(rows):
    a=np.asarray(rows,float)
    if a.ndim!=1 or len(a)<2 or not np.isfinite(a).all(): raise ValueError('Finite replicate statistics required')
    return dict(mean=float(a.mean()),mcse=float(a.std(ddof=1)/np.sqrt(len(a))),
        positive_fraction=float(np.mean(a>0)))


def run():
    rows=[trial(np.random.SeedSequence([PROTOCOL['seed'],i])) for i in range(PROTOCOL['trials'])]
    summaries={}
    for case in PROTOCOL['cases']:
        summaries[case]={k:summarize([r[case][k] for r in rows]) for k in rows[0][case]}
        summaries[case]['paired_same_minus_independent']=summarize([
            r[case]['oracle_noise_same_front']-r[case]['oracle_noise_independent_front'] for r in rows])
    return dict(status='synthetic_shared_field_control_only',physical_acceptance='incomplete',
        protocol=PROTOCOL, numpy_version=np.__version__, trials=rows, summaries=summaries,
        sources={n:digest(Path(__file__).with_name(n)) for n in ('control_front_shared_noise.py','shock_analysis.py','diagnose_plateau_front_balance.py')},
        limitations=['Fixed illustrative Gaussian AR(1) laws, not fitted to physical banks; no physical data read.',
            'Oracle residual subtracts the exact shifted profile; fixed-profile residual is a sensitivity check, not the physical fitted mean.',
            'Independent measurement replicate shares true shifts but not noise; this is not automatically available in physical runs.',
            'Monte Carlo standard errors describe independent synthetic trials, not physical confidence coverage.',
            'No selection among noise cases, significance gate, correction fit or physical acceptance claim.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    r=run()
    with a.output.open('x') as f: json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
