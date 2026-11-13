"""Exploratory training-only plateau mean and window-index drift model.

Frozen model: intercept + spatial linear/quadratic terms + linear window index.
No per-run intercepts, held-out fitting, adaptive region selection or noise-law
qualification. Leave-group-out residuals retain offsets in excluded runs.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from diagnose_physical_plateaus import group_summary
from shock_empirical_fluctuations import load_bank
from shock_analysis import monotone,crossing
from shock_validation import digest


def design(x,windows=64):
    x=np.asarray(x,float)
    if x.ndim!=1 or len(x)<4 or not np.isfinite(x).all() or np.any(np.diff(x)<=0) or windows!=64:
        raise ValueError('Increasing finite spatial nodes and 64 ordered windows required')
    s=2*(x-x[0])/(x[-1]-x[0])-1; t=np.linspace(-1,1,windows)
    return np.stack([np.ones((windows,len(x))),np.broadcast_to(s,(windows,len(x))),
                     np.broadcast_to(s*s,(windows,len(x))),np.broadcast_to(t[:,None],(windows,len(x)))],axis=-1)


def coefficients(raw,basis):
    y=np.asarray(raw,float)
    if y.ndim!=3 or y.shape[1:]!=basis.shape[:2] or len(y)<1 or not np.isfinite(y).all():
        raise ValueError('Finite independent run blocks matching the design required')
    result,_,rank,_=np.linalg.lstsq(basis.reshape(-1,4),y.mean(axis=0).ravel(),rcond=None)
    if rank!=4:raise ValueError('Unidentifiable mean model')
    return result


def residual_summary(residual,sources):
    series=residual.mean(axis=2);half=residual.shape[2]//2
    out=dict(mean=group_summary(series.mean(axis=1),sources),
             time_change=group_summary(series[:,32:].mean(axis=1)-series[:,:32].mean(axis=1),sources),
             spatial_change=group_summary(residual[:,:,half:].mean(axis=(1,2))-residual[:,:,:half].mean(axis=(1,2)),sources),
             per_run_mean=series.mean(axis=1).tolist(),node_mse=float(np.mean(residual**2)),phases={})
    for name,sl in [('full',slice(None)),('first',slice(0,32)),('last',slice(32,64))]:
        a=series[:,sl]
        out['phases'][name]=dict(second_moment=float(np.mean(a*a)),lag_one_product=float(np.mean(a[:,:-1]*a[:,1:])))
    return out


def fit_region(training,heldout,x,train_sources,held_sources):
    training=np.asarray(training,float);heldout=np.asarray(heldout,float)
    if heldout.ndim!=3 or heldout.shape[1:]!=(64,len(x)) or len(heldout)<1 or not np.isfinite(heldout).all():raise ValueError('Finite held-out blocks matching the design required')
    if len({s['seed'] for s in train_sources})!=len(train_sources) or len({s['seed'] for s in held_sources})!=len(held_sources):raise ValueError('Distinct run identities required')
    train_groups=sorted({s['group'] for s in train_sources})
    if len(train_groups)<3 or len(train_sources)!=len(training) or len(held_sources)!=len(heldout):raise ValueError('Matching source identities and at least three training groups required')
    if set(train_groups)&{s['group'] for s in held_sources} or {s['seed'] for s in train_sources}&{s['seed'] for s in held_sources}:raise ValueError('Training/held-out overlap')
    basis=design(x);beta=coefficients(training,basis);prediction=basis@beta
    # All historical groups have three runs; reject unequal groups instead of
    # silently changing the replicate weighting used by the jackknife.
    counts=[sum(s['group']==g for s in train_sources) for g in train_groups]
    if len(set(counts))!=1:raise ValueError('Equal-size training groups required')
    crossfit=np.empty_like(training);loo=[]
    for group in train_groups:
        exclude=np.array([s['group']==group for s in train_sources])
        b=coefficients(training[~exclude],basis);loo.append(b)
        crossfit[exclude]=training[exclude]-basis@b
    loo=np.array(loo);centered=loo-loo.mean(axis=0)
    covariance=(len(loo)-1)/len(loo)*(centered.T@centered)
    mean_se=np.sqrt(np.maximum(0,np.einsum('tni,ij,tnj->tn',basis,covariance,basis)))
    baseline=float(training.mean());full=training-prediction;held=heldout-prediction
    summaries=dict(training_in_sample=residual_summary(full,train_sources),
                   training_leave_group_out=residual_summary(crossfit,train_sources),
                   heldout=residual_summary(held,held_sources),
                   baseline_training=residual_summary(training-baseline,train_sources),
                   baseline_heldout=residual_summary(heldout-baseline,held_sources))
    # Subtracting a shared prediction must preserve every pair of run offsets.
    original_delta=heldout.mean()-training.mean();residual_delta=held.mean()-full.mean()
    if not np.isclose(original_delta,residual_delta,atol=1e-13,rtol=1e-10):raise ValueError('Mean transfer was altered')
    contrasts=dict(mean=basis.mean(axis=(0,1)),
                   time_change=basis[32:].mean(axis=(0,1))-basis[:32].mean(axis=(0,1)),
                   spatial_change=basis[:,len(x)//2:].mean(axis=(0,1))-basis[:,:len(x)//2].mean(axis=(0,1)))
    transfer_uncertainty={}
    for name,contrast in contrasts.items():
        statistic=summaries['heldout'][name];variance=float(contrast@covariance@contrast)
        held_se=statistic['group_mcse']
        combined=float(np.sqrt(held_se**2+max(0.,variance))) if held_se is not None else None
        transfer_uncertainty[name]=dict(mean=statistic['mean'],training_prediction_se=float(np.sqrt(max(0.,variance))),
            heldout_group_mcse=held_se,combined_descriptive_se=combined,
            descriptive_interval=[statistic['mean']-1.96*combined,statistic['mean']+1.96*combined] if combined is not None else None)
    ratios={}
    for phase in ('full','first','last'):
        denominator=summaries['training_leave_group_out']['phases'][phase]['second_moment']
        ratios[phase]=summaries['heldout']['phases'][phase]['second_moment']/denominator if denominator>0 else None
    return dict(coefficients=beta.tolist(),coefficient_order=['intercept','space','space_squared','window_index'],
                leave_group_out_coefficients=loo.tolist(),jackknife_coefficient_covariance=covariance.tolist(),
                predicted_mean_se_range=[float(mean_se.min()),float(mean_se.max())],
                prediction=prediction.tolist(),training_baseline=baseline,summaries=summaries,
                transfer_with_training_mean_uncertainty=transfer_uncertainty,
                heldout_minus_training_mean=float(residual_delta),heldout_to_crossfit_second_moment_ratio=ratios)


def diagnose(manifest_path):
    manifest_path=Path(manifest_path);mh=digest(manifest_path);result={};hashes={}
    for code in ('mc3d','sparta'):
        a=load_bank(manifest_path,code+'-training');b=load_bank(manifest_path,code+'-heldout')
        if not np.array_equal(a['x'],b['x']):raise ValueError('Different grids')
        train=a['fluctuations']+a['mean_profile'];held=b['fluctuations']+b['mean_profile'];x=a['x']
        center=crossing(x,monotone(train.mean(axis=(0,1))),.5);regions={}
        for name,lo,hi in [('upstream',-12.,-8.),('downstream',8.,12.)]:
            if center+lo<x[0]+1 or center+hi>x[-1]-1:raise ValueError('Fixed region lacks clearance')
            mask=(x>=center+lo)&(x<=center+hi)
            regions[name]=fit_region(train[:,:,mask],held[:,:,mask],x[mask],a['sources'],b['sources'])
            regions[name]['bounds']=[float(x[mask][0]),float(x[mask][-1])]
        result[code]=dict(training_center=center,regions=regions)
        for bank in (a,b):hashes[bank['name']]=dict(sha256=bank['bank_sha256'],sources=bank['sources'])
    manifest=json.loads(manifest_path.read_text())
    if digest(manifest_path)!=mh or any(digest(manifest_path.parent/manifest['banks'][n]['file'])!=v['sha256'] for n,v in hashes.items()):raise ValueError('Evidence changed')
    return dict(status='exploratory_mean_transfer_only',physical_acceptance='incomplete',manifest_sha256=mh,banks=hashes,codes=result,
                sources={n:digest(Path(__file__).with_name(n)) for n in ('fit_physical_plateau_mean.py','diagnose_physical_plateaus.py','shock_empirical_fluctuations.py','shock_analysis.py')},
                limitations=['Previously inspected historical splits, not independent qualification.',
                             'Fixed quadratic spatial mean and linear window-index drift; no held-out tuning or per-run centering.',
                             'Leave-group-out coefficient fits exclude each group; regions remain fixed from the full training sample.',
                             'Jackknife uncertainty uses eight training groups and is diagnostic, not a calibrated confidence band.',
                             'Residual second moments retain model error and run offsets; they are not a separable noise covariance.',
                             'Window index is not a fitted physical time law; stationarity and noise coverage remain unqualified.'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=diagnose(args.manifest)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    for code,c in result['codes'].items():
        for name,r in c['regions'].items():
            s=r['summaries'];print(code,name,'beta',r['coefficients'],'held-out time',s['heldout']['time_change']['mean'],'space',s['heldout']['spatial_change']['mean'],'moment ratios',r['heldout_to_crossfit_second_moment_ratio'])
