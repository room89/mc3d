"""Exploratory plateau transfer checks on historical physical profile banks.

Reconstruct original means; use fixed training-selected regions for held-out
runs. Time windows and spatial nodes are not independent replicates. Report
covariance of spatial averages, NOT a fitted separable temporal noise factor.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from shock_empirical_fluctuations import load_bank
from shock_analysis import monotone,crossing
from shock_validation import digest


def group_summary(values,sources):
    groups=sorted({s['group'] for s in sources})
    means=[float(np.mean([v for v,s in zip(values,sources) if s['group']==g])) for g in groups]
    mean=float(np.mean(means));se=float(np.std(means,ddof=1)/np.sqrt(len(means))) if len(means)>1 else None
    return dict(mean=mean,group_count=len(means),group_values=means,group_mcse=se,
                descriptive_interval=[mean-1.96*se,mean+1.96*se] if se is not None else None)


def region_statistics(raw,x,mask,baseline,sources):
    values=raw[:,:,mask];series=values.mean(axis=2);residual=series-baseline
    half=values.shape[2]//2
    phases={}
    for name,sl in [('full',slice(None)),('first',slice(0,32)),('last',slice(32,64))]:
        a=residual[:,sl]
        phases[name]=dict(mean_density=float(series[:,sl].mean()),
                          second_moment_about_training_mean=float(np.mean(a*a)),
                          lag_one_product_about_training_mean=float(np.mean(a[:,:-1]*a[:,1:])))
    return dict(nodes=int(mask.sum()),bounds=[float(x[mask][0]),float(x[mask][-1])],
                mean=group_summary(series.mean(axis=1),sources),
                time_change=group_summary(series[:,32:].mean(axis=1)-series[:,:32].mean(axis=1),sources),
                spatial_change=group_summary(values[:,:,half:].mean(axis=(1,2))-values[:,:,:half].mean(axis=(1,2)),sources),
                per_run_mean=series.mean(axis=1).tolist(),phases=phases)


def analyze(training,heldout):
    x=training['x']
    if not np.array_equal(x,heldout['x']):raise ValueError('Different spatial grids')
    a=training['fluctuations']+training['mean_profile'];b=heldout['fluctuations']+heldout['mean_profile']
    if a.shape[1]!=64 or b.shape[1]!=64:raise ValueError('Exactly 64 ordered windows required')
    if {s['seed'] for s in training['sources']}&{s['seed'] for s in heldout['sources']}:raise ValueError('Overlapping source seeds')
    if {s['group'] for s in training['sources']}&{s['group'] for s in heldout['sources']}:raise ValueError('Overlapping groups')
    center=crossing(x,monotone(a.mean(axis=(0,1))),.5)
    regions={};masks={}
    for name,lo,hi in [('upstream',-12.,-8.),('downstream',8.,12.)]:
        bounds=[center+lo,center+hi]
        if bounds[0]<x[0]+1 or bounds[1]>x[-1]-1:raise ValueError('Fixed region lacks boundary clearance; do not move it')
        mask=(x>=bounds[0])&(x<=bounds[1]);masks[name]=mask
        if mask.sum()<4:raise ValueError('Insufficient plateau nodes')
        baseline=float(a[:,:,mask].mean())
        train=region_statistics(a,x,mask,baseline,training['sources'])
        held=region_statistics(b,x,mask,baseline,heldout['sources'])
        se=float(np.hypot(train['mean']['group_mcse'],held['mean']['group_mcse']))
        delta=held['mean']['mean']-train['mean']['mean']
        ratios={}
        for phase in ('full','first','last'):
            denominator=train['phases'][phase]['second_moment_about_training_mean']
            ratios[phase]=held['phases'][phase]['second_moment_about_training_mean']/denominator if denominator>0 else None
        regions[name]=dict(training_baseline=baseline,training=train,heldout=held,
                           heldout_minus_training_mean=delta,transfer_mean_descriptive_interval=[delta-1.96*se,delta+1.96*se],
                           heldout_training_second_moment_ratio=ratios)
    motion={}
    for name,raw in [('training',a),('heldout',b)]:
        failures=overlaps=0;fronts=[]
        for run in raw:
            for row in run:
                try:
                    fitted=monotone(row);lo=crossing(x,fitted,.1);hi=crossing(x,fitted,.9)
                    if not np.isfinite([lo,hi]).all() or hi<=lo:raise ValueError('Invalid transition')
                    fronts.append([lo,hi])
                    overlaps+=int(x[masks['upstream']][-1]>=lo or x[masks['downstream']][0]<=hi)
                except ValueError:failures+=1
        motion[name]=dict(windows=int(raw.shape[0]*raw.shape[1]),unresolved=failures,region_overlaps=overlaps,
                          transition_envelope=[float(np.min(fronts)),float(np.max(fronts))] if fronts else None)
    return dict(training_center=center,regions=regions,motion=motion)


def diagnose(manifest_path):
    manifest_path=Path(manifest_path);mh=digest(manifest_path)
    result={};banks={}
    for code in ('mc3d','sparta'):
        train=load_bank(manifest_path,code+'-training');held=load_bank(manifest_path,code+'-heldout')
        result[code]=analyze(train,held)
        for bank in (train,held):banks[bank['name']]=dict(sha256=bank['bank_sha256'],sources=bank['sources'])
    manifest=json.loads(manifest_path.read_text())
    if digest(manifest_path)!=mh or any(digest(manifest_path.parent/manifest['banks'][name]['file'])!=bank['sha256'] for name,bank in banks.items()):raise ValueError('Bank evidence changed during diagnosis')
    return dict(status='exploratory_physical_plateau_diagnostic',physical_acceptance='incomplete',
                manifest_sha256=mh,banks=banks,codes=result,
                sources={n:digest(Path(__file__).with_name(n)) for n in ('diagnose_physical_plateaus.py','shock_empirical_fluctuations.py','shock_analysis.py')},
                limitations=['Historical fixed-grid reconstructed profiles; original raw solver files are not re-extracted.',
                             'Existing training/held-out split was previously inspected; no fresh independent validation.',
                             'Fixed training center plus [-12,-8] and [8,12] regions; held-out data never select boundaries.',
                             'Group-level normal intervals are descriptive with only 8/4 groups and no multiplicity adjustment.',
                             'Second moments retain mean drift and run offsets: they are not pure measurement-noise variance.',
                             'Spatial-average covariance is not a separable temporal factor and must not be plugged into the profile fitter.',
                             'No plateau pass, noise-model qualification or physical acceptance is inferred.'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=diagnose(args.manifest)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    for code,r in result['codes'].items():
        print(code,'motion',r['motion'])
        for name,v in r['regions'].items():print(name,'mean transfer',v['heldout_minus_training_mean'],'moment ratios',v['heldout_training_second_moment_ratio'],'time changes',v['training']['time_change']['mean'],v['heldout']['time_change']['mean'])
