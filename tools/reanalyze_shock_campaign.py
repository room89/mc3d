"""Read-only independent-group pilot analysis; never declares physical acceptance."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
from analyze_shock_validation import profile
from audit_shock_validation import gate
from reanalyze_shock_mean import sparta_density,compare
from shock_alignment import initial
from shock_analysis import monotone,crossing
from shock_ensemble_estimator import group_point,estimate,T95,GROUPS
from shock_independent_campaign import verify_baseline
from shock_validation import digest,save
from shock_validation import configuration
from physical_config import prepare


def merged_pairs(trials):
    result=[]
    for run in trials:
        if not run or len(run)%2:raise ValueError('Nonempty even window count required')
        rows=[]
        for a,b in zip(run[::2],run[1::2]):
            x,y,_=initial(a);other,z,_=initial(b)
            if len(x)!=len(other) or not np.allclose(x,other,rtol=1e-14,atol=1e-12):
                raise ValueError('Paired windows have different physical grids')
            mean=(y+z)/2;front=crossing(x,monotone(mean),.5)
            rows.append(dict(front=front,x_aligned=(x-front).tolist(),density_normalized=mean.tolist()))
        result.append(rows)
    return result


def group_edges(trials):
    values=[];boundaries=[]
    for run in trials:
        rows=[]
        for w in run:
            y=np.asarray(w['density_normalized']);n=max(1,len(y)//10)
            rows.append([float(y[:n].mean()),float(y[-n:].mean()-1)])
            boundaries.append(min(-w['x_aligned'][0],w['x_aligned'][-1]))
        values.append(rows)
    a=np.asarray(values)
    return dict(full=a.mean((0,1)).tolist(),first=a[:,:32].mean((0,1)).tolist(),
        last=a[:,32:].mean((0,1)).tolist(),boundary_min=float(min(boundaries)),
        window_outliers=int(np.sum(np.max(abs(a),axis=2)>.01)),window_count=int(a.shape[0]*a.shape[1]))


def plateau_checks(groups):
    if len(groups)!=GROUPS:raise ValueError('Twelve independent groups required')
    summaries={}
    for phase in ('full','first','last'):
        a=np.asarray([g[phase] for g in groups]);mean=a.mean(0);u=T95*a.std(0,ddof=1)/np.sqrt(GROUPS)
        summaries[phase]=[dict(mean=float(m),radius=float(r)) for m,r in zip(mean,u)]
    return dict(method='Nominal Student-t intervals on 12 independent group mean edge offsets',
        mean=[gate(v['mean'],v['radius'],.01) for v in summaries['full']],
        half_change=[gate(a['mean']-b['mean'],a['radius']+b['radius'],.01)
                     for a,b in zip(summaries['first'],summaries['last'])],
        boundary_min=min(g['boundary_min'] for g in groups),
        boundaries_pass=all(g['boundary_min']>=10 for g in groups),
        window_outliers=sum(g['window_outliers'] for g in groups),
        window_count=sum(g['window_count'] for g in groups))


def check_coverage(path,windows):
    report=json.loads(path.read_text());p=report['protocol']
    if not all(p.get(k) for k in ('center_profile','width_band','interpolation_sensitivity')) or p.get('windows')!=windows:
        raise ValueError('Coverage protocol does not match estimator')
    for name in ('shock_ensemble_estimator.py','shock_alignment.py','shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py'):
        if p['sources'][name]!=digest(Path(__file__).with_name(name)):
            raise ValueError('Coverage uses a different estimator')
    required={(k,t,s,n) for k in ('logistic','gompertz') for t in (0.,.8) for s in (0.,.2) for n in (.005,.015,.03)}
    cases=report['cases']
    if len(cases)!=24 or {(c['shape'],c['phi'],c['spatial'],c['noise']) for c in cases}!=required or any(c.get('trial_count',0)<128 or len(c.get('runs',[]))!=c['trial_count'] for c in cases):
        raise ValueError('Complete coverage matrix required')
    return dict(sha256=digest(path),status=report['status'],windows=windows,
        minimum_total_metric_coverage=np.min([c['with_processing_coverage'] for c in cases],axis=0).tolist(),
        minimum_density_coverage=min(c['density_coverage'] for c in cases),
        rejected_trials=sum(len(c['errors']) for c in cases),
        limitation='Completion and matching code are checked; this does not certify nominal coverage')


def verify_mc_inputs(folder,variant,seed):
    manifest=json.loads((folder/'output/manifest.json').read_text())
    if manifest['seed']!=seed:raise ValueError('Actual solver seed differs')
    for role,file in (('config','solver.json'),('physical_input','input_si.json'),('units','units.json')):
        entries=[e for e in manifest['inputs'] if e['role']==role]
        if len(entries)!=1 or entries[0]['sha256']!=digest(folder/file):raise ValueError('Solver input hash differs')
    data,reference=configuration(variant,seed);expected,units=prepare(data)
    actual=json.loads((folder/'solver.json').read_text());actual.pop('units_file')
    if actual!=expected:raise ValueError('Normalized solver configuration differs')
    actual_units=json.loads((folder/'units.json').read_text());actual_units.pop('source_file')
    if actual_units!=units:raise ValueError('Normalization scales differ')
    if json.loads((folder/'reference.json').read_text())!=reference:raise ValueError('Physical reference differs')


def process_group(task):
    name,folder,seeds,variants,binary_hash,sparta_hash,parts,protocol=task
    folder=Path(folder);raw=(folder/'validation.json').read_bytes();report=json.loads(raw)
    if list(report.get('seeds',[123,456,789]))!=seeds:raise ValueError('Group seeds changed')
    verify_baseline(folder,report,variants,seeds,binary_hash,sparta_hash)
    hashes={'validation.json':digest(folder/'validation.json')}
    # Hash the verified manifests and normalization metadata, not only derived profiles.
    for v in variants:
        for s in seeds:
            verify_mc_inputs(folder/f'{v}-{s}',v,s)
            for file in ('input_si.json','solver.json','units.json','reference.json','output/manifest.json'):
                relative=f'{v}-{s}/{file}';hashes[relative]=digest(folder/relative)
            if v in ('C','D','E'):
                for file,sha in report['sparta_runs'][f'{v}-{s}']['hashes'].items():hashes[f'sparta-{v}-{s}/{file}']=sha
    output={}
    for v in variants:
        for code in ('mc3d','sparta') if v in ('C','D','E') else ('mc3d',):
            key=v if code=='mc3d' else 'sparta-'+v
            cache=Path(parts)/f'{name}-{key}.json'
            if cache.exists():
                saved=json.loads(cache.read_text())
                if saved['protocol']!=protocol or saved['hashes']!=hashes:raise ValueError('Group checkpoint changed')
                output[key]=saved['result'];continue
            # Rebuild MC3D density windows directly from the verified raw statistics.
            trials=[profile(folder/f'{v}-{s}') for s in seeds] if code=='mc3d' else [sparta_density(folder/f'sparta-{v}-{s}',{}) for s in seeds]
            phases=dict(full=trials,first=[r[:32] for r in trials],last=[r[32:] for r in trials],doubled=merged_pairs(trials))
            result={p:group_point(t,center_profile=True,interpolation_sensitivity=True) for p,t in phases.items()}
            result['edges']=group_edges(trials)
            save(cache,dict(protocol=protocol,hashes=hashes,result=result));output[key]=result
    if (folder/'validation.json').read_bytes()!=raw:raise ValueError('Source report changed during analysis')
    return name,output,hashes


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('campaign',type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--coverage',required=True,type=Path);p.add_argument('--half-coverage',required=True,type=Path)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.jobs<1:p.error('Positive job count required')
    if a.output.exists():raise FileExistsError(a.output)
    source=a.campaign/'campaign.json';campaign=json.loads(source.read_text())
    if campaign['status']!='execution_complete':raise ValueError('Completed execution campaign required')
    groups=campaign['groups'];cp=campaign['protocol'];variants=cp['variants']
    if set(groups)!={f'g{i:02d}' for i in range(GROUPS)} or any(groups[f'g{i:02d}']['seeds']!=cp['seeds'][i] for i in range(GROUPS)):
        raise ValueError('Group identities changed')
    seeds=[s for g in groups.values() for s in g['seeds']]
    folders=[(a.campaign/g['path']).resolve() for g in groups.values()]
    if len(seeds)!=36 or len(set(seeds))!=36 or len(set(folders))!=12:raise ValueError('Groups are not independent')
    if digest(a.campaign/groups['g00']['path']/'validation.json')!=cp['baseline_report_sha256']:
        raise ValueError('Original baseline report changed')
    coverage=[check_coverage(a.coverage,64),check_coverage(a.half_coverage,32)]
    names=('reanalyze_shock_campaign.py','shock_ensemble_estimator.py','shock_alignment.py',
        'shock_crossfit_gradient.py','shock_mean_profile.py','shock_analysis.py','analyze_shock_validation.py','reanalyze_shock_mean.py',
        'shock_independent_campaign.py','shock_validation.py','normal_shock.py','physical_config.py','argon_vhs.py')
    protocol=dict(campaign_sha256=digest(source),sources={n:digest(Path(__file__).with_name(n)) for n in names},coverage=coverage,
        bootstrap_repeats=1024,resampling_seed=20261109)
    parts=a.output.with_suffix('.parts')
    if parts.exists() and not a.resume:raise FileExistsError(parts)
    parts.mkdir(parents=True,exist_ok=True);results={};input_hashes={}
    tasks=[(name,str((a.campaign/g['path']).resolve()),g['seeds'],variants,cp['binary_sha256'],cp['sparta_sha256'],str(parts),protocol) for name,g in sorted(groups.items())]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for name,result,hashes in pool.map(process_group,tasks):
            results[name]=result;input_hashes[name]=hashes;print(name,'analyzed',flush=True)
    summaries={};halves={};windows={};plateaus={}
    keys=list(variants)+['sparta-'+v for v in variants if v in ('C','D','E')]
    for key in keys:
        phases={phase:estimate([results[g][key][phase] for g in sorted(groups)],repeats=protocol['bootstrap_repeats'],seed=protocol['resampling_seed'],
            center_profile=True,width_band=True,interpolation_sensitivity=True) for phase in ('full','first','last','doubled')}
        summaries[key]=phases['full'];halves[key]=compare(phases['first'],phases['last'])
        windows[key]=compare(phases['full'],phases['doubled'])
        plateaus[key]=plateau_checks([results[g][key]['edges'] for g in sorted(groups)])
    pairs=[('A','B'),('B','C'),('D','C'),('E','C'),('F','C'),('G','C'),('sparta-D','sparta-C'),('sparta-E','sparta-C'),('C','sparta-C')]
    comparisons={a+' vs '+b:compare(summaries[a],summaries[b]) for a,b in pairs if a in summaries and b in summaries}
    if digest(source)!=protocol['campaign_sha256']:raise ValueError('Campaign record changed during analysis')
    if any(digest(Path(__file__).with_name(n))!=sha for n,sha in protocol['sources'].items()):
        raise ValueError('Analysis source code changed during execution')
    save(a.output,dict(status='candidate_pilot_diagnostics_only',physical_acceptance='incomplete',protocol=protocol,
        variants=variants,input_hashes=input_hashes,summaries=summaries,comparisons=comparisons,
        half_comparisons=halves,window_length_comparisons=windows,plateaus=plateaus))


if __name__=='__main__':main()
