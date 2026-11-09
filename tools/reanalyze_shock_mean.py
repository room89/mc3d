"""Separate, experimental mean-profile reanalysis; original report is read-only."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from shock_analysis import monotone,crossing
from shock_mean_profile import GRID,aligned,estimate,measure
from audit_shock_validation import gate


def sparta_density(folder,hashes):
    refpath=folder/'reference.json';raw=refpath.read_bytes();hashes[str(refpath.relative_to(folder.parent))]=hashlib.sha256(raw).hexdigest();ref=json.loads(raw)
    if ref['quick']:raise ValueError('Production data required')
    n1=ref['reference']['left']['number_density_m3'];n2=ref['reference']['right']['number_density_m3'];lam=ref['reference']['mean_free_path_m'];step=ref['steps_per_window'];windows=[]
    for wi in range(64):
        path=folder/f'profile.{(33+wi)*step}';raw=path.read_bytes();hashes[str(path.relative_to(folder.parent))]=hashlib.sha256(raw).hexdigest();lines=raw.decode().splitlines()
        header=next(i for i,line in enumerate(lines) if line.startswith('ITEM: CELLS '));fields=lines[header].split()[2:];groups={}
        for line in lines[header+1:]:
            if not line.strip():continue
            row=dict(zip(fields,map(float,line.split())))
            groups.setdefault(row['xlo'],[]).append(row)
        if len(groups)!=ref['source']['numerics']['cells_x']:raise ValueError('Incomplete grid')
        x=[];y=[]
        for lo,rows in sorted(groups.items()):
            volumes=[(r['xhi']-r['xlo'])*(r['yhi']-r['ylo'])*(r['zhi']-r['zlo']) for r in rows]
            density=sum(v*r['f_avg[1]'] for v,r in zip(volumes,rows))/sum(volumes)
            x.append((lo+rows[0]['xhi'])/(2*lam));y.append((density-n1)/(n2-n1))
        front=crossing(x,monotone(y),.5)
        windows.append(dict(front=front,x_aligned=[v-front for v in x],density_normalized=y))
    return windows


def compare(a,b):
    result={}
    for key in ('width_10_90','gradient_thickness'):
        u=a[key]['statistical_radius']+b[key]['statistical_radius']
        proc=a[key]['processing_sensitivity']+b[key]['processing_sensitivity']
        result[key]=gate(a[key]['mean']-b[key]['mean'],u+proc,.03*b[key]['mean'])
        result[key].update(statistical_radius=u,processing_sensitivity=proc)
    diff=np.asarray(a['density'])-b['density']
    u=a['density_simultaneous_radius']+b['density_simultaneous_radius']
    result['density']=gate(float(abs(diff).max()),u,.01)
    return result


def lag1(values):
    values=np.asarray(values);a=values[:-1];b=values[1:]
    if np.std(a)==0 or np.std(b)==0:return None
    return float(np.corrcoef(a,b)[0,1])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('series',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--alignment',choices=['crossing50','template','template_refined'],default='crossing50');p.add_argument('--jobs',type=int,default=3);p.add_argument('--bootstrap-repeats',type=int);p.add_argument('--resume',action='store_true');a=p.parse_args()
    repeats=a.bootstrap_repeats if a.bootstrap_repeats is not None else (128 if a.alignment=='template_refined' else 512)
    if a.output.exists():raise FileExistsError(a.output)
    raw=(a.series/'validation.json').read_bytes();report=json.loads(raw)
    if any(r['status']!='passed' for k in ('runs','sparta_runs') for r in report[k].values()):raise ValueError('Completed matrix required')
    result=dict(numpy_version=np.__version__,status='incomplete_experimental_estimator',source_report_sha256=hashlib.sha256(raw).hexdigest(),input_hashes={},method=dict(alignment=a.alignment,grid_step=.05,fit_width=1.,sensitivity_width=.5,blocks=[8,16],bootstrap_repeats=repeats,rng_seed=20261109 if a.alignment=='template_refined' else 20261107,confidence='Nominal 95% resampling envelopes; coverage not established for this physical series'),summaries={},half_comparisons={},lag1={},comparisons={})
    checkpoint=a.output.with_suffix('.partial.json')
    sources={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
             ('shock_alignment.py','shock_crossfit_gradient.py','shock_refined_estimator.py','shock_mean_profile.py')}
    result['algorithm_sources']=sources
    if checkpoint.exists():
        if not a.resume:raise FileExistsError(checkpoint)
        saved=json.loads(checkpoint.read_text())
        if saved['source_report_sha256']!=result['source_report_sha256'] or saved['algorithm_sources']!=sources or saved['method']['alignment']!=a.alignment or saved['method']['bootstrap_repeats']!=repeats:
            raise ValueError('Checkpoint protocol or source mismatch')
        for name,digest in saved['input_hashes'].items():
            if hashlib.sha256((a.series/name).read_bytes()).hexdigest()!=digest:raise ValueError('Checkpoint input changed')
        result=saved
    for name in list('ABCDEFG')+['sparta-C','sparta-D','sparta-E']:
        if name in result['summaries']:continue
        if name.startswith('sparta'):
            trials=[sparta_density(a.series/f'{name}-{seed}',result['input_hashes']) for seed in (123,456,789)]
        else:
            path=a.series/f'profiles-{name}.json';data=path.read_bytes();result['input_hashes'][path.name]=hashlib.sha256(data).hexdigest();trials=json.loads(data)
        if a.alignment in ('template','template_refined'):
            from shock_alignment import align
            data,diagnostics=align(trials,refine_passes=int(a.alignment=='template_refined'))
            result.setdefault('alignment_diagnostics',{})[name]=diagnostics
        else:data=aligned(trials)
        if a.alignment=='template_refined':
            from shock_refined_estimator import estimate as refined_estimate
            summary=refined_estimate(trials,repeats=repeats,jobs=a.jobs)
            first=refined_estimate([run[:32] for run in trials],repeats=repeats,jobs=a.jobs)
            last=refined_estimate([run[32:] for run in trials],repeats=repeats,jobs=a.jobs)
            result['method'].update(fit_width=1.5,sensitivity_widths=[1.,2.],refine_passes=1,
                resampling='Full raw-window estimator rebuilt; fixed original seed identities',
                confidence='Nominal 95% block envelope with seed-t heuristic; coverage unverified')
            result['method'].pop('sensitivity_width',None)
        else:
            summary=estimate(data,repeats=repeats)
            first=estimate(data[:,:32],repeats=min(256,repeats));last=estimate(data[:,32:],repeats=min(256,repeats))
        result['summaries'][name]=summary
        result['half_comparisons'][name]=compare(first,last)
        result['lag1'][name]=[dict(zip(('width_10_90','gradient_thickness'),[lag1(v) for v in measure(run).T])) for run in data]
        a.output.parent.mkdir(parents=True,exist_ok=True)
        temporary=checkpoint.with_suffix('.tmp')
        temporary.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temporary.replace(checkpoint)
        print(name,{k:summary[k] for k in ('width_10_90','gradient_thickness')},flush=True)
    for left,right in [('A','B'),('B','C'),('D','C'),('E','C'),('F','C'),('G','C'),('sparta-D','sparta-C'),('sparta-E','sparta-C'),('C','sparta-C')]:
        result['comparisons'][left+' vs '+right]=compare(result['summaries'][left],result['summaries'][right])
    # Keep committed diagnostics compact; curves are reproducible from hashed inputs.
    for summary in result['summaries'].values():
        summary['density_max_pointwise_radius']=max(summary.pop('density_pointwise_radius'))
        summary.pop('density')
    if (a.series/'validation.json').read_bytes()!=raw:raise RuntimeError('Source report changed')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
