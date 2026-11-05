"""Conservative numerical gates for stage 9b; experimental status stays separate."""
import argparse,csv,json,math,statistics
from pathlib import Path
from shock_analysis import pooled,metrics,interval,equivalence
from shock_validation import save,SEEDS

def profile(folder,stride=1):
    cfg=json.loads((folder/'solver.json').read_text());ref=json.loads((folder/'reference.json').read_text());units=json.loads((folder/'units.json').read_text())
    out=folder/'output';manifest=json.loads((out/'manifest.json').read_text())
    if manifest['status']!='completed':raise ValueError('Incomplete run')
    nx,ny,nz=(cfg['cells_'+a] for a in 'xyz');dx=cfg['lx']/nx
    weight=cfg['lx']*cfg['ly']*cfg['lz']/(nx*ny*nz*cfg['particles_per_cell'])
    lam=ref['mean_free_path_m']/units['scales']['length_m'];ratio=ref['right']['number_density_m3']/ref['left']['number_density_m3']
    windows=[];pending=[]
    for wi,bounds in enumerate(cfg['statistics_windows']):
        with (out/f'cell_statistics_{wi}.csv').open() as f:rows=list(csv.DictReader(f,delimiter=';'))
        if len(rows)!=nx*ny*nz or any(r['complete']!='1' for r in rows):raise ValueError('Missing cells or incomplete window')
        if {int(r['cell']) for r in rows}!=set(range(nx*ny*nz)):raise ValueError('Invalid cell indices')
        pending.extend(rows)
        if (wi+1)%stride:continue
        rows=pending;pending=[]
        bounds=[cfg['statistics_windows'][wi-stride+1][0],bounds[1]]
        groups=[[] for _ in range(nx)]
        for r in rows:groups[int(r['cell'])//(ny*nz)].append(r)
        moments=[pooled(g,weight,bounds[1]-bounds[0],sum(float(r['Vgas']) for r in g)/stride) for g in groups]
        x=[(cfg.get('apex_x',-cfg['lx']/2)+(i+.5)*dx)/lam for i in range(nx)]
        rho=[m['density'] for m in moments];result=metrics(x,rho,1,ratio,1)
        result['temperature']=[m['temperature'] for m in moments];result['velocity']=[m['velocity'][0] for m in moments]
        nedge=max(1,nx//10)
        result['plateau_error']=max(abs(statistics.mean(rho[:nedge])-1),abs(statistics.mean(rho[-nedge:])-ratio))/(ratio-1)
        result['boundary_distance']=min(result['front']-x[0],x[-1]-result['front'])
        result['window']=wi;windows.append(result)
    if len(windows)!=64//stride:raise ValueError('Expected 64 production windows')
    summary=json.loads((out/'conservation_summary.json').read_text())
    es=sum(abs(summary[k]['energy']) for k in ('initial','incoming','outgoing'))+1e-30
    ps=sum(summary[k]['absolute_momentum'] for k in ('initial','incoming','outgoing'))+1e-30
    with (out/'conservation.csv').open() as f:
        for r in csv.DictReader(f,delimiter=';'):
            if int(r['Rcount'])!=0:raise ValueError('Count balance failed')
            for k in ('RE','Rp_x','Rp_y','Rp_z'):
                if not math.isfinite(float(r[k])) or abs(float(r[k]))>1e-10*(es if k=='RE' else ps):raise ValueError('Balance failed')
    return windows

def interpolate(x,y,point):
    for i in range(len(x)-1):
        if x[i]<=point<=x[i+1]:return y[i]+(point-x[i])/(x[i+1]-x[i])*(y[i+1]-y[i])
    raise ValueError('Comparison grid outside profile')

def summarize(trials):
    result={};count=len(trials[0]);factor=64//count
    for key in ('width_10_90','gradient_thickness'):
        means=[statistics.mean(w[key] for w in run) for run in trials]
        mean,u=interval(means)
        for block in (8//factor,16//factor):
            block_u=[interval([statistics.mean(w[key] for w in run[i:i+block]) for i in range(0,count,block)])[1] for run in trials]
            u=max(u,max(block_u))
        result[key]=dict(mean=mean,uncertainty=u)
    grid=[-10+i*.1 for i in range(201)]
    profiles=[[statistics.mean(interpolate(w['x_aligned'],w['density_normalized'],x) for w in run) for x in grid] for run in trials]
    result['grid']=grid;result['density']=[];result['density_uncertainty']=[]
    for i in range(len(grid)):
        m,u=interval([p[i] for p in profiles])
        for run in trials:
            values=[interpolate(w['x_aligned'],w['density_normalized'],grid[i]) for w in run]
            for block in (8//factor,16//factor):
                u=max(u,interval([statistics.mean(values[j:j+block]) for j in range(0,count,block)])[1])
        result['density'].append(m);result['density_uncertainty'].append(u)
    result['stationarity']='passed'
    for run in trials:
        if max(w['plateau_error'] for w in run)>.01 or min(w['boundary_distance'] for w in run)<10:result['stationarity']='incomplete'
        for k in ('width_10_90','gradient_thickness'):
            first=[statistics.mean(w[k] for w in run[i:i+8//factor]) for i in range(0,count//2,8//factor)]
            last=[statistics.mean(w[k] for w in run[i:i+8//factor]) for i in range(count//2,count,8//factor)]
            m1,u1=interval(first);m2,u2=interval(last)
            if not equivalence(m1-m2,u1+u2,.03*result[k]['mean']):result['stationarity']='incomplete'
    result['processing_uncertainty']=max(w['processing_uncertainty'] for run in trials for w in run)
    return result

def compare(a,b):
    checks={}
    checks['density']=all(equivalence(x-y,u+v,.01) for x,y,u,v in zip(a['density'],b['density'],a['density_uncertainty'],b['density_uncertainty']))
    for k in ('width_10_90','gradient_thickness'):
        extra=(a['processing_uncertainty']+b['processing_uncertainty']) if k=='gradient_thickness' else 0
        checks[k]=equivalence(a[k]['mean']-b[k]['mean'],a[k]['uncertainty']+b[k]['uncertainty']+extra,.03*b[k]['mean'])
    return dict(status='passed' if all(checks.values()) else 'incomplete',checks=checks)

def sparta_profile(folder):
    ref=json.loads((folder/'reference.json').read_text())
    if ref['quick']:raise ValueError('Smoke run is not production data')
    n1=ref['reference']['left']['number_density_m3'];n2=ref['reference']['right']['number_density_m3'];lam=ref['reference']['mean_free_path_m'];step=ref['steps_per_window'];windows=[]
    for wi in range(64):
        lines=(folder/f'profile.{(33+wi)*step}').read_text().splitlines()
        header=next(i for i,line in enumerate(lines) if line.startswith('ITEM: CELLS '));fields=lines[header].split()[2:]
        rows=[dict(zip(fields,map(float,line.split()))) for line in lines[header+1:] if line.strip()]
        nx=ref['source']['numerics']['cells_x'];groups={}
        for r in rows:groups.setdefault(r['xlo'],[]).append(r)
        if len(groups)!=nx:raise ValueError('Incomplete SPARTA grid')
        x=[];density=[];temp=[];velocity=[]
        for lo,g in sorted(groups.items()):
            vols=[(r['xhi']-r['xlo'])*(r['yhi']-r['ylo'])*(r['zhi']-r['zlo']) for r in g]
            counts=[v*r['f_avg[1]'] for v,r in zip(vols,g)];count=sum(counts)
            if count<=0:raise ValueError('Empty SPARTA bin')
            u=[sum(n*r[f'f_avg[{i}]'] for n,r in zip(counts,g))/count for i in (2,3,4)]
            thermal=sum(n*(r['f_avg[8]']+sum((r[f'f_avg[{i+2}]']-u[i])**2 for i in range(3))*6.633521472521616e-26/(3*1.380649e-23)) for n,r in zip(counts,g))/count
            x.append((lo+g[0]['xhi'])/(2*lam));density.append(count/sum(vols));velocity.append(u[0]);temp.append(thermal)
        result=metrics(x,density,n1,n2,1);edge=max(1,nx//10)
        result.update(temperature=temp,velocity=velocity,plateau_error=max(abs(statistics.mean(density[:edge])-n1),abs(statistics.mean(density[-edge:])-n2))/(n2-n1),boundary_distance=min(result['front']-x[0],x[-1]-result['front']),window=wi)
        windows.append(result)
    return windows

def analyze(root):
    report=json.loads((root/'validation.json').read_text());summaries={};errors={}
    for v in 'ABCDEFG':
        try:
            trials=[profile(root/f'{v}-{s}') for s in SEEDS];summaries[v]=summarize(trials)
            save(root/f'profiles-{v}.json',trials)
        except (OSError,ValueError,KeyError) as e:errors[v]=str(e)
    report['analysis_errors']=errors;report['summaries']=summaries
    comparisons={}
    for left,right in [('A','B'),('B','C'),('D','C'),('E','C'),('F','C'),('G','C')]:
        if left in summaries and right in summaries:comparisons[left+'-'+right]=compare(summaries[left],summaries[right])
    report['stationarity']={'status':'passed' if len(summaries)==7 and all(s['stationarity']=='passed' for s in summaries.values()) else 'incomplete'}
    window_checks={}
    for v in summaries:
        try:window_checks[v]=compare(summaries[v],summarize([profile(root/f'{v}-{seed}',2) for seed in SEEDS]))
        except (OSError,ValueError,KeyError) as e:window_checks[v]=dict(status='incomplete',reason=str(e))
    report['window_sensitivity']={'status':'passed' if len(window_checks)==7 and all(c['status']=='passed' for c in window_checks.values()) else 'incomplete','comparisons':window_checks}
    passed=len(comparisons)==6 and all(c['status']=='passed' for c in comparisons.values()) and report['stationarity']['status']=='passed' and report['window_sensitivity']['status']=='passed'
    report['convergence']={'status':'passed' if passed else 'incomplete','comparisons':comparisons}
    sparta={}
    for v in ('C','D','E'):
        try:sparta[v]=summarize([sparta_profile(root/f'sparta-{v}-{seed}') for seed in SEEDS])
        except (OSError,ValueError,KeyError) as e:errors['sparta-'+v]=str(e)
    cross={}
    if len(sparta)==3 and 'C' in summaries:
        cross['grid']=compare(sparta['D'],sparta['C']);cross['time']=compare(sparta['E'],sparta['C']);cross['mc3d']=compare(summaries['C'],sparta['C'])
        # Also require differences to lie within assessed uncertainty.
        left,right=summaries['C'],sparta['C']
        within=all(abs(a-b)<=u+v for a,b,u,v in zip(left['density'],right['density'],left['density_uncertainty'],right['density_uncertainty']))
        for k in ('width_10_90','gradient_thickness'):
            within &= abs(left[k]['mean']-right[k]['mean'])<=left[k]['uncertainty']+right[k]['uncertainty']+left['processing_uncertainty']+right['processing_uncertainty']
        report['cross_code']={'status':'passed' if passed and within and all(c['status']=='passed' for c in cross.values()) and all(s['stationarity']=='passed' for s in sparta.values()) else 'incomplete','comparisons':cross}
    else:report['cross_code']={'status':'incomplete','reason':'Missing independent production results'}
    report['sparta_summaries']=sparta
    save(root/'validation.json',report)
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();analyze(a.directory)
