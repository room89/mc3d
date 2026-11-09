"""Planar shock engineering runs; no convergence or experimental acceptance claim."""
import argparse,csv,json,math,subprocess,sys,tempfile
from pathlib import Path
from normal_shock import setup
from physical_config import prepare

def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter=';'))

def inspect(out,config,reference):
    manifest=json.loads((out/'manifest.json').read_text())
    if manifest['status']!='completed':raise RuntimeError('Incomplete solver run')
    profiles=[];fronts=[]
    left=reference['left']['number_density_m3'];ratio=reference['right']['number_density_m3']/left
    dx=config['lx']/config['cells_x'];target=(1+ratio)/2
    for window in range(10):
        rows=table(out/f'cell_statistics_{window}.csv')
        if not rows or any(r['complete']!='1' for r in rows):raise RuntimeError('Incomplete window')
        profile=[]
        for r in rows:
            index=int(r['cell']);x=-config['lx']/2+(index+.5)*dx
            row=dict(window=window,x=x,density=float(r['density']))
            for key in ('ux','uy','uz','Tstat','Pxx','Pyy','Pzz','qx'):row[key]=float(r[key]) if r[key] else None
            profile.append(row)
        crossings=[]
        for a,b in zip(profile,profile[1:]):
            if a['density']<target<=b['density']:crossings.append(a['x']+(target-a['density'])*(b['x']-a['x'])/(b['density']-a['density']))
        front=crossings[0] if len(crossings)==1 else None
        for row in profile:row['x_aligned']=row['x']-front if front is not None else None
        profiles.extend(profile);fronts.append(dict(window=window,position=front,status='valid' if len(crossings)==1 else 'missing' if not crossings else 'ambiguous',crossings=crossings))
    summary=json.loads((out/'conservation_summary.json').read_text())
    scale=sum(abs(summary[k]['energy']) for k in ('initial','incoming','outgoing'))+1e-30
    momentum=sum(summary[k]['absolute_momentum'] for k in ('initial','incoming','outgoing'))+1e-30
    for row in table(out/'conservation.csv'):
        if int(row['Rcount'])!=0:raise ValueError('Particle balance failed')
        for key in ('RE','Rp_x','Rp_y','Rp_z'):
            value=float(row[key]);limit=1e-10*(scale if key=='RE' else momentum)
            if not math.isfinite(value) or abs(value)>limit:raise ValueError('Balance failed: '+key)
    with (out/'shock_profile.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
    (out/'shock_front.json').write_text(json.dumps(fronts,indent=2))
    return dict(status='passed',fronts=fronts,interpretation=__doc__)

def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--mode',choices=['quick','pilot','validation'],default='quick');p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--sparta-bin',type=Path);p.add_argument('--resume',action='store_true');p.add_argument('--prepare-only',action='store_true');p.add_argument('--load-mode',choices=['full','balanced','quiet'],default='balanced');a=p.parse_args()
    if a.mode=='validation':
        command=[sys.executable,str(Path(__file__).with_name('shock_validation.py')),'--binary',str(a.binary),'--output-dir',str(a.output_dir),'--load-mode',a.load_mode]
        if a.sparta_bin:command+=['--sparta-bin',str(a.sparta_bin)]
        if a.resume:command+=['--resume']
        if a.prepare_only:command+=['--prepare-only']
        return subprocess.call(command)
    a.output_dir.mkdir(parents=True,exist_ok=False);results=[];repo=Path(__file__).resolve().parents[1]
    for seed in (123,456,789):
        data,ref=setup(a.mode=='quick',seed);cfg,units=prepare(data);prepared=a.output_dir/f'config-{seed}';prepared.mkdir()
        cfg['units_file']=str((prepared/'units.json').resolve())
        units['source_file']=str((prepared/'input_si.json').resolve())
        for name,value in [('solver.json',cfg),('input_si.json',data),('units.json',units),('reference.json',ref)]:
            (prepared/name).write_text(json.dumps(value,indent=2))
        out=a.output_dir/f'run-{seed}';control=a.output_dir/f'control-{seed}.json'
        subprocess.run([sys.executable,str(repo/'tools/control_simulation.py'),'.',a.load_mode,'--file',str(control)],check=True)
        try:
            with (a.output_dir/f'run-{seed}.log').open('w') as log:
                r=subprocess.run([str(a.binary.resolve()),'--config',str(prepared/'solver.json'),'--output-dir',str(out),'--runtime-control',str(control)],stdout=log,stderr=subprocess.STDOUT,timeout=28800)
            result=inspect(out,cfg,ref) if r.returncode==0 else dict(status='incomplete' if r.returncode==130 else 'error')
        except subprocess.TimeoutExpired:result=dict(status='incomplete')
        except RuntimeError as e:result=dict(status='incomplete',reason=str(e))
        except ValueError as e:result=dict(status='failed',reason=str(e))
        except (KeyError,OSError) as e:result=dict(status='error',reason=str(e))
        result['seed']=seed;results.append(result);print(seed,result['status'],flush=True)
        (a.output_dir/'report.json').write_text(json.dumps(results,indent=2))
    return max({'passed':0,'failed':1,'incomplete':2,'error':3}[r['status']] for r in results)
if __name__=='__main__':raise SystemExit(main())
