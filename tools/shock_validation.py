"""Reproducible stage-9b matrix. Execution is not physical acceptance."""
import argparse,hashlib,json,math,subprocess,sys,time,os
from pathlib import Path
from normal_shock import setup
from physical_config import prepare
SEEDS=(123,456,789)
VARIANTS={'A':(400,160),'B':(800,80),'C':(800,160),'D':(400,320),'E':(800,160),'F':(1600,160),'G':(800,160)}

def configuration(variant,seed):
    data,ref=setup(False,seed);nx,np=VARIANTS[variant]
    data['numerics'].update(cells_x=nx,particles_per_cell=np)
    tau=ref['transit_time_s'];data['end_time_s']=12*tau
    data['statistics_windows_s']=[[(4+i/8)*tau,(4+(i+1)/8)*tau] for i in range(64)]
    if variant=='E':data['numerics'].update(cu=.25,collision_courant=.05)
    if variant=='F':data['domain_m'][0]*=2
    if variant=='G':
        data['domain_m'][1]*=2;data['domain_m'][2]*=2
        data['numerics'].update(cells_y=2,cells_z=2)
    return data,ref

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(path)

def verify_run(out,expected_binary):
    m=json.loads((out/'manifest.json').read_text())
    if m['status']!='completed':return False
    if m['binary_sha256']!=expected_binary:raise ValueError('Binary mismatch')
    for entry in m['outputs']:
        if digest(out/entry['path'])!=entry['sha256']:raise ValueError('Output hash mismatch: '+entry['path'])
    return True

def run_limited(command,log,timeout,memory_mb,cwd=None):
    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=cwd)
    deadline=time.monotonic()+timeout
    try:
        while process.poll() is None:
            if time.monotonic()>=deadline:raise subprocess.TimeoutExpired(command,timeout)
            usage=subprocess.run(['ps','-o','rss=','-p',str(process.pid)],capture_output=True,text=True)
            if usage.returncode and process.poll() is None:raise RuntimeError('RSS monitoring unavailable')
            if usage.stdout.strip() and int(usage.stdout.strip())>memory_mb*1024:raise RuntimeError('Memory limit exceeded')
            try:process.wait(timeout=1)
            except subprocess.TimeoutExpired:pass
        return process
    finally:
        if process.poll() is None:process.kill();process.wait()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--memory-limit-mb',type=int,default=8192);p.add_argument('--sparta-bin',type=Path);p.add_argument('--resume',action='store_true');p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--budget-hours',type=float,default=24);p.add_argument('--load-mode',choices=['full','balanced','quiet'],default='balanced')
    a=p.parse_args()
    if a.memory_limit_mb<=0:raise ValueError('Memory limit must be positive')
    if not 0<a.budget_hours<=24:raise ValueError('Budget must be in (0,24] hours')
    root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=a.resume)
    binary=a.binary.resolve();binary_hash=digest(binary);report_file=root/'validation.json'
    report=json.loads(report_file.read_text()) if a.resume and report_file.exists() else {'schema_version':1,'binary_sha256':binary_hash,'spent_seconds':0,'runs':{},'stationarity':{'status':'incomplete'},'convergence':{'status':'incomplete'},'cross_code':{'status':'incomplete'},'experiment':{'status':'incomplete','reason':'Primary experimental data not yet verified'}}
    if report['binary_sha256']!=binary_hash:raise ValueError('Cannot resume with changed binary')
    if report.get('owner_pid'):
        try:os.kill(report['owner_pid'],0)
        except ProcessLookupError:pass
        else:raise ValueError('Another validation runner may still be active')
    report['owner_pid']=os.getpid()
    for entry in report['runs'].values():
        if entry.get('reason')=='Running':
            report['spent_seconds']+=min(28800,max(0,time.time()-entry['started_at']))
            entry.update(status='incomplete',reason='Interrupted runner; output retained')
    if a.sparta_bin and not a.prepare_only and not report.get('sparta_preflight') and report['spent_seconds']<a.budget_hours*3600:
        preflight=root/'sparta-preflight'
        if preflight.exists():raise ValueError('Incomplete preflight retained; inspect before retry')
        began=time.monotonic()
        subprocess.run([sys.executable,str(Path(__file__).with_name('sparta_preflight.py')),'--binary',str(a.sparta_bin.resolve()),'--output-dir',str(preflight)],check=True,timeout=min(1200,max(1,a.budget_hours*3600-report['spent_seconds'])))
        report['spent_seconds']+=time.monotonic()-began
        report['sparta_preflight']=json.loads((preflight/'report.json').read_text())
    save(report_file,report)
    for variant in VARIANTS:
        for seed in SEEDS:
            key=f'{variant}-{seed}';data,ref=configuration(variant,seed)
            folder=root/key;folder.mkdir(exist_ok=True)
            input_file=folder/'input_si.json'
            if input_file.exists() and json.loads(input_file.read_text())!=data:raise ValueError('Configuration mismatch')
            save(input_file,data);cfg,units=prepare(data)
            cfg['units_file']=str(folder/'units.json');units['source_file']=str(input_file)
            save(folder/'solver.json',cfg);save(folder/'units.json',units);save(folder/'reference.json',ref)
            item=report['runs'].setdefault(key,{'status':'incomplete','config_sha256':digest(folder/'solver.json')})
            out=folder/'output'
            if item['status']=='passed':
                if not verify_run(out,binary_hash):raise ValueError('Completed result is incomplete')
                continue
            if out.exists():
                item['reason']='Existing unfinished output retained; not restarted';save(report_file,report);continue
            remaining=a.budget_hours*3600-report['spent_seconds']
            if a.prepare_only or remaining<=0:
                item['reason']='Prepared only' if a.prepare_only else 'Budget exhausted';save(report_file,report);continue
            # Conservative projection from completed runs, scaled by particle count and axial cells.
            previous=[v for v in report['runs'].values() if v.get('seconds') and v['status']=='passed']
            scale=cfg['cells_x']**2*cfg['cells_y']*cfg['cells_z']*cfg['particles_per_cell']
            estimate=max((v['seconds']/v['cost_scale']*scale for v in previous),default=0)*1.5
            if estimate>min(remaining,28800):
                item['reason']='Projected cost exceeds remaining/per-run budget';save(report_file,report);continue
            control=folder/'control.json'
            subprocess.run([sys.executable,str(Path(__file__).with_name('control_simulation.py')),'.',a.load_mode,'--file',str(control)],check=True)
            start=time.monotonic();item.update(status='incomplete',reason='Running',started_at=time.time());save(report_file,report)
            try:
                with (folder/'run.log').open('w') as log:
                    result=run_limited([str(binary),'--config',str(folder/'solver.json'),'--output-dir',str(out),'--runtime-control',str(control)],log,min(28800,remaining),a.memory_limit_mb)
                item['status']='passed' if result.returncode==0 and verify_run(out,binary_hash) else 'incomplete' if result.returncode==130 else 'error'
                item['reason']='Execution only; physical acceptance is separate'
            except subprocess.TimeoutExpired:item.update(status='incomplete',reason='Time budget exceeded')
            except (RuntimeError,OSError) as e:item.update(status='error',reason=str(e))
            finally:
                item['seconds']=time.monotonic()-start;item['cost_scale']=scale
                report['spent_seconds']+=item['seconds'];save(report_file,report)
            print(key,item['status'],flush=True)
    if a.sparta_bin:
        from sparta_shock import export
        sparta=a.sparta_bin.resolve();sparta_hash=digest(sparta)
        if report.get('sparta_sha256',sparta_hash)!=sparta_hash:raise ValueError('SPARTA binary mismatch')
        report['sparta_sha256']=sparta_hash
        independent=report.setdefault('sparta_runs',{})
        for variant in ('C','D','E'):
            for seed in SEEDS:
                key=f'{variant}-{seed}';folder=root/('sparta-'+key)
                entry=independent.setdefault(key,{'status':'incomplete'})
                if entry['status']=='passed':
                    for name,sha in entry['hashes'].items():
                        if digest(folder/name)!=sha:raise ValueError('SPARTA output hash mismatch')
                    continue
                if folder.exists():entry['reason']='Unfinished output retained';continue
                remaining=a.budget_hours*3600-report['spent_seconds']
                if a.prepare_only or remaining<=0:entry['reason']='Not executed within budget';continue
                export(folder,variant,seed)
                start=time.monotonic()
                try:
                    with (folder/'run.log').open('w') as log:
                        result=run_limited([str(sparta),'-in',str(folder/'in.shock')],log,min(remaining,28800),a.memory_limit_mb,cwd=folder)
                    reference=json.loads((folder/'reference.json').read_text());expected=96*reference['steps_per_window']
                    complete=(folder/f'profile.{expected}').exists()
                    entry['status']='passed' if result.returncode==0 and complete else 'error'
                    entry['hashes']={p.name:digest(p) for p in folder.iterdir() if p.is_file()}
                except subprocess.TimeoutExpired:entry.update(status='incomplete',reason='Time budget exceeded')
                except (RuntimeError,OSError) as e:entry.update(status='error',reason=str(e))
                finally:
                    entry['seconds']=time.monotonic()-start;report['spent_seconds']+=entry['seconds'];save(report_file,report)
    from analyze_shock_validation import analyze
    save(report_file,report)
    report=analyze(root)
    report['execution']={'status':'passed' if all(v['status']=='passed' for v in report['runs'].values()) else 'incomplete'};save(report_file,report)
    report.pop('owner_pid',None);save(report_file,report)
    if any(r['status']=='error' for r in list(report['runs'].values())+list(report.get('sparta_runs',{}).values())):return 3
    statuses=[report[k]['status'] for k in ('execution','stationarity','convergence','cross_code','experiment')]
    return 0 if all(s=='passed' for s in statuses) else 1 if 'failed' in statuses else 2
if __name__=='__main__':raise SystemExit(main())
