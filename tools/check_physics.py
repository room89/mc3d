"""Run reproducible analytical acceptance scenarios."""
import argparse,csv,json,math,subprocess,time,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
from analytical_reference import cube_flux

def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter=';'))

def evaluate(out,name,full):
    manifest=json.loads((out/'manifest.json').read_text())
    if manifest['status']!='completed':return 'incomplete',{}
    for entry in manifest['outputs']:
        if hashlib.sha256((out/entry['path']).read_bytes()).hexdigest()!=entry['sha256']:
            raise ValueError('Output hash mismatch: '+entry['path'])
    summary=json.loads((out/'conservation_summary.json').read_text())
    cells=rows(out/'cell_statistics_0.csv')
    if not cells or any(r['complete']!='1' for r in cells):return 'incomplete',{}
    checks={}
    def check(key,value,expected,tolerance):
        checks[key]=dict(value=value,expected=expected,tolerance=tolerance,passed=math.isfinite(value) and abs(value-expected)<=tolerance)
    scale=sum(abs(summary[k]['energy']) for k in ('initial','incoming','outgoing'))+abs(summary['body_heat'])+1e-30
    pscale=sum(summary[k]['absolute_momentum'] for k in ('initial','incoming','outgoing'))+1e-30
    ledger=rows(out/'conservation.csv')
    for key in ('Rcount','RE','Rwall','Rp_x','Rp_y','Rp_z'):
        check(key,max(abs(float(r[key])) for r in ledger),0,0 if key=='Rcount' else 1e-10*(pscale if key.startswith('Rp') else scale))
    if name=='cube':
        surface=rows(out/'surface_statistics_0.csv')
        for r in surface:check('wall_'+r['polygon'],float(r['wall_temperature_at_output']),1,1e-12)
        if full:
            force,heat=cube_flux()
            check('Fx',sum(float(r['Fx']) for r in surface),force[0],.01*force[0])
            check('Qdot',sum(float(r['power']) for r in surface),heat,.01*heat)
    if name=='uniform' and full:
        a0=sum(float(r['A0']) for r in cells)
        u=[sum(float(r['A0'])*float(r['u'+k]) for r in cells)/a0 for k in 'xyz']
        temp=sum(float(r['A0'])*(3*float(r['Tstat'])+sum((float(r['u'+k])-u[i])**2 for i,k in enumerate('xyz'))) for r in cells)/(3*a0)
        density=sum(float(r['density'])*float(r['Vgas']) for r in cells)/sum(float(r['Vgas']) for r in cells)
        for key,value,expected,tol in [('density',density,1,.01),('ux',u[0],5*math.sqrt(2),.05*math.sqrt(2)),('uy',u[1],0,.01),('uz',u[2],0,.01),('T',temp,1,.01)]:check(key,value,expected,tol)
    if name=='closed':
        check('input',summary['incoming']['N'],0,0);check('output',summary['outgoing']['N'],0,0)
        check('energy',summary['final']['energy']-summary['initial']['energy']+summary['wall_energy_change'],0,1e-10*scale)
        check('exchange',int(abs(summary['body_heat'])>1e-8),1,0)
        cfg=json.loads((out/'effective_config.json').read_text()) if (out/'effective_config.json').exists() else manifest['effective_config']
        weight=cfg['lx']*cfg['ly']*cfg['lz']/(cfg['cells_x']*cfg['cells_y']*cfg['cells_z']*cfg['particles_per_cell'])
        energies=[]
        for filename in ('initial.dat','final.dat'):
            gas=sum(float(r['N'])*float(r['E'])*weight for r in rows(out/filename))
            wall=sum(float(r['temperature'])*float(r['heat_capacity']) for r in rows(out/(filename+'.surface.csv')))
            energies.append(gas+wall)
        # Legacy gas snapshots print eight significant digits; ledger checks above remain strict.
        check('snapshot_energy',energies[1]-energies[0],0,1e-7*max(energies))
    return ('passed' if all(c['passed'] for c in checks.values()) else 'failed'),checks

def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--mode',choices=['quick','full'],default='quick');p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--load-mode',choices=['full','balanced','quiet'],default='balanced');p.add_argument('--resume',action='store_true');a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];a.binary=a.binary.resolve();a.output_dir.mkdir(parents=True,exist_ok=a.resume)
    report={'mode':a.mode,'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),'runs':[]}
    micro=a.binary.parent/'mc3d_tests'
    micro_xml=a.output_dir/'microphysics.xml'
    subprocess.run([str(micro),'--gtest_filter=HsScatteringTest.*:CollisionFrequencyTest.*:VhsTest.*:StatisticalMomentsTest.*', '--gtest_output=xml:'+str(micro_xml.resolve())],check=True,stdout=(a.output_dir/'microphysics.log').open('w'))
    report['microphysics']={'status':'passed','tests':ET.parse(micro_xml).getroot().attrib}
    for name in ('uniform','cube','closed'):
        for seed in (123,456,789):
            dest=a.output_dir/f'{name}-{seed}';cfg=json.loads((repo/f'configs/checks/analytical_{name}.json').read_text());cfg['seed']=seed
            if a.mode=='quick':cfg.update(cells_x=5,cells_y=5,cells_z=5,particles_per_cell=20,end_time=.1,statistics_windows=[[0,.1]])
            config=a.output_dir/f'{name}-{seed}.json'
            if dest.exists():
                if not a.resume or json.loads(config.read_text())!=cfg:raise ValueError('Existing configuration mismatch')
                manifest=json.loads((dest/'manifest.json').read_text())
                if manifest['binary_sha256']!=report['binary_sha256']:raise ValueError('Existing binary mismatch')
                status,checks=evaluate(dest,name,a.mode=='full')
                report['runs'].append(dict(scenario=name,seed=seed,status=status,checks=checks,reused=True))
                (a.output_dir/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
                continue
            config.write_text(json.dumps(cfg,indent=2))
            control=a.output_dir/f'{name}-{seed}.control'
            subprocess.run(['python3',str(repo/'tools/control_simulation.py'),'.',a.load_mode,'--file',str(control),'--create-parent'],check=True)
            result={'scenario':name,'seed':seed};start=time.monotonic()
            try:
                with (a.output_dir/f'{name}-{seed}.log').open('w') as log:
                    proc=subprocess.run([str(a.binary),'--config',str(config),'--output-dir',str(dest),'--runtime-control',str(control)],stdout=log,stderr=subprocess.STDOUT,timeout=28800)
                if proc.returncode:result['status']='incomplete' if proc.returncode==130 else 'error'
                else:result['status'],result['checks']=evaluate(dest,name,a.mode=='full')
            except subprocess.TimeoutExpired:result['status']='incomplete'
            except (ValueError,KeyError,OSError) as e:result.update(status='error',error=str(e))
            result['seconds']=time.monotonic()-start;report['runs'].append(result)
            (a.output_dir/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False));print(name,seed,result['status'],flush=True)
    (a.output_dir/'report.txt').write_text('\n'.join(f"{r['scenario']} seed={r['seed']}: {r['status']}" for r in report['runs'])+'\n')
    return max({'passed':0,'failed':1,'incomplete':2,'error':3}[r['status']] for r in report['runs'])
if __name__=='__main__':
    try:raise SystemExit(main())
    except KeyboardInterrupt:raise SystemExit(2)
    except (OSError,ValueError,KeyError,subprocess.CalledProcessError) as error:
        print(str(error),flush=True);raise SystemExit(3)
