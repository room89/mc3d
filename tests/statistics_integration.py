"""Conservation, overlapping windows and observational-only statistics."""
import csv,json,math,subprocess,sys,tempfile,signal,time,hashlib
from pathlib import Path
binary=Path(sys.argv[1]).resolve();repo=Path(__file__).resolve().parents[1]
def rows(path):
    with path.open() as f:
        result=list(csv.DictReader(f,delimiter=';'))
    assert all(None not in r and None not in r.values() for r in result),path
    return result
with tempfile.TemporaryDirectory(prefix='mc3d-statistics-') as temporary:
    root=Path(temporary)
    for name in ['closed','open','free','mirror','fixed_wall','dynamic_wall','dynamic_wall_parallel']:
        model=name.removesuffix('_parallel')
        config=repo/'configs/checks'/('statistics_'+model+'.json')
        data=json.loads(config.read_text());out=root/name
        if name.endswith('_parallel'):
            data['thread_pool_size']=4
            config=root/'parallel.json';config.write_text(json.dumps(data))
        subprocess.run([str(binary),'--config',str(config),'--output-dir',str(out)],check=True,capture_output=True)
        summary=json.loads((out/'conservation_summary.json').read_text())
        manifest=json.loads((out/'manifest.json').read_text())
        entry=next(e for e in manifest['outputs'] if e['path']=='conservation_summary.json')
        assert entry['sha256']==hashlib.sha256((out/'conservation_summary.json').read_bytes()).hexdigest()
        if name.endswith('_parallel'):
            subprocess.run([sys.executable,str(repo/'tools/replay_run.py'),str(out/'manifest.json'),'--output-dir',str(root/'replay')],check=True,capture_output=True)
        assert summary['residual_count']==0
        if name=='closed':assert math.isclose(summary['initial']['mass'],1,rel_tol=1e-12)
        scale=summary['initial']['energy']+abs(summary['incoming']['energy'])+abs(summary['outgoing']['energy'])+abs(summary['body_heat'])+1e-30
        momentum_scale=summary['initial']['absolute_momentum']+summary['incoming']['absolute_momentum']+summary['outgoing']['absolute_momentum']+1e-30
        for row in rows(out/'conservation.csv'):
            assert int(row['Rcount'])==0
            for key in ['RE','Rwall']:assert abs(float(row[key]))<=1e-10*scale,(name,key,row)
            for key in ['Rp_x','Rp_y','Rp_z']:assert abs(float(row[key]))<=1e-10*momentum_scale,(name,key,row)
        for i,(start,end) in enumerate(data['statistics_windows']):
            table=rows(out/f'cell_statistics_{i}.csv')
            assert table
            if name=='closed':assert math.isclose(sum(float(r['mean_N']) for r in table),summary['initial']['N'],rel_tol=1e-12)
            for row in table:
                assert math.isclose(float(row['observed_time']),end-start,rel_tol=1e-12)
                if row['status']=='valid':assert float(row['Tstat'])>=-1e-12
                else:assert row['Tstat']==''
            surface=rows(out/f'surface_statistics_{i}.csv')
            if model=='fixed_wall':
                assert all(float(row['wall_temperature_at_output'])==data['wall_temperature'] for row in surface)
        if model in ('fixed_wall','dynamic_wall'):assert abs(summary['body_heat'])>1e-8
        if model=='mirror':
            assert summary['incoming']['N']==summary['outgoing']['N']==0
            assert abs(summary['mirror_heat'])<1e-12
            assert sum(abs(x) for x in summary['mirror_impulse'])>0
        disabled=data|{'statistics_enabled':False,'balance_enabled':False}
        disabled_file=root/(name+'-off.json');disabled_file.write_text(json.dumps(disabled))
        off=root/(name+'-off')
        subprocess.run([str(binary),'--config',str(disabled_file),'--output-dir',str(off)],check=True,capture_output=True)
        assert not list(off.glob('*statistics*')) and not list(off.glob('conservation*'))
        for f in off.iterdir():
            if f.suffix in ('.dat','.net','.csv'):
                assert f.read_bytes()==(out/f.name).read_bytes(),(name,f.name)
    # Output cadence must not change any accumulated moment, balance or final state.
    # Include overlapping windows and output-time temperatures of a dynamic wall.
    for model in ('closed', 'dynamic_wall'):
        source=json.loads((repo/f'configs/checks/statistics_{model}.json').read_text())
        outputs=[]
        for label,interval in [('dense',0),('sparse',source['end_time'])]:
            config=root/f'{model}-{label}.json'
            config.write_text(json.dumps(source|{'snapshot_interval':interval}))
            out=root/f'{model}-{label}'
            subprocess.run([str(binary),'--config',str(config),'--output-dir',str(out)],check=True,capture_output=True)
            outputs.append(out)
            # The first window closes before the final step and must stay untouched.
            assert (out/'cell_statistics_0.csv').stat().st_mtime_ns < (out/'conservation.csv').stat().st_mtime_ns
        for pattern in ('*statistics_*.csv','conservation*','final.dat*'):
            for path in outputs[0].glob(pattern):
                assert path.read_bytes()==(outputs[1]/path.name).read_bytes(),(model,path.name)
    empty=root/'zero.json';empty.write_text(json.dumps(data|{'end_time':0,'statistics_windows':[]}))
    zero=root/'zero'
    subprocess.run([str(binary),'--config',str(empty),'--output-dir',str(zero)],check=True,capture_output=True)
    assert rows(zero/'conservation.csv')==[]
    assert json.loads((zero/'conservation_summary.json').read_text())['steps']==0
    assert all(r['status']=='unobserved' for r in rows(zero/'cell_statistics_0.csv'))
    partial=json.loads((repo/'configs/checks/statistics_closed.json').read_text())
    partial.update(end_time=1000,statistics_windows=[[0,500],[500,1000]])
    partial_file=root/'partial.json';partial_file.write_text(json.dumps(partial));partial_out=root/'partial'
    process=subprocess.Popen([str(binary),'--config',str(partial_file),'--output-dir',str(partial_out)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
    try:
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            ledger=partial_out/'conservation.csv'
            if ledger.exists() and len(ledger.read_text().splitlines())>2:break
            assert process.poll() is None
            time.sleep(.01)
        else:raise AssertionError('No completed step before interruption')
        process.send_signal(signal.SIGINT)
        _,error=process.communicate(timeout=15)
        assert process.returncode==130,error.decode()
        manifest=json.loads((partial_out/'manifest.json').read_text())
        assert manifest['status']=='interrupted'
        snapshot=rows(partial_out/'cell_statistics_0.csv')
        assert snapshot and all(r['complete']=='0' for r in snapshot)
        assert all(0<float(r['observed_time'])<500 for r in snapshot)
    finally:
        if process.poll() is None:process.kill();process.wait()
    failure=root/'write-failure';failure.mkdir();(failure/'conservation.csv.tmp').mkdir()
    r=subprocess.run([str(binary),'--config',str(repo/'configs/checks/statistics_closed.json')],cwd=failure,capture_output=True)
    assert r.returncode!=0 and b'Cannot write statistics' in r.stderr
    reserved=root/'reserved.json';reserved.write_text(json.dumps(data|{'final_snapshot':'conservation_summary.json'}))
    r=subprocess.run([str(binary),'--config',str(reserved),'--output-dir',str(root/'reserved-out')],capture_output=True)
    assert r.returncode!=0 and b'reserved' in r.stderr
    # Invalid windows fail at parsing rather than partially running a simulation.
    for i,windows in enumerate([[[.1,.01]],[[-.1,.01]],[[0,1]],[[0]]]):
        bad=root/f'bad{i}.json';bad.write_text(json.dumps(data|{'statistics_windows':windows}))
        r=subprocess.run([str(binary),'--config',str(bad),'--output-dir',str(root/f'bad-out{i}')],capture_output=True)
        assert r.returncode!=0,windows
print('Statistics integration passed')
