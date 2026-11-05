"""Executable-level manifest tests; invoked by CTest with the built binary."""
import hashlib
import csv
import math
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from physical_config import prepare
from replay_run import replay
binary=Path(sys.argv[1]).resolve()

def run(config,root,name,extra=()):
    path=root/(name+'.json');path.write_text(json.dumps(config))
    output=root/name
    result=subprocess.run([str(binary),'--config',str(path),'--output-dir',str(output),*extra],capture_output=True)
    return result,output

with tempfile.TemporaryDirectory(prefix='mc3d manifest ') as tmp:
    root=Path(tmp)
    base=dict(seed=123,thread_pool_size=1,cells_x=2,cells_y=2,cells_z=2,particles_per_cell=10,end_time=.001,final_snapshot='final.dat',geometry_type='none')
    result,out=run(base,root,'original');assert result.returncode==0,result.stderr.decode()
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['status']=='completed' and manifest['seed']==123 and manifest['steps']>0
    assert manifest['build']['revision'] and manifest['build']['compiler']
    assert manifest['binary_sha256']==hashlib.sha256(binary.read_bytes()).hexdigest()
    assert replay(out/'manifest.json',root/'repeat')==0
    repeat=json.loads((root/'repeat'/'manifest.json').read_text())
    first={x['path']:x['sha256'] for x in manifest['outputs'] if x['kind']=='physical'}
    second={x['path']:x['sha256'] for x in repeat['outputs'] if x['kind']=='physical'}
    assert first==second,(first,second)
    altered=dict(manifest);altered['binary_sha256']='0'*64
    alternative=out/'changed-binary.json';alternative.write_text(json.dumps(altered))
    try:replay(alternative,root/'wrong-binary');raise AssertionError('Binary mismatch accepted')
    except ValueError:pass
    assert replay(alternative,root/'allowed-binary',allow_different=True)==0
    allowed=json.loads((root/'allowed-binary'/'manifest.json').read_text())
    assert allowed['provenance']['allow_different_binary'] is True
    result,_=run(base,root,'original');assert result.returncode!=0
    result,out=run({**base,'final_snapshot':'missing/final.dat'},root,'failure')
    assert result.returncode!=0 and json.loads((out/'manifest.json').read_text())['status']=='failed'
    result,_=run({**base,'final_snapshot':'../escape.dat'},root,'escape');assert result.returncode!=0
    result,_=run({**base,'final_snapshot':'manifest.json'},root,'reserved');assert result.returncode!=0
    result,out=run({k:v for k,v in base.items() if k!='seed'},root,'random')
    assert result.returncode==0 and isinstance(json.loads((out/'manifest.json').read_text())['seed'],int)
    # The unit table must reject a physical/numerical CLI override.
    si=dict(temperature_k=300,number_density_m3=1e20,speed_m_s=500,domain_m=[.01]*3,end_time_s=1e-7,numerics=dict(cells_x=2,cells_y=2,cells_z=2,particles_per_cell=10,seed=123,thread_pool_size=1))
    config,units=prepare(si);up=root/'units.json';up.write_text(json.dumps(units));config['units_file']=str(up)
    result,out=run(config,root,'si');assert result.returncode==0,result.stderr.decode()
    result,out=run(config,root,'mismatch',('--kn','1'));assert result.returncode!=0
    assert json.loads((out/'manifest.json').read_text())['status']=='failed'
    corrupted=json.loads(up.read_text());corrupted['scales']['force_n']*=2
    up.write_text(json.dumps(corrupted))
    result,out=run(config,root,'bad-scales');assert result.returncode!=0
    up.write_text(json.dumps(units))
    # Same physical flow with power-of-two alternative normalization.
    si['gas']={'model':'none'};si['outputs']={'final_snapshot':'final.dat'}
    physical=[]
    for index,scales in enumerate(({},dict(length_unit_m=.02,temperature_unit_k=1200))):
        config,units=prepare({**si,**scales});up=root/f'eq-units-{index}.json';up.write_text(json.dumps(units));config['units_file']=str(up)
        result,eqout=run(config,root,f'equivalent-{index}');assert result.returncode==0,result.stderr.decode()
        with (eqout/'final.dat').open() as f:rows=list(csv.DictReader(f,delimiter=';'))
        multipliers=dict(x=units['scales']['length_m'],y=units['scales']['length_m'],z=units['scales']['length_m'],T=units['scales']['temperature_k'],vx=units['scales']['velocity_m_s'],vy=units['scales']['velocity_m_s'],vz=units['scales']['velocity_m_s'],E=units['scales']['particle_energy_j'],ro=units['scales']['number_density_m3'],N=1)
        physical.append([{k:float(v)*multipliers[k] for k,v in row.items()} for row in rows])
    for first,second in zip(*physical):
        for key in first:assert math.isclose(first[key],second[key],rel_tol=3e-7,abs_tol=1e-30),(key,first[key],second[key])
    # Graceful signal handling, without sleeps in the simulation thread.
    path=root/'long.json';path.write_text(json.dumps({**base,'end_time':10000}))
    with (root/'signal.log').open('w') as log:
        process=subprocess.Popen([str(binary),'--config',str(path),'--output-dir',str(root/'signal')],stdout=log,stderr=log)
        deadline=time.monotonic()+20
        while time.monotonic()<deadline:
            try:
                if json.loads((root/'signal'/'manifest.json').read_text()).get('steps',0)>0:break
            except (FileNotFoundError,json.JSONDecodeError):pass
            time.sleep(.01)
        process.send_signal(signal.SIGTERM);process.wait(timeout=20)
    assert json.loads((root/'signal'/'manifest.json').read_text())['status']=='interrupted'
    # Corruption must be detected before starting a replay.
    archived=out/'manifest.json'
    origin=root/'original'/'manifest.json';data=json.loads(origin.read_text())
    (origin.parent/data['inputs'][0]['path']).write_text('corrupt')
    try:replay(origin,root/'corrupt');raise AssertionError('corruption accepted')
    except ValueError:pass
print('Manifest integration checks passed')
