"""Shock configuration, compatibility, replay and quick balance checks."""
import copy,json,subprocess,sys,tempfile,math,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from normal_shock import setup
from physical_config import prepare
from check_shock import inspect
from replay_run import replay
binary=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix='mc3d-shock-') as d:
    root=Path(d)
    def run(c,name):
        f=root/(name+'.json');f.write_text(json.dumps(c));out=root/name
        r=subprocess.run([str(binary),'--config',str(f),'--output-dir',str(out)],capture_output=True)
        return r,out
    for seed in (123,456,789):
        data,ref=setup(True,seed);c,_=prepare(data)
        r,out=run(c,str(seed));assert r.returncode==0,r.stderr.decode()
        assert inspect(out,c,ref)['status']=='passed'
        with (out/'initial.dat').open() as f:initial=list(csv.DictReader(f,delimiter=';'))
        for row in initial:
            state=c['initial_step']['left' if float(row['x'])<0 else 'right']
            expected=c['particles_per_cell']*state['density_ratio']
            assert int(row['N']) in (math.floor(expected),math.ceil(expected)),row
        for side in ('left','right'):
            selected=[r for r in initial if (float(r['x'])<0)==(side=='left')]
            n=sum(int(r['N']) for r in selected);state=c['initial_step'][side]
            for axis in 'xyz':
                u=sum(int(r['N'])*float(r['v'+axis]) for r in selected)/n
                assert abs(u-state['velocity']['xyz'.index(axis)])<6*math.sqrt(state['temperature']/n)

        if seed==123:assert replay(out/'manifest.json',root/'replay')==0
    for i,patch in enumerate([{'boundary_x_pos':'mirror'},{'initial_step':dict(c['initial_step'],x=.1234567)},{'initial_step':dict(c['initial_step'],x=100)},{'reservoirs':{'unknown':c['reservoirs']['x_neg']}},{'reservoirs':{'x_neg':{'density_ratio':1}}},{'reservoirs':{'x_neg':dict(c['reservoirs']['x_neg'],temperature=0)}},{'reservoirs':{'x_neg':dict(c['reservoirs']['x_neg'],density_ratio=-1)}},{'reservoirs':{'x_neg':dict(c['reservoirs']['x_neg'],velocity=[1,2])}}]):
        r,_=run(c|patch,'invalid'+str(i));assert r.returncode!=0,patch
    base=dict(seed=123,thread_pool_size=1,lx=1,ly=1,lz=1,cells_x=4,cells_y=1,cells_z=1,particles_per_cell=20,temperature=1,s=2,end_time=.02,geometry_type='none',final_snapshot='final.dat')
    r,old=run(base,'legacy');assert r.returncode==0
    state=dict(density_ratio=1,temperature=1,velocity=[2*math.sqrt(2),0,0])
    r,new=run(base|{'reservoirs':{f:state for f in ('x_neg','x_pos','y_neg','y_pos','z_neg','z_pos')}},'explicit');assert r.returncode==0
    assert (old/'final.dat').read_bytes()==(new/'final.dat').read_bytes()
print('Shock integration passed')
