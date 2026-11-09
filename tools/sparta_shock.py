"""Export the same SI gas and reservoirs to pinned SPARTA, VHS via alpha=1."""
import argparse,json,math
from pathlib import Path
from shock_validation import configuration,save
from argon_vhs import fit,MASS,KB

def export(folder,variant,seed,quick=False):
    data,ref=configuration(variant,seed);gas,_=fit();num=data['numerics'];L=data['domain_m'];n=data['number_density_m3']
    folder.mkdir(parents=True,exist_ok=False)
    weight=n*math.prod(L)/(num['cells_x']*num['cells_y']*num['cells_z']*num['particles_per_cell'])
    states=[data['reservoirs'][f] for f in ('x_neg','x_pos')]
    minimum_cell=min(L[i]/num[f'cells_{a}'] for i,a in enumerate('xyz'))
    dt=.5*minimum_cell/max(abs(s['velocity_m_s'][0])+math.sqrt(2*KB*s['temperature_k']/MASS) for s in states)
    from sparta_preflight import rate
    dt=min(dt,.1/max(s['number_density_m3']*rate(s['temperature_k'],gas) for s in states))
    if variant=='E':dt/=2
    tau=ref['transit_time_s'];perwindow=max(1,math.ceil(tau/8/dt));dt=tau/(8*perwindow)
    (folder/'argon.species').write_text(f'Ar 39.948 {MASS:.17g} 0 0 0 0 0 1 0\n')
    (folder/'argon.vss').write_text(f"Ar {gas['vhs_diameter_m']:.17g} {gas['vhs_omega']:.17g} 300 1\n")
    lines=[f'seed {seed}','dimension 3','boundary o p p',f'create_box {-L[0]/2} {L[0]/2} {-L[1]/2} {L[1]/2} {-L[2]/2} {L[2]/2}',f"create_grid {num['cells_x']} {num['cells_y']} {num['cells_z']}",'species argon.species Ar',f'global fnum {weight:.17g}']
    for name,s in zip(('up','down'),states):lines.append(f"mixture {name} Ar nrho {s['number_density_m3']:.17g} temp {s['temperature_k']:.17g} vstream "+' '.join(map(str,s['velocity_m_s'])))
    for name,lo,hi in [('up',-L[0]/2,0),('down',0,L[0]/2)]:
        lines += [f'region {name} block {lo} {hi} {-L[1]/2} {L[1]/2} {-L[2]/2} {L[2]/2}',f'create_particles {name} n 0 region {name}']
    lines+=['collide vss all argon.vss','fix inlet emit/face up xlo','fix outlet emit/face down xhi',f'timestep {dt:.17g}','stats 100','stats_style step cpu np ncoll',f'run {10 if quick else 32*perwindow}','compute moments grid all all nrho u v w usq vsq wsq','compute thermal thermal/grid all all temp',f'fix avg ave/grid all 1 {10 if quick else perwindow} {10 if quick else perwindow} c_moments[*] c_thermal[*]',f'dump profiles grid all {10 if quick else perwindow} profile.* id xlo xhi ylo yhi zlo zhi f_avg[*]','dump_modify profiles format float %.17g',f'run {20 if quick else 64*perwindow}']
    (folder/'in.shock').write_text('\n'.join(lines)+'\n');save(folder/'reference.json',dict(source=data,reference=ref,fnum=weight,timestep=dt,steps_per_window=perwindow,quick=quick))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--variant',choices=['C','D','E'],default='C');p.add_argument('--seed',type=int,default=123);p.add_argument('--quick',action='store_true');a=p.parse_args();export(a.output_dir,a.variant,a.seed,a.quick)
