"""Engineering normal shock setup; not Alsmeyer experimental conditions."""
import argparse,json,math
from pathlib import Path
from physical_config import prepare
from argon_vhs import MASS,KB,fit

def downstream(mach,n,t):
    if not all(math.isfinite(x) and x>0 for x in (mach,n,t)) or mach<=1:raise ValueError('Require M>1, n>0, T>0')
    gamma=5/3;r=(gamma+1)*mach**2/((gamma-1)*mach**2+2)
    pressure=(2*gamma*mach**2-(gamma-1))/(gamma+1)
    u=mach*math.sqrt(gamma*KB*t/MASS)
    return dict(number_density_m3=n*r,temperature_k=t*pressure/r,velocity_m_s=[u/r,0,0])

def setup(quick=False,seed=123):
    n=1e20;t=300.;u=3.8*math.sqrt(5/3*KB*t/MASS)
    gas,_=fit();w=gas['vhs_omega'];d=gas['vhs_diameter_m']
    power=2-2*w
    kernel=math.pi*d*d/math.gamma(2.5-w)*(4*KB*t/MASS)**.5*math.gamma((power+3)/2)/math.gamma(1.5)
    mean_free_path=math.sqrt(8*KB*t/(math.pi*MASS))/(n*kernel)
    length=40*mean_free_path;tau=length/u;end=(.05 if quick else 2)*tau
    left=dict(number_density_m3=n,temperature_k=t,velocity_m_s=[u,0,0]);right=downstream(3.8,n,t)
    for f in (lambda s:s['number_density_m3']*s['velocity_m_s'][0],lambda s:s['number_density_m3']*(MASS*s['velocity_m_s'][0]**2+KB*s['temperature_k']),lambda s:s['number_density_m3']*s['velocity_m_s'][0]*(.5*MASS*s['velocity_m_s'][0]**2+2.5*KB*s['temperature_k'])):
        assert math.isclose(f(left),f(right),rel_tol=1e-12)
    data=dict(temperature_k=t,number_density_m3=n,speed_m_s=u,domain_m=[length,.1*mean_free_path,.1*mean_free_path],end_time_s=end,snapshot_interval_s=end/10,gas={'model':'vhs'},geometry={'type':'none'},reservoirs={'x_neg':left,'x_pos':right},initial_step={'x_m':0,'left':left,'right':right},numerics=dict(cells_x=40 if quick else 400,cells_y=1,cells_z=1,particles_per_cell=20 if quick else 80,cu=.5,collision_courant=.1,thread_pool_size=4,seed=seed,boundary_x_neg='hyperfree',boundary_x_pos='hyperfree',boundary_y_neg='periodic',boundary_y_pos='periodic',boundary_z_neg='periodic',boundary_z_pos='periodic'),outputs=dict(precompute_snapshot='initial.dat',final_snapshot='final.dat',final_cell_snapshot='',write_default_snapshot_before_compute=False,write_default_snapshot_after_compute=False,write_speed_before=False,write_speed_after=False,write_times=False),statistics_windows_s=[[i*end/10,(i+1)*end/10] for i in range(10)])
    return data,dict(description=__doc__,mean_free_path_m=mean_free_path,transit_time_s=tau,left=left,right=right)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--quick',action='store_true');p.add_argument('--seed',type=int,default=123);a=p.parse_args()
    data,reference=setup(a.quick,a.seed);config,units=prepare(data);a.output_dir.mkdir(parents=True,exist_ok=False)
    config['units_file']=str((a.output_dir/'units.json').resolve())
    units['source_file']=str((a.output_dir/'input_si.json').resolve())
    for name,value in [('input_si.json',data),('solver.json',config),('units.json',units),('reference.json',reference)]:
        (a.output_dir/name).write_text(json.dumps(value,indent=2)+'\n')
