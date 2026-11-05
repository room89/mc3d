"""Prepare normalized solver input from SI. All dimensional fields have SI suffixes."""
import argparse
import json
import math
from pathlib import Path
from argon_vhs import fit, MASS, KB


def positive(x, name):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0:
        raise ValueError(f'{name}: expected finite positive number')
    return x


def prepare(data):
    def finite(value):
        if isinstance(value, float) and not math.isfinite(value):raise ValueError('Nonfinite input')
        if isinstance(value,dict):
            for x in value.values():finite(x)
        if isinstance(value,list):
            for x in value:finite(x)
    finite(data)
    allowed={'temperature_k','number_density_m3','speed_m_s','flow_angle_rad','domain_m','apex_m','length_unit_m','temperature_unit_k','wall_temperature_k','wall_heat_capacity_j_k','force_reference_area_m2','end_time_s','snapshot_interval_s','gas','geometry','numerics','outputs','provenance'}
    if set(data)-allowed: raise ValueError(f'Unknown SI fields: {set(data)-allowed}')
    t=positive(data['temperature_k'],'temperature_k'); n=positive(data['number_density_m3'],'number_density_m3')
    domain=data['domain_m']
    if len(domain)!=3: raise ValueError('domain_m requires three lengths')
    for x in domain: positive(x,'domain_m')
    l0=positive(data.get('length_unit_m',domain[0]),'length_unit_m');t0=positive(data.get('temperature_unit_k',t),'temperature_unit_k')
    gas=data.get('gas',{'model':'vhs'})
    if set(gas)-{'model','mass_kg','diameter_m','reference_temperature_k','omega'}: raise ValueError('Unknown gas field')
    model=gas.get('model','vhs');profile,_=fit();m=positive(gas.get('mass_kg',MASS),'mass_kg')
    if model not in ('hs','vhs','none'): raise ValueError('Unknown collision model')
    v0=math.sqrt(KB*t0/m);time0=l0/v0;rho=n*m
    units={'length_m':l0,'temperature_k':t0,'velocity_m_s':v0,'time_s':time0,'number_density_m3':n,'mass_density_kg_m3':rho,'force_n':rho*v0*v0*l0*l0,'energy_j':rho*v0*v0*l0**3,'heat_capacity_j_k':rho*v0*v0*l0**3/t0,'particle_energy_j':m*v0*v0,'pressure_pa':rho*v0*v0,'area_m2':l0*l0}
    c=dict(zip(('lx','ly','lz'),(x/l0 for x in domain)))
    speed=data['speed_m_s']
    if not math.isfinite(speed) or speed<0: raise ValueError('speed_m_s must be finite and nonnegative')
    c.update(collision_model=model,temperature=t/t0,wall_temperature=positive(data.get('wall_temperature_k',t),'wall_temperature_k')/t0,s=speed/math.sqrt(2*KB*t/m),alpha=data.get('flow_angle_rad',0),end_time=positive(data['end_time_s'],'end_time_s')/time0)
    if model=='none':c['kn']=1
    else:
        d=positive(gas.get('diameter_m',profile['vhs_diameter_m']) if model=='vhs' else gas['diameter_m'],'diameter_m')
        w=gas.get('omega',profile['vhs_omega']) if model=='vhs' else .5
        tr=positive(gas.get('reference_temperature_k',300),'reference_temperature_k')
        if not math.isfinite(w) or not .5<=w<=1: raise ValueError('omega must be in [0.5,1]')
        power=2-2*w
        mean_k=math.pi*d*d/math.gamma(2.5-w)*(4*KB*tr/m)**(w-.5)*(4*KB*tr/m)**(power/2)*math.gamma((power+3)/2)/math.gamma(1.5)
        c['kn']=math.sqrt(8*KB*tr/(math.pi*m))/(n*mean_k*domain[0])
        if model=='vhs':c.update(molecular_mass_kg=m,vhs_diameter_m=d,vhs_omega=w,vhs_reference_temperature_k=tr,temperature_unit_k=t0)
    for key,target,scale in [('wall_heat_capacity_j_k','wall_heat_capacity',units['heat_capacity_j_k']),('force_reference_area_m2','force_reference_area',l0*l0),('snapshot_interval_s','snapshot_interval',time0)]:
        if key in data:
            c[target]=positive(data[key],key)/scale
    if speed==0 and 'force_reference_area' in c: raise ValueError('Force coefficients require nonzero flow speed')
    numerics=data.get('numerics',{})
    if set(numerics)-{'cells_x','cells_y','cells_z','particles_per_cell','cu','seed','thread_pool_size','runtime_control','boundary_x_neg','boundary_x_pos','boundary_y_neg','boundary_y_pos','boundary_z_neg','boundary_z_pos'}: raise ValueError('Unknown numerical field')
    for k in ('cells_x','cells_y','cells_z','particles_per_cell','thread_pool_size'):
        value=numerics.get(k,dict(cells_x=20,cells_y=20,cells_z=3,particles_per_cell=50,thread_pool_size=4)[k])
        if type(value) is not int or value<=0:raise ValueError(k+': expected positive integer')
        c[k]=value
    if 'cu' in numerics:positive(numerics['cu'],'cu')
    if 'seed' in numerics and (type(numerics['seed']) is not int or not 0<=numerics['seed']<=4294967295):raise ValueError('seed: expected uint32')
    for key,value in numerics.items():
        if key.startswith('boundary_') and value not in ('none','mirror','periodic','free','hyperfree'):raise ValueError('Unknown boundary type')
    c.update(numerics)
    if 'apex_m' in data:
        if len(data['apex_m'])!=3:raise ValueError('apex_m requires three coordinates')
        for axis,x in zip('xyz',data['apex_m']):c['apex_'+axis]=x/l0
    outputs=data.get('outputs',{})
    if set(outputs)-{'final_snapshot','precompute_snapshot','final_cell_snapshot','snapshots_binary','write_default_snapshot_before_compute','write_default_snapshot_after_compute','write_speed_before','write_speed_after','write_times'}: raise ValueError('Unknown output field')
    c.update(outputs)
    geo=data.get('geometry',{'type':'none'})
    kind=geo.get('type','none')
    lengths={'wedge':{'x','width','length'},'cube':{'x','width','length','height'},'pyramid':{'x','width','length','height'},'cylinder':{'x','radius','length'},'none':set(),'stl':set()}
    if kind not in ('none','cube','wedge'):raise ValueError('Gas volume currently supports only built-in cube and wedge; STL, cylinder and pyramid are unsupported')
    c['geometry_type']='none' if kind=='stl' else kind
    for k,value in geo.items():
        if k=='type':continue
        if k.endswith('_m') and k[:-2] in lengths[kind]:
            if k!='x_m':positive(value,k)
            c[f'geometry_{kind}_{k[:-2]}']=value/l0
        elif k=='angle_rad' and kind=='wedge':c['geometry_wedge_alpha']=value
        elif k=='segments' and kind=='cylinder':
            if type(value) is not int or value<3:raise ValueError('Cylinder segments must be an integer >=3')
            c['geometry_cylinder_segments']=value
        elif k=='file' and kind=='stl':c['geometry_file']=value
        elif k=='file_unit_m' and kind=='stl':positive(value,k)
        elif k=='scale':c['geometry_scale']=positive(value,k)
        elif k=='fragment_length_m':c['geometry_fragment_length']=positive(value,k)/l0
        elif k=='move_m':
            if len(value)!=3:raise ValueError('move_m requires three coordinates')
            for axis,x in zip('xyz',value):c['geometry_move_'+axis]=x/l0
        elif k in ('fix_polygons','reverse_normals'):c['geometry_'+k]=bool(value)
        else:raise ValueError(f'Unknown geometry field: {k}')
    if kind=='stl':c['geometry_scale']=positive(geo['file_unit_m'],'file_unit_m')*geo.get('scale',1)/l0
    else:
        for key in lengths[kind]:
            if key+'_m' not in geo:raise ValueError(f'Missing geometry field: {key}_m')
    if kind in ('cube','wedge'):
        x=c[f'geometry_{kind}_x'];length=c[f'geometry_{kind}_length'];width=c[f'geometry_{kind}_width']
        if kind=='wedge':
            angle=c.get('geometry_wedge_alpha',math.pi/18)
            if not 0<angle<math.pi/2:raise ValueError('Wedge angle must be in (0,pi/2)')
            half_height=length*math.sin(angle)
        else:half_height=c['geometry_cube_height']/2
        scale=c.get('geometry_scale',1)
        low=[x,-half_height,-width/2];high=[x+length,half_height,width/2]
        for i,axis in enumerate('xyz'):
            minimum=c.get('apex_'+axis,-c['l'+axis]/2);maximum=minimum+c['l'+axis]
            move=c.get('geometry_move_'+axis,0)
            if not minimum<low[i]*scale+move<high[i]*scale+move<maximum:
                raise ValueError('Body must lie strictly inside domain; clipped inflow faces are unsupported')
    count=c.get('particles_per_cell',50)*c.get('cells_x',20)*c.get('cells_y',20)*c.get('cells_z',3)
    units['molecules_per_model_particle']=n*math.prod(domain)/positive(count,'reference count')
    for key,value in units.items():positive(value,key)
    reference_gas={'model':model,'mass_kg':m}
    if model!='none':reference_gas.update(diameter_m=d,omega=w,reference_temperature_k=tr)
    return c,{'reference_gas':reference_gas,'schema_version':1,'scales':units,'normalized_parameters':c.copy(),'provenance':data.get('provenance'), 'kn_used':model!='none'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('--output-dir',required=True,type=Path);args=parser.parse_args()
    data=json.loads(args.input.read_text());config,units=prepare(data)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):raise ValueError('Preparation directory must be empty')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if config.get('runtime_control'):config['runtime_control']=str((args.input.resolve().parent/config['runtime_control']).resolve());units['normalized_parameters']['runtime_control']=config['runtime_control']
    if config.get('geometry_file'):config['geometry_file']=str((args.input.resolve().parent/config['geometry_file']).resolve());units['normalized_parameters']['geometry_file']=config['geometry_file']
    units['source_file']=str((args.output_dir/'input_si.json').resolve())
    if isinstance(units.get('provenance'),dict) and 'files' in units['provenance']:
        units['provenance']['files']=[str((args.input.resolve().parent/p).resolve()) for p in units['provenance']['files']]
    config['units_file']=str((args.output_dir/'units.json').resolve())
    for name,obj in [('solver.json',config),('units.json',units),('input_si.json',data)]:
        (args.output_dir/name).write_text(json.dumps(obj,indent=2,allow_nan=False,ensure_ascii=False)+'\n')
