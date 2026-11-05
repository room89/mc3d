import copy
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from physical_config import prepare
from argon_vhs import MASS,KB

class UnitsTest(unittest.TestCase):
    def case(self):return dict(temperature_k=300,number_density_m3=1e20,speed_m_s=500,domain_m=[.01,.02,.03],end_time_s=1e-5,gas=dict(model='hs',diameter_m=4e-10))
    def test_scales_and_hs_rate(self):
        d=self.case();c,u=prepare(d);s=u['scales']
        self.assertAlmostEqual(c['kn']*d['domain_m'][0],1/(math.sqrt(2)*d['number_density_m3']*math.pi*(4e-10)**2),places=14)
        self.assertAlmostEqual(c['s']*math.sqrt(2*KB*300/MASS),500,places=12)
        self.assertAlmostEqual(c['end_time']*s['time_s'],1e-5,places=15)
        self.assertAlmostEqual(s['heat_capacity_j_k']*300,s['energy_j'],places=20)
        self.assertAlmostEqual(s['particle_energy_j'],KB*300,places=30)
        for scale in s.values():self.assertAlmostEqual((1.234/scale)*scale,1.234,places=12)
    def test_equivalent_scales_and_weight(self):
        d=self.case();a,u=prepare(d);d.update(length_unit_m=.02,temperature_unit_k=1200);b,v=prepare(d)
        for key,scale in [('lx','length_m'),('temperature','temperature_k'),('end_time','time_s')]:self.assertAlmostEqual(a[key]*u['scales'][scale],b[key]*v['scales'][scale],places=12)
        self.assertEqual(a['kn'],b['kn']);self.assertEqual(a['s'],b['s'])
        d['numerics']={'particles_per_cell':100};_,w=prepare(d)
        self.assertAlmostEqual(w['scales']['molecules_per_model_particle']*2,u['scales']['molecules_per_model_particle'])
    def test_invalid(self):
        for field,value in [('kn',.1),('temperature_k',float('nan')),('speed_m_s',-1),('domain_m',[1,0,1]),('numerics',{'cells_x':0}),('gas',{'model':'unknown'})]:
            d=self.case();d[field]=value
            with self.subTest(field=field),self.assertRaises((ValueError,KeyError)):prepare(d)
    def test_geometry(self):
        for kind,params in [('none',{}),('cube',dict(x_m=0,width_m=.001,length_m=.001,height_m=.001)),('wedge',dict(x_m=0,width_m=.001,length_m=.001)),('pyramid',dict(x_m=0,width_m=.001,length_m=.001,height_m=.001)),('cylinder',dict(x_m=0,radius_m=.001,length_m=.001)),('stl',dict(file='body.stl',file_unit_m=.001))]:
            d=self.case();d['geometry']={'type':kind,**params};c,_=prepare(d)
            self.assertIn('geometry_type',c)
        d['geometry']={'type':'stl','file':'body.stl'}
        with self.assertRaises(KeyError):prepare(d)
    def test_free_molecular_force_and_wall_energy_scales(self):
        # A cold incident beam on a stationary absorbing plane transfers m*U
        # per molecule and n*U*A molecules per second.
        d=self.case();area=1e-6;capacity=2e-9;d.update(force_reference_area_m2=area,wall_heat_capacity_j_k=capacity)
        c,u=prepare(d);s=u['scales'];velocity=d['speed_m_s']/s['velocity_m_s']
        physical_force=velocity**2*c['force_reference_area']*s['force_n']
        expected=d['number_density_m3']*MASS*d['speed_m_s']**2*area
        self.assertAlmostEqual(physical_force/expected,1,places=12)
        delta_temperature=17
        normalized_heat=c['wall_heat_capacity']*delta_temperature/s['temperature_k']
        self.assertAlmostEqual(normalized_heat*s['energy_j']/(capacity*delta_temperature),1,places=12)
