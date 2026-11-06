import sys,math,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from normal_shock import downstream,setup,MASS,KB
from physical_config import prepare
class NormalShockTest(unittest.TestCase):
    def test_limit(self):
        s=downstream(1+1e-8,1,300)
        self.assertAlmostEqual(s['number_density_m3'],1,places=6)
        self.assertAlmostEqual(s['temperature_k'],300,places=4)
    def test_invalid(self):
        for m,n,t in [(1,1,300),(0,1,300),(3,0,300),(3,1,-1),(float('nan'),1,300)]:
            with self.assertRaises(ValueError):downstream(m,n,t)
    def test_generator_and_units(self):
        data,r=setup(True);c,u=prepare(data)
        self.assertEqual(c['reservoirs']['x_neg']['density_ratio'],1)
        self.assertGreater(c['reservoirs']['x_pos']['density_ratio'],1)
        self.assertAlmostEqual(c['initial_step']['right']['temperature'],data['initial_step']['right']['temperature_k']/300)
        self.assertEqual(c['boundary_y_neg'],'periodic')
        self.assertAlmostEqual(c['snapshot_interval'],c['end_time']/10)
        self.assertEqual(u['reference_gas']['omega'],c['vhs_omega'])
