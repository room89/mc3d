import copy
import json
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from prepare_alsmeyer_matrix import reference_configuration
from physical_config import prepare
from export_alsmeyer_sparta import render


class AlsmeyerSpartaTest(unittest.TestCase):
    def case(self, variant='C'):
        evidence = json.loads((Path(__file__).resolve().parents[2]/'validation/stage9b/experiment.json').read_text())
        d, r = reference_configuration(variant, 30001, evidence)
        _, u = prepare(d)
        return d, r, u

    def test_exported_geometry_gas_reservoirs_and_measurement_schedule(self):
        d, r, u = self.case(); files, meta = render(d, r, u['reference_gas'])
        lines = files['in.shock'].splitlines()
        box = [float(x) for x in next(x for x in lines if x.startswith('create_box ')).split()[1:]]
        self.assertAlmostEqual(box[1]-box[0], .04392)
        self.assertEqual(next(x for x in lines if x.startswith('create_grid ')), 'create_grid 800 1 1')
        for name, state in zip(('up', 'down'), d['reservoirs'].values()):
            row = next(x for x in lines if x.startswith('mixture '+name+' ')).split()
            self.assertEqual(float(row[row.index('nrho')+1]), state['number_density_m3'])
            self.assertEqual(float(row[row.index('temp')+1]), state['temperature_k'])
        self.assertEqual(float(files['argon.vss'].split()[1]), u['reference_gas']['diameter_m'])
        self.assertEqual(float(files['argon.vss'].split()[2]), u['reference_gas']['omega'])
        self.assertEqual(meta['fnum'], u['scales']['molecules_per_model_particle'])
        self.assertAlmostEqual(meta['burn_steps']*meta['timestep'], d['statistics_windows_s'][0][0])
        self.assertAlmostEqual((meta['burn_steps']+meta['measurement_steps'])*meta['timestep'], d['end_time_s'])
        self.assertEqual(meta['measurement_steps']//meta['dump_period_steps'], 64)
        self.assertLessEqual(meta['timestep'], min(meta['flight_dt_limit'], meta['collision_dt_limit']))

    def test_half_step_and_quick_are_not_confused(self):
        d, r, u = self.case(); _, normal = render(d, r, u['reference_gas'])
        e, er, eu = self.case('E'); _, half = render(e, er, eu['reference_gas'])
        self.assertAlmostEqual(half['flight_dt_limit'], normal['flight_dt_limit']/2)
        self.assertAlmostEqual(half['collision_dt_limit'], normal['collision_dt_limit']/2)
        _, quick = render(d, r, u['reference_gas'], quick=True)
        self.assertEqual((quick['burn_steps'], quick['measurement_steps']), (10, 20))
        self.assertEqual(quick['physical_acceptance'], 'incomplete')
        self.assertEqual(quick['timestep'], normal['timestep'])

    def test_mismatched_states_and_windows_are_rejected(self):
        d, r, u = self.case(); bad = copy.deepcopy(d)
        bad['initial_step']['left']['temperature_k'] += 1
        with self.assertRaisesRegex(ValueError, 'states differ'): render(bad, r, u['reference_gas'])
        bad = copy.deepcopy(d); bad['statistics_windows_s'][0][0] = 0.
        with self.assertRaisesRegex(ValueError, '64 windows'): render(bad, r, u['reference_gas'])
        gas = dict(u['reference_gas'], reference_temperature_k=301)
        with self.assertRaisesRegex(ValueError, '300 K'): render(d, r, gas)


if __name__ == '__main__': unittest.main()
