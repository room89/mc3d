import json
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from prepare_alsmeyer_matrix import reference_configuration
from physical_config import prepare
from argon_vhs import KB


class AlsmeyerMatrixTest(unittest.TestCase):
    def evidence(self):
        return json.loads((Path(__file__).resolve().parents[2]/'validation/stage9b/experiment.json').read_text())

    def test_primary_domain_reservoirs_and_windows(self):
        evidence = self.evidence(); data, ref = reference_configuration('C', 30001, evidence)
        self.assertAlmostEqual(data['domain_m'][0], .04392)
        self.assertAlmostEqual(data['number_density_m3']*KB*300, .05*101325/760)
        self.assertEqual(data['reservoirs']['x_neg'], data['initial_step']['left'])
        self.assertEqual(data['reservoirs']['x_pos'], data['initial_step']['right'])
        self.assertEqual(ref['left'], data['reservoirs']['x_neg'])
        self.assertEqual(ref['right'], data['reservoirs']['x_pos'])
        tau = ref['transit_time_s']
        self.assertEqual(len(data['statistics_windows_s']), 64)
        self.assertAlmostEqual(data['statistics_windows_s'][0][0]/tau, 4)
        self.assertAlmostEqual(data['statistics_windows_s'][-1][1]/tau, 12)
        cfg, units = prepare(data)
        self.assertAlmostEqual(cfg['kn']*data['domain_m'][0], ref['mean_free_path_m'])
        self.assertEqual(cfg['cells_x'], 800)
        self.assertEqual(cfg['particles_per_cell'], 160)

    def test_convergence_variants_change_only_intended_axes(self):
        evidence = self.evidence(); c, cr = reference_configuration('C', 30001, evidence)
        f, fr = reference_configuration('F', 30001, evidence)
        self.assertEqual(f['domain_m'][0], 2*c['domain_m'][0])
        self.assertEqual(f['numerics']['cells_x'], 2*c['numerics']['cells_x'])
        self.assertEqual(fr['transit_time_s'], cr['transit_time_s'])
        self.assertEqual(f['statistics_windows_s'], c['statistics_windows_s'])
        e, _ = reference_configuration('E', 30001, evidence)
        self.assertEqual(e['numerics']['collision_courant'], c['numerics']['collision_courant']/2)
        self.assertEqual(e['numerics']['cu'], c['numerics']['cu']/2)
        a, _ = reference_configuration('A', 30001, evidence)
        self.assertEqual(a['numerics']['cells_x'], 400)
        self.assertEqual(a['domain_m'], c['domain_m'])


if __name__ == '__main__': unittest.main()
