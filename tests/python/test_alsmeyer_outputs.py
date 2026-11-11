import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from shock_validation import digest
from verify_alsmeyer_outputs import checked_files, mc3d_outputs, sparta_dump, sparta_outputs


class OutputVerificationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def save(self, name, obj):
        path = self.root/name; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(obj)); return path

    def fixture_mc(self):
        cfg = dict(seed=123, cells_x=1, cells_y=1, cells_z=1, end_time=96,
                   statistics_windows=[[32+i, 33+i] for i in range(64)])
        for name, data in [('solver.json', cfg), ('units.json', {}), ('input_si.json', {})]:
            self.save(name, data)
        inputs = []
        for role, name in [('config', 'solver.json'), ('units', 'units.json'), ('physical_input', 'input_si.json')]:
            p = self.root/name; q = self.root/'output/inputs'/name
            q.parent.mkdir(parents=True, exist_ok=True); q.write_bytes(p.read_bytes())
            inputs.append(dict(role=role, path='inputs/'+name, sha256=digest(q)))
        outputs = []
        for i in range(64):
            p = self.root/f'output/cell_statistics_{i}.csv'
            p.write_text('window;start;end;observed_time;complete;cell;status;Vgas;A0;density;Tstat\n'+
                         f'{i};{32+i};{33+i};1;1;0;valid;1;1;1;1\n')
            outputs.append(dict(path=p.name, sha256=digest(p)))
        manifest = dict(status='completed', error='', binary_sha256='binary', seed=123,
                        simulation_time=96, effective_config=cfg, units={}, inputs=inputs, outputs=outputs)
        self.save('output/manifest.json', manifest); return manifest

    def test_mc_complete_does_not_mean_physical_pass(self):
        self.fixture_mc(); result = mc3d_outputs(self.root, 'binary')
        self.assertEqual(result['windows'], 64); self.assertEqual(result['physical_acceptance'], 'incomplete')

    def test_mc_rehashed_incomplete_window_and_changed_effective_config_rejected(self):
        manifest = self.fixture_mc(); p = self.root/'output/cell_statistics_63.csv'
        p.write_text(p.read_text().replace(';1;0;valid;', ';0;0;valid;'))
        manifest['outputs'][-1]['sha256'] = digest(p); self.save('output/manifest.json', manifest)
        with self.assertRaisesRegex(ValueError, 'Incomplete or invalid'):
            mc3d_outputs(self.root, 'binary')
        manifest = self.fixture_mc(); manifest['effective_config']['seed'] = 124
        self.save('output/manifest.json', manifest)
        with self.assertRaisesRegex(ValueError, 'effective parameter'):
            mc3d_outputs(self.root, 'binary')

    def fixture_sparta(self):
        ref = dict(quick=False, steps_per_window=1, burn_steps=32, measurement_steps=64, dump_period_steps=1,
                   source=dict(numerics=dict(cells_x=2, cells_y=1, cells_z=1), domain_m=[2, 2, 2]))
        self.save('reference.json', ref)
        for name in ('in.shock', 'argon.species', 'argon.vss'): (self.root/name).write_text('fixture')
        (self.root/'run.log').write_text('fixture execution')
        outputs = []
        fields = ' '.join(f'f_avg[{i}]' for i in range(1, 9))
        for step in range(33, 97):
            p = self.root/f'profile.{step}'
            p.write_text(f'ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF CELLS\n2\nITEM: BOX BOUNDS oo pp pp\n-1 1\n-1 1\n-1 1\nITEM: CELLS id xlo xhi ylo yhi zlo zhi {fields}\n'+
                         '1 -1 0 -1 1 -1 1 1 0 0 0 1 1 1 1\n2 0 1 -1 1 -1 1 1 0 0 0 1 1 1 1\n')
            outputs.append(dict(path=p.name, sha256=digest(p)))
        execution = dict(returncode=0, binary_sha256='binary', outputs=outputs, log_sha256=digest(self.root/'run.log'),
                         input_hashes={n:digest(self.root/n) for n in ('in.shock','argon.species','argon.vss','reference.json')})
        return ref, execution

    def test_sparta_missing_window_and_nonzero_exit_rejected(self):
        ref, execution = self.fixture_sparta()
        self.assertEqual(sparta_outputs(self.root, execution, 'binary')['physical_acceptance'], 'incomplete')
        changed = copy.deepcopy(execution); changed['outputs'].pop()
        with self.assertRaisesRegex(ValueError, 'measurement missing'):
            sparta_outputs(self.root, changed, 'binary')
        changed = copy.deepcopy(execution); changed['returncode'] = 1
        with self.assertRaisesRegex(ValueError, 'Successful production'):
            sparta_outputs(self.root, changed, 'binary')

    def test_sparta_rehashed_duplicate_coordinates_and_wrong_timestep_rejected(self):
        ref, execution = self.fixture_sparta(); p = self.root/'profile.96'
        p.write_text(p.read_text().replace('2 0 1 ', '2 -1 0 '))
        execution['outputs'][-1]['sha256'] = digest(p)
        with self.assertRaisesRegex(ValueError, 'Duplicate SPARTA cell'):
            sparta_outputs(self.root, execution, 'binary')
        with self.assertRaisesRegex(ValueError, 'header, time'):
            sparta_dump(p, ref, 95)

    def test_hash_and_path_escape_are_rejected(self):
        (self.root/'artifact').write_text('new data')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            checked_files(self.root, [dict(path='artifact', sha256='old')])
        with self.assertRaisesRegex(ValueError, 'Unsafe'):
            checked_files(self.root, [dict(path='../artifact', sha256='anything')])


if __name__ == '__main__': unittest.main()
