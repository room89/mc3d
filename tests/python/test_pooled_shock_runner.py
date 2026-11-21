from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from run_pooled_shock import child_status

class PooledRunnerTest(unittest.TestCase):
    def test_execution_requires_both_solvers(self):
        plan=dict(binary_sha256='a',sparta_sha256='b')
        r=dict(plan,seeds=[1,2,3],variants=['C'],runs={f'C-{s}':{'status':'passed'} for s in (1,2,3)})
        self.assertEqual(child_status(r,[1,2,3],plan,2,True),'incomplete')
        r['sparta_runs']=r['runs'].copy()
        self.assertEqual(child_status(r,[1,2,3],plan,2,True),'execution_complete')
        self.assertEqual(child_status(r,[1,2,3],plan,3,True),'incomplete')

    def test_changed_identities_are_rejected(self):
        plan=dict(binary_sha256='a',sparta_sha256='b')
        r=dict(plan,seeds=[1,2,3],variants=['C'])
        with self.assertRaises(ValueError):child_status(r,[4,5,6],plan,2,False)
        with self.assertRaises(ValueError):child_status(dict(r,binary_sha256='other'),[1,2,3],plan,2,False)
        self.assertEqual(child_status(r,[1,2,3],plan,2,False),'prepared')
        self.assertEqual(child_status(r,[1,2,3],plan,3,False),'incomplete')
