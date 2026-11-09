from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_independent_campaign import seed_groups,execution_complete


class IndependentCampaignTest(unittest.TestCase):
    def test_groups_have_no_shared_seed(self):
        groups=seed_groups([123,456,789],10000)
        self.assertEqual(len(groups),12)
        self.assertEqual(len(set(s for g in groups for s in g)),36)
        self.assertEqual(groups[0],[123,456,789])
        with self.assertRaises(ValueError):seed_groups([10003,456,789],10000)
        with self.assertRaises(ValueError):seed_groups([123,456,789],2**31)

    def test_mc_completion_does_not_substitute_for_sparta(self):
        report={'runs':{f'C-{s}':dict(status='passed') for s in (123,456,789)}}
        self.assertFalse(execution_complete(report,['C'],[123,456,789]))
        report['sparta_runs']=report['runs'].copy()
        self.assertTrue(execution_complete(report,['C'],[123,456,789]))
        self.assertFalse(execution_complete(report,['C','D'],[123,456,789]))
        self.assertFalse(execution_complete({},[],[]))


if __name__=='__main__':unittest.main()
