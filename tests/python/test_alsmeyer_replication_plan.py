from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from plan_alsmeyer_replication import assignments


class ReferenceReplicationPlanTest(unittest.TestCase):
    def test_existing_subgroup_reused_exactly_once(self):
        groups=assignments([30001,30002,30003]);flat=[s for g in groups for triple in g for s in triple]
        self.assertEqual(len(groups),12);self.assertTrue(all(len(g)==6 for g in groups))
        self.assertTrue(all(len(t)==3 for g in groups for t in g))
        self.assertEqual(groups[0][0],[30001,30002,30003])
        self.assertEqual(len(flat),216);self.assertEqual(len(set(flat)),216)
        self.assertEqual(groups[-1][-1],[40210,40211,40212])

    def test_overlap_and_invalid_seed_ranges_rejected(self):
        for existing,start in [([30001,30001,30003],40000),([30001,30002,30003],30001),
                               ([30001,30002,30003],2**31-100),([0,2,3],40000),([True,2,3],40000)]:
            with self.assertRaises(ValueError):assignments(existing,start)

    def test_three_seeds_are_required_in_existing_subgroup(self):
        with self.assertRaises(ValueError):assignments([1,2])


if __name__=='__main__':unittest.main()
