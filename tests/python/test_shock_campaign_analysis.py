import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from reanalyze_shock_campaign import merged_pairs,plateau_checks,main
from shock_alignment import initial


class CampaignAnalysisTest(unittest.TestCase):
    def test_pairing_uses_physical_coordinates_before_realignment(self):
        x=np.linspace(-10,10,201);a=1/(1+np.exp(-(x-.2)));b=1/(1+np.exp(-(x+.3)))
        rows=[dict(x_aligned=(x-f).tolist(),front=f,density_normalized=y.tolist()) for y,f in ((a,.2),(b,-.3))]
        merged=merged_pairs([rows])[0][0];xx,yy,_=initial(merged)
        np.testing.assert_allclose(xx,x,atol=1e-14)
        np.testing.assert_allclose(yy,(a+b)/2,atol=1e-14)
        with self.assertRaises(ValueError):merged_pairs([rows[:1]])
        rows[1]['x_aligned']=(x+.4).tolist()
        with self.assertRaisesRegex(ValueError,'different physical grids'):merged_pairs([rows])

    def test_plateau_outliers_remain_visible_when_means_pass(self):
        groups=[dict(full=[0.,0.],first=[0.,0.],last=[0.,0.],boundary_min=12.,window_outliers=2,window_count=192) for _ in range(12)]
        r=plateau_checks(groups)
        self.assertEqual(r['window_outliers'],24)
        self.assertTrue(all(g['passed'] for g in r['mean']))
        self.assertTrue(r['boundaries_pass'])
        groups[0]['boundary_min']=9.
        self.assertFalse(plateau_checks(groups)['boundaries_pass'])

    def test_preparation_is_not_analyzable_as_completed_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'campaign.json').write_text(json.dumps({'status':'prepared'}))
            output=root/'result.json'
            with patch.object(sys,'argv',['analysis',str(root),'--output',str(output),'--coverage','unused','--half-coverage','unused']):
                with self.assertRaisesRegex(ValueError,'Completed execution'):main()
            self.assertFalse(output.exists())


if __name__=='__main__':unittest.main()
