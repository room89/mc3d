import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
try:
 import numpy as np
 from shock_alignment import align,template
except ImportError:np=None
@unittest.skipIf(np is None,'Requires numpy')
class AlignmentTest(unittest.TestCase):
 def fixture(self):
  x=np.linspace(-14,14,561)
  return [[dict(x_aligned=x.tolist(),density_normalized=(1/(1+np.exp(-(x-d)))).tolist()) for d in [-.7,0,.8]] for _ in range(3)]
 def test_translation(self):
  data,details=align(self.fixture())
  self.assertLess(np.max(abs(data[0,0]-data[0,2])),1e-4)
  self.assertAlmostEqual(details[0]['windows'][2]['front']-details[0]['windows'][0]['front'],1.5,places=3)
 def test_independent_template(self):
  trials=self.fixture();before=template(trials,0)
  trials[0][0]['density_normalized']=[0.]*561
  after=template(trials,0)
  np.testing.assert_array_equal(before,after)
 def test_invalid(self):
  trials=self.fixture();trials[0][0]['density_normalized']=[float('nan')]*561
  with self.assertRaises(ValueError):align(trials)
 def test_global_translation(self):
  trials=self.fixture();a,_=align(trials)
  for run in trials:
   for w in run:w['x_aligned']=[x+3.7 for x in w['x_aligned']]
  b,_=align(trials)
  np.testing.assert_allclose(a,b,atol=1e-4)
 def test_missing_coverage(self):
  trials=self.fixture()
  for run in trials:
   for w in run:
    w['x_aligned']=w['x_aligned'][220:341];w['density_normalized']=w['density_normalized'][220:341]
  with self.assertRaises(ValueError):align(trials)

@unittest.skipIf(np is None,'Requires numpy')
class RefinedAlignmentTest(unittest.TestCase):
 fixture=AlignmentTest.fixture
 def test_refined_template_independence(self):
  trials=self.fixture();before=template(trials,0,refine_passes=1)
  trials[0][0]['density_normalized']=[float('nan')]*561
  after=template(trials,0,refine_passes=1)
  np.testing.assert_array_equal(before,after)
 def test_refined_translation(self):
  data,details=align(self.fixture(),refine_passes=1)
  self.assertLess(np.max(abs(data[0,0]-data[0,2])),1e-4)
  self.assertEqual(details[0]['refine_passes'],1)
 def test_raw_coordinates(self):
  trials=self.fixture()
  for run in trials:
   for w in run:w['x']=w.pop('x_aligned')
  data,_=align(trials,refine_passes=1)
  self.assertTrue(np.isfinite(data).all())
 def test_unsupported_refinement(self):
  with self.assertRaises(ValueError):align(self.fixture(),refine_passes=2)
