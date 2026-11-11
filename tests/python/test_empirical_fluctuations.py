from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from shock_empirical_fluctuations import load_bank,fixture
from shock_validation import digest,save
from validate_shock_estimator import shape


class EmpiricalFluctuationTest(unittest.TestCase):
    def create(self,root):
        x=np.linspace(-14,14,71);mean=shape(x,'logistic')
        # Nonzero run means deliberately retained; ensemble mean exactly zero.
        offsets=np.linspace(-.02,.02,12)[:,None,None]
        noise=np.broadcast_to(offsets,(12,64,len(x))).copy()
        noise+=.003*np.sin(np.arange(64)*2*np.pi/64)[None,:,None]*np.cos(x)[None,None,:]
        noise-=noise.mean((0,1))
        file=root/'heldout.npz';np.savez(file,x=x,mean_profile=mean,fluctuations=noise)
        groups=[f'g{i:02d}' for i in range(8,12)];sources=[dict(group=g,seed=100+i*3+j) for i,g in enumerate(groups) for j in range(3)]
        manifest=root/'manifest.json';save(manifest,dict(status='conditional_empirical_banks_prepared',physical_acceptance='incomplete',
            protocol=dict(heldout_groups=groups),banks={'mc3d-heldout':dict(file=file.name,sha256=digest(file),shape=list(noise.shape),sources=sources,mean_profile=mean.tolist())}))
        return manifest,noise

    def test_fixture_retains_noise_offsets_and_native_time_order_without_clipping(self):
        with tempfile.TemporaryDirectory() as d:
            manifest,original=self.create(Path(d));bank=load_bank(manifest,'mc3d-heldout')
            trials,meta=fixture('gompertz',169100001,bank,windows=64,motion_amplitude=.05)
            for run,choice in zip(trials,meta['choices']):
                source=choice['source_index'];indices=choice['window_indices']
                residual=np.asarray([w['density_normalized'] for w in run])-np.asarray([shape(bank['x']-s,'gompertz') for s in choice['true_shifts']])
                np.testing.assert_allclose(residual,original[source,indices],rtol=1e-12,atol=1e-14)
                np.testing.assert_allclose(residual.mean(0),original[source].mean(0),rtol=1e-12,atol=1e-14)
                self.assertGreater(abs(residual.mean()),.001)
                self.assertEqual(indices,[(choice['start_window']+i)%64 for i in range(64)])
            self.assertEqual(fixture('gompertz',169100001,bank,64,.05),(trials,meta))
            half,details=fixture('logistic',169100002,bank,32,0.)
            self.assertTrue(all(len(run)==32 for run in half));self.assertTrue(all(set(c['true_shifts'])=={0.} for c in details['choices']))
            self.assertFalse(bank['fluctuations'].flags.writeable)

    def test_corruption_and_rehashed_invalid_dimensions_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest,noise=self.create(root);file=root/'heldout.npz';file.write_bytes(file.read_bytes()+b'changed')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):load_bank(manifest,'mc3d-heldout')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest,noise=self.create(root);file=root/'heldout.npz';np.savez(file,x=np.linspace(-14,14,71),mean_profile=np.zeros(71),fluctuations=noise[:,:63])
            data=json.loads(manifest.read_text());data['banks']['mc3d-heldout']['sha256']=digest(file);save(manifest,data)
            with self.assertRaisesRegex(ValueError,'dimensions'):load_bank(manifest,'mc3d-heldout')

    def test_duplicate_physical_sources_and_noncentered_bank_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest,noise=self.create(root);data=json.loads(manifest.read_text());sources=data['banks']['mc3d-heldout']['sources'];sources[1]['seed']=sources[0]['seed'];save(manifest,data)
            with self.assertRaisesRegex(ValueError,'distinct seeds'):load_bank(manifest,'mc3d-heldout')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest,noise=self.create(root);file=root/'heldout.npz';x=np.linspace(-14,14,71);np.savez(file,x=x,mean_profile=shape(x,'logistic'),fluctuations=noise+.01)
            data=json.loads(manifest.read_text());data['banks']['mc3d-heldout']['sha256']=digest(file);save(manifest,data)
            with self.assertRaisesRegex(ValueError,'ensemble centered'):load_bank(manifest,'mc3d-heldout')

if __name__=='__main__':unittest.main()
