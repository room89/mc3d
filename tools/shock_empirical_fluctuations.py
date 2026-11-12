"""Conditional synthetic stress fixtures from verified historical fluctuation banks.

A bank contains historical fluctuations, front motion and run offsets. It is not
an isolated measurement-noise model. Reusing it does not create new independent
physical runs. Only the added synthetic clean profile defines fixture truth.
"""
import json
import math
from pathlib import Path
import numpy as np
from shock_validation import digest
from validate_shock_estimator import shape


def load_bank(manifest_path, name):
    manifest_path=Path(manifest_path);before=digest(manifest_path)
    manifest=json.loads(manifest_path.read_text())
    if manifest.get('status')!='conditional_empirical_banks_prepared' or manifest.get('physical_acceptance')!='incomplete':
        raise ValueError('Conditional historical bank manifest required')
    if name not in ('mc3d-training','mc3d-heldout','sparta-training','sparta-heldout'):
        raise ValueError('Unknown bank name')
    entry=manifest['banks'][name];filename=entry['file']
    if Path(filename).name!=filename or Path(filename).suffix!='.npz':
        raise ValueError('Local NPZ bank file required')
    path=manifest_path.parent/filename
    if digest(path)!=entry['sha256']:raise ValueError('Bank hash mismatch')
    with np.load(path,allow_pickle=False) as data:
        if set(data.files)!={'x','mean_profile','fluctuations'}:raise ValueError('Unexpected bank arrays')
        x=np.array(data['x'],dtype=float);mean=np.array(data['mean_profile'],dtype=float);noise=np.array(data['fluctuations'],dtype=float)
    count=24 if name.endswith('-training') else 12
    if x.ndim!=1 or len(x)<4 or not np.isfinite(x).all() or np.any(np.diff(x)<=0):
        raise ValueError('Finite increasing physical grid required')
    if mean.shape!=x.shape or noise.shape!=(count,64,len(x)) or not np.isfinite(mean).all() or not np.isfinite(noise).all():
        raise ValueError('Invalid finite bank dimensions')
    if entry['shape']!=list(noise.shape) or len(entry['sources'])!=count:
        raise ValueError('Bank metadata differs from arrays')
    split='training' if count==24 else 'heldout';groups=manifest['protocol'][split+'_groups'];sources=entry['sources']
    if len(groups)!=(8 if count==24 else 4) or len(set(groups))!=len(groups):raise ValueError('Invalid bank group split')
    if (any(type(s['seed'])!=int or s['seed']<=0 for s in sources) or len({s['seed'] for s in sources})!=count or
        set(s['group'] for s in sources)!=set(groups) or any(sum(s['group']==g for s in sources)!=3 for g in groups)):
        raise ValueError('Require distinct seeds and three sources per group')
    if not np.allclose(noise.mean((0,1)),0.,rtol=0,atol=1e-12):
        raise ValueError('Bank is not ensemble centered')
    if not np.allclose(mean,entry['mean_profile'],rtol=1e-12,atol=1e-14):raise ValueError('Bank mean metadata differs')
    if digest(manifest_path)!=before or digest(path)!=entry['sha256']:raise ValueError('Bank changed while loading')
    for array in (x,mean,noise):array.setflags(write=False)
    return dict(name=name,x=x,mean_profile=mean,fluctuations=noise,sources=sources,
                manifest_sha256=before,bank_sha256=entry['sha256'])


def fixture(kind,seed,bank,windows=64,motion_amplitude=.05):
    if kind not in ('logistic','gompertz') or windows not in (32,64) or type(seed)!=int or seed<=0:
        raise ValueError('Supported shape, seed and 32/64 windows required')
    if not math.isfinite(motion_amplitude) or motion_amplitude<0:
        raise ValueError('Finite nonnegative front motion required')
    rng=np.random.default_rng(seed);x=bank['x'];noise=bank['fluctuations'];trials=[];choices=[]
    for _ in range(3):
        source=int(rng.integers(len(noise)));start=int(rng.integers(64));indices=(start+np.arange(windows))%64
        shifts=rng.uniform(-motion_amplitude,motion_amplitude,size=windows);run=[]
        for index,shift in zip(indices,shifts):
            y=shape(x-shift,kind)+noise[source,index]
            run.append(dict(x=x.tolist(),density_normalized=y.tolist()))
        trials.append(run);choices.append(dict(source_index=source,source=bank['sources'][source],start_window=start,
            window_indices=indices.tolist(),true_shifts=shifts.tolist()))
    return trials,dict(kind=kind,seed=seed,windows=windows,motion_amplitude=motion_amplitude,
        bank=bank['name'],manifest_sha256=bank['manifest_sha256'],bank_sha256=bank['bank_sha256'],choices=choices,
        limitation='Conditional empirical stress only; cyclic ordering has a wrap boundary and per-run mean fluctuations are retained.')
