"""Historical training-only precision diagnostic; never a physical acceptance gate."""
import argparse
import json
from pathlib import Path
import numpy as np
from shock_empirical_fluctuations import load_bank
from screen_simultaneous_band import critical_value
from shock_validation import digest
from shock_analysis import crossing, monotone


def grouped(raw, sources):
    a = np.asarray(raw, float)
    if a.ndim != 3 or len(sources) != len(a) or not np.isfinite(a).all():
        raise ValueError('Finite runs/windows/nodes and matching sources required')
    groups = sorted({s['group'] for s in sources})
    if len(groups) != 8 or any(sum(s['group'] == g for s in sources) != 3 for g in groups):
        raise ValueError('Eight complete three-run groups required')
    return np.array([a[[s['group'] == g for s in sources]].mean(axis=(0, 1)) for g in groups])


def band(group_profiles, critical):
    return group_profiles.mean(axis=0), critical * group_profiles.std(axis=0, ddof=1) / np.sqrt(8)



def temporal_diagnostic(raw, sources, x):
    """Descriptive half-profile changes; no stationarity test or fitted correction."""
    if raw.shape[1] != 64:
        raise ValueError('Historical 64-window series required')
    first = grouped(raw[:, :32], sources)
    second = grouped(raw[:, 32:], sources)
    front = lambda y: float(crossing(x, monotone(y), .5))
    mean = grouped(raw, sources)
    return dict(mean_front=front(mean.mean(axis=0)),
        group_fronts=[front(y) for y in mean],
        first_half_front=front(first.mean(axis=0)),
        second_half_front=front(second.mean(axis=0)),
        half_profile_change=(second.mean(axis=0)-first.mean(axis=0)).tolist(),
        group_half_profile_changes=(second-first).tolist(),
        max_abs_half_profile_change=float(abs(second.mean(axis=0)-first.mean(axis=0)).max()))


def audit(root):
    manifest = Path(root)/'output/shock-physical-fluctuation-banks/manifest.json'
    before = digest(manifest)
    banks = {c: load_bank(manifest, c+'-training') for c in ('mc3d', 'sparta')}
    x = banks['mc3d']['x']
    if not np.array_equal(x, banks['sparta']['x']):
        raise ValueError('Identical fixed physical grids required; no interpolation')
    # Two code bands jointly at nominal 95%; sum their widths for a difference.
    # This does not assume independence between codes or pair their seeds.
    critical = critical_value(len(x), .025)
    results = {}
    for code, b in banks.items():
        raw = b['fluctuations'] + b['mean_profile']
        mean, half = band(grouped(raw, b['sources']), critical)
        results[code] = dict(mean=mean.tolist(), half_width=half.tolist(),
            group_profiles=grouped(raw,b['sources']).tolist(), sources=b['sources'],
            temporal=temporal_diagnostic(raw,b['sources'],x))
    difference = np.array(results['mc3d']['mean'])-results['sparta']['mean']
    uncertainty = np.array(results['mc3d']['half_width'])+results['sparta']['half_width']
    if digest(manifest) != before:
        raise ValueError('Manifest changed')
    for code, b in banks.items():
        if load_bank(manifest, code+'-training')['bank_sha256'] != b['bank_sha256']:
            raise ValueError('Bank changed')
    return dict(status='historical_precision_diagnostic',physical_acceptance='incomplete',
        manifest_sha256=before,bank_sha256={c:b['bank_sha256'] for c,b in banks.items()},
        source_sha256=digest(Path(__file__)),nodes=len(x),x=x.tolist(),critical=critical,
        nominal_joint_alpha=.05,codes=results,difference=difference.tolist(),
        difference_half_width=uncertainty.tolist(),
        max_abs_difference=float(abs(difference).max()),max_half_width=float(uncertainty.max()),
        max_abs_difference_plus_half_width=float((abs(difference)+uncertainty).max()),
        precision_target=.005,nodes_exceeding_precision=int(np.sum(uncertainty>.005)),
        limitations=['Historical training data, not fresh qualification; no held-out data read.',
            '560 source nodes instead of 113 synthetic nodes; family error reallocated to both codes and all nodes.',
            'Nominal Student bands require approximately Gaussian independent group means; physical coverage is unqualified.',
            'Time-averaged laboratory-frame density, not intrinsic shock shape, thickness or proof of stationarity.',
            'No fitted mean, alignment, rescaling, per-run centering, pairing of solvers or physical parameter changes.',
            'No extrapolated campaign size: group variance and physical laws need prospective verification.'])


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=audit(a.root)
    with a.output.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:r[k] for k in ('nodes','max_abs_difference','max_half_width','max_abs_difference_plus_half_width','nodes_exceeding_precision')}))
