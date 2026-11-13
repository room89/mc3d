"""Read-only physical-time and group diagnostic; no new mean model is fitted."""
import argparse
import json
from pathlib import Path
import numpy as np
from fit_physical_plateau_mean import diagnose as mean_diagnose
from shock_empirical_fluctuations import load_bank
from shock_validation import digest


def validate_bounds(bounds):
    a = np.asarray(bounds, float)
    if a.shape != (64, 2) or not np.isfinite(a).all() or np.any(a[:, 1] <= a[:, 0]):
        raise ValueError('64 finite positive-duration windows required')
    if not np.allclose(a[1:, 0], a[:-1, 1], rtol=1e-12, atol=1e-15):
        raise ValueError('Contiguous ordered windows required')
    return a


def summarize(series, bounds, sources):
    y = np.asarray(series, float)
    b = validate_bounds(bounds)
    if y.shape != (len(sources), 64) or not len(y) or not np.isfinite(y).all():
        raise ValueError('Finite run series with matching identities required')
    t = b.mean(axis=1)
    groups = {}
    for group in sorted({s['group'] for s in sources}):
        v = y[[s['group'] == group for s in sources]]
        groups[group] = dict(seeds=[s['seed'] for s in sources if s['group'] == group],
            mean=float(v.mean()), window_means=v.mean(axis=0).tolist(),
            eight_window_means=v.mean(axis=0).reshape(8, 8).mean(axis=1).tolist(),
            half_change=float(v[:, 32:].mean()-v[:, :32].mean()))
    changes = np.array([g['half_change'] for g in groups.values()])
    return dict(groups=groups, eight_window_means=y.mean(axis=0).reshape(8, 8).mean(axis=1).tolist(),
        block_midpoints_s=t.reshape(8, 8).mean(axis=1).tolist(),
        half_change=float(y[:, 32:].mean()-y[:, :32].mean()),
        positive_group_half_changes=int(np.sum(changes > 0)), group_count=len(groups),
        group_half_change_range=[float(changes.min()), float(changes.max())],
        leave_one_group_out_half_changes={g: float(np.delete(changes, i).mean())
            for i, g in enumerate(groups)} if len(groups) > 1 else {},
        normalized_time_vs_index_max_error=float(np.max(abs(2*(t-t[0])/(t[-1]-t[0])-1-np.linspace(-1, 1, 64)))))


def diagnose(root):
    root = Path(root).resolve()
    manifest = root/'output/shock-physical-fluctuation-banks/manifest.json'
    bank_manifest = json.loads(manifest.read_text())
    campaign_path = root/'output/shock-independent-c/campaign.json'
    campaign = json.loads(campaign_path.read_text())
    expected = bank_manifest['protocol']['previous_protocol']['campaign_sha256']
    if digest(campaign_path) != expected:
        raise ValueError('Campaign changed')
    evidence = {str(campaign_path.relative_to(root)): expected}
    mean = mean_diagnose(manifest)
    result = {}
    def read(path, expected_hash):
        relative = str(path.resolve().relative_to(root))
        if digest(path) != expected_hash:
            raise ValueError('Historical evidence changed: '+relative)
        evidence[relative] = expected_hash
        return path.read_text()
    for code in ('mc3d', 'sparta'):
        result[code] = {}
        for split in ('training', 'heldout'):
            bank = load_bank(manifest, code+'-'+split)
            all_bounds = []
            for source in bank['sources']:
                group, seed = source['group'], source['seed']
                checkpoint = root/'output/shock-physical-noise-diagnostics'/f'{group}.json'
                checkpoint_hash = bank_manifest['previous_checkpoint_sha256'][group]
                hashes = json.loads(read(checkpoint, checkpoint_hash))['result']['input_sha256']
                folder = (campaign_path.parent/campaign['groups'][group]['path']).resolve()
                def historical(relative):
                    path = folder/relative
                    return read(path, hashes[str(path.relative_to(root))])
                if code == 'mc3d':
                    cfg = json.loads(historical(f'C-{seed}/solver.json'))
                    units = json.loads(historical(f'C-{seed}/units.json'))
                    run = json.loads(historical(f'C-{seed}/output/manifest.json'))
                    if run['status'] != 'completed' or cfg['statistics_windows'] != run['effective_config']['statistics_windows']:
                        raise ValueError('Incomplete or mismatched MC3D schedule')
                    bounds = np.asarray(cfg['statistics_windows'])*units['scales']['time_s']
                else:
                    ref = json.loads(historical(f'sparta-C-{seed}/reference.json'))
                    step, dt = ref['steps_per_window'], ref['timestep']
                    if ref['quick'] or step <= 0 or dt <= 0:
                        raise ValueError('Production SPARTA schedule required')
                    deck = historical(f'sparta-C-{seed}/in.shock')
                    required = [f'fix avg ave/grid all 1 {step} {step} c_moments[*] c_thermal[*]', f'run {32*step}', f'run {64*step}']
                    timesteps = [float(line.split()[1]) for line in deck.splitlines() if line.startswith('timestep ')]
                    if timesteps != [dt] or any(line not in deck.splitlines() for line in required):
                        raise ValueError('Unexpected SPARTA averaging schedule')
                    bounds = np.array([[(32+i)*step*dt, (33+i)*step*dt] for i in range(64)])
                    for i in range(64):
                        text = historical(f'sparta-C-{seed}/profile.{(33+i)*step}')
                        lines = text.splitlines()
                        if lines[0] != 'ITEM: TIMESTEP' or int(lines[1]) != (33+i)*step:
                            raise ValueError('SPARTA dump time mismatch')
                all_bounds.append(validate_bounds(bounds))
            if not all(np.allclose(b, all_bounds[0], rtol=1e-12, atol=1e-15) for b in all_bounds):
                raise ValueError('Run schedules differ; pooling is not supported')
            bounds = all_bounds[0]
            region = mean['codes'][code]['regions']['downstream']
            lo, hi = region['bounds']
            mask = (bank['x'] >= lo) & (bank['x'] <= hi)
            raw = (bank['fluctuations']+bank['mean_profile'])[:, :, mask]
            residual = raw-np.asarray(region['prediction'])
            result[code][split] = dict(bounds_s=bounds.tolist(),
                raw=summarize(raw.mean(axis=2), bounds, bank['sources']),
                residual=summarize(residual.mean(axis=2), bounds, bank['sources']))
    if any(digest(root/name) != sha for name, sha in evidence.items()):
        raise ValueError('Evidence changed while reading')
    if digest(manifest) != mean['manifest_sha256']:
        raise ValueError('Bank manifest changed')
    for name, bank in mean['banks'].items():
        if digest(manifest.parent/bank_manifest['banks'][name]['file']) != bank['sha256']:
            raise ValueError('Bank changed')
    return dict(status='exploratory_time_group_diagnostic', physical_acceptance='incomplete',
        manifest_sha256=mean['manifest_sha256'], banks=mean['banks'], input_sha256=evidence,
        source_sha256=digest(Path(__file__)), mean_sources=mean['sources'], codes=result,
        limitations=['Historical inspected data, not fresh validation; no model selection or fitting to held-out runs.',
            'MC3D bounds are configured integration windows converted to SI; SPARTA bounds are averaging cycles ending at verified dump steps.',
            'SPARTA samples discrete step endpoints; cycle midpoints are labels, not continuous-time measurement equivalence.',
            'Residuals use the full training mean model, so training residuals are in-sample.',
            'Group trends and eight-window averages are descriptive; stationarity and noise coverage remain unqualified.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = diagnose(args.root)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
