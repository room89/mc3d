"""Training-only downstream/front and global particle-balance diagnostic."""
import argparse
import csv
import io
import json
from pathlib import Path
import numpy as np
from fit_physical_plateau_mean import design, coefficients
from shock_empirical_fluctuations import load_bank
from shock_analysis import crossing, monotone
from shock_validation import digest


def correlation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape or a.ndim != 1 or not np.isfinite([a, b]).all():
        raise ValueError('Matching finite series required')
    a = a-a.mean(); b = b-b.mean()
    denominator = np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/denominator) if denominator > 0 else None


def associations(y, z, sources):
    groups = {}
    for group in sorted({s['group'] for s in sources}):
        mask = [s['group'] == group for s in sources]
        a, b = y[mask].mean(axis=0), z[mask].mean(axis=0)
        groups[group] = dict(correlation=correlation(a, b),
            target_half_change=float(a[32:].mean()-a[:32].mean()),
            predictor_half_change=float(b[32:].mean()-b[:32].mean()))
    return dict(group_curves=groups, ensemble_correlation=correlation(y.mean(axis=0), z.mean(axis=0)),
        per_run_correlations=[correlation(a, b) for a, b in zip(y, z)])


def particle_balance(text, bounds):
    """Integrate pre-step population; assign entire net flux by step midpoint."""
    rows = list(csv.DictReader(io.StringIO(text), delimiter=';'))
    if not rows:
        raise ValueError('Empty conservation history')
    a = np.array([[float(r[k]) for k in ('time', 'end_time', 'N_before', 'N_after', 'N_in', 'N_out', 'Rcount')] for r in rows])
    bounds = np.asarray(bounds, float)
    if bounds.shape != (64, 2) or not np.isfinite(bounds).all() or np.any(bounds[:,1] <= bounds[:,0]) or not np.isfinite(a).all():
        raise ValueError('Invalid history or windows')
    if not np.allclose(bounds[1:,0], bounds[:-1,1], rtol=1e-12, atol=1e-14):
        raise ValueError('Noncontiguous windows')
    if np.any(a[:,1] <= a[:,0]) or not np.allclose(a[1:,0], a[:-1,1], rtol=1e-12, atol=1e-14):
        raise ValueError('Noncontiguous steps')
    if np.any(a[:,2:] != np.floor(a[:,2:])) or np.any(a[:,2:6] < 0):
        raise ValueError('Invalid particle counts')
    if np.any(a[:,3]-a[:,2] != a[:,4]-a[:,5]) or np.any(a[:,6] != 0) or np.any(a[1:,2] != a[:-1,3]):
        raise ValueError('Particle balance failure')
    if a[0,0] > bounds[0,0] or a[-1,1] < bounds[-1,1]-1e-12:
        raise ValueError('Incomplete time coverage')
    population, flux = [], []
    midpoint = a[:,:2].mean(axis=1)
    for lo, hi in bounds:
        weights = np.maximum(0, np.minimum(a[:,1], hi)-np.maximum(a[:,0], lo))
        if not np.isclose(weights.sum(), hi-lo, rtol=1e-10, atol=1e-14):
            raise ValueError('Incomplete window coverage')
        population.append(float(weights@a[:,2]/weights.sum()))
        mask = (midpoint >= lo) & (midpoint < hi)
        flux.append(float(np.sum(a[mask,4]-a[mask,5])))
    return dict(population=population, net_particles_by_step_midpoint=flux,
        checked_steps=len(rows), max_step_over_window=float(np.max(a[:,1]-a[:,0])/np.min(np.diff(bounds,axis=1))))


def diagnose(root):
    root = Path(root).resolve()
    manifest = root/'output/shock-physical-fluctuation-banks/manifest.json'
    mh = digest(manifest); bm = json.loads(manifest.read_text())
    cp = root/'output/shock-independent-c/campaign.json'
    expected = bm['protocol']['previous_protocol']['campaign_sha256']
    evidence = {}; result = {}; bank_hashes = {}
    def read(path, sha):
        if digest(path) != sha:
            raise ValueError('Historical evidence changed')
        evidence[str(path.resolve().relative_to(root))] = sha
        return path.read_text()
    campaign = json.loads(read(cp, expected))
    for code in ('mc3d', 'sparta'):
        bank = load_bank(manifest, code+'-training'); bank_hashes[bank['name']] = bank['bank_sha256']
        raw = bank['fluctuations']+bank['mean_profile']; x = bank['x']; sources = bank['sources']
        if len(sources) != 24 or {s['group'] for s in sources} != {f'g{i:02d}' for i in range(8)}:
            raise ValueError('Historical training split required')
        center = crossing(x, monotone(raw.mean(axis=(0,1))), .5)
        mask = (x >= center+8) & (x <= center+12)
        basis = design(x[mask]); prediction = basis@coefficients(raw[:,:,mask], basis)
        y = (raw[:,:,mask]-prediction).mean(axis=2)
        front = np.array([[crossing(x, monotone(w), .5) for w in run] for run in raw])
        entry = dict(sources=sources, residual_series=y.tolist(), front_lambda=front.tolist(),
            front_association=associations(y, front, sources),
            downstream_bounds_lambda=[float(x[mask][0]),float(x[mask][-1])])
        if not np.isfinite(front).all():
            raise ValueError('Unresolved front')
        if code == 'mc3d':
            balances = []
            for source in sources:
                group, seed = source['group'], source['seed']
                checkpoint = root/'output/shock-physical-noise-diagnostics'/f'{group}.json'
                hashes = json.loads(read(checkpoint,bm['previous_checkpoint_sha256'][group]))['result']['input_sha256']
                folder = (cp.parent/campaign['groups'][group]['path']/f'C-{seed}').resolve()
                mp = folder/'output/manifest.json'
                run = json.loads(read(mp,hashes[str(mp.relative_to(root))]))
                if run['status'] != 'completed':
                    raise ValueError('Incomplete run')
                outputs = {v['path']:v['sha256'] for v in run['outputs']}
                balances.append(particle_balance(read(folder/'output/conservation.csv',outputs['conservation.csv']),run['effective_config']['statistics_windows']))
            population = np.array([b['population'] for b in balances]); flux = np.array([b['net_particles_by_step_midpoint'] for b in balances])
            entry.update(particle_balances=balances, population_association=associations(y,population,sources),
                net_flux_association=associations(y,flux,sources), front_population_association=associations(front,population,sources))
        else:
            entry['balance_status'] = 'not_audited_no_equivalent_particle_flux_series_in_this_diagnostic'
        result[code] = entry
    if digest(manifest) != mh or any(digest(root/n) != sha for n,sha in evidence.items()):
        raise ValueError('Evidence changed during analysis')
    if any(digest(manifest.parent/bm['banks'][n]['file']) != sha for n,sha in bank_hashes.items()):
        raise ValueError('Bank changed')
    return dict(status='exploratory_training_front_particle_balance', physical_acceptance='incomplete',
        manifest_sha256=mh, bank_sha256=bank_hashes, input_sha256=evidence, codes=result,
        sources={n:digest(Path(__file__).with_name(n)) for n in ('diagnose_plateau_front_balance.py','fit_physical_plateau_mean.py','shock_analysis.py','shock_empirical_fluctuations.py')},
        limitations=['Training runs only; no held-out series read and no new predictor fitted.',
            'Front is the isotonic 50 percent crossing on the bank fixed grid; it shares density noise with the plateau.',
            'Correlations are contemporaneous descriptive associations with autocorrelated series; no p-values or causality claim.',
            'Population is duration-weighted pre-step count; whole-step net particle flux is assigned by step midpoint, not exact fractional flux.',
            'Global particle conservation is audited for MC3D only; energy and momentum are outside this diagnostic.',
            'Full-training mean subtraction is in-sample; no stationary noise or physical acceptance qualification.'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    report = diagnose(args.root)
    with args.output.open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
