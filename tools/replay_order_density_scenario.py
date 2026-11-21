"""Replay all 128 empirical trials with profiles, without adding coverage trials."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
import calibrate_empirical_order as driver
from diagnose_order_density_coverage import diagnose
from shock_alignment import GRID
from shock_validation import digest
from validate_shock_estimator import shape


def save_new(path, data):
    temporary = path.with_suffix('.tmp')
    with temporary.open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')
    if path.exists():
        raise FileExistsError(path)
    temporary.rename(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('audit', 'parts', 'manifest', 'verified', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--bank', default='mc3d-heldout')
    parser.add_argument('--shape', default='logistic')
    parser.add_argument('--motion', type=float, default=0.)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    evidence = diagnose(args.audit, args.parts, args.bank, args.shape, args.motion)
    protocol = evidence['protocol']
    source_dir = Path(__file__).parent
    def verify_inputs():
        if any(digest(source_dir/name) != value for name, value in protocol['sources'].items()):
            raise ValueError('Frozen source changed')
        if driver.verify_banks(args.manifest, args.verified) != protocol['bank_evidence']:
            raise ValueError('Bank evidence changed')
    verify_inputs()
    args.output.mkdir(parents=True, exist_ok=args.resume)
    plan = dict(status='diagnostic_replay_only', physical_acceptance='incomplete',
                scenario=evidence['scenario'], protocol=protocol, trials=128,
                script_sha256=digest(Path(__file__)),
                checkpoints={r['checkpoint']:r['sha256'] for r in evidence['trials']},
                limitations=['Replays are not new independent coverage observations.',
                             'All 128 original trials retained; no correction fitted.'])
    plan_path = args.output/'plan.json'
    if plan_path.exists():
        if not args.resume or json.loads(plan_path.read_text()) != plan:
            raise ValueError('Replay plan changed')
    else:
        save_new(plan_path, plan)
    original = driver.estimate
    captured = {}
    def capture(groups, **kwargs):
        result = original(groups, **kwargs)
        captured.update(groups=groups, estimate=result)
        return result
    driver.estimate = capture
    center = -math.log(math.log(2)) if args.shape == 'gompertz' else 0.
    exact = shape(GRID+center, args.shape)
    completed = {}
    for row in evidence['trials']:
        checkpoint = args.parts/row['checkpoint']
        if digest(checkpoint) != row['sha256']:
            raise ValueError('Checkpoint changed before replay')
        saved = json.loads(checkpoint.read_text())
        target = args.output/row['checkpoint']
        if target.exists():
            report = json.loads(target.read_text())
            if (report['checkpoint_sha256'] != row['sha256'] or
                report['plan_sha256'] != digest(plan_path) or
                report['scalar_result'] != saved['result'] or
                report['spec'] != saved['spec']):
                raise ValueError('Existing replay identity mismatch')
        else:
            captured.clear()
            scalar = driver.trial(saved['spec'], protocol, args.manifest)
            if scalar != saved['result']:
                raise ValueError('Replay differs from original scalar result')
            diagnostics = {}
            for branch in ('primary_estimate', 'alternative_estimate'):
                estimate = captured['estimate'][branch]
                error = np.asarray(estimate['density'])-exact
                index = int(np.argmax(abs(error)))
                radius = estimate['density_statistical_radius']+estimate['density_processing_sensitivity']
                diagnostics[branch] = dict(signed_error=error.tolist(),
                    x_at_maximum=float(GRID[index]), signed_error_at_maximum=float(error[index]),
                    maximum_absolute_error=float(abs(error[index])), own_radius=radius,
                    own_band_covers=bool(abs(error[index]) <= radius))
            verify_inputs()
            if digest(checkpoint) != row['sha256']:
                raise ValueError('Checkpoint changed during replay')
            report = dict(status='diagnostic_replay_only', checkpoint_sha256=row['sha256'],
                          plan_sha256=digest(plan_path), spec=saved['spec'], scalar_result=scalar,
                          exact_scalar_replay=True, grid=GRID.tolist(), truth=exact.tolist(),
                          diagnostics=diagnostics, **captured)
            save_new(target, report)
        completed[target.name] = digest(target)
        print('Replayed', len(completed), 'of 128; exact original scalars', flush=True)
    verify_inputs()
    if any(digest(args.parts/name) != value for name, value in plan['checkpoints'].items()):
        raise ValueError('Original evidence changed')
    if any(digest(args.output/name) != value for name, value in completed.items()):
        raise ValueError('Replay evidence changed')
    save_new(args.output/'complete.json', dict(status='all_128_replays_verified',
             physical_acceptance='incomplete', plan_sha256=digest(plan_path), replays=completed))


if __name__ == '__main__':
    main()
