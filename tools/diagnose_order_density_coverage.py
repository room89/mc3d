"""Read-only, all-trial density diagnostics for a completed empirical scenario."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def diagnose(audit_path, parts, bank, shape, motion):
    raw = audit_path.read_bytes()
    audit = json.loads(raw)
    series = parts.name.removesuffix('.parts')
    hashes = audit['empirical'][series]['checkpoint_sha256']
    rows = []
    protocol = None
    for name, expected in hashes.items():
        path = parts / name
        content = path.read_bytes()
        if sha(content) != expected:
            raise ValueError(f'Checkpoint changed: {name}')
        data = json.loads(content)
        if data['spec'][:-1] != [bank, shape, motion]:
            continue
        if protocol is None:
            protocol = data['protocol']
        if data['protocol'] != protocol:
            raise ValueError('Mixed protocols')
        r = data['result']
        if 'error' in r:
            raise ValueError('Rejected trial must not be omitted from diagnosis')
        error, radius, statistical, processing, difference = (
            r[k] for k in ('density_error', 'density_radius',
                          'density_statistical_radius', 'density_processing_sensitivity',
                          'order_density_difference'))
        if not all(math.isfinite(v) and v >= 0 for v in
                   (error, radius, statistical, processing, difference)) or radius <= 0:
            raise ValueError('Invalid density diagnostics')
        if not math.isclose(radius, statistical + processing, rel_tol=1e-12):
            raise ValueError('Radius decomposition mismatch')
        if r['density_covers'] != (error <= radius) or processing < difference:
            raise ValueError('Coverage or envelope mismatch')
        if r['repetition'] != data['spec'][-1]:
            raise ValueError('Trial identity mismatch')
        rows.append(dict(repetition=r['repetition'], checkpoint=name, sha256=expected,
                         error=error, radius=radius, statistical=statistical,
                         processing=processing, order_difference=difference,
                         error_radius_ratio=error/radius,
                         uncovered_excess=max(0., error-radius), covers=error <= radius,
                         alternative_error_lower_bound=max(0., error-difference),
                         alternative_error_upper_bound=error+difference))
    rows.sort(key=lambda r: r['repetition'])
    if len(rows) != 128 or [r['repetition'] for r in rows] != list(range(128)):
        raise ValueError('Exactly the complete fixed 128-trial scenario is required')
    summary = {}
    for label, subset in (('all', rows), ('covered', [r for r in rows if r['covers']]),
                          ('missed', [r for r in rows if not r['covers']])):
        summary[label] = dict(count=len(subset), fields={
            key: dict(min=min(values), median=statistics.median(values),
                      mean=statistics.mean(values), max=max(values))
            for key in ('error', 'radius', 'statistical', 'processing', 'order_difference',
                        'error_radius_ratio', 'uncovered_excess')
            if (values := [r[key] for r in subset])})
    for r in rows:
        if sha((parts/r['checkpoint']).read_bytes()) != r['sha256']:
            raise ValueError('Checkpoint changed during diagnosis')
    if audit_path.read_bytes() != raw:
        raise ValueError('Audit changed during diagnosis')
    return dict(status='diagnostic_not_qualification', physical_acceptance='incomplete',
                scenario=[bank, shape, motion], windows=protocol['windows'],
                audit_sha256=sha(raw), script_sha256=sha(Path(__file__).read_bytes()),
                protocol=protocol, summary=summary, trials=rows, limitations=[
                    'Conditional on a finite empirical bank; all 128 trials retained.',
                    'Covered/missed summaries are descriptive partitions, not independent tests.',
                    'Checkpoints retain supremum errors, not signed profiles or error locations.',
                    'Alternative error bounds follow the triangle inequality and are not confidence bounds.',
                    'No correction, radius inflation, tolerance change, or qualification fitted here.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--parts', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bank', default='mc3d-heldout')
    parser.add_argument('--shape', default='logistic')
    parser.add_argument('--motion', type=float, default=0.)
    args = parser.parse_args()
    result = diagnose(args.audit, args.parts, args.bank, args.shape, args.motion)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    main()
