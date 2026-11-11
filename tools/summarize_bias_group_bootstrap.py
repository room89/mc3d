#!/usr/bin/env python3
"""Summarize a completed exploratory source-group bootstrap, never acceptance."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from shock_validation import digest


def curve_statistics(curves, baseline):
    curves = np.asarray(curves, dtype=float)
    baseline = np.asarray(baseline, dtype=float)
    if curves.shape != (64, 401) or baseline.shape != (401,):
        raise ValueError('Require 64 complete curves on 401 nodes')
    if not np.isfinite(curves).all() or not np.isfinite(baseline).all():
        raise ValueError('Nonfinite bootstrap curves')
    deviations = curves - baseline
    norms = np.max(np.abs(deviations), axis=1)
    return dict(
        maximum_absolute_bootstrap_mean_shift=float(np.max(np.abs(deviations.mean(0)))),
        maximum_pointwise_standard_deviation=float(np.max(curves.std(0, ddof=1))),
        empirical_95_percent_supnorm=float(np.quantile(norms, .95, method='higher')),
        maximum_supnorm=float(norms.max()),
    )


def summarize(directory, prior_path):
    directory, prior_path = Path(directory), Path(prior_path)
    complete_path, plan_path = directory/'complete.json', directory/'plan.json'
    if not complete_path.exists():
        raise ValueError('Incomplete diagnostic: complete.json is required')
    plan_hash, complete_hash, prior_hash = map(digest, (plan_path, complete_path, prior_path))
    plan, complete, prior = [json.loads(p.read_text()) for p in (plan_path, complete_path, prior_path)]
    if (plan['scheduled_replicates'], plan['replicates_per_case'], complete['completed_replicates']) != (256, 64, 256):
        raise ValueError('Require the complete fixed 256-replicate protocol')
    if complete['plan_sha256'] != plan_hash or plan['prior_evidence_sha256'] != prior_hash:
        raise ValueError('Evidence hash mismatch')
    if (plan['seed_start'], plan['shape'], plan['added_motion'], plan['status']) != (190000000, 'logistic', 0., 'exploratory_source_group_bootstrap'):
        raise ValueError('Unexpected diagnostic protocol')
    if any(prior['bank_hashes'].get(name) != sha for name, sha in plan['bank_hashes'].items()):
        raise ValueError('Prior bank evidence mismatch')
    if complete['failed_replicates'] != 0:
        raise ValueError('Failed replicates present; no successful-only summary allowed')
    if plan['banks'] != ['mc3d-training', 'mc3d-heldout'] or plan['modes'] != ['windowwise', 'mean32']:
        raise ValueError('Unexpected scenario set')
    expected = {f'{name}-{mode}-{i}.json' for name, mode in itertools.product(plan['banks'], plan['modes']) for i in range(64)}
    if set(complete['checkpoint_sha256']) != expected:
        raise ValueError('Missing or extra checkpoints in completion record')
    arrays, baselines, statistics = {}, {}, {}
    for case_index, (name, mode) in enumerate(itertools.product(plan['banks'], plan['modes'])):
        key = f'{name}-{mode}'
        rows = []
        for i in range(64):
            filename = f'{key}-{i}.json'; path = directory/filename
            if digest(path) != complete['checkpoint_sha256'][filename]:
                raise ValueError('Checkpoint hash mismatch')
            row = json.loads(path.read_text())
            if (row['bank'], row['mode'], row['replicate'], row['seed'], row['plan_sha256'], row['status']) != (name, mode, i, plan['seed_start']+case_index*10000+i, plan_hash, 'ok'):
                raise ValueError('Checkpoint identity or status mismatch')
            rows.append(row['curve'])
        arrays[key] = np.asarray(rows)
        baselines[key] = np.asarray(prior['cases'][f'{name}-logistic-{mode}-nonlinear']['error'])
        statistics[key] = curve_statistics(arrays[key], baselines[key])
    comparisons = {}
    for mode in plan['modes']:
        train, held = f'mc3d-training-{mode}', f'mc3d-heldout-{mode}'
        observed = baselines[held]-baselines[train]
        deviations = (arrays[held]-baselines[held])-(arrays[train]-baselines[train])
        stats = curve_statistics(deviations, np.zeros(401))
        stats['observed_training_heldout_supnorm'] = float(np.max(np.abs(observed)))
        stats['observed_over_empirical_95_percent_supnorm'] = stats['observed_training_heldout_supnorm']/stats['empirical_95_percent_supnorm'] if stats['empirical_95_percent_supnorm'] > 0 else None
        comparisons[mode] = stats
    for filename, sha in complete['checkpoint_sha256'].items():
        if digest(directory/filename) != sha:
            raise ValueError('Checkpoint changed during summary')
    if tuple(map(digest, (plan_path, complete_path, prior_path))) != (plan_hash, complete_hash, prior_hash):
        raise ValueError('Evidence changed during summary')
    return dict(status='complete_exploratory_bootstrap_summary', physical_acceptance='incomplete',
                completed_replicates=256, failed_replicates=0, cases=statistics, comparisons=comparisons,
                plan_sha256=plan_hash, complete_sha256=complete_hash, prior_sha256=prior_hash,
                limitations=[
                    'Empirical 95th percentiles use the higher order statistic of 64 resamples; not validated confidence bands.',
                    'Training and heldout resamples use independent seeds; differences are paired by replicate index.',
                    'Observed discrepancy is descriptive; no p-value, coverage claim, or tuned acceptance threshold.',
                    'Only 8 training and 4 heldout original groups; finite-bank and exact-template limitations remain.',
                    'No qualification on new physical ensembles or change to the production estimator.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = summarize(args.directory, args.prior)
        with args.output.open('x') as handle:
            json.dump(result, handle, indent=2); handle.write('\n')
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, f'{error}\n')


if __name__ == '__main__':
    main()
