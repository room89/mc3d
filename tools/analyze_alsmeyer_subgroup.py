"""Build subgroup statistics from verified reference profiles; no acceptance verdict."""
import argparse
import json
from pathlib import Path
import numpy as np
from extract_alsmeyer_profiles import SOURCES as EXTRACTION_SOURCES, canonical_windows
from reanalyze_shock_campaign import merged_pairs, group_edges
from shock_ensemble_estimator import group_point
from shock_high_order_estimator import METHOD
from shock_validation import digest, save

KEYS = tuple('mc3d-'+v for v in 'ABCDEFG') + tuple('sparta-'+v for v in 'CDE')
SOURCES = tuple(sorted(set(EXTRACTION_SOURCES) | {
    'analyze_alsmeyer_subgroup.py', 'reanalyze_shock_campaign.py',
    'shock_ensemble_estimator.py', 'shock_high_order_estimator.py',
    'shock_alignment.py', 'shock_crossfit_gradient.py', 'shock_mean_profile.py'}))


def prepare_runs(runs, seeds):
    if len(seeds)!=3 or len(set(seeds))!=3 or any(type(s)!=int or not 0<s<2**31 for s in seeds):
        raise ValueError('Three distinct positive seed identities required')
    if len(runs)!=3 or sorted(r['seed'] for r in runs)!=sorted(seeds):
        raise ValueError('Subgroup seed identities differ')
    trials=[]; shared=None
    for run in sorted(runs,key=lambda r:r['seed']):
        times=np.asarray(run['statistics_windows_s'],float)
        scale=run['collision_mean_free_path_m']
        if times.shape!=(64,2) or not np.isfinite(times).all() or not np.isfinite(scale) or scale<=0:
            raise ValueError('Finite physical windows and collision length required')
        duration=times[:,1]-times[:,0]
        if (times[0,0]<0 or np.any(duration<=0) or
                not np.allclose(duration,duration[0],rtol=1e-10,atol=0) or
                not np.allclose(times[1:,0],times[:-1,1],rtol=1e-12,atol=0)):
            raise ValueError('Contiguous equal-duration windows required for pair averaging')
        windows=run['windows']
        if [w['window'] for w in windows]!=list(range(64)):
            raise ValueError('Missing, repeated or reordered windows')
        rows=canonical_windows(windows,len(windows[0]['x_aligned']))
        grid=np.asarray(rows[0]['x_aligned'])+rows[0]['front']
        if shared is None:shared=(times,scale,grid)
        elif (not np.array_equal(times,shared[0]) or scale!=shared[1] or
              grid.shape!=shared[2].shape or not np.allclose(grid,shared[2],rtol=1e-12,atol=1e-11)):
            raise ValueError('Subgroup physical sampling or grid differs between seeds')
        trials.append(rows)
    return trials, shared[1]


def subgroup(runs,seeds):
    trials,scale=prepare_runs(runs,seeds)
    phases=dict(full=trials,first=[r[:32] for r in trials],last=[r[32:] for r in trials],
                doubled=merged_pairs(trials))
    result={p:group_point(rows,center_profile=True,interpolation_sensitivity=True)
            for p,rows in phases.items()}
    return dict(phases=result,edges=group_edges(trials),seeds=sorted(seeds),
                collision_mean_free_path_m=scale)


def read_verified_profiles(folder,verified):
    """Load the complete verified matrix before any estimator-specific work."""
    folder,verified=map(Path,(folder,verified))
    report_path=folder/'extraction.json'; snapshots={report_path:digest(report_path),verified:digest(verified)}
    report=json.loads(report_path.read_text());proof=json.loads(verified.read_text())
    if (report['status']!='verified_profiles_extracted' or report['case_count']!=30 or report['dataset_count']!=10 or
            report['verification_sha256']!=snapshots[verified] or proof['status']!='matrix_outputs_verified' or
            proof['tool_sha256']!=digest(Path(__file__).with_name('verify_alsmeyer_outputs.py')) or
            report['input_verification']!=proof['inputs']):
        raise ValueError('Matching complete extraction and output verification required')
    expected_sources={n:digest(Path(__file__).with_name(n)) for n in EXTRACTION_SOURCES}
    if report['sources']!=expected_sources:raise ValueError('Extraction source version changed')
    if set(report['artifacts'])!={k+'.json' for k in KEYS}:raise ValueError('Complete ten-dataset matrix required')
    datasets={};seeds=report['input_verification']['seeds']
    for key in KEYS:
        path=folder/(key+'.json');snapshots[path]=digest(path)
        if snapshots[path]!=report['artifacts'][path.name]:raise ValueError('Extracted profile changed: '+key)
        data=json.loads(path.read_text())
        if data['status']!='profiles_only':raise ValueError('Unexpected profile artifact')
        for run in data['runs']:
            identity=key.replace('-', '/',1)+'-'+str(run['seed'])
            record=proof['cases'][identity]
            expected=record['manifest_sha256'] if key.startswith('mc3d-') else record['execution_sha256']
            if run['source_sha256']!=expected:raise ValueError('Profile source differs from verified case')
        prepare_runs(data['runs'],seeds)
        datasets[key]=data['runs']
    return datasets,seeds,report,snapshots


def analyze(folder,verified,output):
    folder,verified,output=map(Path,(folder,verified,output))
    if output.exists():raise FileExistsError(output)
    sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES}
    datasets,seeds,report,snapshots=read_verified_profiles(folder,verified)
    report_path=folder/'extraction.json'
    results={key:subgroup(runs,seeds) for key,runs in datasets.items()}
    if any(digest(p)!=h for p,h in snapshots.items()):raise ValueError('Profile inputs changed during analysis')
    if sources!={n:digest(Path(__file__).with_name(n)) for n in SOURCES}:raise ValueError('Analysis source changed')
    result=dict(status='subgroup_statistics_only',physical_acceptance='incomplete',method=METHOD,
                sources=sources,extraction_sha256=snapshots[report_path],verification_sha256=snapshots[verified],
                artifacts=report['artifacts'],datasets=results,
                limitations=['One subgroup cannot supply independent-group uncertainty or a convergence verdict.',
                             'No window is discarded; unequal sampling is rejected before pair averaging.',
                             'Method qualification and raw-output verification remain separate requirements.'])
    output.parent.mkdir(parents=True,exist_ok=True);save(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('profiles',type=Path);p.add_argument('--verified',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(analyze(a.profiles,a.verified,a.output)['status'])


if __name__=='__main__':main()
