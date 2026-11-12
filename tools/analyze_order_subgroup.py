"""Build both alignment orders from verified reference profiles; no acceptance."""
import argparse
from pathlib import Path
from analyze_alsmeyer_subgroup import SOURCES as BASE_SOURCES, prepare_runs, read_verified_profiles
from reanalyze_shock_campaign import merged_pairs, group_edges
from shock_order_sensitivity import METHOD, subgroup_point
from shock_validation import digest, save

SOURCES=tuple(sorted(set(BASE_SOURCES)|{'analyze_order_subgroup.py','shock_order_sensitivity.py'}))


def subgroup(runs,seeds):
    trials,scale=prepare_runs(runs,seeds)
    phases=dict(full=trials,first=[r[:32] for r in trials],last=[r[32:] for r in trials],doubled=merged_pairs(trials))
    points={phase:subgroup_point(rows) for phase,rows in phases.items()}
    return dict(phases=points,edges=group_edges(trials),seeds=sorted(seeds),collision_mean_free_path_m=scale,
                source_windows_per_seed=dict(full=64,first=32,last=32,doubled=64),
                effective_windows_per_seed={phase:len(rows[0]) for phase,rows in phases.items()})


def analyze(folder,verified,output):
    folder,verified,output=map(Path,(folder,verified,output))
    if output.exists():raise FileExistsError(output)
    sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES}
    datasets,seeds,report,snapshots=read_verified_profiles(folder,verified)
    results={key:subgroup(runs,seeds) for key,runs in datasets.items()}
    if any(digest(path)!=sha for path,sha in snapshots.items()):raise ValueError('Profile inputs changed during analysis')
    if sources!={n:digest(Path(__file__).with_name(n)) for n in SOURCES}:raise ValueError('Analysis source changed')
    result=dict(status='order_subgroup_statistics_only',physical_acceptance='incomplete',method=METHOD,
        sources=sources,extraction_sha256=snapshots[folder/'extraction.json'],verification_sha256=snapshots[verified],
        artifacts=report['artifacts'],datasets=results,
        limitations=['One subgroup cannot supply twelve-group uncertainty, method qualification or physical acceptance.',
          'Both orders use the same retained windows; they are processing alternatives, not independent physical observations.',
          'Doubled phase retains all 64 source windows as 32 adjacent equal-duration pairs.',
          'Production qualification and numerical/experimental gates remain separate requirements.'])
    output.parent.mkdir(parents=True,exist_ok=True);save(output,result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('profiles',type=Path)
    parser.add_argument('--verified',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();print(analyze(args.profiles,args.verified,args.output)['status'])

if __name__=='__main__':main()
