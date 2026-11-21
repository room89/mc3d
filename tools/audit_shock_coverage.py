"""Exact binomial diagnostics for completed synthetic coverage reports.

No acceptance verdict: non-rejection does not prove nominal coverage. Scenario
random streams may be shared, so counts are never pooled across scenarios.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
from shock_validation import digest,save


@lru_cache(maxsize=None)
def log_coefficients(n):
    return [math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1) for k in range(n+1)]


def tail(n,k,p,upper=True):
    if p==0:return float(k<=0) if upper else 1.
    if p==1:return 1. if upper else float(k>=n)
    indices=range(k,n+1) if upper else range(k+1)
    coeff=log_coefficients(n)
    return min(1.,math.fsum(math.exp(coeff[i]+i*math.log(p)+(n-i)*math.log1p(-p)) for i in indices))


def exact_interval(k,n,alpha=.05):
    if type(k)!=int or type(n)!=int or not 0<=k<=n or n<1 or not 0<alpha<1:
        raise ValueError('Invalid binomial counts or alpha')
    def solve(upper):
        lo,hi=0.,1.
        for _ in range(60):
            mid=(lo+hi)/2;prob=tail(n,k,mid,upper)
            if (prob<alpha/2)==upper:lo=mid
            else:hi=mid
        return (lo+hi)/2
    return [0. if k==0 else solve(True),1. if k==n else solve(False)]


def case_counts(case):
    rows=case['runs'];n=case['trial_count']
    if n<1 or len(rows)!=n or len({r['repetition'] for r in rows})!=n:
        raise ValueError('Missing or repeated trials')
    counts=[0,0,0]
    for row in rows:
        if 'error' in row:continue
        flags=row['with_processing_covers']+[row['density_covers']]
        if len(flags)!=3 or any(type(v)!=bool for v in flags):raise ValueError('Boolean coverage outcomes required')
        counts=[k+int(v) for k,v in zip(counts,flags)]
    return n,counts


def audit(reports):
    total_cases=sum(len(r['cases']) for r in reports)
    if not total_cases:raise ValueError('No cases')
    family_tests=3*total_cases;cases=[]
    for ri,report in enumerate(reports):
        for c in report['cases']:
            n,counts=case_counts(c);metrics={}
            for name,k in zip(('width_10_90','gradient_thickness','density'),counts):
                metrics[name]=dict(hits=k,trials=n,coverage=k/n,exact_95_interval=exact_interval(k,n),
                    familywise_95_interval=exact_interval(k,n,.05/family_tests),
                    one_sided_undercoverage_p=tail(n,k,.95,False))
            cases.append(dict(report_index=ri,shape=c['shape'],phi=c['phi'],spatial=c['spatial'],noise=c['noise'],metrics=metrics))
    return dict(status='diagnostic_only_no_acceptance_verdict',nominal_coverage=.95,
        family_tests=family_tests,cases=cases,
        method='Clopper-Pearson exact two-sided intervals; Bonferroni family intervals across all supplied scenario/metric pairs; exact one-sided binomial test against p=0.95.',
        limitations=['Independent Bernoulli trial assumption within each scenario; no independence assumed between scenarios.',
            'No counts pooled across scenarios, because scenario random streams may be shared.',
            'Non-rejection is not proof of nominal coverage; confidence intervals may be too wide to qualify a method.',
            'Choice of report family and inspection after model selection must be disclosed; this is not a retrospective acceptance rule.',
            'Physical precision and coverage tolerances are unchanged.'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('reports',nargs='+',type=Path)
    p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    reports=[json.loads(p.read_text()) for p in a.reports]
    result=audit(reports);result['report_hashes']=[digest(p) for p in a.reports]
    result['tool_sha256']=digest(Path(__file__));save(a.output,result)
    print('Diagnostic audit saved; no acceptance verdict')


if __name__=='__main__':main()
