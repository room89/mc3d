"""Prepare or execute the additional C subgroups in a verified pooled plan.

Preparation is the default. Execution requires --execute and completed matching
coverage artifacts; their adequacy must be assessed separately before invocation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from threading import Event
from reanalyze_pooled_shock import validate_plan,check_coverage
from reanalyze_shock_campaign import verify_mc_inputs
from shock_independent_campaign import verify_baseline,execution_complete
from shock_validation import digest,save


def child_status(report,seeds,plan,returncode,execute):
    if report.get('seeds')!=seeds or report.get('variants')!=['C']:
        raise ValueError('Child identities differ')
    if report.get('binary_sha256')!=plan['binary_sha256'] or report.get('sparta_sha256')!=plan['sparta_sha256']:
        raise ValueError('Child binaries differ')
    if execution_complete(report,['C'],seeds) and returncode in (0,1,2):return 'execution_complete'
    if not execute and returncode==2:return 'prepared'
    return 'incomplete'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('--execute',action='store_true')
    p.add_argument('--coverage',type=Path);p.add_argument('--half-coverage',type=Path)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--load-mode',choices=['full','balanced','quiet'],default='balanced')
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    if not 1<=a.jobs<=6:p.error('Choose 1–6 concurrent subgroup jobs')
    if a.execute and (not a.coverage or not a.half_coverage):p.error('Execution requires completed coverage artifacts')
    root=a.plan.resolve().parent;plan_hash=digest(a.plan);plan=json.loads(a.plan.read_text());groups=validate_plan(plan,root)
    baseline=root/plan['baseline_campaign']
    if digest(baseline)!=plan['baseline_campaign_sha256']:raise ValueError('Baseline campaign changed')
    binary=(root/plan['binary']).resolve();sparta=(root/plan['sparta_binary']).resolve()
    if digest(binary)!=plan['binary_sha256'] or digest(sparta)!=plan['sparta_sha256']:raise ValueError('Binary changed')
    for g in groups:
        item=g['subgroups'][0];folder=root/item['path'];report=json.loads((folder/'validation.json').read_text())
        if digest(folder/'validation.json')!=item['validation_sha256']:raise ValueError('Existing report changed')
        verify_baseline(folder,report,['C'],item['seeds'],plan['binary_sha256'],plan['sparta_sha256'])
        for s in item['seeds']:verify_mc_inputs(folder/f'C-{s}','C',s)
    # Only planned descendants are writable; existing baseline data are read-only.
    tasks=[]
    for g in groups:
        for i,item in enumerate(g['subgroups'][1:],1):
            name=g['name']+f'-s{i}';folder=(root/item['path']).resolve()
            if folder!=root/name:raise ValueError('Unexpected new subgroup destination')
            tasks.append((name,folder,item['seeds']))
    qualification=[check_coverage(a.coverage,64),check_coverage(a.half_coverage,32)] if a.execute else None
    names=('run_pooled_shock.py','shock_validation.py','normal_shock.py','physical_config.py','argon_vhs.py',
           'sparta_shock.py','sparta_preflight.py','control_simulation.py','analyze_shock_validation.py','shock_analysis.py')
    protocol=dict(plan_sha256=plan_hash,sources={n:digest(Path(__file__).with_name(n)) for n in names})
    manifest=root/'execution.json'
    # Advisory lock closes the check-then-start race; kernel releases it on exit.
    lock=(root/'execution.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise ValueError('Another pooled runner is active')
    state=json.loads(manifest.read_text()) if manifest.exists() else dict(protocol=protocol,subgroups={})
    if manifest.exists() and not a.resume:raise ValueError('Use --resume for existing execution manifest')
    if state['protocol']!=protocol:raise ValueError('Execution protocol changed')
    state.update(status='running' if a.execute else 'preparing',coverage=qualification,jobs=a.jobs,load_mode=a.load_mode)
    save(manifest,state);stop=Event()
    def interrupt(signum,frame):raise KeyboardInterrupt('Pooled campaign interrupted')
    signal.signal(signal.SIGTERM,interrupt)
    def run(task):
        name,folder,seeds=task
        command=[sys.executable,str(Path(__file__).with_name('shock_validation.py')),'--binary',str(binary),
            '--sparta-bin',str(sparta),'--output-dir',str(folder),'--no-time-budget','--load-mode',a.load_mode,
            '--variants','C','--seeds',*map(str,seeds)]
        if folder.exists():command.append('--resume')
        if not a.execute:command.append('--prepare-only')
        with (root/(name+'.log')).open('a') as log:
            if stop.is_set():raise RuntimeError('Campaign stopped')
            proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                while proc.poll() is None:
                    if stop.is_set():raise RuntimeError('Campaign stopped')
                    try:proc.wait(timeout=1)
                    except subprocess.TimeoutExpired:pass
            finally:
                if proc.poll() is None:
                    try:os.killpg(proc.pid,signal.SIGTERM)
                    except ProcessLookupError:pass
                    try:proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        try:os.killpg(proc.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                        proc.wait()
        report=json.loads((folder/'validation.json').read_text())
        status=child_status(report,seeds,plan,proc.returncode,a.execute)
        return name,dict(status=status,seeds=seeds,path=name,returncode=proc.returncode)
    pool=ThreadPoolExecutor(max_workers=a.jobs)
    try:
        futures=[pool.submit(run,t) for t in tasks]
        for future in as_completed(futures):
            name,item=future.result();state['subgroups'][name]=item;save(manifest,state)
            print(name,item['status'],flush=True)
        if digest(a.plan)!=plan_hash or digest(baseline)!=plan['baseline_campaign_sha256']:
            raise ValueError('Plan or baseline changed')
        if any(digest(Path(__file__).with_name(n))!=h for n,h in protocol['sources'].items()):
            raise ValueError('Execution source changed')
        statuses=[v['status'] for v in state['subgroups'].values()]
        state['status']='execution_complete' if len(statuses)==60 and all(s=='execution_complete' for s in statuses) else (
            'prepared' if not a.execute and len(statuses)==60 and all(s in ('prepared','execution_complete') for s in statuses) else 'incomplete')
    except BaseException:
        stop.set();state['status']='interrupted_or_error';raise
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
        state['physical_acceptance']='incomplete; execution is not physical acceptance';save(manifest,state);lock.close()
    return 0 if state['status']=='execution_complete' else 2


if __name__=='__main__':raise SystemExit(main())
