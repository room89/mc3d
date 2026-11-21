"""Run 12 independent three-seed groups, reusing one verified baseline group.

Execution completion is not physical acceptance. Defaults select only case C
for a precision pilot. All paths in the campaign record are relative to it.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import signal
from threading import Event
from shock_validation import SEEDS,VARIANTS,configuration,digest,save,verify_run


def seed_groups(baseline,seed_start):
    groups=[list(baseline)]+[[seed_start+3*i+j for j in range(3)] for i in range(1,12)]
    flat=[s for group in groups for s in group]
    if len(baseline)!=3 or len(set(flat))!=36 or any(not isinstance(s,int) or isinstance(s,bool) or s<=0 or s>=2**31 for s in flat):
        raise ValueError('Campaign needs 36 distinct positive 31-bit seeds')
    return groups


def execution_complete(report,variants,seeds):
    if not variants or len(seeds)!=3 or len(set(seeds))!=3:return False
    return all(report.get('runs',{}).get(f'{v}-{s}',{}).get('status')=='passed'
               for v in variants for s in seeds) and all(
        report.get('sparta_runs',{}).get(f'{v}-{s}',{}).get('status')=='passed'
        for v in variants if v in ('C','D','E') for s in seeds)


def verify_baseline(folder,report,variants,seeds,binary_hash,sparta_hash):
    if report['binary_sha256']!=binary_hash or report.get('sparta_sha256')!=sparta_hash:
        raise ValueError('Baseline binary hashes do not match requested binaries')
    if not execution_complete(report,variants,seeds):raise ValueError('Selected baseline cases are incomplete')
    for v in variants:
        for s in seeds:
            run=folder/f'{v}-{s}'
            if json.loads((run/'input_si.json').read_text())!=configuration(v,s)[0]:
                raise ValueError('Baseline physical configuration differs')
            if not verify_run(run/'output',binary_hash):raise ValueError('Baseline manifest incomplete')
            if v in ('C','D','E'):
                entry=report['sparta_runs'][f'{v}-{s}']
                if not entry.get('hashes'):raise ValueError('Missing baseline SPARTA provenance')
                for name,sha in entry['hashes'].items():
                    if digest(folder/f'sparta-{v}-{s}'/name)!=sha:raise ValueError('Baseline SPARTA output changed')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','binary','sparta-bin','output-dir'):p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--variants',nargs='+',choices=list(VARIANTS),default=['C'])
    p.add_argument('--seed-start',type=int,default=10000);p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--load-mode',choices=['full','balanced','quiet'],default='balanced')
    p.add_argument('--prepare-only',action='store_true');p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if a.jobs<1 or len(set(a.variants))!=len(a.variants):p.error('Positive job count and distinct variants required')
    variants=[v for v in VARIANTS if v in a.variants]
    root=a.output_dir.resolve();baseline=a.baseline.resolve()
    if baseline==root or baseline.is_relative_to(root):raise ValueError('Baseline must remain outside the new campaign')
    raw=(baseline/'validation.json').read_bytes();report=json.loads(raw)
    seeds=report.get('seeds',SEEDS);groups=seed_groups(seeds,a.seed_start)
    binary=a.binary.resolve();sparta=a.sparta_bin.resolve()
    binary_hash=digest(binary);sparta_hash=digest(sparta)
    verify_baseline(baseline,report,variants,seeds,binary_hash,sparta_hash)
    from hashlib import sha256
    sources=('shock_independent_campaign.py','shock_validation.py','normal_shock.py','physical_config.py',
             'argon_vhs.py','sparta_shock.py','sparta_preflight.py','control_simulation.py',
             'analyze_shock_validation.py','shock_analysis.py')
    protocol=dict(variants=variants,seeds=groups,baseline=os.path.relpath(baseline,root),
        baseline_report_sha256=sha256(raw).hexdigest(),binary=os.path.relpath(binary,root),
        binary_sha256=binary_hash,sparta_binary=os.path.relpath(sparta,root),sparta_sha256=sparta_hash,
        sources={n:digest(Path(__file__).with_name(n)) for n in sources})
    root.mkdir(parents=True,exist_ok=a.resume);manifest=root/'campaign.json'
    state=json.loads(manifest.read_text()) if manifest.exists() else dict(protocol=protocol,groups={})
    if state['protocol']!=protocol:raise ValueError('Campaign protocol changed; cannot resume')
    if state.get('owner_pid'):
        try:os.kill(state['owner_pid'],0)
        except ProcessLookupError:pass
        else:raise ValueError('Another campaign runner may still be active')
    state.update(owner_pid=os.getpid(),status='preparing' if a.prepare_only else 'running',jobs=a.jobs,load_mode=a.load_mode)
    state['groups']['g00']=dict(path=protocol['baseline'],seeds=groups[0],status='execution_complete')
    save(manifest,state)
    stop=Event()

    def interrupt(signum,frame):
        raise KeyboardInterrupt('Campaign interrupted')
    signal.signal(signal.SIGTERM,interrupt)

    def run_group(i):
        name=f'g{i:02d}';folder=root/name
        command=[sys.executable,str(Path(__file__).with_name('shock_validation.py')),
            '--binary',str(binary),'--sparta-bin',str(sparta),'--output-dir',str(folder),
            '--no-time-budget','--load-mode',a.load_mode,'--variants',*variants,
            '--seeds',*map(str,groups[i])]
        if folder.exists():command.append('--resume')
        if a.prepare_only:command.append('--prepare-only')
        with (root/(name+'.log')).open('a') as log:
            if stop.is_set():raise RuntimeError('Campaign stopped')
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                while process.poll() is None:
                    if stop.is_set():raise RuntimeError('Campaign stopped')
                    try:process.wait(timeout=1)
                    except subprocess.TimeoutExpired:pass
            finally:
                if process.poll() is None:
                    try:os.killpg(process.pid,signal.SIGTERM)
                    except ProcessLookupError:pass
                    try:process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        try:os.killpg(process.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                        process.wait()
        code=process.returncode
        data=json.loads((folder/'validation.json').read_text()) if (folder/'validation.json').exists() else {}
        if data and (data.get('seeds')!=groups[i] or data.get('variants')!=variants):
            raise ValueError('Child group identities changed')
        if data and (data.get('binary_sha256')!=binary_hash or data.get('sparta_sha256')!=sparta_hash):
            raise ValueError('Child binary identities changed')
        complete=execution_complete(data,variants,groups[i])
        status='execution_complete' if complete and code in (0,1,2) else 'prepared' if a.prepare_only and code==2 else 'incomplete'
        return name,dict(path=name,seeds=groups[i],status=status,returncode=code)

    pool=ThreadPoolExecutor(max_workers=min(a.jobs,11))
    try:
        futures=[pool.submit(run_group,i) for i in range(1,12)]
        for future in as_completed(futures):
            name,item=future.result();state['groups'][name]=item;save(manifest,state)
            print(name,item['status'],flush=True)
        if (baseline/'validation.json').read_bytes()!=raw:raise ValueError('Baseline report changed during campaign')
        if any(digest(Path(__file__).with_name(n))!=sha for n,sha in protocol['sources'].items()):
            raise ValueError('Campaign source code changed during execution')
        complete=len(state['groups'])==12 and all(g['status']=='execution_complete' for g in state['groups'].values())
        prepared=len(state['groups'])==12 and all(g['status'] in ('prepared','execution_complete') for g in state['groups'].values())
        state['status']='execution_complete' if complete else 'prepared' if a.prepare_only and prepared else 'incomplete'
        state['physical_acceptance']='incomplete; aggregate uncertainty and physical gates are separate'
    except BaseException:
        stop.set()
        state['status']='interrupted_or_error'
        raise
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
        state.pop('owner_pid',None);save(manifest,state)
    return 0 if state['status']=='execution_complete' else 2


if __name__=='__main__':raise SystemExit(main())
