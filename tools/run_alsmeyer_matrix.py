"""Prepare or execute the verified primary-reference matrix with pinned binaries.

Execution is separate from method qualification and physical acceptance. Resume
reuses only verified completed cases; incomplete attempts are never overwritten.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
from threading import Event
import time
from control_simulation import set_mode
from reanalyze_pooled_shock import check_coverage
from shock_validation import digest, save
from verify_alsmeyer_matrix import verify
from verify_alsmeyer_outputs import mc3d_outputs, sparta_outputs

SOURCES = ('run_alsmeyer_matrix.py', 'reanalyze_pooled_shock.py', 'verify_alsmeyer_outputs.py', 'verify_alsmeyer_matrix.py',
           'prepare_alsmeyer_matrix.py', 'export_alsmeyer_sparta.py', 'control_simulation.py',
           'physical_config.py', 'shock_validation.py', 'normal_shock.py', 'argon_vhs.py', 'sparta_preflight.py')


def read(path):
    return json.loads(path.read_text())


def input_hashes(kind, folder):
    names = ('solver.json', 'input_si.json', 'units.json', 'reference.json') if kind == 'mc3d' else (
        'in.shock', 'argon.species', 'argon.vss', 'reference.json')
    return {name: digest(folder/name) for name in names}


def reuse(kind, folder, binary_hash, inputs):
    record = folder/'execution.json'
    if not record.exists():
        # Do not adopt unrecorded runs or mix fresh output with a previous attempt.
        if (folder/'output').exists() or (folder/'run.log').exists() or any(folder.glob('profile.*')):
            raise ValueError('Unrecorded output exists: '+folder.name)
        return None
    execution = read(record)
    if execution.get('status') != 'outputs_verified' or execution.get('returncode') != 0:
        raise ValueError('Incomplete attempt must be inspected; refusing overwrite: '+folder.name)
    if execution['binary_sha256'] != binary_hash or execution['input_hashes'] != inputs:
        raise ValueError('Completed execution provenance changed: '+folder.name)
    if execution['log_sha256'] != digest(folder/'run.log'):
        raise ValueError('Completed log changed: '+folder.name)
    result = mc3d_outputs(folder, binary_hash) if kind == 'mc3d' else sparta_outputs(folder, execution, binary_hash)
    if result != execution['verification']:
        raise ValueError('Completed output identity changed: '+folder.name)
    return dict(**result, execution_sha256=digest(record), reused=True)


def run_case(kind, folder, binary, binary_hash, expected_inputs, load_mode, stop):
    if stop.is_set(): raise RuntimeError('Campaign stopped')
    if input_hashes(kind, folder) != expected_inputs or digest(binary) != binary_hash:
        raise ValueError('Input or binary changed before execution')
    previous = reuse(kind, folder, binary_hash, expected_inputs)
    if previous is not None: return previous
    if kind == 'mc3d':
        control = folder/'load.control'; set_mode(control, load_mode)
        command = [str(binary), '--config', str(folder/'solver.json'), '--output-dir', str(folder/'output'),
                   '--runtime-control', str(control)]
    else:
        command = [str(binary), '-in', 'in.shock']
    record = folder/'execution.json'
    execution = dict(status='starting', physical_acceptance='incomplete', binary_sha256=binary_hash,
                     input_hashes=expected_inputs, command=command, started_unix=time.time())
    save(record, execution)
    try:
        with (folder/'run.log').open('x') as log:
            proc = subprocess.Popen(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                                    start_new_session=True)
            try:
                execution.update(status='running', pid=proc.pid); save(record, execution)
                while proc.poll() is None:
                    if stop.is_set(): raise RuntimeError('Campaign stopped')
                    try: proc.wait(timeout=1)
                    except subprocess.TimeoutExpired: pass
            finally:
                if proc.poll() is None:
                    try: os.killpg(proc.pid, signal.SIGTERM)
                    except ProcessLookupError: pass
                    try: proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        try: os.killpg(proc.pid, signal.SIGKILL)
                        except ProcessLookupError: pass
                        proc.wait()
                execution['returncode'] = proc.returncode
        execution['log_sha256'] = digest(folder/'run.log')
        if proc.returncode != 0: raise ValueError('Solver exited with code '+str(proc.returncode))
        if input_hashes(kind, folder) != expected_inputs or digest(binary) != binary_hash:
            raise ValueError('Input or binary changed during execution')
        if kind == 'sparta':
            execution['outputs'] = [dict(path=p.name, sha256=digest(p)) for p in sorted(folder.glob('profile.*'))]
            result = sparta_outputs(folder, execution, binary_hash)
        else:
            result = mc3d_outputs(folder, binary_hash)
        execution.update(status='outputs_verified', verification=result)
    except BaseException as error:
        execution.update(status='interrupted_or_error', error=str(error)); raise
    finally:
        execution['finished_unix'] = time.time(); save(record, execution)
    return dict(**result, execution_sha256=digest(record), reused=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mc3d', type=Path, required=True); p.add_argument('--sparta', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True); p.add_argument('--sparta-bin', type=Path, required=True)
    p.add_argument('--execute', action='store_true'); p.add_argument('--resume', action='store_true')
    p.add_argument('--coverage', type=Path); p.add_argument('--half-coverage', type=Path)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--load-mode', choices=('full','balanced','quiet'), default='full')
    a = p.parse_args()
    if not 1 <= a.jobs <= 3: p.error('Choose 1–3 concurrent solver jobs')
    if a.execute and not (a.coverage and a.half_coverage): p.error('Completed matching coverage reports required')
    roots = dict(mc3d=a.mc3d.resolve(), sparta=a.sparta.resolve())
    binaries = dict(mc3d=a.binary.resolve(), sparta=a.sparta_bin.resolve())
    with ExitStack() as stack:
        for root in sorted(set(roots.values())):
            lock = stack.enter_context((root/'reference-execution.lock').open('a'))
            try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError: raise ValueError('Another reference runner is active')
        proof = verify(roots['mc3d'], roots['sparta'], a.evidence)
        hashes = {k:digest(v) for k,v in binaries.items()}
        protocol = dict(inputs=proof, binary_sha256=hashes, binaries={k:str(v) for k,v in binaries.items()},
                        roots={k:str(v) for k,v in roots.items()},
                        sources={n:digest(Path(__file__).with_name(n)) for n in SOURCES})
        manifest = roots['mc3d']/'execution.json'
        if manifest.exists() and not a.resume: raise ValueError('Use --resume for an existing campaign')
        state = read(manifest) if manifest.exists() else dict(protocol=protocol, cases={})
        if state['protocol'] != protocol: raise ValueError('Campaign protocol changed')
        tasks = []
        for kind, filename in [('mc3d','preparation.json'), ('sparta','export.json')]:
            for case in read(roots[kind]/filename)['cases']:
                folder = roots[kind]/case['folder']; inputs = input_hashes(kind, folder)
                # Check all existing output before launching any solver.
                reuse(kind, folder, hashes[kind], inputs)
                tasks.append((kind, folder, binaries[kind], hashes[kind], inputs))
        coverage = [check_coverage(a.coverage,64), check_coverage(a.half_coverage,32)] if a.execute else None
        if state.get('coverage') is not None and coverage != state['coverage']:
            raise ValueError('Coverage evidence changed; do not replace a running campaign protocol')
        state.update(status='running' if a.execute else 'prepared', jobs=a.jobs, load_mode=a.load_mode,
                     coverage=coverage, physical_acceptance='incomplete', planned_cases=len(tasks))
        save(manifest,state)
        if not a.execute:
            print('Prepared 30-case execution protocol; no solver launched.'); return 0
        stop = Event(); pool = ThreadPoolExecutor(max_workers=a.jobs)
        def interrupt(signum, frame): raise KeyboardInterrupt('Reference campaign interrupted')
        old_handler = signal.signal(signal.SIGTERM, interrupt)
        try:
            futures = {pool.submit(run_case,*task,a.load_mode,stop):task for task in tasks}
            for future in as_completed(futures):
                kind, folder, *_ = futures[future]; result = future.result()
                state['cases'][kind+'/'+folder.name] = result; save(manifest,state)
                print(kind,folder.name,result['status'],flush=True)
            if proof != verify(roots['mc3d'],roots['sparta'],a.evidence): raise ValueError('Prepared matrix changed')
            if any(digest(binaries[k]) != v for k,v in hashes.items()): raise ValueError('Binary changed')
            if any(digest(Path(__file__).with_name(n)) != h for n,h in protocol['sources'].items()):
                raise ValueError('Execution source changed')
            state['status'] = 'matrix_outputs_verified'
        except BaseException:
            stop.set(); state['status'] = 'interrupted_or_error'; raise
        finally:
            pool.shutdown(wait=True,cancel_futures=True); save(manifest,state)
            signal.signal(signal.SIGTERM,old_handler)
    return 0


if __name__ == '__main__': raise SystemExit(main())
