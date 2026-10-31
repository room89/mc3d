#!/usr/bin/env python3
"""Run MC3D with sampled process/system memory diagnostics (macOS and Linux)."""
from __future__ import annotations

import argparse
import bisect
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import statistics
import subprocess
import sys
import threading
import time

PROCESS_KEYS = ('rss_bytes', 'cpu_seconds', 'minor_faults', 'major_faults')
SYSTEM_KEYS = ('available_bytes', 'free_bytes', 'swap_used_bytes', 'swap_in_bytes',
               'swap_out_bytes', 'compressed_bytes', 'compression_bytes',
               'decompression_bytes', 'page_in_bytes', 'page_out_bytes',
               'psi_some_seconds', 'psi_full_seconds')
COUNTERS = ('process.cpu_seconds', 'process.minor_faults', 'process.major_faults',
            'system.swap_in_bytes', 'system.swap_out_bytes', 'system.compression_bytes',
            'system.decompression_bytes', 'system.page_in_bytes', 'system.page_out_bytes',
            'system.psi_some_seconds', 'system.psi_full_seconds')


def json_line(stream, value):
    stream.write(json.dumps(value, allow_nan=False) + '\n')
    stream.flush()


def cpu_time(value):
    """ps TIME: [[days-]hours:]minutes:seconds, including fractional seconds."""
    days, separator, clock = value.partition('-')
    if not separator:
        clock, days = days, '0'
    fields = [float(x) for x in clock.split(':')]
    if len(fields) not in (2, 3):
        raise ValueError('Invalid ps CPU time')
    seconds = 0.0
    for field in fields:
        seconds = seconds * 60 + field
    return int(days) * 86400 + seconds


def parse_ps(text):
    rss, cpu = text.split()
    return {'rss_bytes': int(rss) * 1024, 'cpu_seconds': cpu_time(cpu)}


def parse_vm_stat(text):
    page_size = int(re.search(r'page size of (\d+) bytes', text)[1])
    fields = {key.strip(): int(value) for key, value in
              re.findall(r'^([^:\n]+):\s*(\d+)\.?\s*$', text, re.M)}
    mapping = {'free_bytes': 'Pages free', 'compressed_bytes': 'Pages occupied by compressor',
               'compression_bytes': 'Compressions', 'decompression_bytes': 'Decompressions',
               'swap_in_bytes': 'Swapins', 'swap_out_bytes': 'Swapouts',
               'page_in_bytes': 'Pageins', 'page_out_bytes': 'Pageouts'}
    return {key: fields[source] * page_size for key, source in mapping.items() if source in fields}


def parse_swap(text):
    match = re.search(r'used\s*=\s*([\d.]+)([KMGT]?)', text)
    if not match:
        raise ValueError('Missing used swap')
    return {'swap_used_bytes': int(float(match[1]) * 1024 ** (' KMGT'.index(match[2]) if match[2] else 0))}


def parse_linux_stat(text, ticks, page_size):
    # comm can contain spaces and parentheses; fields after its final ')' start at field 3.
    fields = text[text.rindex(')') + 2:].split()
    return {'rss_bytes': int(fields[21]) * page_size,
            'cpu_seconds': (int(fields[11]) + int(fields[12])) / ticks,
            'minor_faults': int(fields[7]), 'major_faults': int(fields[9])}


def parse_meminfo(text):
    fields = {key: int(value) * 1024 for key, value in
              re.findall(r'^(\w+):\s*(\d+)\s+kB', text, re.M)}
    result = {}
    for target, source in [('available_bytes', 'MemAvailable'), ('free_bytes', 'MemFree')]:
        if source in fields:
            result[target] = fields[source]
    if 'SwapTotal' in fields and 'SwapFree' in fields:
        result['swap_used_bytes'] = fields['SwapTotal'] - fields['SwapFree']
    return result


def parse_linux_vmstat(text, page_size):
    fields = dict(line.split() for line in text.splitlines() if line.strip())
    return {key: int(fields[source]) * page_size for key, source in
            [('swap_in_bytes', 'pswpin'), ('swap_out_bytes', 'pswpout')] if source in fields}


def parse_psi(text):
    return {'psi_' + scope + '_seconds': int(total) / 1_000_000 for scope, total in
            re.findall(r'^(some|full).*?\btotal=(\d+)', text, re.M)}


class Sampler:
    def __init__(self, pid):
        self.pid = pid
        self.platform = platform.system()
        self.page_size = os.sysconf('SC_PAGE_SIZE')
        self.ticks = os.sysconf('SC_CLK_TCK')

    def sample(self):
        result = {'process': dict.fromkeys(PROCESS_KEYS), 'system': dict.fromkeys(SYSTEM_KEYS),
                  'errors': {}}

        def read(scope, name, function):
            try:
                result[scope].update(function())
            except (OSError, ValueError, IndexError, TypeError, subprocess.SubprocessError) as exc:
                result['errors'][name] = str(exc)

        def command(args):
            return subprocess.run(args, check=True, capture_output=True, text=True,
                                  timeout=1, env={**os.environ, 'LC_ALL': 'C'}).stdout

        if self.platform == 'Linux':
            read('process', 'proc_stat', lambda: parse_linux_stat(
                Path(f'/proc/{self.pid}/stat').read_text(), self.ticks, self.page_size))
            read('system', 'meminfo', lambda: parse_meminfo(Path('/proc/meminfo').read_text()))
            read('system', 'vmstat', lambda: parse_linux_vmstat(Path('/proc/vmstat').read_text(), self.page_size))
            read('system', 'psi', lambda: parse_psi(Path('/proc/pressure/memory').read_text()))
        elif self.platform == 'Darwin':
            read('process', 'ps', lambda: parse_ps(command(['/bin/ps', '-p', str(self.pid), '-o', 'rss=,time='])))
            read('system', 'vm_stat', lambda: parse_vm_stat(command(['/usr/bin/vm_stat'])))
            read('system', 'swapusage', lambda: parse_swap(command(['/usr/sbin/sysctl', 'vm.swapusage'])))
        else:
            raise RuntimeError('Only macOS and Linux are supported')
        result['unavailable'] = {f'{scope}.{key}': 'unsupported or unavailable; see errors'
                                 for scope in ('process', 'system')
                                 for key, value in result[scope].items() if value is None}
        return result


def counter_deltas(previous, current):
    values, resets = {}, []
    for name in COUNTERS:
        scope, key = name.split('.')
        before = previous[scope].get(key) if previous else None
        after = current[scope].get(key)
        values[name] = None
        if before is not None and after is not None:
            if after >= before:
                values[name] = after - before
            else:
                resets.append(name)
    return values, resets


def parse_event(line, timestamp):
    match = re.search(r'text=(.*)', line)
    if not match:
        return None
    text = match[1]
    event = {'t': timestamp}
    if match := re.fullmatch(r'Time step:\s*(\d+)\s*', text):
        event.update(kind='step', step=int(match[1]))
    elif match := re.fullmatch(r'stopwatch:([^;]+);total_time:([\d.eE+-]+)\s*', text):
        duration = float(match[2])
        if not math.isfinite(duration) or duration < 0:
            return None
        event.update(kind='phase', phase=match[1], duration_s=duration,
                     estimated_start=timestamp-duration)
    elif match := re.fullmatch(r'particle_buffer_size:\s*(\d+)\s*', text):
        event.update(kind='particles', count=int(match[1]))
    else:
        return None
    return event


def read_output(pipe, raw, events, start, errors):
    pending = b''
    discarding = False
    try:
        while chunk := pipe.read1(65536):
            raw.write(chunk)
            # Raw bytes are preserved; oversized individual lines are not parsed.
            for part in chunk.splitlines(keepends=True):
                complete = part.endswith(b'\n')
                if not discarding:
                    pending += part
                    if len(pending) > 65536:
                        pending, discarding = b'', True
                if complete:
                    if not discarding:
                        event = parse_event(pending.decode('utf-8', errors='replace').rstrip(), time.monotonic()-start)
                        if event:
                            json_line(events, event)
                    pending, discarding = b'', False
        # An unterminated final line is retained in run.log, not treated as a complete event.
        raw.flush()
    except Exception as exc:
        errors.append(str(exc))


def records(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def summarize(directory, interval, observation_end=None):
    samples = list(records(directory/'samples.jsonl'))
    timestamps = [s['t'] for s in samples]
    steps, initial = [], []
    current = None
    for event in records(directory/'events.jsonl'):
        if event['kind'] == 'step':
            if current is not None and not current['complete']:
                current['end'] = event['t']
            current = {'step': event['step'], 'start': event['t'], 'complete': False, 'phases': []}
            steps.append(current)
        elif current is None:
            initial.append(event)
        elif event['kind'] == 'phase':
            current['phases'].append(event)
            if event['phase'] == 'ClusterStep':
                current.update(complete=True, end=event['t'], duration_s=event['duration_s'])
        else:
            current['outgoing_particles'] = event['count']
    for step in steps:
        end = step.get('end', max(step['start'], observation_end if observation_end is not None else (timestamps[-1] if timestamps else step['start'])))
        lo = max(0, bisect.bisect_right(timestamps, step['start'])-1)
        hi = min(len(samples), bisect.bisect_left(timestamps, end)+1)
        window = samples[lo:hi]
        gaps = [b['t']-a['t'] for a,b in zip(window, window[1:])]
        deltas, resets = counter_deltas(window[0], window[-1]) if len(window)>1 else (dict.fromkeys(COUNTERS), [])
        resets = sorted(set(resets).union(*(s.get('counter_resets', []) for s in window[1:])))
        for key in resets:
            deltas[key] = None
        rss = [s['process']['rss_bytes'] for s in window if s['process']['rss_bytes'] is not None]
        span = window[-1]['t']-window[0]['t'] if len(window)>1 else 0
        cpu = deltas['process.cpu_seconds']
        system_ranges = {}
        for key in SYSTEM_KEYS:
            values = [s['system'][key] for s in window if s['system'].get(key) is not None]
            system_ranges[key] = {'min': min(values, default=None), 'max': max(values, default=None)}
        step['metrics'] = {'cpu_percent': cpu/span*100 if cpu is not None and span>0 else None,
                           'system_ranges': system_ranges,'sample_count': len(window), 'window_start': window[0]['t'] if window else None,
                           'window_end': window[-1]['t'] if window else None,
                           'max_observed_rss_bytes': max(rss, default=None), 'counter_deltas': deltas,
                           'counter_resets': resets, 'max_gap_s': max(gaps, default=None),
                           'coarse': end-step['start'] < interval,
                           'incomplete_coverage': not window or window[0]['t']>step['start'] or window[-1]['t']<end,
                           'sampling_gap': any(g > interval*1.5 for g in gaps)}
    durations = sorted(s['duration_s'] for s in steps if s['complete'])
    rss = [s['process']['rss_bytes'] for s in samples if s['process']['rss_bytes'] is not None]
    result = {'initialization_duration_s': steps[0]['start'] if steps else None,
              'completed_steps': len(durations), 'incomplete_steps': [s['step'] for s in steps if not s['complete']],
              'median_step_s': statistics.median(durations) if durations else None,
              'p95_step_s': durations[math.ceil(.95*len(durations))-1] if durations else None,
              'max_step_s': max(durations, default=None), 'max_observed_rss_bytes': max(rss, default=None),
              'initialization_events': initial, 'steps': steps,
              'slowest_steps': sorted((s for s in steps if s['complete']), key=lambda s:s['duration_s'], reverse=True)[:10],
              'limitations': 'Sampled RSS is not peak RSS. Windows include bracketing samples and may overlap. '
              'Log receipt times include pipe delivery delay. System swap is not process swap; correlation does not prove causality.'}
    (directory/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    lines = ['# MC3D memory diagnostics', '', result['limitations'], '',
             f"Time before first step (including process startup): {result['initialization_duration_s']} s",
             f"Completed steps: {len(durations)}; incomplete: {result['incomplete_steps']}",
             f"Step seconds: median={result['median_step_s']}, p95={result['p95_step_s']}, max={result['max_step_s']}",
             f"Maximum observed RSS bytes: {result['max_observed_rss_bytes']}", '',
             '| Step | Seconds | RSS bytes | CPU seconds Δ | Major faults Δ | System swap in/out bytes Δ | Coverage |',
             '|---|---:|---:|---:|---:|---|---|']
    for step in result['slowest_steps']:
        m=step['metrics']; d=m['counter_deltas']
        flags=', '.join(k for k in ('coarse','incomplete_coverage','sampling_gap') if m[k]) or 'bracketed'
        lines.append(f"| {step['step']} | {step['duration_s']:.6f} | {m['max_observed_rss_bytes']} | {d['process.cpu_seconds']} | {d['process.major_faults']} | {d['system.swap_in_bytes']} / {d['system.swap_out_bytes']} | {flags} |")
    for step in result['slowest_steps']:
        lines += ['', f"Step {step['step']} phases: " + ', '.join(f"{p['phase']}={p['duration_s']:.6f}s" for p in step['phases'])]
    lines += ['', 'Missing values are null/None, not zero. See samples.jsonl for collector errors, unavailable metrics, compression, PSI and sample durations.']
    (directory/'summary.md').write_text('\n'.join(lines)+'\n')
    return result


def stop_process(process, sig):
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
    # The leader can exit before descendants; do not leave its group running.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run(command, directory, interval):
    directory.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    metadata = {'schema_version': 1, 'command': command, 'cwd': os.getcwd(),
                'platform': platform.platform(), 'interval_s': interval, 'started_unix_s': time.time(),
                'pid': None, 'exit_code': None, 'collector_error': None, 'interrupted': False}
    process = None
    reader = None
    errors = []
    status = 1
    try:
        with (directory/'run.log').open('wb') as raw, (directory/'events.jsonl').open('w') as events, (directory/'samples.jsonl').open('w') as metrics:
            try:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, start_new_session=True)
                metadata['pid'] = process.pid
                (directory/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
                reader = threading.Thread(target=read_output, args=(process.stdout, raw, events, start, errors), daemon=True)
                reader.start()
                sampler = Sampler(process.pid)
                previous = None
                deadline = start
                while process.poll() is None:
                    before = time.monotonic()-start
                    sample = sampler.sample()
                    after = time.monotonic()-start
                    sample.update(t=(before+after)/2, collection_duration_s=after-before)
                    sample['deltas'], sample['counter_resets'] = counter_deltas(previous, sample)
                    json_line(metrics, sample)
                    previous = sample
                    metadata['availability'] = {scope: {k: v is not None or metadata.get('availability', {}).get(scope, {}).get(k, False) for k,v in sample[scope].items()} for scope in ('process','system')}
                    if errors:
                        raise RuntimeError('Log reader: '+errors[0])
                    deadline = max(deadline+interval, time.monotonic())
                    try:
                        process.wait(timeout=max(0, deadline-time.monotonic()))
                    except subprocess.TimeoutExpired:
                        pass
                status = process.wait()
            except KeyboardInterrupt:
                metadata['interrupted'] = True
                if process:
                    stop_process(process, signal.SIGINT)
                status = 130
            except Exception as exc:
                metadata['collector_error'] = str(exc)
                if process:
                    stop_process(process, signal.SIGTERM)
                status = 1
            finally:
                if reader:
                    reader.join(timeout=2)
                    if reader.is_alive():
                        stop_process(process, signal.SIGTERM)
                        reader.join(timeout=2)
                    if reader.is_alive():
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        reader.join(timeout=2)
                    if reader.is_alive() or errors:
                        metadata['collector_error'] = 'Log reader did not finish' if reader.is_alive() else errors[0]
                        status = 1
                if process:
                    metadata['child_returncode'] = process.poll()
                    if not reader or not reader.is_alive():
                        process.stdout.close()
        summarize(directory, interval, time.monotonic()-start)
    except Exception as exc:
        metadata['collector_error'] = str(exc)
        status = 1
    finally:
        metadata['exit_code'] = status if status >= 0 else 128-status
        metadata['elapsed_s'] = time.monotonic()-start
        (directory/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata['exit_code']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--interval', type=float, default=1.)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not math.isfinite(args.interval) or args.interval <= 0:
        parser.error('Provide a command after -- and a finite positive interval')
    if args.out.exists():
        parser.error('--out must be a new directory')
    if platform.system() not in ('Darwin', 'Linux'):
        parser.error('Only macOS and Linux are supported')
    return run(command, args.out, args.interval)


if __name__ == '__main__':
    sys.exit(main())
