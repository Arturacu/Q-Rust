#!/usr/bin/env python3
"""Preserve existing schedules and insert balanced contemporaneous old-Q-Rust controls.

Public adapters and independent oracles remain unchanged. The primary result
files retain the existing external/matched schemas and can be checked by those
harnesses. Separate control files and a checksum-bound dispatch log retain the
extra requests. No verification runs in the timing process.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import select
import subprocess
import sys
import time
from types import SimpleNamespace

def sha(value):
    if isinstance(value, Path):
        value = value.read_bytes()
    if isinstance(value, str):
        value = value.encode()
    return hashlib.sha256(value).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False)


def dump_new(path, obj):
    with path.open('x') as stream:
        stream.write(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n')


def modules(repo):
    sys.path.insert(0, str(repo / 'experiments'))
    return importlib.import_module('external_comparison'), importlib.import_module('benchpress_matched')


def baseline_provenance(external_baseline, matched_baseline, old_binary):
    """Bind an immediate old control to both supplied baseline protocols.

    Historical compiler identity comes from the timed protocols, not the
    current repository HEAD (which can include later reports or harness work).
    """
    paths = {'external': Path(external_baseline) / 'protocol.json',
             'matched': Path(matched_baseline) / 'protocol.json'}
    protocols = {suite: json.loads(path.read_text()) for suite, path in paths.items()}
    checksum_path = Path(external_baseline) / 'protocol.sha256'
    if sha(paths['external']) != checksum_path.read_text().strip():
        raise ValueError('External baseline protocol checksum changed')
    actual = sha(Path(old_binary))
    hashes = [protocols['external']['qrust_binary_sha256'], protocols['matched']['binary_sha256']]
    if any(expected != actual for expected in hashes):
        raise ValueError('Old binary does not match both supplied baseline protocols')
    commits = [protocols[suite]['compiler_commit'] for suite in ('external', 'matched')]
    if not commits[0] or commits[0] != commits[1]:
        raise ValueError('Baseline compiler commits disagree')
    return {'old_binary_sha256': actual, 'old_compiler_commit': commits[0],
            'baseline_protocol_sha256': {suite: sha(path) for suite, path in paths.items()},
            'old_compiler_provenance': protocols['external']['compiler_provenance']}


def expanded_schedule(schedule, session):
    """Exact original projection, with old/new order balanced across sessions."""
    expanded = []
    for index, item in enumerate(schedule):
        primary = dict(item, primary_schedule_index=index, variant='primary')
        if item['tool'] != 'qrust':
            expanded.append(primary)
            continue
        # The hash depends only on case/trial, not the shuffled session order.
        # Flipping parity by session yields one old-first and one new-first
        # observation for every corresponding seed and warmup across sessions.
        parity = int(sha(item['case_id'] + ':' + str(item['trial']))[:8], 16)
        control = dict(primary, variant='control')
        expanded.extend([control, primary] if (parity + session) % 2 == 0 else [primary, control])
    return expanded


def prepare(args):
    ext, matched = modules(args.repo)
    baseline = baseline_provenance(args.external_baseline, args.matched_baseline, args.old_binary)
    if args.binary == args.old_binary or sha(args.binary) == baseline['old_binary_sha256']:
        raise ValueError('Candidate and old binary must differ')
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip()
    # The source state must be reviewable and committed before timing.
    dirty = subprocess.check_output(['git', 'diff', 'HEAD', '--', '*.rs', '**/Cargo.toml', '**/Cargo.lock', 'Cargo.toml', 'Cargo.lock'], cwd=args.repo)
    untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=args.repo).decode().split('\0')
    if dirty or any(name.endswith('.rs') or Path(name).name in {'Cargo.toml', 'Cargo.lock'} for name in untracked if name):
        raise ValueError('Compiler source must be committed before freezing')
    args.output.mkdir(parents=True, exist_ok=False)
    ext.prepare(SimpleNamespace(corpus=args.corpus, output=args.output / 'external',
        binary=args.binary, compiler_commit=head, trials=7, sessions=2, timeout=60.0,
        seed=20261018, machine_note='Sequential paired followup. Same primary schedules; old Q-Rust runs adjacent to new Q-Rust with balanced order. No concurrent builds, simulation or benchmark sessions.'))
    matched.prepare(SimpleNamespace(repo=args.repo, binary=args.binary,
        output=args.output / 'matched', cases=args.cases))
    for suite, old_path in [('external', args.external_baseline), ('matched', args.matched_baseline)]:
        old = json.loads((old_path / 'protocol.json').read_text())
        current = json.loads((args.output / suite / 'protocol.json').read_text())
        for field in ('cases', 'schedules', 'versions', 'thread_environment', 'trials', 'sessions', 'tools', 'timeout_seconds'):
            if current[field] != old[field]:
                raise ValueError(f'{suite} changed baseline contract field {field}')
        # Every existing adapter/oracle/scheduler hash must still match.
        source_hashes = {(Path(filename).name if Path(filename).is_absolute() else filename): checksum
                         for filename, checksum in current['source_sha256'].items()}
        for filename, expected in old['source_sha256'].items():
            key = Path(filename).name if Path(filename).is_absolute() else filename
            if source_hashes.get(key) != expected:
                raise ValueError('Existing harness/adapter/oracle changed: ' + filename)
    protocols = {suite: json.loads((args.output / suite / 'protocol.json').read_text())
                 for suite in ('external', 'matched')}
    protocol = {
        'schema': 1, 'created_utc': ext.utc_now(), 'compiler_commit': head,
        'candidate_binary': str(args.binary), 'candidate_binary_sha256': sha(args.binary),
        'old_binary': str(args.old_binary), **baseline,
        'repo': str(args.repo), 'controller_sha256': sha(Path(__file__).resolve()),
        'analysis_source_sha256': sha(Path(__file__).with_name('paired_analysis.py')),
        'primary_protocol_sha256': {suite: sha(args.output / suite / 'protocol.json') for suite in protocols},
        'primary_scheduled_counts': {'external': {'timed': 8400, 'warmup': 1200}, 'matched': {'timed': 192, 'warmup': 64}},
        'control_scheduled_counts': {'external': {'timed': 1400, 'warmup': 200}, 'matched': {'timed': 48, 'warmup': 16}},
        'schedules': {suite: {str(session): expanded_schedule(p['schedules'][str(session)], session)
                             for session in (1, 2)} for suite, p in protocols.items()},
        'pairing': 'For each original Q-Rust slot, old and candidate requests run consecutively; SHA256(case_id:trial) parity chooses order and the second session reverses it. Original primary schedules remain exact projections. Repeated seeds describe timing and seed sensitivity, not additional independent circuits.',
        'timing': 'Each request retains its own warm in-process compile_ns and complete request_wall_ns. No old request runs inside a candidate timer. Startup, other requests, output writes and verification are excluded from compiler timers. Workers persist, all SDK imports precede requests, and requests are sequential.',
        'resources': 'External retains its existing 60-second request policy. Matched retains sampled per-process 2-GiB and aggregate 3-GiB caps approximately every 100 ms and 30-second request timeouts. Matched aggregate RSS now includes the additional resident old-Q-Rust worker; sampled RSS is not exact peak memory and does not include descendants. Limits can overshoot between samples.',
        'restart_policy': 'External uses its existing replacement warmup after any worker loss; matched uses its existing failed-initial-warmup block and replacement-warmup policy. Old control failure does not block the candidate or alter SDK eligibility. Failed/unsupported outcomes remain recorded.',
        'eligibility': 'Every paired per-case/session comparison requires all scheduled old and new timed trials plus every initial/replacement warmup to pass the original independent semantic checker. No successful subset replaces an ineligible cell. Matched unsupported non-Clifford output remains unverified, not incorrect.',
        'native_benchpress': 'Existing native workout evidence remains unchanged. No native Q-Rust gym exists and no native workout is rerun in this followup.',
    }
    dump_new(args.output / 'paired_protocol.json', protocol)
    (args.output / 'paired_protocol.sha256').write_text(sha(args.output / 'paired_protocol.json') + '\n')
    print(json.dumps({'paired_protocol_sha256': sha(args.output / 'paired_protocol.json'), 'compiler_commit': head,
                      'scheduled_timed': 10040, 'scheduled_warmup': 1480}), flush=True)


def verify(args):
    ext, matched = modules(args.repo)
    path = args.output / 'paired_protocol.json'
    if sha(path) != (args.output / 'paired_protocol.sha256').read_text().strip():
        raise ValueError('Paired protocol changed')
    pair = json.loads(path.read_text())
    for path, expected in [(args.binary, pair['candidate_binary_sha256']),
                           (args.old_binary, pair['old_binary_sha256']),
                           (Path(__file__).resolve(), pair['controller_sha256']),
                           (Path(__file__).with_name('paired_analysis.py'), pair['analysis_source_sha256'])]:
        if sha(path) != expected:
            raise ValueError('Frozen binary/controller changed: ' + str(path))
    protocols = {}
    for suite in ('external', 'matched'):
        directory = args.output / suite
        if sha(directory / 'protocol.json') != pair['primary_protocol_sha256'][suite]:
            raise ValueError('Primary protocol changed')
        if suite == 'external':
            p, _ = ext.load_protocol(directory)
            ext.check_runtime(p, args.binary)
        else:
            p = matched.verify_protocol(SimpleNamespace(output=directory, binary=args.binary))
        for session in (1, 2):
            if pair['schedules'][suite][str(session)] != expanded_schedule(p['schedules'][str(session)], session):
                raise ValueError('Expanded schedule changed')
        protocols[suite] = p
    return pair, protocols


def time_session(args):
    ext, matched = modules(args.repo)
    pair, protocols = verify(args)
    protocol, suite, session = protocols[args.suite], args.suite, args.session
    output = args.output / suite
    lock = args.output / 'active.lock'
    descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(descriptor)
    registry, workers, startup_errors, blocked, handles = [], {}, {}, {}, []
    started = time.monotonic()
    completed = False
    dispatch_count = 0
    primary_count = control_count = 0
    control_dir = output / f'control_session_{session}'
    primary_dir = output / f'session_{session:02d}' if suite == 'external' else output

    class ObservedWorker(ext.ProfileWorker):
        def __init__(self, recipe, binary, directory, label):
            self.label, self.peak_rss, self.peak_all_rss = label, 0, 0
            registry.append(self)
            super().__init__(recipe, binary, directory)

        def receive(self, timeout):
            if suite == 'external':
                return super().receive(timeout)
            import psutil
            deadline = time.monotonic() + timeout
            while True:
                values = {}
                for worker in registry:
                    process = getattr(worker, 'process', None)
                    if process is None:
                        continue
                    try:
                        values[worker.label] = psutil.Process(process.pid).memory_info().rss
                    except psutil.NoSuchProcess:
                        values[worker.label] = 0
                    except psutil.AccessDenied as error:
                        raise RuntimeError('Cannot enforce sampled worker RSS limit') from error
                own, total = values.get(self.label, 0), sum(values.values())
                self.peak_rss, self.peak_all_rss = max(self.peak_rss, own), max(self.peak_all_rss, total)
                if own > matched.PER_WORKER_RSS_LIMIT or total > matched.ALL_WORKER_RSS_LIMIT:
                    raise MemoryError('Sampled worker RSS limit exceeded')
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f'{self.label} exceeded {timeout} seconds')
                if select.select([self.process.stdout], [], [], min(0.1, remaining))[0]:
                    line = self.process.stdout.readline()
                    if not line:
                        raise RuntimeError(f'Worker exited: {self.process.poll()}')
                    return json.loads(line)

    def open_rows(path):
        stream = path.open('x')
        handles.append(stream)
        return stream

    def measured(worker, request):
        worker.peak_rss = worker.peak_all_rss = 0
        result = ext.timed_request(worker, request, protocol['timeout_seconds'])
        if suite == 'matched':
            if result.get('error', '').startswith('MemoryError:'):
                result['status'] = 'memory_limit'
            result.update(sampled_worker_rss_bytes=worker.peak_rss, sampled_all_worker_rss_bytes=worker.peak_all_rss)
        if result.get('status') == 'compiled':
            result['output_sha256'] = sha(result['out_qasm'])
        return result

    try:
        control_dir.mkdir(exist_ok=False)
        if suite == 'external':
            primary_dir.mkdir(exist_ok=False)
            (primary_dir / 'protocol.json').write_bytes((output / 'protocol.json').read_bytes())
            info = {'session': session, 'protocol_sha256': sha(output / 'protocol.json'),
                    'started_utc': ext.utc_now(), 'finished_utc': None, 'completed': False}
            ext.dump(primary_dir / 'session.json', info)
            streams = {'primary': {kind: open_rows(primary_dir / f'{kind}.jsonl') for kind in ('timed', 'warmup')},
                       'control': {kind: open_rows(control_dir / f'{kind}.jsonl') for kind in ('timed', 'warmup')}}
        else:
            streams = {'primary': {'all': open_rows(output / f'session_{session}.jsonl')},
                       'control': {'all': open_rows(control_dir / 'rows.jsonl')}}
        dispatch_path = output / f'dispatch_session_{session}.jsonl'
        dispatch_stream = open_rows(dispatch_path)
        for recipe in protocol['tools']:
            label = 'primary:' + recipe
            try:
                workers[label] = ObservedWorker(recipe, args.binary, primary_dir, label)
            except Exception as error:
                startup_errors[label] = f'{type(error).__name__}: {error}'
        try:
            workers['control:qrust'] = ObservedWorker('qrust', args.old_binary, control_dir, 'control:qrust')
        except Exception as error:
            startup_errors['control:qrust'] = f'{type(error).__name__}: {error}'
        if suite == 'external':
            ext.dump(primary_dir / 'startup.json', {'errors': {key.split(':')[1]: value for key, value in startup_errors.items() if key.startswith('primary:')}})
        dump_new(control_dir / 'startup.json', {'errors': startup_errors})
        cases = {case['id']: case for case in protocol['cases']}

        def record(row, variant, phase, expanded_index):
            nonlocal dispatch_count
            kind = 'all' if suite == 'matched' else ('timed' if phase == 'timed' else 'warmup')
            stream = streams[variant][kind]
            line = canonical(row)
            stream.write(line + '\n')
            stream.flush()
            dispatch_stream.write(canonical({'dispatch_index': dispatch_count,
                'expanded_schedule_index': expanded_index, 'variant': variant,
                'primary_schedule_index': expanded[expanded_index]['primary_schedule_index'],
                'phase': phase, 'case_id': row['case_id'], 'tool': row['tool'], 'trial': row['trial'],
                'session': session, 'status': row['status'], 'record_sha256': sha(line)}) + '\n')
            dispatch_stream.flush()
            dispatch_count += 1

        expanded = pair['schedules'][suite][str(session)]
        for expanded_index, planned in enumerate(expanded):
            variant = planned['variant']
            index = planned['primary_schedule_index']
            entry = protocol['schedules'][str(session)][index]
            case, recipe, trial = cases[entry['case_id']], entry['tool'], entry['trial']
            phase = 'warmup' if trial == -1 else 'timed'
            label, cell = variant + ':' + recipe, (variant, case['id'], recipe)
            request = {key: case[key] for key in ('qasm', 'n_logical', 'physical_width', 'edges', 'topology')}
            request['seed'] = entry['seed']
            base = dict(entry, session=session)
            if suite == 'matched':
                base.update(phase=phase, schedule_index=index)
            worker = workers.get(label)
            if worker is None:
                result = {'status': 'startup_error' if suite == 'matched' else 'error', 'error': startup_errors[label],
                          'request_wall_ns': None}
                if suite == 'external':
                    result['phase'] = 'worker_startup'
            elif suite == 'matched' and cell in blocked:
                result = {'status': 'blocked_by_warmup', 'warmup_status': blocked[cell]}
            else:
                replacement_failed = False
                if trial >= 0 and worker.process is None:
                    warmed = measured(worker, request)
                    replacement = dict(base, **warmed)
                    if suite == 'external':
                        replacement['reason'] = 'worker_restart'
                    else:
                        replacement['phase'] = 'replacement_warmup'
                    record(replacement, variant, 'replacement_warmup', expanded_index)
                    if warmed['status'] != 'compiled':
                        replacement_failed = True
                        if suite == 'matched':
                            blocked[cell] = warmed['status']
                            result = {'status': 'blocked_by_warmup', 'warmup_status': blocked[cell]}
                        else:
                            result = dict(warmed, phase='replacement_warmup',
                                error='Replacement warmup failed; measured request not attempted: ' + warmed.get('error', 'unknown failure'))
                if not replacement_failed:
                    result = measured(worker, request)
            if suite == 'matched' and phase == 'warmup' and result['status'] != 'compiled':
                blocked[cell] = result['status']
            row = {**result, **base}
            if suite == 'external':
                row['id'] = f"{case['id']}__{recipe}__{trial}"
                if row.get('status') == 'compiled':
                    try:
                        row['artifact_sha256'] = ext.artifact_hash(case, row)
                    except (KeyError, TypeError) as error:
                        row.update(status='error', error=f'Incomplete worker output: {error}')
            record(row, variant, phase, expanded_index)
            if variant == 'primary':
                primary_count += 1
            else:
                control_count += 1
            if (expanded_index + 1) % (56 if suite == 'external' else 20) == 0:
                print(json.dumps({'suite': suite, 'session': session, 'scheduled_done': expanded_index + 1,
                    'scheduled_total': len(expanded), 'last_case': case['id'],
                    'wall_seconds': round(time.monotonic() - started, 2)}), flush=True)
        completed = True
    finally:
        for worker in registry:
            worker.close()
        for stream in handles:
            stream.close()
        lock.unlink()
    if not completed:
        raise RuntimeError('Incomplete session')
    if suite == 'external':
        info.update(finished_utc=ext.utc_now(), completed=True,
            timed_rows=len(protocol['cases']) * len(protocol['tools']) * protocol['trials'],
            timed_sha256=sha(primary_dir / 'timed.jsonl'), warmup_sha256=sha(primary_dir / 'warmup.jsonl'))
        ext.dump(primary_dir / 'session.json', info)
    else:
        dump_new(output / f'session_{session}_receipt.json', {
            'completed': True, 'session': session, 'wall_seconds': time.monotonic() - started,
            'protocol_sha256': sha(output / 'protocol.json'),
            'rows_sha256': sha(output / f'session_{session}.jsonl'), 'scheduled_requests': primary_count})
    bound_paths = {str(Path(stream.name).relative_to(output)): sha(Path(stream.name)) for stream in handles}
    dump_new(output / f'paired_session_{session}_receipt.json', {
        'completed': True, 'suite': suite, 'session': session,
        'paired_protocol_sha256': sha(args.output / 'paired_protocol.json'),
        'primary_protocol_sha256': sha(output / 'protocol.json'),
        'primary_scheduled': primary_count, 'control_scheduled': control_count,
        'dispatch_records': dispatch_count, 'files_sha256': bound_paths,
        'wall_seconds': time.monotonic() - started})
    print(json.dumps({'suite': suite, 'session': session, 'completed': True,
                      'wall_seconds': round(time.monotonic() - started, 2)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'verify', 'time'])
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--old-binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--corpus', type=Path)
    parser.add_argument('--cases', type=Path)
    parser.add_argument('--external-baseline', type=Path)
    parser.add_argument('--matched-baseline', type=Path)
    parser.add_argument('--suite', choices=['external', 'matched'])
    parser.add_argument('--session', type=int, choices=[1, 2])
    args = parser.parse_args()
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.resolve())
    if args.command == 'prepare' and any(getattr(args, name) is None for name in
            ('corpus', 'cases', 'external_baseline', 'matched_baseline')):
        parser.error('prepare requires corpus, cases and both baseline directories')
    if args.command == 'time' and (args.suite is None or args.session is None):
        parser.error('time requires suite and session')
    if args.command == 'prepare':
        prepare(args)
    elif args.command == 'verify':
        verify(args)
        print('Frozen contract verified')
    else:
        time_session(args)


if __name__ == '__main__':
    main()
