#!/usr/bin/env python3
"""Bounded matched extension using pinned Benchpress inputs and coupling maps.

This is a separate study, not a native Benchpress workout or its overall score.
Preparation, timing and independent verification are deliberately separate.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import select
import statistics
import subprocess
import sys
import time

import psutil


TOOLS = ['qrust', 'qiskit', 'cirq', 'tket']
SEED = 20261018
TRIALS = 3
SESSIONS = 2
TIMEOUT = 30
PER_WORKER_RSS_LIMIT = 2 * 1024 ** 3
ALL_WORKER_RSS_LIMIT = 3 * 1024 ** 3
PACKAGES = ['qiskit', 'cirq-core', 'pytket', 'numpy', 'scipy', 'networkx', 'psutil']


def digest(data):
    return hashlib.sha256(data.encode() if isinstance(data, str) else data).hexdigest()


def dump_new(path, obj):
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def versions():
    return {name: importlib.metadata.version(name) for name in PACKAGES}


def prepare(args):
    from comparison_study import THREAD_ENV
    args.output.mkdir(parents=True, exist_ok=False)
    cases = json.loads(args.cases.read_text())
    if len(cases) != 8 or len({case['id'] for case in cases}) != 8:
        raise ValueError('Expected eight prespecified unique pilot cases')
    for case in cases:
        if digest(case['qasm']) != case['qasm_sha256']:
            raise ValueError('Case checksum mismatch')
    files = [args.repo / 'experiments' / name for name in (
        'comparison_study.py', 'comparison_qiskit_tket.py', 'comparison_cirq.py', 'comparison_validation.py')]
    files += [Path(__file__).resolve(), Path(__file__).with_name('scalable_clifford.py')]
    schedules = {}
    for session in range(1, SESSIONS + 1):
        rng = random.Random(SEED + session * 1000)
        order = [case['id'] for case in cases]
        rng.shuffle(order)
        rows = []
        for case_id in order:
            for trial in range(-1, TRIALS):
                tools = list(TOOLS)
                rng.shuffle(tools)
                rows.extend({'case_id': case_id, 'tool': tool, 'trial': trial,
                             'seed': SEED + max(0, trial)} for tool in tools)
        schedules[str(session)] = rows
    protocol = {
        'schema': 1, 'study': 'matched Benchpress-derived Clifford pilot',
        'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'compiler_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip(),
        'binary_sha256': digest(args.binary.read_bytes()),
        'source_sha256': {str(path): digest(path.read_bytes()) for path in files},
        'cases_sha256': digest(args.cases.read_bytes()), 'cases': cases,
        'upstream_protocol_sha256': digest(args.cases.with_name('protocol.json').read_bytes()),
        'versions': versions(), 'python': sys.version, 'platform': platform.platform(),
        'machine': platform.machine(), 'memory_bytes': psutil.virtual_memory().total,
        'thread_environment': THREAD_ENV, 'tools': TOOLS, 'sessions': SESSIONS,
        'trials': TRIALS, 'warmups': 1, 'timeout_seconds': TIMEOUT,
        'worker_startup_timeout_seconds': 120,
        'per_worker_rss_limit_bytes': PER_WORKER_RSS_LIMIT,
        'all_workers_rss_limit_bytes': ALL_WORKER_RSS_LIMIT,
        'schedules': schedules,
        'target': 'Same unitary input, U3/CX, bidirectional square/heavy-hex edges, clean spare wires and free initial/final maps; no calibration optimization.',
        'timing': 'Frozen existing public adapters; warm in-process compilation, separately recorded parse and request-wall time. Workers persist. Three measured trials after warmup. No checker/profiling/native campaign concurrently.',
        'resource_observer': 'Controller samples known worker-process RSS approximately every 100 ms and before responses; child processes are not included. Single-process/thread settings are recorded. Sampled high-water observations are not exact peak allocation. Monitoring is outside the compiler timer but consumes some CPU. Limits can be overshot between samples. The aggregate limit includes idle resident workers and may be affected by request order.',
        'restart_policy': 'Failed initial warmup blocks that recipe/case for this session. After a timed worker failure, one replacement warmup is required before another timed request; its failure blocks remaining trials. Every scheduled trial retains a status.',
        'validation': 'Independent map-aware Clifford Choi-isometry comparison with exact signed-Pauli arithmetic after numerically bounded 1q Clifford recognition. All logical inputs including zero-initialized spare wires; non-Clifford output fails closed. Basis/edge checks separate. Verification only after both sessions.',
        'scope': 'Feasibility pilot of four selected Clifford workloads, not representative general-program performance, a holdout after future tuning, or native Benchpress scores. Unequal optimization effort remains. Ordinary OS load/CPU frequency uncontrolled.',
    }
    dump_new(args.output / 'protocol.json', protocol)
    print(json.dumps({'protocol_sha256': digest((args.output / 'protocol.json').read_bytes()),
                      'cases': len(cases), 'timed_requests': len(cases) * len(TOOLS) * TRIALS * SESSIONS}))


def verify_protocol(args):
    p = json.loads((args.output / 'protocol.json').read_text())
    if p['binary_sha256'] != digest(args.binary.read_bytes()) or p['versions'] != versions():
        raise ValueError('Frozen binary/package versions changed')
    for name, expected in p['source_sha256'].items():
        if digest(Path(name).read_bytes()) != expected:
            raise ValueError('Frozen source changed: ' + name)
    return p


def time_session(args):
    from comparison_study import Worker
    protocol = verify_protocol(args)
    rows_path = args.output / f'session_{args.session}.jsonl'
    if rows_path.exists():
        raise FileExistsError(rows_path)
    registry = []

    class ObservedWorker(Worker):
        def __init__(self, tool):
            self.peak_rss = self.peak_all_rss = 0
            registry.append(self)
            super().__init__(tool, args.binary, args.output)

        def receive(self, timeout):
            deadline = time.monotonic() + timeout
            while True:
                values = {}
                for item in registry:
                    process = getattr(item, 'process', None)
                    if process is None:
                        continue
                    try:
                        parent = psutil.Process(process.pid)
                        values[item.tool] = parent.memory_info().rss
                    except psutil.NoSuchProcess:
                        values[item.tool] = 0
                    except psutil.AccessDenied as error:
                        raise RuntimeError('Cannot enforce sampled worker RSS limit') from error
                own, total = values.get(self.tool, 0), sum(values.values())
                self.peak_rss, self.peak_all_rss = max(self.peak_rss, own), max(self.peak_all_rss, total)
                if own > PER_WORKER_RSS_LIMIT or total > ALL_WORKER_RSS_LIMIT:
                    raise MemoryError('Sampled worker RSS limit exceeded')
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(f'{self.tool} exceeded {timeout} seconds')
                if select.select([self.process.stdout], [], [], min(0.1, remaining))[0]:
                    line = self.process.stdout.readline()
                    if not line:
                        raise RuntimeError(f'Worker exited: {self.process.poll()}')
                    return json.loads(line)

        def measured_call(self, request):
            self.peak_rss = self.peak_all_rss = 0
            start = time.perf_counter_ns()
            result = super().call(request, TIMEOUT)
            if result.get('error', '').startswith('MemoryError:'):
                result['status'] = 'memory_limit'
            result.update(request_wall_ns=time.perf_counter_ns() - start,
                          sampled_worker_rss_bytes=self.peak_rss,
                          sampled_all_worker_rss_bytes=self.peak_all_rss)
            if 'out_qasm' in result:
                result['output_sha256'] = digest(result['out_qasm'])
            return result

    cases = {case['id']: case for case in protocol['cases']}
    workers, startup_errors, blocked = {}, {}, {}
    started = time.monotonic()
    try:
        for tool in TOOLS:
            try:
                workers[tool] = ObservedWorker(tool)
            except Exception as error:
                startup_errors[tool] = f'{type(error).__name__}: {error}'
                for item in registry:
                    if item.tool == tool:
                        item.close()
        with rows_path.open('x') as stream:
            def record(item):
                stream.write(json.dumps(item, separators=(',', ':'), allow_nan=False) + '\n')
                stream.flush()
            for index, scheduled in enumerate(protocol['schedules'][str(args.session)]):
                case, tool, trial = cases[scheduled['case_id']], scheduled['tool'], scheduled['trial']
                cell = (case['id'], tool)
                phase = 'warmup' if trial == -1 else 'timed'
                request = {key: case[key] for key in ('qasm', 'n_logical', 'physical_width', 'edges', 'topology')}
                request['seed'] = scheduled['seed']
                base = dict(scheduled, session=args.session, phase=phase, schedule_index=index)
                if tool in startup_errors:
                    result = {'status': 'startup_error', 'error': startup_errors[tool]}
                elif cell in blocked:
                    result = {'status': 'blocked_by_warmup', 'warmup_status': blocked[cell]}
                else:
                    if trial >= 0 and workers[tool].process is None:
                        replacement = workers[tool].measured_call(request)
                        record(dict(base, phase='replacement_warmup', **replacement))
                        if replacement['status'] != 'compiled':
                            blocked[cell] = replacement['status']
                    result = ({'status': 'blocked_by_warmup', 'warmup_status': blocked[cell]}
                              if cell in blocked else workers[tool].measured_call(request))
                if phase == 'warmup' and result['status'] != 'compiled':
                    blocked[cell] = result['status']
                record(dict(base, **result))
                if (index + 1) % 16 == 0:
                    print(json.dumps({'session': args.session, 'scheduled_done': index + 1,
                                      'last_case': case['id'], 'seconds': round(time.monotonic() - started, 2)}), flush=True)
    finally:
        for item in registry:
            item.close()
    dump_new(args.output / f'session_{args.session}_receipt.json', {
        'completed': True, 'session': args.session, 'wall_seconds': time.monotonic() - started,
        'protocol_sha256': digest((args.output / 'protocol.json').read_bytes()),
        'rows_sha256': digest(rows_path.read_bytes()), 'scheduled_requests': len(protocol['schedules'][str(args.session)])})


def validate(args):
    from comparison_validation import common_metrics, parse_qasm
    from scalable_clifford import validate_clifford
    import scalable_clifford
    if Path(scalable_clifford.__file__).resolve() != Path(__file__).with_name('scalable_clifford.py').resolve():
        raise ValueError('Imported checker differs from the frozen sibling source')
    protocol = verify_protocol(args)
    cases = {case['id']: case for case in protocol['cases']}
    cache, summaries, total = {}, {}, Counter()
    for session in range(1, SESSIONS + 1):
        receipt = json.loads((args.output / f'session_{session}_receipt.json').read_text())
        if not receipt['completed'] or receipt['protocol_sha256'] != digest((args.output / 'protocol.json').read_bytes()):
            raise ValueError('Incomplete/unbound timing session')
        rows_path = args.output / f'session_{session}.jsonl'
        if receipt['rows_sha256'] != digest(rows_path.read_bytes()):
            raise ValueError('Timing rows changed')
    for session in range(1, SESSIONS + 1):
        raw = [json.loads(line) for line in (args.output / f'session_{session}.jsonl').read_text().splitlines()]
        normal = [row for row in raw if row['phase'] != 'replacement_warmup']
        observed = [{k: row[k] for k in ('case_id', 'tool', 'trial', 'seed')} for row in normal]
        if observed != protocol['schedules'][str(session)]:
            raise ValueError('Scheduled outcomes missing or duplicated')
        with (args.output / f'validated_session_{session}.jsonl').open('x') as out:
            for row in raw:
                case = cases[row['case_id']]
                if row['status'] == 'compiled':
                    if row['output_sha256'] != digest(row['out_qasm']):
                        raise ValueError('Output hash mismatch')
                    key = digest(json.dumps([case['id'], row['output_sha256'], row['initial_layout'], row['final_layout']]))
                    if key not in cache:
                        check = validate_clifford(case['qasm'], row['out_qasm'], case['n_logical'], case['physical_width'],
                                                  row['initial_layout'], row['final_layout'], case['edges'])
                        try:
                            width, gates = parse_qasm(row['out_qasm'])
                            check['metrics'] = common_metrics(width, gates)
                        except (ValueError, TypeError, SyntaxError, ZeroDivisionError, OverflowError) as error:
                            check.update(passed=False, metrics_error=str(error))
                        cache[key] = check
                    row['validation'] = cache[key]
                    row['validation_cache_key'] = key
                status = ('passed' if row.get('validation', {}).get('passed') else
                          'unverified_output' if row['status'] == 'compiled' else row['status'])
                total[f"{row['phase']}:{row['tool']}:{status}"] += 1
                out.write(json.dumps(row, separators=(',', ':')) + '\n')
        cells = []
        for case_id in cases:
            for tool in TOOLS:
                measured = [row for row in raw if row['phase'] == 'timed' and row['case_id'] == case_id and row['tool'] == tool]
                warmups = [row for row in raw if row['phase'] != 'timed' and row['case_id'] == case_id and row['tool'] == tool]
                eligible = (len(measured) == TRIALS and len(warmups) >= 1
                            and all(row.get('validation', {}).get('passed') for row in [*measured, *warmups]))
                cell = {'case_id': case_id, 'tool': tool, 'eligible': eligible,
                        'status_counts': dict(Counter(row['status'] for row in measured)),
                        'validation_errors': sorted({row['validation'].get('error', '') for row in measured if 'validation' in row and not row['validation']['passed']})}
                if eligible:
                    fields = ['compile_ns', 'parse_ns', 'request_wall_ns', 'sampled_worker_rss_bytes']
                    metrics = {name: [row[name] for row in measured] for name in fields}
                    for name in ['cx_count', 'gate_count', 'depth', 'two_q_depth']:
                        metrics[name] = [row['validation']['metrics'][name] for row in measured]
                    cell['metrics'] = {name: {'median': statistics.median(values), 'min': min(values), 'max': max(values)}
                                       for name, values in metrics.items()}
                cells.append(cell)
        summaries[str(session)] = cells
    summary = {'schema': 1, 'protocol_sha256': digest((args.output / 'protocol.json').read_bytes()),
               'counts': dict(total), 'unique_artifacts': len(cache), 'sessions': summaries,
               'validation_limits': protocol['validation'], 'scope': protocol['scope']}
    dump_new(args.output / 'summary.json', summary)
    print(json.dumps({'counts': dict(total), 'unique_artifacts': len(cache)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'time', 'validate'])
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', type=Path)
    parser.add_argument('--session', type=int, choices=[1, 2])
    args = parser.parse_args()
    args.repo, args.binary, args.output = args.repo.resolve(), args.binary.resolve(), args.output.resolve()
    sys.path.insert(0, str(args.repo / 'experiments'))
    if args.command == 'prepare':
        if args.cases is None:
            parser.error('--cases required')
        prepare(args)
    elif args.command == 'time':
        if args.session is None:
            parser.error('--session required')
        time_session(args)
    else:
        validate(args)


if __name__ == '__main__':
    main()
