#!/usr/bin/env python3
"""Replay archived matched-pilot validation without a compiler binary or SDKs.

The archive keeps the original absolute-path protocol unchanged. Resolve its
source hashes by unique filename in the explicitly supplied archived harness,
then bind every output to the original session receipt and planned schedule.
Unverified outputs remain unverified; they are not discarded or called wrong.
"""
from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(evidence, harness):
    from comparison_validation import common_metrics, parse_qasm
    from scalable_clifford import validate_clifford
    import scalable_clifford
    import comparison_validation
    if Path(scalable_clifford.__file__).resolve() != (harness / 'experiments/scalable_clifford.py').resolve():
        raise ValueError('Checker import differs from the supplied archived harness')
    if Path(comparison_validation.__file__).resolve() != (harness / 'experiments/comparison_validation.py').resolve():
        raise ValueError('Analytic oracle import differs from the supplied archived harness')
    directory = evidence / 'matched'
    protocol_path = directory / 'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    names = [Path(name).name for name in protocol['source_sha256']]
    if len(names) != len(set(names)):
        raise ValueError('Ambiguous archived source filenames')
    for name, expected in protocol['source_sha256'].items():
        if sha(harness / 'experiments' / Path(name).name) != expected:
            raise ValueError('Archived harness source changed: ' + name)
    cases = {case['id']: case for case in protocol['cases']}
    if len(cases) != len(protocol['cases']):
        raise ValueError('Repeated case identifier')
    for case in cases.values():
        if hashlib.sha256(case['qasm'].encode()).hexdigest() != case['qasm_sha256']:
            raise ValueError('Input checksum differs')
    summary = json.loads((directory / 'summary.json').read_text())
    counts, cache, rows_checked = Counter(), {}, 0
    if summary['protocol_sha256'] != sha(protocol_path):
        raise ValueError('Summary protocol binding differs')
    if set(summary['sessions']) != {str(s) for s in range(1, protocol['sessions'] + 1)}:
        raise ValueError('Summary session coverage differs')
    for session in range(1, protocol['sessions'] + 1):
        receipt = json.loads((directory / f'session_{session}_receipt.json').read_text())
        raw_path = directory / f'session_{session}.jsonl'
        if (not receipt['completed'] or receipt['protocol_sha256'] != sha(protocol_path)
                or receipt['rows_sha256'] != sha(raw_path)):
            raise ValueError('Session receipt differs')
        raw = [json.loads(line) for line in raw_path.read_text().splitlines()]
        saved = [json.loads(line) for line in (directory / f'validated_session_{session}.jsonl').read_text().splitlines()]
        if len(raw) != len(saved):
            raise ValueError('Validation rows missing or duplicated')
        observed = [{k: row[k] for k in ('case_id', 'tool', 'trial', 'seed')}
                    for row in raw if row['phase'] != 'replacement_warmup']
        if observed != protocol['schedules'][str(session)]:
            raise ValueError('Schedule differs')
        for row in raw:
            index = row['schedule_index']
            if type(index) is not int or not 0 <= index < len(observed):
                raise ValueError('Invalid schedule index')
            expected = protocol['schedules'][str(session)][index]
            if row['session'] != session or any(row[k] != v for k, v in expected.items()):
                raise ValueError('Schedule metadata differs')
            if row['phase'] == 'replacement_warmup':
                if row['trial'] < 0:
                    raise ValueError('Replacement warmup has invalid trial')
            elif row['phase'] != ('warmup' if row['trial'] == -1 else 'timed'):
                raise ValueError('Phase differs from schedule')
        if receipt['scheduled_requests'] != len(observed):
            raise ValueError('Scheduled request count differs')
        for row, validated in zip(raw, saved):
            original = {k: v for k, v in validated.items()
                        if k not in ('validation', 'validation_cache_key')}
            if row != original:
                raise ValueError('Validated artifact differs from raw outcome')
            if row['status'] == 'compiled':
                case = cases[row['case_id']]
                if hashlib.sha256(row['out_qasm'].encode()).hexdigest() != row['output_sha256']:
                    raise ValueError('Output checksum differs')
                key = hashlib.sha256(json.dumps([case['id'], row['output_sha256'], row['initial_layout'], row['final_layout']]).encode()).hexdigest()
                if key not in cache:
                    check = validate_clifford(case['qasm'], row['out_qasm'], case['n_logical'], case['physical_width'], row['initial_layout'], row['final_layout'], case['edges'])
                    try:
                        width, gates = parse_qasm(row['out_qasm'])
                        check['metrics'] = common_metrics(width, gates)
                    except (ValueError, TypeError, SyntaxError, ZeroDivisionError, OverflowError) as error:
                        check.update(passed=False, metrics_error=str(error))
                    cache[key] = check
                if validated.get('validation_cache_key') != key or validated.get('validation') != cache[key]:
                    raise ValueError('Independent validation replay differs')
            elif 'validation' in validated or 'validation_cache_key' in validated:
                raise ValueError('Non-compiled outcome carries fabricated validation')
            status = ('passed' if validated.get('validation', {}).get('passed') else
                      'unverified_output' if row['status'] == 'compiled' else row['status'])
            counts[f"{row['phase']}:{row['tool']}:{status}"] += 1
            rows_checked += 1
        cells = []
        for case_id in cases:
            for tool in protocol['tools']:
                measured = [row for row in saved if row['phase'] == 'timed' and row['case_id'] == case_id and row['tool'] == tool]
                warmups = [row for row in saved if row['phase'] != 'timed' and row['case_id'] == case_id and row['tool'] == tool]
                eligible = (len(measured) == protocol['trials'] and len(warmups) >= 1
                            and all(row.get('validation', {}).get('passed') for row in [*measured, *warmups]))
                cell = {'case_id': case_id, 'tool': tool, 'eligible': eligible,
                        'status_counts': dict(Counter(row['status'] for row in measured)),
                        'validation_errors': sorted({row['validation'].get('error', '') for row in measured if 'validation' in row and not row['validation']['passed']})}
                if eligible:
                    metrics = {name: [row[name] for row in measured] for name in
                               ('compile_ns', 'parse_ns', 'request_wall_ns', 'sampled_worker_rss_bytes')}
                    for name in ('cx_count', 'gate_count', 'depth', 'two_q_depth'):
                        metrics[name] = [row['validation']['metrics'][name] for row in measured]
                    cell['metrics'] = {name: {'median': statistics.median(values), 'min': min(values), 'max': max(values)}
                                       for name, values in metrics.items()}
                cells.append(cell)
        if cells != summary['sessions'][str(session)]:
            raise ValueError('Summary cell aggregates differ')
    if dict(counts) != summary['counts'] or len(cache) != summary['unique_artifacts']:
        raise ValueError('Summary coverage differs')
    if summary['scope'] != protocol['scope'] or summary['validation_limits'] != protocol['validation']:
        raise ValueError('Summary scope or validation limits differ')
    return {'passed': True, 'rows_checked': rows_checked, 'unique_artifacts': len(cache),
            'counts': dict(counts), 'protocol_sha256': sha(protocol_path),
            'scope': 'Receipt, schedule, input/output hashes, complete validation, common metrics, eligibility and all per-cell aggregates replay; no compilation or timing rerun.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    harness = args.harness.resolve()
    sys.path.insert(0, str(harness / 'experiments'))
    result = audit(args.evidence.resolve(), harness)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
