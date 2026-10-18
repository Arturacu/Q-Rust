#!/usr/bin/env python3
"""Check full dispatch coverage, independently validate outputs, and compare pairs.

Runs only after all four timing sessions. Numerical tests are the frozen
external oracle and the frozen map-aware Clifford oracle; their limitations
are preserved. Every warmup and timed output is checked, including controls.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import statistics
from types import SimpleNamespace

import paired_campaign as campaign


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def phase_of(row, suite):
    if suite == 'matched':
        return row['phase']
    return ('replacement_warmup' if row.get('reason') == 'worker_restart' else
            'warmup' if row['trial'] == -1 else 'timed')


def classify_evidence(result, suite):
    """Interpret checker evidence separately from legacy outcome labels.

    The original matched harness calls every non-passing compiled output
    unverified. Preserve that label for comparability, but never conceal a
    recognized Clifford counterexample or an explicit target-contract failure.
    """
    if result.get('passed') is True:
        return 'verified'
    error = result.get('error', '')
    if result.get('output_width_error') or result.get('metrics_error'):
        return 'structural_contract_failure'
    if suite == 'external':
        if result.get('evidence_class') in ('full_isometry', 'sampled_isometry'):
            return 'semantic_counterexample'
        if result.get('evidence_class') == 'unchecked_width':
            return 'unsupported_width'
        if result.get('evidence_class') == 'invalid':
            return 'structural_contract_failure'
        return 'unclassified_checker_failure'
    if result.get('evidence_class') == 'clifford_choi_isometry' and result.get('tableau_equal') is False:
        return 'semantic_counterexample'
    if error == 'Recognition residual budget exceeded':
        return 'recognition_inconclusive'
    if error.startswith(('Non-Clifford local gate:', 'Unsupported/non-Clifford gate:',
                         'Unsupported symbolic Clifford:')):
        return 'unsupported_recognition'
    if error.startswith(('Invalid widths', 'Invalid logical-to-physical map', 'QASM width mismatch',
                         'Invalid coupling edge', 'Output outside declared basis:', 'Output violates coupling graph')):
        return 'structural_contract_failure'
    return 'unclassified_checker_failure'


def check_dispatch(pair, protocol, directory, suite, session):
    """Bind every raw row to exactly one scheduled or replacement dispatch."""
    receipt = json.loads((directory / f'paired_session_{session}_receipt.json').read_text())
    require(receipt['completed'] and receipt['session'] == session and receipt['suite'] == suite, 'Incomplete or wrong session')
    require(receipt['paired_protocol_sha256'] == campaign.sha(directory.parent / 'paired_protocol.json'), 'Unbound paired session')
    require(receipt['primary_protocol_sha256'] == campaign.sha(directory / 'protocol.json'), 'Unbound primary session')
    for name, expected in receipt['files_sha256'].items():
        path = directory / name
        require(path.resolve().is_relative_to(directory.resolve()), 'Unsafe receipt path')
        require(campaign.sha(path) == expected, 'Changed recorded file: ' + name)
    if suite == 'external':
        primary_dir = directory / f'session_{session:02d}'
        control_dir = directory / f'control_session_{session}'
        recorded = {variant: [row for kind in ('timed', 'warmup') for row in rows(folder / f'{kind}.jsonl')]
                    for variant, folder in [('primary', primary_dir), ('control', control_dir)]}
        info = json.loads((primary_dir / 'session.json').read_text())
        require(info['completed'] and info['protocol_sha256'] == campaign.sha(directory / 'protocol.json'), 'Primary external receipt mismatch')
        for kind in ('timed', 'warmup'):
            require(info[f'{kind}_sha256'] == campaign.sha(primary_dir / f'{kind}.jsonl'), 'Primary external row hash mismatch')
    else:
        recorded = {'primary': rows(directory / f'session_{session}.jsonl'),
                    'control': rows(directory / f'control_session_{session}/rows.jsonl')}
        info = json.loads((directory / f'session_{session}_receipt.json').read_text())
        require(info['completed'] and info['protocol_sha256'] == campaign.sha(directory / 'protocol.json'), 'Primary matched receipt mismatch')
        require(info['rows_sha256'] == campaign.sha(directory / f'session_{session}.jsonl'), 'Primary matched row hash mismatch')
    indexed = {}
    expected_hashes = Counter()
    for variant, observations in recorded.items():
        for row in observations:
            identity = variant, campaign.sha(campaign.canonical(row))
            indexed[identity] = row
            expected_hashes[identity] += 1
    dispatch = rows(directory / f'dispatch_session_{session}.jsonl')
    require(receipt['dispatch_records'] == len(dispatch), 'Dispatch count mismatch')
    planned = pair['schedules'][suite][str(session)]
    scheduled_index, seen, ordered = 0, Counter(), []
    for index, record in enumerate(dispatch):
        require(record['dispatch_index'] == index and record['session'] == session, 'Dispatch identity mismatch')
        require(scheduled_index < len(planned), 'Extra dispatch record')
        slot = planned[scheduled_index]
        require(record['expanded_schedule_index'] == scheduled_index, 'Expanded dispatch order mismatch')
        for field in ('variant', 'primary_schedule_index', 'case_id', 'tool', 'trial'):
            require(record[field] == slot[field], 'Dispatch slot mismatch: ' + field)
        key = record['variant'], record['record_sha256']
        require(key in indexed, 'Dispatch references missing raw record')
        seen[key] += 1
        row = indexed[key]
        require(all(row[field] == slot[field] for field in ('case_id', 'tool', 'trial', 'seed')), 'Raw row differs from scheduled request')
        require(row['session'] == session and phase_of(row, suite) == record['phase'] and row['status'] == record['status'], 'Raw phase/status mismatch')
        if record['phase'] == 'replacement_warmup':
            require(slot['trial'] >= 0, 'Replacement before initial warmup')
            require(index + 1 < len(dispatch) and dispatch[index + 1]['expanded_schedule_index'] == scheduled_index
                    and dispatch[index + 1]['phase'] == 'timed', 'Replacement must immediately precede its timed status')
        else:
            require(record['phase'] == ('warmup' if slot['trial'] == -1 else 'timed'), 'Scheduled phase mismatch')
            scheduled_index += 1
        ordered.append((record['variant'], row))
    require(scheduled_index == len(planned) and seen == expected_hashes, 'Missing/extra/duplicated dispatch or raw record')
    for variant in ('primary', 'control'):
        count = sum(record['variant'] == variant and record['phase'] != 'replacement_warmup' for record in dispatch)
        require(count == receipt[variant + '_scheduled'], 'Scheduled receipt count mismatch')
    return ordered


def checked_artifact(case, row, suite, cache, validate_output, validate_clifford, parse_qasm, common_metrics):
    if row['status'] != 'compiled':
        return row['status'], None, None
    require(campaign.sha(row['out_qasm']) == row['output_sha256'], 'Output QASM hash mismatch')
    bound = {'case': case, 'output_sha256': row['output_sha256'], 'initial_layout': row['initial_layout'],
             'final_layout': row['final_layout'], 'output_qubits': row['output_qubits']}
    key = campaign.sha(campaign.canonical(bound))
    if key not in cache:
        if suite == 'external':
            edges = None if case['topology'] == 'all_to_all' else [pair for a, b in case['edges'] for pair in ((a, b), (b, a))]
            result = validate_output(case['source_qasm'], row['out_qasm'], case['n_logical'],
                case['physical_width'], edges, row['initial_layout'], row['final_layout'])
        else:
            result = validate_clifford(case['qasm'], row['out_qasm'], case['n_logical'],
                case['physical_width'], row['initial_layout'], row['final_layout'], case['edges'])
            try:
                width, gates = parse_qasm(row['out_qasm'])
                result['metrics'] = common_metrics(width, gates)
            except (ValueError, TypeError, SyntaxError, ZeroDivisionError, OverflowError) as error:
                result.update(passed=False, metrics_error=str(error))
        if row['output_qubits'] != case['physical_width']:
            result = dict(result, passed=False, output_width_error='Worker output width differs from target')
        result = dict(result, failure_classification=classify_evidence(result, suite))
        cache[key] = result
    result = cache[key]
    status = ('passed' if result['passed'] is True else 'failed' if suite == 'external' else 'unverified_output')
    return status, key, result


def summarize_cells(protocol, checked):
    grouped = defaultdict(list)
    for row in checked:
        grouped[row['variant'], row['case_id'], row['tool']].append(row)
    cells = []
    for variant, tools in [('primary', protocol['tools']), ('control', ['qrust'])]:
        for case in protocol['cases']:
            for tool in tools:
                group = grouped[variant, case['id'], tool]
                timed = [row for row in group if row['observed_phase'] == 'timed']
                warmed = [row for row in group if row['observed_phase'] != 'timed']
                eligible = (len(timed) == protocol['trials'] and len(warmed) >= 1 and
                            all(row['outcome'] == 'passed' for row in timed + warmed))
                cell = {'variant': variant, 'case_id': case['id'], 'tool': tool, 'eligible': eligible,
                        'timed_outcomes': dict(Counter(row['outcome'] for row in timed)),
                        'warmup_outcomes': dict(Counter(row['outcome'] for row in warmed)),
                        'evidence_classifications': dict(Counter(row.get('validation', {}).get('failure_classification', 'not_compiled') for row in group))}
                if eligible:
                    values = {name: [row[name] for row in timed] for name in ('compile_ns', 'parse_ns', 'request_wall_ns')}
                    values.update({name: [row['validation']['metrics'][name] for row in timed]
                                   for name in ('cx_count', 'gate_count', 'depth', 'two_q_depth')})
                    if all('sampled_worker_rss_bytes' in row for row in timed):
                        values['sampled_worker_rss_bytes'] = [row['sampled_worker_rss_bytes'] for row in timed]
                    cell['metrics'] = {name: {'median': statistics.median(data), 'min': min(data), 'max': max(data)}
                                       for name, data in values.items()}
                cells.append(cell)
    return cells


def compare_pairs(protocol, cells):
    indexed = {(cell['variant'], cell['case_id'], cell['tool']): cell for cell in cells}
    pairs = []
    for case in protocol['cases']:
        old = indexed['control', case['id'], 'qrust']
        new = indexed['primary', case['id'], 'qrust']
        pair = {'case_id': case['id'], 'eligible': old['eligible'] and new['eligible'],
                'old_eligible': old['eligible'], 'new_eligible': new['eligible']}
        if pair['eligible']:
            pair['metrics'] = {}
            for name in new['metrics']:
                left, right = old['metrics'][name]['median'], new['metrics'][name]['median']
                pair['metrics'][name] = {'old_median': left, 'new_median': right, 'new_minus_old': right - left,
                                        'new_over_old': right / left if left else None}
        pairs.append(pair)
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--official', action='store_true', help='Also run original primary validators before the full independent replay')
    args = parser.parse_args()
    args.repo, args.campaign, args.output = args.repo.resolve(), args.campaign.resolve(), args.output.resolve()
    ext, matched = campaign.modules(args.repo)
    pair_path = args.campaign / 'paired_protocol.json'
    require(campaign.sha(pair_path) == (args.campaign / 'paired_protocol.sha256').read_text().strip(), 'Paired protocol checksum mismatch')
    pair = json.loads(pair_path.read_text())
    require(campaign.sha(Path(__file__).resolve()) == pair['analysis_source_sha256'], 'Frozen analysis source changed')
    require(campaign.sha(Path(campaign.__file__).resolve()) == pair['controller_sha256'], 'Frozen controller source changed')
    require(not (args.campaign / 'active.lock').exists(), 'Timing session is active')
    protocols, loaded = {}, {}
    for suite in ('external', 'matched'):
        directory = args.campaign / suite
        require(campaign.sha(directory / 'protocol.json') == pair['primary_protocol_sha256'][suite], 'Primary protocol changed')
        protocol = json.loads((directory / 'protocol.json').read_text())
        # Validate the complete schedule/hash envelope for every session BEFORE
        # importing either numerical oracle or simulating any output.
        for filename, expected in protocol['source_sha256'].items():
            path = (args.repo / 'experiments' / Path(filename).name if Path(filename).is_absolute()
                    else args.repo / filename)
            require(campaign.sha(path) == expected, 'Frozen adapter/oracle changed: ' + filename)
        for session in (1, 2):
            require(pair['schedules'][suite][str(session)] == campaign.expanded_schedule(protocol['schedules'][str(session)], session), 'Expanded schedule differs from primary projection')
            loaded[suite, session] = check_dispatch(pair, protocol, directory, suite, session)
        protocols[suite] = protocol
    if args.official:
        ext.validate(SimpleNamespace(study=args.campaign / 'external'))
        matched.validate(SimpleNamespace(output=args.campaign / 'matched', binary=Path(pair['candidate_binary'])))
    from comparison_validation import validate_output, parse_qasm, common_metrics
    from scalable_clifford import validate_clifford
    args.output.mkdir(parents=True, exist_ok=False)
    results = {}
    for suite, protocol in protocols.items():
        cache, counts, classifications, sessions = {}, Counter(), Counter(), {}
        case_by_id = {case['id']: case for case in protocol['cases']}
        for session in (1, 2):
            checked = []
            path = args.output / f'{suite}_session_{session}.jsonl'
            with path.open('x') as stream:
                for index, (variant, raw) in enumerate(loaded[suite, session]):
                    case = case_by_id[raw['case_id']]
                    if raw['status'] == 'compiled':
                        require(type(raw['compile_ns']) is int and raw['compile_ns'] > 0 and
                                type(raw['parse_ns']) is int and raw['parse_ns'] >= 0 and
                                raw['compile_ns'] + raw['parse_ns'] <= raw['request_wall_ns'], 'Invalid timer components')
                        if suite == 'external' and phase_of(raw, suite) != 'replacement_warmup':
                            require(ext.artifact_hash(case, raw) == raw['artifact_sha256'], 'External artifact/layout hash mismatch')
                    outcome, key, result = checked_artifact(case, raw, suite, cache, validate_output,
                                                          validate_clifford, parse_qasm, common_metrics)
                    row = dict(raw, variant=variant, outcome=outcome, observed_phase=phase_of(raw, suite))
                    if result is not None:
                        row.update(validation=result, validation_cache_key=key)
                    stream.write(campaign.canonical(row) + '\n')
                    checked.append(row)
                    counts[f"{variant}:{row['observed_phase']}:{raw['tool']}:{outcome}"] += 1
                    classifications[f"{variant}:{raw['tool']}:{result['failure_classification'] if result is not None else 'not_compiled'}"] += 1
                    if (index + 1) % 200 == 0:
                        print(json.dumps({'suite': suite, 'session': session, 'checked': index + 1, 'unique_artifacts': len(cache)}), flush=True)
            cells = summarize_cells(protocol, checked)
            pairs = compare_pairs(protocol, cells)
            sessions[str(session)] = {'validated_rows_sha256': campaign.sha(path), 'cells': cells, 'old_new_pairs': pairs,
                                     'eligible_pairs': sum(pair['eligible'] for pair in pairs)}
        results[suite] = {'counts': dict(counts), 'evidence_classifications': dict(classifications), 'unique_artifacts': len(cache), 'sessions': sessions,
                         'validation_contract': protocol['validation']}
    summary = {'schema': 1, 'paired_protocol_sha256': campaign.sha(pair_path),
               'analysis_source_sha256': campaign.sha(Path(__file__).resolve()),
               'suites': results, 'scope': 'Per-session old/new Q-Rust comparisons with complete-cell and all-warmup eligibility. New SDK comparisons retain the exact primary schedule and recipes. Repeated trials/seeds are not additional circuits; numerical evidence is not a theorem. Structural counts are not device execution fidelity.'}
    campaign.dump_new(args.output / 'summary.json', summary)
    print(json.dumps({'summary_sha256': campaign.sha(args.output / 'summary.json'),
                     'suites': {name: {'unique_artifacts': item['unique_artifacts'], 'counts': item['counts']} for name, item in results.items()}}), flush=True)


if __name__ == '__main__':
    main()
