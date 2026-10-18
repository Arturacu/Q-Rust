#!/usr/bin/env python3
"""Untimed old/new Q-Rust smoke over the frozen 100 external and 8 matched cases.

Worker durations are retained solely as diagnostics. This has no warmups or
repetitions and is not part of the formal performance campaign. Validation
starts only after every compiler worker has closed.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

import paired_campaign as campaign
from paired_analysis import checked_artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--old-binary', type=Path, required=True)
    parser.add_argument('--external-baseline', type=Path, required=True)
    parser.add_argument('--matched-baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.resolve())
    ext, _ = campaign.modules(args.repo)
    baseline_provenance = campaign.baseline_provenance(args.external_baseline, args.matched_baseline, args.old_binary)
    protocols = {}
    for suite in ('external', 'matched'):
        baseline = getattr(args, suite + '_baseline')
        protocol = json.loads((baseline / 'protocol.json').read_text())
        for filename, expected in protocol['source_sha256'].items():
            path = (args.repo / 'experiments' / Path(filename).name if Path(filename).is_absolute()
                    else args.repo / filename)
            if campaign.sha(path) != expected:
                raise ValueError('Frozen harness/adapter/oracle changed: ' + filename)
        protocols[suite] = protocol
    args.output.mkdir(parents=True, exist_ok=False)
    campaign.dump_new(args.output / 'protocol.json', {
        'schema': 1, 'scope': 'Untimed correctness/structural-quality smoke. No warmups/repetitions. Diagnostic durations are not performance evidence.',
        'binary_sha256': campaign.sha(args.binary), 'old_binary_sha256': campaign.sha(args.old_binary),
        'baseline_provenance': baseline_provenance,
        'source_provenance': ext.source_provenance(args.repo),
        'driver_sha256': campaign.sha(Path(__file__).resolve()),
        'analysis_source_sha256': campaign.sha(Path(__file__).with_name('paired_analysis.py')),
        'controller_source_sha256': campaign.sha(Path(campaign.__file__).resolve()),
        'baseline_protocol_sha256': {suite: campaign.sha(getattr(args, suite + '_baseline') / 'protocol.json') for suite in protocols},
        'cases': {suite: protocol['cases'] for suite, protocol in protocols.items()},
        'timeout_seconds': {suite: protocol['timeout_seconds'] for suite, protocol in protocols.items()},
        'seed': 20261018})
    workers, captured = {}, []
    started = time.monotonic()
    try:
        for variant, binary in [('old', args.old_binary), ('new', args.binary)]:
            directory = args.output / variant
            directory.mkdir()
            workers[variant] = ext.ProfileWorker('qrust', binary, directory)
        with (args.output / 'raw.jsonl').open('x') as stream:
            for suite, protocol in protocols.items():
                for index, case in enumerate(protocol['cases']):
                    request = {key: case[key] for key in ('qasm', 'n_logical', 'physical_width', 'edges', 'topology')}
                    request['seed'] = 20261018
                    for variant, worker in workers.items():
                        result = ext.timed_request(worker, request, protocol['timeout_seconds'])
                        row = dict(result, suite=suite, variant=variant, case_id=case['id'], seed=20261018)
                        if row['status'] == 'compiled':
                            row['output_sha256'] = campaign.sha(row['out_qasm'])
                        stream.write(campaign.canonical(row) + '\n')
                        stream.flush()
                        captured.append(row)
                    if (index + 1) % 10 == 0 or suite == 'matched':
                        print(json.dumps({'phase': 'untimed_compilation', 'suite': suite,
                            'cases_done': index + 1, 'elapsed_seconds': round(time.monotonic() - started, 2)}), flush=True)
    finally:
        for worker in workers.values():
            worker.close()
    from comparison_validation import validate_output, parse_qasm, common_metrics
    from scalable_clifford import validate_clifford
    lookup = {(suite, case['id']): case for suite, protocol in protocols.items() for case in protocol['cases']}
    cache, counts, classifications, checked = {}, Counter(), Counter(), {}
    with (args.output / 'validated.jsonl').open('x') as stream:
        for index, row in enumerate(captured):
            case = lookup[row['suite'], row['case_id']]
            status, key, result = checked_artifact(case, row, row['suite'], cache,
                validate_output, validate_clifford, parse_qasm, common_metrics)
            value = dict(row, outcome=status)
            if result is not None:
                value.update(validation=result, validation_cache_key=key)
            counts[f"{row['suite']}:{row['variant']}:{status}"] += 1
            classifications[f"{row['suite']}:{row['variant']}:{result['failure_classification'] if result is not None else 'not_compiled'}"] += 1
            checked[row['suite'], row['case_id'], row['variant']] = value
            stream.write(campaign.canonical(value) + '\n')
            if (index + 1) % 20 == 0:
                print(json.dumps({'phase': 'untimed_validation', 'checked': index + 1,
                    'total': len(captured), 'unique_artifacts': len(cache)}), flush=True)
    comparisons = []
    for suite, case_id in lookup:
        old, new = checked[suite, case_id, 'old'], checked[suite, case_id, 'new']
        item = {'suite': suite, 'case_id': case_id, 'old_outcome': old['outcome'], 'new_outcome': new['outcome'],
                'old_evidence_classification': old.get('validation', {}).get('failure_classification', 'not_compiled'),
                'new_evidence_classification': new.get('validation', {}).get('failure_classification', 'not_compiled'),
                'semantically_eligible': old['outcome'] == new['outcome'] == 'passed'}
        if all('metrics' in value.get('validation', {}) for value in (old, new)):
            item['metrics'] = {name: {'old': old['validation']['metrics'][name],
                'new': new['validation']['metrics'][name],
                'delta': new['validation']['metrics'][name] - old['validation']['metrics'][name]}
                for name in ('cx_count', 'gate_count', 'depth', 'two_q_depth')}
            item['structural_regression'] = any(item['metrics'][name]['delta'] > 0
                                               for name in ('cx_count', 'depth', 'two_q_depth'))
        comparisons.append(item)
    summary = {'schema': 1, 'protocol_sha256': campaign.sha(args.output / 'protocol.json'),
        'raw_sha256': campaign.sha(args.output / 'raw.jsonl'), 'validated_sha256': campaign.sha(args.output / 'validated.jsonl'),
        'scope': 'Untimed diagnostic; no speed claims. Structural deltas on unverified outputs are diagnostic only. Large-output non-Clifford recognition failure is unverified, not proven incorrect.',
        'counts': dict(counts), 'evidence_classifications': dict(classifications), 'unique_artifacts': len(cache), 'comparisons': comparisons,
        'structural_regressions': [item for item in comparisons if item.get('structural_regression')],
        'new_definite_failures': [item for item in comparisons if item['new_evidence_classification'] in ('semantic_counterexample', 'structural_contract_failure')],
        'new_not_passed': [item for item in comparisons if item['new_outcome'] != 'passed']}
    campaign.dump_new(args.output / 'summary.json', summary)
    print(json.dumps({'counts': dict(counts), 'structural_regressions': len(summary['structural_regressions']),
        'new_not_passed': len(summary['new_not_passed']), 'unique_artifacts': len(cache)}), flush=True)


if __name__ == '__main__':
    main()
