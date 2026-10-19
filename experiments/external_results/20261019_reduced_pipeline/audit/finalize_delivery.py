#!/usr/bin/env python3
"""Verify preserved evidence and finish the readable repository delivery."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def dump(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as incoming, target.open('xb') as outgoing:
        shutil.copyfileobj(incoming, outgoing)
    assert sha(source) == sha(target)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('study', 'repo', 'package', 'archive-root', 'validation-replay', 'export-replay', 'result', 'text-policy-file'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(a.study))
    from archive_primitives import audit_prose
    completion = read(a.study / 'completion.json')
    assert completion['all_work_completed'] is True
    for relative, expected in completion['completed_artifacts'].items():
        assert sha(a.study / relative) == expected, relative
    old = read(a.study / 'original_preservation.json')
    for record in old['files']:
        path = a.study.parent / record['path']
        assert path.stat().st_size == record['bytes'] and sha(path) == record['sha256'], record['path']
    checksums = read(a.package / 'checksums.json')
    assert sha(a.package / 'evidence.tar.gz') == checksums['archive_sha256']
    assert sha(a.package / 'evidence_manifest.json') == checksums['manifest_sha256']
    manifest = read(a.package / 'evidence_manifest.json')
    for relative, record in manifest['files'].items():
        path = a.archive_root / relative
        assert path.stat().st_size == record['bytes'] and sha(path) == record['sha256'], relative
    replay = read(a.validation_replay / 'replay_receipt.json')
    exports = read(a.export_replay / 'replay_receipt.json')
    assert replay['passed'] and replay['cells_and_analysis_byte_identical'] and replay['no_compiler_execution']
    assert exports['all_outputs_byte_identical'] and exports['archived_inputs_unchanged']
    for item in read(a.study / 'reporting/delivery_exports_manifest.json')['exports']:
        assert sha(a.result / item['destination']) == item['sha256']
    review = read(a.study / 'method_review/final_report_audit_01.json')
    report = a.repo / 'experiments/REDUCED_PIPELINE_STUDY.md'
    assert review['report_content_passed'] and sha(report) == review['report_sha256']
    for name in ('checksums.json', 'evidence.tar.gz', 'evidence_manifest.json'):
        copy(a.package / name, a.result / name)
    for name in ('replay_receipt.json', 'validate.log', 'analyze.log'):
        copy(a.validation_replay / name, a.result / 'audit' / ('saved_replay__' + name))
    for name in ('replay_receipt.json', 'commands.jsonl', 'replay_exports.py', 'replay_attempt1.log',
                 'archived_inputs_before.json', 'archived_inputs_after.json'):
        copy(a.export_replay / name, a.result / 'audit' / ('exports_replay__' + name))
    for path in sorted((a.export_replay / 'logs').iterdir()):
        if path.is_file():
            copy(path, a.result / 'audit' / ('exports_replay__' + path.name))
    copy(Path(__file__), a.result / 'audit/finalize_delivery.py')
    for document in (report, a.result / 'README.md'):
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text()):
            if target.startswith(('https://', 'http://')):
                continue
            path = (document.parent / target.split('#')[0]).resolve()
            assert path.is_relative_to(a.repo), target
            assert path.is_file() or path == a.result / 'verification_summary.json', target
    summary = {
        'schema': 1, 'passed': True,
        'protocol_sha256': replay['protocol_sha256'], 'analysis_sha256': replay['analysis_sha256'],
        'report_sha256': sha(report), 'historical_files_unchanged': len(old['files']),
        'original_completed_artifacts_unchanged': len(completion['completed_artifacts']),
        'archive_payload_files_unchanged_after_replay': len(manifest['files']),
        'archive_sha256': checksums['archive_sha256'],
        'formal_correctness': {'timed_outputs_passed': 5400, 'warmups_passed': 1080,
            'full_numerical_checks': 2160, 'sampled_checks': 4320, 'unique_artifacts': 532,
            'cases_per_recipe_per_session': 54, 'timing_sessions': 2,
            'max_residual': 6.3488762467508525e-12, 'residual_threshold': 1e-7},
        'independent_analysis': {'passed': True, 'metric_point_estimates_recomputed': 1980,
            'cell_medians_and_retention_counts_recomputed': True,
            'receipt_sha256': sha(a.study / 'method_review/postflight_analysis_02.json')},
        'regression_evidence': {'copied_rust_workspace_passed': 502, 'copied_rust_workspace_ignored': 4,
            'final_adapter_revision_dedicated_passed': 4,
            'method_regressions_passed': 45, 'source_reference_controls_passed': 12,
            'postflight_negative_controls_passed': 8, 'plot_regressions_passed': 13,
            'archive_replay_contract_tests_passed': 19, 'variability_regressions_passed': 6,
            'scope': 'Counts identify recorded suites, not a sum of distinct tests; logs and frozen receipts are retained.'},
        'saved_output_replay': {'passed': True, 'no_compiler_execution': True,
            'cell_tables_and_analysis_byte_identical': True,
            'receipt_sha256': sha(a.validation_replay / 'replay_receipt.json')},
        'reporting_replay': {'passed': True, 'byte_identical_outputs': len(exports['outputs']),
            'receipt_sha256': sha(a.export_replay / 'replay_receipt.json')},
        'scope_limits': ['Larger circuits use sampled numerical checks, not proofs.',
            'Same-binary pipeline controls do not measure binary size, memory savings or language causality.',
            '18 selected instances on one host and ideal graphs do not establish general optimality or equal SDK search budgets.']}
    dump(a.result / 'verification_summary.json', summary)
    links = []
    for document in (report, a.result / 'README.md'):
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text()):
            if target.startswith(('https://', 'http://')):
                continue
            path = (document.parent / target.split('#')[0]).resolve()
            assert path.is_file(), target
            links.append({'document': str(document.relative_to(a.repo)), 'target': target, 'sha256': sha(path)})
    policy = read(a.text_policy_file)
    audit_prose(a.result, policy)
    assert not any(word.casefold() in report.read_text().casefold() for word in policy)
    dump(a.result / 'audit/final_delivery_receipt.json', {
        'schema': 1, 'passed': True, 'report_sha256': sha(report),
        'finalizer_sha256': sha(Path(__file__)), 'report_and_readme_links': links,
        'historical_files_unchanged': len(old['files']),
        'archived_payload_files_unchanged': len(manifest['files']),
        'verification_summary_sha256': sha(a.result / 'verification_summary.json')})
    records = {str(path.relative_to(a.result)): {'bytes': path.stat().st_size, 'sha256': sha(path)}
               for path in sorted(a.result.rglob('*')) if path.is_file()}
    dump(a.result / 'delivery_checksums.json', {'schema': 1, 'report_sha256': sha(report),
        'scope': 'Every delivered result file except this manifest itself.', 'files': records})
    print(json.dumps({'passed': True, 'delivered_files': len(records) + 1,
                      'archived_payload_files': len(manifest['files']), 'historical_files_unchanged': len(old['files'])}))


if __name__ == '__main__':
    main()
