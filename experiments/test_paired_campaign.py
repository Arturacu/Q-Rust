"""Controller regressions use tiny fake workers, never timed SDK compilations."""
import contextlib
import io
import json
import os
from pathlib import Path
import random
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import paired_campaign as campaign

REPO = Path(os.environ.get('QRUST_TEST_REPO', Path(__file__).resolve().parent.parent))
EXT, MATCHED = campaign.modules(REPO)


class BaselineProvenanceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.old_binary = self.root / 'old_binary'
        self.old_binary.write_bytes(b'previously timed compiler')
        self.directories = {suite: self.root / suite for suite in ('external', 'matched')}
        for directory in self.directories.values():
            directory.mkdir()
        self.protocols = {
            'external': {'qrust_binary_sha256': campaign.sha(self.old_binary),
                         'compiler_commit': 'previous-timed-source',
                         'compiler_provenance': {'git_head': 'previous-timed-source'}},
            'matched': {'binary_sha256': campaign.sha(self.old_binary),
                        'compiler_commit': 'previous-timed-source'},
        }
        self.write_protocols()

    def write_protocols(self):
        for suite, protocol in self.protocols.items():
            (self.directories[suite] / 'protocol.json').write_text(json.dumps(protocol))
        (self.directories['external'] / 'protocol.sha256').write_text(
            campaign.sha(self.directories['external'] / 'protocol.json') + '\n')

    def provenance(self):
        return campaign.baseline_provenance(self.directories['external'], self.directories['matched'],
                                            self.old_binary)

    def test_identity_comes_from_both_timed_protocols(self):
        actual = self.provenance()
        self.assertEqual(actual['old_compiler_commit'], 'previous-timed-source')
        self.assertEqual(actual['old_binary_sha256'], campaign.sha(self.old_binary))
        self.assertEqual(actual['old_compiler_provenance'], self.protocols['external']['compiler_provenance'])
        self.assertEqual(actual['baseline_protocol_sha256'],
                         {suite: campaign.sha(path / 'protocol.json') for suite, path in self.directories.items()})

    def test_rejects_different_binary(self):
        self.old_binary.write_bytes(b'not the timed compiler')
        with self.assertRaisesRegex(ValueError, 'both supplied baseline protocols'):
            self.provenance()

    def test_rejects_one_protocol_binding_another_binary(self):
        for suite, field in [('external', 'qrust_binary_sha256'), ('matched', 'binary_sha256')]:
            with self.subTest(suite=suite):
                old = self.protocols[suite][field]
                self.protocols[suite][field] = '0' * 64
                self.write_protocols()
                with self.assertRaisesRegex(ValueError, 'both supplied baseline protocols'):
                    self.provenance()
                self.protocols[suite][field] = old

    def test_rejects_missing_or_disagreeing_timed_commit(self):
        for commit in ('another-source', ''):
            with self.subTest(commit=commit):
                self.protocols['matched']['compiler_commit'] = commit
                self.write_protocols()
                with self.assertRaisesRegex(ValueError, 'compiler commits disagree'):
                    self.provenance()

    def test_rejects_external_protocol_with_changed_checksum(self):
        (self.directories['external'] / 'protocol.json').write_text(json.dumps({
            **self.protocols['external'], 'unbound_change': True}))
        with self.assertRaisesRegex(ValueError, 'protocol checksum changed'):
            self.provenance()


class ScheduleTests(unittest.TestCase):
    def test_original_projection_and_pair_balance(self):
        for count, trials, tools in [(100, 7, ['qrust', 'qiskit', 'cirq', 'tket', 'cirq_one_pass', 'tket_light']),
                                    (8, 3, ['qrust', 'qiskit', 'cirq', 'tket'])]:
            first = {}
            for session in (1, 2):
                rng = random.Random(901 + session)
                ids = [f'case_{index}' for index in range(count)]
                rng.shuffle(ids)
                schedule = []
                for case_id in ids:
                    for trial in range(-1, trials):
                        shuffled = list(tools)
                        rng.shuffle(shuffled)
                        schedule.extend({'case_id': case_id, 'tool': tool, 'trial': trial, 'seed': 100 + max(0, trial)}
                                        for tool in shuffled)
                expanded = campaign.expanded_schedule(schedule, session)
                self.assertEqual(len(expanded), count * (trials + 1) * (len(tools) + 1))
                primary = [{key: value for key, value in item.items()
                            if key not in ('variant', 'primary_schedule_index')}
                           for item in expanded if item['variant'] == 'primary']
                self.assertEqual(primary, schedule)
                for index, item in enumerate(expanded):
                    if item['variant'] != 'control':
                        continue
                    key = item['case_id'], item['trial']
                    before = (index + 1 < len(expanded) and
                              expanded[index + 1]['primary_schedule_index'] == item['primary_schedule_index'])
                    companion = expanded[index + 1] if before else expanded[index - 1]
                    self.assertEqual(companion['variant'], 'primary')
                    self.assertEqual(companion['primary_schedule_index'], item['primary_schedule_index'])
                    if session == 1:
                        first[key] = before
                    else:
                        self.assertNotEqual(first[key], before)


class FakeWorker:
    outcomes = {}
    calls = []

    def __init__(self, tool, binary, outdir):
        self.tool, self.binary, self.outdir = tool, binary, outdir
        self.process = object()

    def call(self, request, timeout):
        key = self.binary.name
        FakeWorker.calls.append(key)
        outcomes = FakeWorker.outcomes.get(key, [])
        status = outcomes.pop(0) if outcomes else 'compiled'
        if status != 'compiled':
            self.process = None
            return {'status': status, 'error': status}
        self.process = object()
        return {'status': 'compiled', 'compile_ns': 10, 'parse_ns': 1,
                'out_qasm': request['qasm'], 'output_qubits': 1,
                'initial_layout': [0], 'final_layout': [0]}

    def close(self):
        self.process = None


class DispatchTests(unittest.TestCase):
    def run_fake(self, suite, outcomes=None):
        FakeWorker.outcomes = {key: list(value) for key, value in (outcomes or {}).items()}
        FakeWorker.calls = []
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        output = Path(temporary.name)
        directory = output / suite
        directory.mkdir()
        qasm = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[1]; h q[0];'
        case = {'id': 'tiny', 'qasm': qasm, 'source_qasm': qasm, 'n_logical': 1,
                'physical_width': 1, 'edges': [], 'topology': 'all_to_all',
                'source_sha256': campaign.sha(qasm), 'normalized_sha256': campaign.sha(qasm)}
        schedule = [{'case_id': 'tiny', 'tool': 'qrust', 'trial': trial, 'seed': 20 + max(0, trial)}
                    for trial in (-1, 0, 1)]
        protocol = {'tools': ['qrust'], 'cases': [case], 'timeout_seconds': 1, 'trials': 2,
                    'schedules': {'1': schedule}}
        pair = {'schedules': {suite: {'1': campaign.expanded_schedule(schedule, 1)}}}
        campaign.dump_new(output / 'paired_protocol.json', pair)
        campaign.dump_new(directory / 'protocol.json', protocol)
        args = SimpleNamespace(repo=REPO, binary=Path('candidate'), old_binary=Path('old'),
                               output=output, suite=suite, session=1)
        with patch.object(campaign, 'verify', return_value=(pair, {suite: protocol})), \
                patch.object(EXT, 'ProfileWorker', FakeWorker), contextlib.redirect_stdout(io.StringIO()):
            campaign.time_session(args)
        receipt = json.loads((directory / 'paired_session_1_receipt.json').read_text())
        self.assertTrue(receipt['completed'])
        self.assertEqual(receipt['primary_scheduled'], 3)
        self.assertEqual(receipt['control_scheduled'], 3)
        for path, expected in receipt['files_sha256'].items():
            self.assertEqual(campaign.sha(directory / path), expected)
        records = [json.loads(line) for line in (directory / 'dispatch_session_1.jsonl').read_text().splitlines()]
        return directory, records

    def test_external_primary_schema_and_independent_control(self):
        directory, records = self.run_fake('external')
        primary = [json.loads(line) for line in (directory / 'session_01/timed.jsonl').read_text().splitlines()]
        control = [json.loads(line) for line in (directory / 'control_session_1/timed.jsonl').read_text().splitlines()]
        self.assertEqual(len(primary), 2)
        self.assertEqual(len(control), 2)
        self.assertTrue(all(row['compile_ns'] == 10 for row in primary + control))
        self.assertEqual(len(records), 6)
        self.assertEqual([row['variant'] for row in records],
                         ['control' if name == 'old' else 'primary' for name in FakeWorker.calls])

    def test_external_control_restart_does_not_omit_candidate(self):
        directory, records = self.run_fake('external', {'old': ['timeout', 'compiled']})
        self.assertEqual(len(records), 7)
        self.assertEqual(sum(row['phase'] == 'replacement_warmup' for row in records), 1)
        primary = [json.loads(line) for line in (directory / 'session_01/timed.jsonl').read_text().splitlines()]
        self.assertTrue(all(row['status'] == 'compiled' for row in primary))

    def test_matched_failed_control_warmup_blocks_only_control(self):
        directory, records = self.run_fake('matched', {'old': ['timeout']})
        primary = [json.loads(line) for line in (directory / 'session_1.jsonl').read_text().splitlines()]
        control = [json.loads(line) for line in (directory / 'control_session_1/rows.jsonl').read_text().splitlines()]
        self.assertEqual([row['status'] for row in control], ['timeout', 'blocked_by_warmup', 'blocked_by_warmup'])
        self.assertTrue(all(row['status'] == 'compiled' for row in primary))
        self.assertEqual(len(records), 6)

    def test_matched_failed_timed_requires_replacement(self):
        directory, records = self.run_fake('matched', {'candidate': ['compiled', 'timeout', 'compiled']})
        primary = [json.loads(line) for line in (directory / 'session_1.jsonl').read_text().splitlines()]
        self.assertEqual([row['phase'] for row in primary], ['warmup', 'timed', 'replacement_warmup', 'timed'])
        self.assertEqual([row['status'] for row in primary], ['compiled', 'timeout', 'compiled', 'compiled'])
        self.assertEqual(len(records), 7)


if __name__ == '__main__':
    unittest.main()
