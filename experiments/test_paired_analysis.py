"""Fault injection against paired dispatch and complete-cell evidence rules."""
import json
import unittest

import paired_analysis as analysis
import paired_campaign as campaign
import test_paired_campaign as fixtures


class IntegrityTests(unittest.TestCase):
    run_fake = fixtures.DispatchTests.run_fake
    def load(self, suite):
        directory, records = self.run_fake(suite)
        pair = json.loads((directory.parent / 'paired_protocol.json').read_text())
        protocol = json.loads((directory / 'protocol.json').read_text())
        return directory, pair, protocol

    def rebind_file(self, directory, filename):
        path = directory / 'paired_session_1_receipt.json'
        receipt = json.loads(path.read_text())
        receipt['files_sha256'][filename] = campaign.sha(directory / filename)
        path.write_text(json.dumps(receipt))

    def test_external_and_matched_dispatch_envelopes(self):
        for suite in ('external', 'matched'):
            directory, pair, protocol = self.load(suite)
            ordered = analysis.check_dispatch(pair, protocol, directory, suite, 1)
            self.assertEqual(len(ordered), 6)

    def test_raw_mutation_fails_even_after_rebinding_file_receipt(self):
        directory, pair, protocol = self.load('matched')
        path = directory / 'control_session_1/rows.jsonl'
        observations = analysis.rows(path)
        observations[0]['compile_ns'] += 1
        path.write_text(''.join(campaign.canonical(row) + '\n' for row in observations))
        self.rebind_file(directory, 'control_session_1/rows.jsonl')
        with self.assertRaisesRegex(ValueError, 'missing raw record'):
            analysis.check_dispatch(pair, protocol, directory, 'matched', 1)

    def test_dispatch_reorder_fails_even_after_rebinding_receipt(self):
        directory, pair, protocol = self.load('external')
        path = directory / 'dispatch_session_1.jsonl'
        observations = analysis.rows(path)
        observations[0], observations[1] = observations[1], observations[0]
        for index, row in enumerate(observations):
            row['dispatch_index'] = index
        path.write_text(''.join(campaign.canonical(row) + '\n' for row in observations))
        self.rebind_file(directory, 'dispatch_session_1.jsonl')
        with self.assertRaisesRegex(ValueError, 'dispatch order mismatch'):
            analysis.check_dispatch(pair, protocol, directory, 'external', 1)

    def test_missing_warmup_is_ineligible(self):
        protocol = {'trials': 1, 'tools': ['qrust'], 'cases': [{'id': 'tiny'}]}
        observed = [{'variant': variant, 'case_id': 'tiny', 'tool': 'qrust',
                     'observed_phase': 'timed', 'outcome': 'passed'} for variant in ('primary', 'control')]
        cells = analysis.summarize_cells(protocol, observed)
        self.assertTrue(all(not cell['eligible'] for cell in cells))

    def test_unverified_warmup_is_ineligible(self):
        protocol = {'trials': 1, 'tools': ['qrust'], 'cases': [{'id': 'tiny'}]}
        observed = [{'variant': variant, 'case_id': 'tiny', 'tool': 'qrust',
                     'observed_phase': phase, 'outcome': 'unverified_output' if phase == 'warmup' else 'passed'}
                    for variant in ('primary', 'control') for phase in ('warmup', 'timed')]
        cells = analysis.summarize_cells(protocol, observed)
        self.assertTrue(all(not cell['eligible'] for cell in cells))


if __name__ == '__main__':
    unittest.main()
