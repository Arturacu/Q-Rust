"""Fail-closed preflight checks with independent one-qubit mathematics."""
import unittest

import paired_analysis as analysis
import paired_campaign as campaign
import test_paired_campaign  # makes the repository experiments importable
from comparison_validation import common_metrics, parse_qasm, validate_output
from scalable_clifford import validate_clifford


class QualityTests(unittest.TestCase):
    def check(self, emitted, suite='external', width=1, initial_layout=None):
        source = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[1];'
        output = source + emitted
        case = {'id': 'identity', 'qasm': source, 'source_qasm': source, 'n_logical': 1,
                'physical_width': 1, 'edges': [], 'topology': 'all_to_all'}
        row = {'status': 'compiled', 'out_qasm': output, 'output_sha256': campaign.sha(output),
               'initial_layout': [0] if initial_layout is None else initial_layout,
               'final_layout': [0], 'output_qubits': width}
        return analysis.checked_artifact(case, row, suite, {}, validate_output, validate_clifford,
                                         parse_qasm, common_metrics)

    def test_general_oracle_catches_phase_error_invisible_on_zero_input(self):
        status, _, result = self.check('u3(0,0,pi) q[0];')
        self.assertEqual(status, 'failed')
        self.assertFalse(result['passed'])
        self.assertEqual(result['evidence_class'], 'full_isometry')
        self.assertEqual(result['failure_classification'], 'semantic_counterexample')

    def test_identity_output_passes(self):
        status, _, result = self.check('u3(0,0,0) q[0];')
        self.assertEqual(status, 'passed')
        self.assertTrue(result['passed'])

    def test_nonclifford_output_remains_unverified(self):
        status, _, result = self.check('u3(0,0,pi/4) q[0];', suite='matched')
        self.assertEqual(status, 'unverified_output')
        self.assertFalse(result['passed'])
        self.assertIn('Non-Clifford', result['error'])
        self.assertEqual(result['failure_classification'], 'unsupported_recognition')

    def test_inconsistent_reported_width_cannot_pass(self):
        status, _, result = self.check('u3(0,0,0) q[0];', width=2)
        self.assertEqual(status, 'failed')
        self.assertIn('output_width_error', result)
        self.assertEqual(result['failure_classification'], 'structural_contract_failure')

    def test_recognized_clifford_counterexample_is_not_hidden_by_legacy_label(self):
        status, _, result = self.check('u3(0,0,pi) q[0];', suite='matched')
        self.assertEqual(status, 'unverified_output')  # original outcome contract
        self.assertIs(result['tableau_equal'], False)
        self.assertEqual(result['failure_classification'], 'semantic_counterexample')

    def test_bad_map_is_a_structural_failure(self):
        status, _, result = self.check('u3(0,0,0) q[0];', suite='matched', initial_layout=[1])
        self.assertEqual(status, 'unverified_output')
        self.assertEqual(result['failure_classification'], 'structural_contract_failure')

    def test_illegal_coupling_is_a_structural_failure(self):
        source = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[2];'
        output = source + 'cx q[0],q[1];'
        case = {'id': 'illegal_edge', 'qasm': source, 'source_qasm': source, 'n_logical': 2,
                'physical_width': 2, 'edges': [], 'topology': 'linear'}
        row = {'status': 'compiled', 'out_qasm': output, 'output_sha256': campaign.sha(output),
               'initial_layout': [0, 1], 'final_layout': [0, 1], 'output_qubits': 2}
        status, _, result = analysis.checked_artifact(case, row, 'matched', {}, validate_output,
            validate_clifford, parse_qasm, common_metrics)
        self.assertEqual(status, 'unverified_output')
        self.assertEqual(result['failure_classification'], 'structural_contract_failure')


if __name__ == '__main__':
    unittest.main()
