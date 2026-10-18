"""Regression tests for external-source eligibility and semantic normalization."""
import hashlib
import unittest

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import Gate

from external_corpus import (canonical_source, digest, git_blob_digest,
                             normalization_evidence, qasmbench_inventory,
                             unitary_prefix)
from bench_baselines import normalize_qasm


class ExternalCorpusTests(unittest.TestCase):
    def test_inventory_preserves_all_declared_scope_exclusions(self):
        paths = [('small/a_n5/a_n5.qasm', 250000),
                 ('small/a_n5/a_n5_transpiled.qasm', 10),
                 ('medium/a_n15/a_n15.qasm', 10),
                 ('small/b_n5/b_n5.qasm', 250001)]
        tree = {'tree': [{'path': path, 'size': size, 'sha': 'sha'} for path, size in paths]}
        rows = qasmbench_inventory(tree)
        self.assertEqual(len(rows), 4)
        self.assertEqual(sum(row['selection_status'] == 'candidate' for row in rows), 1)
        self.assertEqual({row['exclusion_reason'] for row in rows},
                         {None, 'upstream_compiler_transpiled_variant', 'width_above_14',
                          'raw_source_above_250000_bytes'})

    def test_blob_digest_binds_exact_bytes(self):
        data = b'abc\r\n'
        self.assertEqual(git_blob_digest(data), hashlib.sha1(b'blob 5\0abc\r\n').hexdigest())
        self.assertNotEqual(git_blob_digest(data), git_blob_digest(b'abc\n'))

    def test_terminal_measurements_and_barriers_recorded(self):
        circuit = QuantumCircuit(2, 2)
        circuit.h(0)
        circuit.barrier()
        circuit.cx(0, 1)
        circuit.measure([0, 1], [0, 1])
        prefix, removed = unitary_prefix(circuit)
        self.assertEqual(removed, {'terminal_measurements': 2, 'barriers': 1})
        self.assertEqual(prefix.count_ops(), {'h': 1, 'cx': 1})

    def test_mid_circuit_measurement_is_not_silently_stripped(self):
        circuit = QuantumCircuit(2, 1)
        circuit.measure(0, 0)
        circuit.x(1)
        with self.assertRaisesRegex(ValueError, 'measurement_precedes'):
            unitary_prefix(circuit)

    def test_reset_is_outside_unitary_track(self):
        circuit = QuantumCircuit(1)
        circuit.reset(0)
        with self.assertRaisesRegex(ValueError, 'nonunitary'):
            unitary_prefix(circuit)

    def test_control_flow_is_outside_unitary_track(self):
        circuit = QuantumCircuit(1, 1)
        with circuit.if_test((circuit.clbits[0], True)):
            circuit.x(0)
        with self.assertRaisesRegex(ValueError, 'nonunitary|dynamic'):
            unitary_prefix(circuit)

    def test_source_defined_standard_name_is_matrix_checked(self):
        gate = Gate('rx', 1, [0.4])
        definition = QuantumCircuit(1)
        definition.ry(0.4, 0)
        gate.definition = definition
        circuit = QuantumCircuit(1)
        circuit.append(gate, [0])
        source, lineage = canonical_source(circuit)
        self.assertIn('ry(0.4)', source)
        self.assertEqual(lineage['definition_expansions'], {'rx': 1})
        self.assertTrue(normalization_evidence(circuit, source, source)['passed'])

    def test_controlled_operation_compares_complete_matrix(self):
        # A phase on the controlled operation's target subspace is not a
        # discardable whole-circuit phase. CRX-like naming must not hide it.
        gate = Gate('crx', 2, [0.4])
        definition = QuantumCircuit(2)
        definition.crx(0.4, 0, 1)
        definition.p(0.2, 0)
        gate.definition = definition
        circuit = QuantumCircuit(2)
        circuit.append(gate, [0, 1])
        source, lineage = canonical_source(circuit)
        self.assertEqual(lineage['definition_expansions'], {'crx': 1})
        self.assertIn('u1(0.2)', source)
        self.assertTrue(normalization_evidence(circuit, source, normalize_qasm(source))['passed'])

    def test_common_expansion_is_checked_against_original(self):
        circuit = QuantumCircuit(2)
        circuit.crz(0.739, 0, 1)
        circuit.rzz(-0.315, 0, 1)
        source, _ = canonical_source(circuit)
        evidence = normalization_evidence(circuit, source, normalize_qasm(source))
        self.assertTrue(evidence['passed'])
        self.assertEqual(evidence['mode'], 'full_unitary_numerical')
        self.assertEqual(evidence['input_count'], 4)

    def test_normalization_negative_control_detects_wrong_unitary(self):
        circuit = QuantumCircuit(1)
        circuit.x(0)
        wrong = 'OPENQASM 2.0;\nqreg q[1];\nz q[0];\n'
        self.assertFalse(normalization_evidence(circuit, wrong, wrong)['passed'])

    def test_matrix_valued_high_level_parameters_are_expanded(self):
        from qiskit.circuit.library import UnitaryGate
        circuit = QuantumCircuit(1)
        circuit.append(UnitaryGate(np.array([[0, 1], [1, 0]], dtype=complex)), [0])
        source, lineage = canonical_source(circuit)
        self.assertEqual(lineage['definition_expansions'], {'unitary': 1})
        self.assertTrue(normalization_evidence(circuit, source, source)['passed'])

    def test_permutation_definition_expands_to_correct_swaps(self):
        from qiskit.circuit.library import PermutationGate
        circuit = QuantumCircuit(3)
        circuit.append(PermutationGate([1, 2, 0]), [0, 1, 2])
        source, lineage = canonical_source(circuit)
        self.assertEqual(lineage['definition_expansions'], {'permutation_to_swaps': 1})
        self.assertTrue(normalization_evidence(circuit, source, source)['passed'])

    def test_gate_cap_does_not_depend_on_compiler_outcome(self):
        circuit = QuantumCircuit(1)
        circuit.x(0)
        circuit.x(0)
        with self.assertRaisesRegex(ValueError, 'expanded_source_above'):
            canonical_source(circuit, max_gates=1)


if __name__ == '__main__':
    unittest.main()
