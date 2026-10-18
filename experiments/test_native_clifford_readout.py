import random
import unittest

import numpy as np

from comparison_validation import Gate, simulate
from native_clifford_readout import (Measurement, Program, compare_programs,
                                     parse_program_qasm, parse_tket_native, readout_signature)


def dense_distribution(program, labels):
    """Independent projective branching, preserving incoherence after readout."""
    vector = np.zeros((1 << program.width, 1), dtype=complex)
    vector[0, 0] = 1
    branches = [(vector, {})]
    indices = np.arange(1 << program.width)
    for op in program.operations:
        if isinstance(op, Measurement):
            expanded = []
            for state, classical in branches:
                for value in [0, 1]:
                    projected = state.copy()
                    projected[((indices >> op.qubit) & 1) != value] = 0
                    if np.linalg.norm(projected) > 1e-14:
                        expanded.append((projected, {**classical, op.classical_label: value}))
            branches = expanded
        else:
            branches = [(simulate(program.width, (op,), state), classical)
                        for state, classical in branches]
    result = np.zeros(1 << len(labels))
    for state, classical in branches:
        result[sum(classical[label] << i for i, label in enumerate(labels))] += float(np.vdot(state, state).real)
    return result


def signature_distribution(program):
    labels = program.measured_labels
    parities, info = readout_signature(program, labels)
    result = np.zeros(1 << len(labels))
    for bits in range(len(result)):
        if all((bits & mask).bit_count() % 2 == parity for mask, parity in parities):
            result[bits] = 2.0 ** -info['uniform_support_log2_size']
    return result


def qasm(width, operations, bits=3):
    return f'OPENQASM 2.0; include "qelib1.inc"; qreg q[{width}]; creg c[{bits}];' + operations


class NativeReadoutChecks(unittest.TestCase):
    def test_random_projective_branching_crosscheck(self):
        rng = random.Random(2026101811)
        for case in range(100):
            width = rng.randrange(1, 5)
            operations = []
            labels = []
            for _ in range(25):
                if len(labels) < 4 and rng.random() < .2:
                    label = f'c[{len(labels)}]'
                    labels.append(label)
                    operations.append(Measurement(rng.randrange(width), label))
                elif width > 1 and rng.random() < .5:
                    operations.append(Gate(rng.choice(['cx', 'cy', 'cz', 'swap']), (), tuple(rng.sample(range(width), 2))))
                else:
                    operations.append(Gate(rng.choice(['h', 's', 'sdg', 'x', 'y', 'z', 'sx']), (), (rng.randrange(width),)))
            if not labels:
                labels.append('c[0]')
                operations.append(Measurement(rng.randrange(width), labels[0]))
            if case % 3 == 0:
                operations.append(Measurement(rng.randrange(width), labels[0]))
            program = Program(width, tuple(operations), tuple(labels))
            np.testing.assert_allclose(signature_distribution(program), dense_distribution(program, program.measured_labels), atol=1e-12)

    def test_joint_correlations_not_only_marginals(self):
        bell = parse_program_qasm(qasm(2, 'h q[0]; cx q[0],q[1]; measure q[0] -> c[0]; measure q[1] -> c[1];'))
        independent = parse_program_qasm(qasm(2, 'h q[0]; h q[1]; measure q[0] -> c[0]; measure q[1] -> c[1];'))
        self.assertFalse(compare_programs(bell, independent)['passed'])
        self.assertEqual(readout_signature(bell, bell.measured_labels)[0], ((3, 0),))

    def test_wire_reuse_after_measurement(self):
        anti = parse_program_qasm(qasm(1, 'h q[0]; measure q[0] -> c[0]; x q[0]; measure q[0] -> c[1];'))
        corr = parse_program_qasm(qasm(1, 'h q[0]; measure q[0] -> c[0]; measure q[0] -> c[1];'))
        self.assertEqual(readout_signature(anti, anti.measured_labels)[0], ((3, 1),))
        self.assertFalse(compare_programs(anti, corr)['passed'])
        np.testing.assert_allclose(signature_distribution(anti), dense_distribution(anti, anti.measured_labels))

    def test_classical_overwrite_keeps_measurement_dephasing(self):
        program = parse_program_qasm(qasm(1, 'h q[0]; measure q[0] -> c[0]; h q[0]; measure q[0] -> c[0];'))
        np.testing.assert_allclose(signature_distribution(program), [0.5, 0.5])
        np.testing.assert_allclose(signature_distribution(program), dense_distribution(program, program.measured_labels))

    def test_partial_measurement_traces_unobserved_quantum_outputs(self):
        source = parse_program_qasm(qasm(2, 'h q[0]; cx q[0],q[1]; measure q[0] -> c[0];'))
        output = parse_program_qasm(qasm(4, 'h q[3]; measure q[3] -> c[0]; x q[0]; h q[1];'))
        result = compare_programs(source, output)
        self.assertTrue(result['passed'], result)
        self.assertFalse(result['quantum_state_equivalence'])
        self.assertFalse(result['arbitrary_input_equivalence'])

    def test_classical_labels_determine_comparison_order(self):
        source = parse_program_qasm(qasm(2, 'x q[0]; measure q[0] -> c[0]; measure q[1] -> c[1];'))
        valid = parse_program_qasm(qasm(3, 'x q[2]; measure q[0] -> c[1]; measure q[2] -> c[0];'))
        invalid = parse_program_qasm(qasm(3, 'x q[2]; measure q[0] -> c[0]; measure q[2] -> c[1];'))
        self.assertTrue(compare_programs(source, valid)['passed'])
        self.assertFalse(compare_programs(source, invalid)['passed'])

    def test_unused_classical_storage_ignored_but_missing_observed_labels_rejected(self):
        source = parse_program_qasm(qasm(1, 'measure q[0] -> c[0];', bits=9))
        output = parse_program_qasm(qasm(2, 'measure q[1] -> c[0];', bits=1))
        self.assertTrue(compare_programs(source, output)['passed'])
        wrong = parse_program_qasm(qasm(1, 'measure q[0] -> c[1];'))
        self.assertFalse(compare_programs(source, wrong)['passed'])

    def test_non_clifford_reset_and_classical_feedback_fail_closed(self):
        source = parse_program_qasm(qasm(1, 'measure q[0] -> c[0];'))
        noncliff = parse_program_qasm(qasm(1, 't q[0]; measure q[0] -> c[0];'))
        self.assertFalse(compare_programs(source, noncliff)['passed'])
        for operations in ['reset q[0]; measure q[0] -> c[0];', 'if(c==0) x q[0]; measure q[0] -> c[0];', 'barrier q[9]; measure q[0] -> c[0];']:
            with self.assertRaises(ValueError):
                parse_program_qasm(qasm(1, operations))

    def test_readout_success_is_not_arbitrary_input_equivalence(self):
        source = parse_program_qasm(qasm(2, 'measure q[0] -> c[0]; measure q[1] -> c[1];'))
        output = parse_program_qasm(qasm(2, 'cx q[0],q[1]; measure q[0] -> c[0]; measure q[1] -> c[1];'))
        result = compare_programs(source, output)
        self.assertTrue(result['passed'])
        self.assertFalse(result['arbitrary_input_equivalence'])

    def test_recognition_error_budget_is_not_reset_by_measurements(self):
        source = parse_program_qasm(qasm(1, 'measure q[0] -> c[0];'))
        output = Program(1, tuple([Gate('rz', (1e-10,), (0,))] * 1500 + [Measurement(0, 'c[0]')]), ('c[0]',))
        result = compare_programs(source, output)
        self.assertFalse(result['passed'])
        self.assertTrue(result['parity_subgroups_equal'])
        self.assertGreater(result['recognition_residual_sum'], result['numerical_recognition_limit'])

    def test_duplicate_observation_labels_fail_closed(self):
        source = parse_program_qasm(qasm(1, 'measure q[0] -> c[0];'))
        with self.assertRaises(ValueError):
            readout_signature(source, ['c[0]', 'c[0]'])

    def test_native_tket_order_and_physical_wire_compression(self):
        native = {'qubits': [['node', [33]], ['node', [101]]],
                  'bits': [['c', [0]], ['c', [1]]], 'implicit_permutation': [[['node', [33]], ['node', [101]]], [['node', [101]], ['node', [33]]]],
                  'commands': [
                      {'op': {'type': 'H'}, 'args': [['node', [33]]]},
                      {'op': {'type': 'Measure'}, 'args': [['node', [33]], ['c', [0]]]},
                      {'op': {'type': 'X'}, 'args': [['node', [33]]]},
                      {'op': {'type': 'Measure'}, 'args': [['node', [33]], ['c', [1]]]},
                  ]}
        expected = parse_program_qasm(qasm(1, 'h q[0]; measure q[0] -> c[0]; x q[0]; measure q[0] -> c[1];'))
        self.assertTrue(compare_programs(expected, parse_tket_native(native))['passed'])


if __name__ == '__main__':
    unittest.main()
