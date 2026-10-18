import math
import random
import unittest

import numpy as np

from comparison_validation import parse_qasm, simulate
from scalable_clifford import State, canonical, multiply, representatives, validate_clifford


def qasm(width, operations):
    return f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{width}];\n' + operations


class CliffordChecks(unittest.TestCase):
    def check(self, source, output, logical, physical=None, initial=None, final=None, **kwargs):
        physical = logical if physical is None else physical
        initial = list(range(logical)) if initial is None else initial
        final = list(range(logical)) if final is None else final
        return validate_clifford(source, output, logical, physical, initial, final,
                                 output_basis=None, **kwargs)

    def test_twenty_four_representatives(self):
        self.assertEqual(len(representatives()), 24)
        for matrix, _ in representatives():
            np.testing.assert_allclose(matrix.conj().T @ matrix, np.eye(2), atol=2e-15)

    def test_signed_stabilizers_against_dense_analytic_evolution(self):
        rng = random.Random(2026101809)
        for width in range(1, 5):
            state = State(width)
            operations = []
            for _ in range(80):
                name = rng.choice(['h', 's'] + (['cx', 'cy', 'cz', 'swap'] if width > 1 else []))
                wires = rng.sample(range(width), 1 if name in {'h', 's'} else 2)
                operations.append(name + ' ' + ','.join(f'q[{q}]' for q in wires) + ';')
                state.apply(name, *wires)
                _, gates = parse_qasm(qasm(width, '\n'.join(operations)))
                zero = np.zeros((1 << width, 1), complex)
                zero[0, 0] = 1
                vector = simulate(width, gates, zero)[:, 0]
                matrices = [np.eye(2), np.array([[0, 1], [1, 0]]),
                            np.array([[1, 0], [0, -1]]), np.array([[0, -1j], [1j, 0]])]
                for x, z, r in state.rows:
                    operator = np.ones((1, 1), complex)
                    for qubit in reversed(range(width)):
                        index = ((x >> qubit) & 1) + 2 * ((z >> qubit) & 1)
                        operator = np.kron(operator, matrices[index])
                    np.testing.assert_allclose((-1) ** r * operator @ vector, vector, atol=1e-13)

    def test_canonical_group_independent_of_generator_basis(self):
        state = State(7)
        for q in range(7):
            state.apply('h', q)
            state.apply('s', q)
        for q in range(6):
            state.apply('cx', q, q + 1)
        expected = canonical(state.rows, 7)
        rows = list(state.rows)
        rng = random.Random(904)
        for _ in range(100):
            a, b = rng.sample(range(7), 2)
            rows[a] = multiply(rows[a], rows[b])
            rng.shuffle(rows)
        self.assertEqual(canonical(rows, 7), expected)

    def test_u3_clifford_and_global_phase(self):
        source = qasm(1, 'h q[0];')
        self.assertTrue(self.check(source, qasm(1, f'u3({math.pi/2},0,{math.pi}) q[0];'), 1)['passed'])
        # XZXZ is -I, a permissible single global phase.
        self.assertTrue(self.check(qasm(1, ''), qasm(1, 'x q[0]; z q[0]; x q[0]; z q[0];'), 1)['passed'])

    def test_relative_phase_error_fails(self):
        result = self.check(qasm(5, ''), qasm(5, 'z q[4];'), 5)
        self.assertFalse(result['passed'])
        self.assertFalse(result['tableau_equal'])

    def test_input_output_permutations(self):
        source = qasm(2, 'h q[0]; cx q[0],q[1];')
        output = qasm(3, 'h q[2]; cx q[2],q[0]; swap q[2],q[1];')
        self.assertTrue(self.check(source, output, 2, 3, [2, 0], [1, 0])['passed'])
        self.assertFalse(self.check(source, output, 2, 3, [2, 0], [0, 1])['passed'])

    def test_unused_wires_and_dirty_ancillas(self):
        source = qasm(2, 'cx q[0],q[1];')
        output = qasm(4, 'cx q[0],q[1]; z q[3];')
        self.assertTrue(self.check(source, output, 2, 4)['passed'])
        dirty = qasm(4, 'cx q[0],q[1]; cx q[1],q[3];')
        self.assertFalse(self.check(source, dirty, 2, 4)['passed'])
        dirty = qasm(4, 'cx q[0],q[1]; x q[3];')
        self.assertFalse(self.check(source, dirty, 2, 4)['passed'])

    def test_temporary_ancilla_use(self):
        source = qasm(2, 'cx q[0],q[1];')
        output = qasm(3, 'cx q[0],q[2]; cx q[2],q[1]; cx q[0],q[2];')
        self.assertTrue(self.check(source, output, 2, 3)['passed'])

    def test_hundred_qubit_permuted_isometry(self):
        source = qasm(100, 'h q[0];' + ''.join(f'cx q[0],q[{q}];' for q in range(1, 100)))
        mapping = list(reversed(range(100)))
        output = qasm(133, 'h q[99];' + ''.join(f'cx q[99],q[{99-q}];' for q in range(1, 100)))
        result = self.check(source, output, 100, 133, mapping, mapping)
        self.assertTrue(result['passed'], result)
        self.assertEqual(result['choi_qubits'], 233)

    def test_non_clifford_and_invalid_inputs_fail_closed(self):
        for operation in ['t q[0];', 'u3(0.1,0,0) q[0];', 'rz(1/0) q[0];', 'rz(1+) q[0];', 'measure q[0] -> c[0];']:
            self.assertFalse(self.check(qasm(1, ''), qasm(1, operation), 1)['passed'])
        self.assertFalse(self.check(qasm(2, ''), qasm(2, ''), 2, initial=[0, 0])['passed'])

    def test_three_cycle_mapping_direction(self):
        source = qasm(3, 'h q[0]; s q[1]; cx q[0],q[2];')
        output = qasm(3, 'h q[0]; s q[1]; cx q[0],q[2]; swap q[0],q[1]; swap q[1],q[2];')
        self.assertTrue(self.check(source, output, 3, final=[2, 0, 1])['passed'])
        self.assertFalse(self.check(source, output, 3, final=[1, 2, 0])['passed'])

    def test_basis_and_coupling_checks(self):
        source = qasm(3, 'cx q[0],q[2];')
        self.assertFalse(self.check(source, source, 3, edges=[(0, 1), (1, 2)])['passed'])
        self.assertFalse(validate_clifford(qasm(1, 'h q[0];'), qasm(1, 'h q[0];'), 1, 1, [0], [0])['passed'])
        reverse = qasm(2, 'cx q[1],q[0];')
        self.assertTrue(self.check(reverse, reverse, 2, edges=[(0, 1)])['passed'])
        self.assertFalse(self.check(reverse, reverse, 2, edges=[(0, 1)], directed_edges=True)['passed'])
        self.assertFalse(self.check(reverse, reverse, 2, edges=[(0, 4)])['passed'])

    def test_accumulated_recognition_error_is_bounded(self):
        # Every local angle is below the recognition threshold, but their sum
        # exceeds the whole-circuit allowance; accepting the tableau alone is wrong.
        output = qasm(1, 'rz(0.0000000001) q[0];' * 1500)
        result = self.check(qasm(1, ''), output, 1)
        self.assertFalse(result['passed'])
        self.assertTrue(result['tableau_equal'])
        self.assertGreater(result['reconstruction_error_sum'], result['total_tolerance'])


if __name__ == '__main__':
    unittest.main()
