"""Independent semantic regressions for the Cirq comparison adapter."""
from __future__ import annotations

import importlib.util
import unittest

from comparison_cirq import SYNTHESIS_ATOL, _serialize, _to_cx, compile_one


HAS_DEPS = all(importlib.util.find_spec(name) is not None
               for name in ("cirq", "qiskit", "numpy"))


def qasm(width, body):
    return f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{width}];\n{body}\n'


@unittest.skipUnless(HAS_DEPS, "Optional comparison dependencies not installed")
class CirqCzClosureTests(unittest.TestCase):
    def assert_lowered_equivalent(self, operation):
        import cirq
        import numpy as np

        original = cirq.Circuit(operation)
        lowered = _to_cx(original, cirq)
        self.assertEqual(sum(op.gate == cirq.CX for op in lowered.all_operations()), 1)
        self.assertTrue(all(len(op.qubits) == 1 or op.gate == cirq.CX
                            for op in lowered.all_operations()))
        expected, actual = cirq.unitary(original), cirq.unitary(lowered)
        phase = np.vdot(expected, actual)
        phase /= abs(phase)
        self.assertLessEqual(np.max(np.abs(actual - phase * expected)), SYNTHESIS_ATOL)
        # Serialization stays a one-for-one representation of the already
        # lowered gates; it must not conceal an extra two-qubit decomposition.
        self.assertEqual(_serialize(lowered, 2, cirq).count("cx q["), 1)

    def test_imported_cu1_odd_pi_is_recognized_semantically(self):
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        for angle in ("pi", "-pi", "3*pi"):
            for operands in ("q[0],q[1]", "q[1],q[0]"):
                with self.subTest(angle=angle, operands=operands):
                    imported = circuit_from_qasm(qasm(2, f"cu1({angle}) {operands};"))
                    imported = imported.transform_qubits(
                        lambda q: cirq.LineQubit(int(q.name.split("_")[1])))
                    operation = next(imported.all_operations())
                    self.assertNotIsInstance(operation.gate, cirq.CZPowGate)
                    self.assertTrue(cirq.CZTargetGateset(allow_partial_czs=False).validate(operation))
                    self.assert_lowered_equivalent(operation)

    def test_global_phase_is_one_scalar_for_entire_two_qubit_unitary(self):
        import cirq
        import numpy as np

        qubits = cirq.LineQubit.range(2)
        cz = np.diag([1, 1, 1, -1]).astype(complex)
        self.assert_lowered_equivalent(cirq.MatrixGate(np.exp(0.37j) * cz).on(*qubits))
        self.assert_lowered_equivalent(cirq.CZPowGate(exponent=-3, global_shift=0.23).on(*qubits))
        # A phase on only one computational-basis state cannot be discarded.
        input_dependent_phase = np.diag([1, 1, np.exp(0.37j), -1])
        with self.assertRaisesRegex(ValueError, "outside full CZ/CX"):
            _to_cx(cirq.Circuit(cirq.MatrixGate(input_dependent_phase).on(*qubits)), cirq)

    def test_fractional_and_other_two_qubit_gates_are_rejected(self):
        import cirq
        import numpy as np

        qubits = cirq.LineQubit.range(2)
        for gate in (cirq.CZ ** 0.5, cirq.CX ** 0.5, cirq.SWAP,
                     cirq.MatrixGate(np.diag([1, 1, 1, 1j]))):
            with self.subTest(gate=gate), self.assertRaises(ValueError):
                _to_cx(cirq.Circuit(gate.on(*qubits)), cirq)

    def test_cz_rounding_uses_matrix_absolute_tolerance_without_relative_slack(self):
        import cirq
        import numpy as np

        qubits = cirq.LineQubit.range(2)
        for phase_error in (0.5 * SYNTHESIS_ATOL, 2.0 * SYNTHESIS_ATOL,
                            0.5e-9 * np.pi):
            # The last exponent error was accepted by math.isclose's default
            # relative tolerance despite exceeding the stated synthesis atol.
            gates = (cirq.CZPowGate(exponent=1 + phase_error / np.pi),
                     cirq.MatrixGate(np.diag([1, 1, 1, -np.exp(1j * phase_error)])))
            for gate in gates:
                with self.subTest(phase_error=phase_error, gate=gate):
                    if phase_error < SYNTHESIS_ATOL:
                        self.assert_lowered_equivalent(gate.on(*qubits))
                    else:
                        with self.assertRaises(ValueError):
                            _to_cx(cirq.Circuit(gate.on(*qubits)), cirq)

    def test_native_cx_is_preserved(self):
        import cirq

        operation = cirq.CX(*reversed(cirq.LineQubit.range(2)))
        self.assertEqual(list(_to_cx(cirq.Circuit(operation), cirq).all_operations()), [operation])


@unittest.skipUnless(HAS_DEPS, "Optional comparison dependencies not installed")
class CirqComparisonTests(unittest.TestCase):
    def assert_equivalent(self, source, logical_width, physical_width, edges, topology):
        import numpy as np
        from qiskit import qasm2
        from qiskit.quantum_info import Statevector

        result = compile_one(source, logical_width, physical_width, edges, topology, seed=42)
        original = qasm2.loads(source, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
        output = qasm2.loads(result["out_qasm"], custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
        initial, final = result["initial_layout"], result["final_layout"]
        self.assertEqual(len(set(initial)), logical_width)
        self.assertEqual(len(set(final)), logical_width)
        self.assertTrue(all(0 <= x < physical_width for x in initial + final))
        self.assertEqual(result["output_qubits"], physical_width)
        self.assertGreater(result["compile_ns"], 0)
        self.assertGreater(result["parse_ns"], 0)
        self.assertFalse(result["seed_used"])
        self.assertTrue(set(output.count_ops()).issubset({"u3", "cx"}))
        if topology != "all_to_all":
            allowed = {tuple(sorted(edge)) for edge in edges}
            for instruction in output.data:
                if instruction.operation.num_qubits == 2:
                    pair = tuple(sorted(output.find_bit(q).index for q in instruction.qubits))
                    self.assertIn(pair, allowed)

        def embed(index, layout):
            return sum(((index >> logical) & 1) << physical
                       for logical, physical in enumerate(layout))

        # Check the entire logical isometry, with spectator qubits initialized
        # to zero. A single global overlap across all columns also catches
        # input-dependent phases that per-basis-state fidelities would miss.
        actual_columns, expected_columns = [], []
        for index in range(1 << logical_width):
            actual = Statevector.from_int(embed(index, initial), 1 << physical_width).evolve(output).data
            logical_expected = Statevector.from_int(index, 1 << logical_width).evolve(original).data
            expected = np.zeros(1 << physical_width, dtype=complex)
            for logical_index, amplitude in enumerate(logical_expected):
                expected[embed(logical_index, final)] = amplitude
            actual_columns.append(actual)
            expected_columns.append(expected)
        overlap = abs(np.vdot(np.asarray(expected_columns), np.asarray(actual_columns))) / (1 << logical_width)
        self.assertAlmostEqual(overlap, 1.0, places=9)
        return result

    def test_all_to_all_ccx_is_timed_and_closed_over_cx(self):
        result = self.assert_equivalent(
            qasm(3, "h q[0];\nry(0.37) q[1];\nccx q[0],q[1],q[2];"),
            3, 3, [(0, 1), (0, 2), (1, 2)], "all_to_all",
        )
        self.assertEqual(result["initial_layout"], [0, 1, 2])
        self.assertEqual(result["final_layout"], [0, 1, 2])
        self.assertIn("cx", result["out_qasm"])

    def test_sparse_routing_composes_initial_and_final_maps(self):
        source = qasm(4, "h q[0];\nry(0.43) q[2];\ncx q[0],q[1];\n"
                      "cx q[1],q[2];\ncx q[2],q[3];\nrz(0.41) q[1];\n"
                      "cx q[3],q[0];\ncx q[0],q[2];")
        result = self.assert_equivalent(source, 4, 4, [(0, 1), (1, 2), (2, 3)], "linear")
        self.assertNotEqual(result["initial_layout"], result["final_layout"])

    def test_spare_physical_wires_and_inactive_logical_wires(self):
        self.assert_equivalent(qasm(3, "h q[2];"), 3, 5,
                               [(0, 1), (1, 2), (1, 3), (3, 4)], "ibm_quito")

    def test_fully_cancelled_circuit_keeps_identity_layout(self):
        result = self.assert_equivalent(qasm(3, "h q[1];\nh q[1];"), 3, 5,
                                        [(0, 1), (1, 2), (1, 3), (3, 4)], "ibm_quito")
        self.assertEqual(result["initial_layout"], [0, 1, 2])
        self.assertEqual(result["final_layout"], [0, 1, 2])

    def test_sparse_high_arity_input(self):
        self.assert_equivalent(qasm(3, "h q[0];\nccx q[0],q[1],q[2];"),
                               3, 3, [(0, 1), (1, 2)], "linear")

    def test_cu1_full_and_fractional_angles_preserve_all_input_states(self):
        for angle in ("pi", "-pi", "3*pi", "pi/2", "0.37", "-0.41"):
            for topology in ("all_to_all", "linear"):
                with self.subTest(angle=angle, topology=topology):
                    self.assert_equivalent(
                        qasm(3, f"cu1({angle}) q[2],q[0];\n"
                             "h q[1];\ncx q[0],q[1];\n"),
                        3, 3, [(0, 1), (1, 2)], topology,
                    )

    def test_rejects_nonunitary_input(self):
        with self.assertRaisesRegex(ValueError, "only unitary"):
            compile_one(qasm(1, "reset q[0];"), 1, 1, [], "all_to_all", 42)


if __name__ == "__main__":
    unittest.main()
