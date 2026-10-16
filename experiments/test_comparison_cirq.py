"""Independent semantic regressions for the Cirq comparison adapter."""
from __future__ import annotations

import importlib.util
import unittest

from comparison_cirq import compile_one


HAS_DEPS = all(importlib.util.find_spec(name) is not None
               for name in ("cirq", "qiskit", "numpy"))


def qasm(width, body):
    return f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{width}];\n{body}\n'


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

    def test_rejects_nonunitary_input(self):
        with self.assertRaisesRegex(ValueError, "only unitary"):
            compile_one(qasm(1, "reset q[0];"), 1, 1, [], "all_to_all", 42)


if __name__ == "__main__":
    unittest.main()
