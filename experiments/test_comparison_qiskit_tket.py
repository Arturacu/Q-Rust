"""Layout and semantic regression checks for the measured compiler adapters."""
import importlib.util
import itertools
import random
import unittest

from comparison_qiskit_tket import compile_one


def qasm(n, body):
    return f'OPENQASM 2.0; include "qelib1.inc"; qreg q[{n}]; {body}'


def encoded_index(index, layout):
    return sum(((index >> logical) & 1) << physical for logical, physical in enumerate(layout))


AVAILABLE = all(importlib.util.find_spec(package) is not None for package in ("qiskit", "pytket", "numpy"))


@unittest.skipUnless(AVAILABLE, "requires the pinned comparison environment")
class AdapterSemantics(unittest.TestCase):
    def assert_equivalent(self, tool, source, n, width, edges, topology):
        import numpy as np
        from qiskit import qasm2
        from qiskit.quantum_info import Operator

        result = compile_one(tool, source, n, width, edges, topology, 42)
        before = Operator(qasm2.loads(source, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)).data
        output = qasm2.loads(result["out_qasm"], custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
        after = Operator(output).data
        inputs = [encoded_index(i, result["initial_layout"]) for i in range(1 << n)]
        outputs = [encoded_index(i, result["final_layout"]) for i in range(1 << n)]
        # Full isometry check includes all output amplitudes; projecting and
        # renormalizing could conceal dirty ancillas or an incorrect layout.
        expected = np.zeros((1 << width, 1 << n), dtype=complex)
        expected[outputs, :] = before
        actual = after[:, inputs]
        overlap = np.vdot(expected, actual)
        self.assertGreater(abs(overlap), 0)
        np.testing.assert_allclose(actual, expected * overlap / abs(overlap), atol=2e-9)
        allowed_edges = {tuple(sorted(e)) for e in edges}
        for instruction in output.data:
            if len(instruction.qubits) == 2:
                self.assertEqual(instruction.operation.name, "cx")
                pair = tuple(sorted(output.find_bit(q).index for q in instruction.qubits))
                self.assertIn(pair, allowed_edges)
            else:
                self.assertEqual(len(instruction.qubits), 1)
        self.assertEqual(output.num_qubits, width)
        self.assertGreater(result["compile_ns"], 0)
        self.assertGreater(result["parse_ns"], 0)
        return result

    def test_three_cycle_implicit_permutation(self):
        source = qasm(3, "h q[0]; cx q[0],q[1]; cx q[1],q[0]; cx q[1],q[2]; cx q[2],q[1];")
        for tool in ("qiskit", "tket"):
            with self.subTest(tool=tool):
                result = self.assert_equivalent(tool, source, 3, 3, list(itertools.combinations(range(3), 2)), "all_to_all")
                if tool == "tket":
                    self.assertEqual(result["initial_layout"], [0, 1, 2])
                    self.assertEqual(result["final_layout"], [1, 2, 0])
                    self.assertEqual(result["out_qasm"].count("cx "), 2)

    def test_sparse_routing_with_extra_idle_physical_qubits(self):
        source = qasm(4, "h q[0]; rx(0.2) q[3]; cx q[0],q[1]; cx q[0],q[2]; cx q[0],q[3]; ry(0.3) q[1]; cx q[1],q[3]; cx q[2],q[3];")
        edges = [(0, 1), (1, 2), (2, 3), (3, 4)]
        for tool in ("qiskit", "tket"):
            with self.subTest(tool=tool):
                self.assert_equivalent(tool, source, 4, 5, edges, "linear")

    def test_idle_logical_wires_and_empty_circuit_are_preserved(self):
        for tool in ("qiskit", "tket"):
            for body in ("", "x q[2];", "cx q[0],q[2];"):
                with self.subTest(tool=tool, body=body):
                    self.assert_equivalent(tool, qasm(3, body), 3, 5, [(0, 1), (1, 2), (2, 3), (3, 4)], "linear")

    def test_toffoli_is_lowered_inside_compilation(self):
        source = qasm(3, "h q[0]; rx(0.4) q[1]; ccx q[0],q[1],q[2];")
        for tool in ("qiskit", "tket"):
            with self.subTest(tool=tool):
                self.assert_equivalent(tool, source, 3, 3, [(0, 1), (1, 2)], "linear")

    def test_tket_cleanup_preserves_sparse_connectivity_and_isometry(self):
        architectures = {
            "linear": [(0, 1), (1, 2), (2, 3), (3, 4)],
            "star": [(0, 1), (0, 2), (0, 3), (0, 4)],
            "ring": [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
            "ibm_quito": [(0, 1), (1, 2), (1, 3), (3, 4)],
        }
        for seed in range(3):
            rng = random.Random(seed)
            body = []
            for _ in range(16):
                a, b = rng.sample(range(4), 2)
                body.append(f"rx({rng.uniform(-3, 3):.17g}) q[{a}]; cx q[{a}],q[{b}];")
            for topology, edges in architectures.items():
                with self.subTest(seed=seed, topology=topology):
                    self.assert_equivalent("tket", qasm(4, " ".join(body)), 4, 5, edges, topology)

    def test_single_qubit_target(self):
        for tool in ("qiskit", "tket"):
            with self.subTest(tool=tool):
                self.assert_equivalent(tool, qasm(1, "rx(0.3) q[0]; rz(-0.7) q[0];"), 1, 1, [], "all_to_all")

    def test_arbitrary_input_contract(self):
        # These gates have observable action away from |0...0>; a compile
        # pipeline that assumes known-zero input is inadmissible here.
        source = qasm(3, "z q[0]; cx q[0],q[2]; rz(0.4) q[1];")
        for tool in ("qiskit", "tket"):
            with self.subTest(tool=tool):
                self.assert_equivalent(tool, source, 3, 3, [(0, 1), (1, 2)], "linear")

    def test_rejects_dynamic_or_wrong_width_input(self):
        for tool in ("qiskit", "tket"):
            for source in (qasm(2, ""), qasm(3, "reset q[0];")):
                with self.subTest(tool=tool, source=source), self.assertRaises(ValueError):
                    compile_one(tool, source, 3, 3, [(0, 1), (1, 2)], "linear", 42)


class AdapterArguments(unittest.TestCase):
    def test_invalid_edges_and_width_are_rejected_before_compilation(self):
        for n, width, edges in ((3, 2, []), (0, 0, []), (2, 2, [(0, 2)]), (2, 2, [(0, 0)])):
            with self.subTest(n=n, width=width, edges=edges), self.assertRaises(ValueError):
                compile_one("qiskit", "", n, width, edges, "linear", 42)

    def test_unknown_tool_is_rejected(self):
        with self.assertRaises(ValueError):
            compile_one("unknown", "", 1, 1, [], "all_to_all", 42)


if __name__ == "__main__":
    unittest.main()
