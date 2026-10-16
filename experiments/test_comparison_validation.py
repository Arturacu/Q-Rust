"""Independent checker regression tests (optional NumPy dependency)."""
import math
import unittest

try:
    import numpy as np
    from comparison_validation import common_metrics, parse_qasm, simulate, validate_output
except ImportError:
    np = None


def qasm(width, operations=""):
    return f'OPENQASM 2.0; include "qelib1.inc"; qreg q[{width}]; {operations}'


@unittest.skipIf(np is None, "optional comparison-checker NumPy dependency not installed")
class ComparisonValidationTests(unittest.TestCase):
    def check(self, source, output, n=1, p=None, **kwargs):
        return validate_output(qasm(n, source), qasm(p or n, output), n, p or n,
                               kwargs.pop("edges", None), **kwargs)

    def test_analytic_hadamard_with_shared_global_phase(self):
        result = self.check("h q[0];", "u3(pi/2,0,pi) q[0];")
        self.assertTrue(result["passed"], result)
        self.assertEqual(result["evidence_class"], "full_isometry")
        self.assertEqual(result["checked_input_count"], 2)

    def test_wrong_input_dependent_phase_is_not_global_phase(self):
        result = self.check("", "u3(0,0,pi) q[0];")
        self.assertFalse(result["passed"])
        self.assertGreater(result["max_input_residual"], 1)

    def test_controlled_full_turn_retains_control_phase(self):
        result = self.check("crz(2*pi) q[0],q[1];", "", n=2)
        self.assertFalse(result["passed"])
        correct = self.check("crz(2*pi) q[0],q[1];", "u3(0,0,pi) q[0];", n=2)
        self.assertTrue(correct["passed"], correct)

    def test_dirty_ancillas_not_projected_and_normalized(self):
        for width, evidence in ((2, "full_isometry"), (7, "sampled_isometry")):
            with self.subTest(width=width):
                result = self.check("", f"u3(pi/2,0,pi) q[{width - 1}];", p=width,
                                    initial_layout=[0], final_layout=[0])
                self.assertFalse(result["passed"])
                self.assertEqual(result["evidence_class"], evidence)
                self.assertAlmostEqual(result["max_ancilla_leakage_probability"], 0.5)
                self.assertAlmostEqual(result["isometry_fidelity"], 0.5)

    def test_entangled_ancilla_is_rejected(self):
        result = self.check("", "cx q[0],q[1];", p=2, initial_layout=[0], final_layout=[0])
        self.assertFalse(result["passed"])
        self.assertAlmostEqual(result["max_ancilla_leakage_probability"], 1)

    def test_valid_routed_permutation_and_false_map(self):
        swaps = "cx q[0],q[1]; cx q[1],q[0]; cx q[0],q[1];"
        result = self.check("", swaps, n=2, initial_layout=[0, 1], final_layout=[1, 0])
        self.assertTrue(result["passed"], result)
        wrong = self.check("", swaps, n=2, initial_layout=[0, 1], final_layout=[0, 1])
        self.assertFalse(wrong["passed"])

    def test_padded_permutation_and_full_layout_are_supported(self):
        for layout in ([2], [2, 0, 1]):
            result = self.check("h q[0];", "u3(pi/2,0,pi) q[2];", p=3,
                                initial_layout=layout, final_layout=layout)
            self.assertTrue(result["passed"], result)

    def test_invalid_layouts_fail_closed(self):
        for layout in ([0, 0], [0, 3], [0], [False, 1], [0, 1, 2, 3]):
            with self.subTest(layout=layout):
                result = self.check("", "", n=2, p=3, initial_layout=layout, final_layout=[0, 1])
                self.assertFalse(result["passed"], result)
                self.assertEqual(result["evidence_class"], "invalid")

    def test_illegal_basis_and_directed_edge_fail(self):
        wrong_basis = self.check("h q[0];", "h q[0];")
        self.assertFalse(wrong_basis["passed"])
        self.assertFalse(wrong_basis["basis_legal"])
        wrong_direction = self.check("cx q[1],q[0];", "cx q[1],q[0];", n=2, edges=[(0, 1)])
        self.assertFalse(wrong_direction["passed"])
        self.assertTrue(wrong_direction["basis_legal"])
        self.assertFalse(wrong_direction["coupling_legal"])

    def test_sampled_inputs_reproducible_and_use_common_phase(self):
        first = self.check("", "", n=7, seed=42)
        second = self.check("", "", n=7, seed=42)
        self.assertTrue(first["passed"])
        self.assertEqual(first["checked_input_count"], 8)
        self.assertEqual(first["max_input_residual"], second["max_input_residual"])
        bad = self.check("", "u3(0,0,pi) q[0];", n=7)
        self.assertFalse(bad["passed"])

    def test_common_metric_definition(self):
        width, gates = parse_qasm(qasm(3, "u3(0,0,pi) q[0]; u3(0,0,pi) q[2]; cx q[0],q[1]; u3(0,0,pi) q[0]; cx q[1],q[2];"))
        self.assertEqual(common_metrics(width, gates), {"gate_count": 5, "cx_count": 2, "depth": 3, "two_q_depth": 2})

    def test_malformed_nonunitary_and_nonfinite_inputs_fail(self):
        for bad in ("reset q[0];", "measure q[0] -> c[0];", "u3(1e999,0,0) q[0];",
                    "cx q[0],q[0];", "u3(0,0) q[0];", "x q[1];"):
            with self.subTest(bad=bad):
                result = self.check(bad, "")
                self.assertFalse(result["passed"])
                self.assertEqual(result["evidence_class"], "invalid")

    def test_analytic_fourier_oracle(self):
        from generate_fixtures import make_qft
        width, gates = parse_qasm(make_qft(3).qasm())
        actual = simulate(width, gates, np.eye(8))
        expected = np.exp(2j * math.pi * np.outer(np.arange(8), np.arange(8)) / 8) / math.sqrt(8)
        self.assertLess(float(np.max(np.abs(actual - expected))), 1e-12)

    def test_standard_gate_matrices_against_optional_qiskit_oracle(self):
        try:
            from qiskit import qasm2
            from qiskit.quantum_info import Operator
        except ImportError:
            self.skipTest("optional second independent Qiskit oracle not installed")
        cases = ["h q[0];", "x q[0];", "y q[0];", "z q[0];", "s q[0];", "sdg q[0];",
                 "t q[0];", "tdg q[0];", "sx q[0];", "sxdg q[0];", "id q[0];",
                 "rx(0.7) q[0];", "ry(0.7) q[0];", "rz(0.7) q[0];",
                 "u1(0.7) q[0];", "u2(0.2,0.4) q[0];", "u3(0.2,0.3,0.4) q[0];",
                 "cx q[1],q[0];", "cy q[0],q[1];", "cz q[0],q[1];", "ch q[0],q[1];",
                 "crx(0.7) q[0],q[1];", "cry(0.7) q[0],q[1];", "crz(0.7) q[0],q[1];",
                 "cu1(0.7) q[0],q[1];", "rzz(0.7) q[0],q[1];", "swap q[0],q[2];", "ccx q[0],q[2],q[1];"]
        for operation in cases:
            with self.subTest(operation=operation):
                text = qasm(3, operation)
                width, gates = parse_qasm(text)
                actual = simulate(width, gates, np.eye(8))
                expected = Operator(qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)).data
                self.assertLess(float(np.max(np.abs(actual - expected))), 1e-12)


if __name__ == "__main__":
    unittest.main()
