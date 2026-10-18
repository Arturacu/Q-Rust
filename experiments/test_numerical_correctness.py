"""Adversarial whole-unitary regressions for the public compiler pipeline.

Run with QRUST_COMPARISON_BINARY pointing at a freshly built comparison_runner.
The expected matrices below are analytic products of 2x2/4x4 operators, not
Q-Rust matrices or its verifier. Every basis column is checked with one common
global phase. These bounded numerical tests are evidence, not a general proof.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import subprocess
import unittest

try:
    import numpy as np
    from comparison_validation import parse_qasm, simulate, validate_output
except ImportError:
    np = None


TOLERANCE = 1e-7
REPETITIONS = 1000
TINY_ANGLE = 9e-10


@dataclass
class NumericalCase:
    name: str
    width: int
    operations: str
    expected: object
    category: str
    topology: str = "all_to_all"

    @property
    def qasm(self):
        return (f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{self.width}];\n'
                + self.operations)

    def request(self):
        edges = ([[i, i + 1] for i in range(self.width - 1)]
                 if self.topology == "line" else
                 [[i, j] for i in range(self.width) for j in range(i + 1, self.width)])
        return {"qasm": self.qasm, "n_logical": self.width,
                "physical_width": self.width, "topology": self.topology, "edges": edges}


def rotation(axis, angle):
    """exp(-i angle sigma_axis / 2), directly from the Pauli identity."""
    pauli = {"x": [[0, 1], [1, 0]], "y": [[0, -1j], [1j, 0]],
             "z": [[1, 0], [0, -1]]}[axis]
    return (math.cos(angle / 2) * np.eye(2, dtype=complex)
            - 1j * math.sin(angle / 2) * np.array(pauli, dtype=complex))


def u_matrix(theta, phi, lam):
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array([[c, -np.exp(1j * lam) * s],
                     [np.exp(1j * phi) * s, np.exp(1j * (phi + lam)) * c]])


def controlled(matrix):
    """Control q0, target q1 in little endian order; inactive control is I."""
    result = np.eye(4, dtype=complex)
    result[np.ix_([1, 3], [1, 3])] = matrix
    return result


def repeat_case(name, width, operations, period_matrix, category, count=REPETITIONS,
                topology="all_to_all"):
    return NumericalCase(name, width, operations * count,
                         np.linalg.matrix_power(period_matrix, count), category, topology)


def numerical_cases():
    cases = []
    # Same-axis products have a closed-form oracle; no long matrix product is
    # needed to establish the original accumulated-angle counterexample.
    for axis in "xyz":
        for sign in (-1, 1):
            angle = sign * TINY_ANGLE
            cases.append(NumericalCase(
                f"tiny_r{axis}_{'positive' if sign > 0 else 'negative'}", 1,
                f"r{axis}({angle:.17g}) q[0];\n" * REPETITIONS,
                rotation(axis, angle * REPETITIONS), "tiny_rotations"))
        angle = 1.1e-9
        cases.append(NumericalCase(f"above_old_threshold_r{axis}", 1,
                                  f"r{axis}({angle:.17g}) q[0];\n" * REPETITIONS,
                                  rotation(axis, angle * REPETITIONS), "tiny_rotations"))

    for name, theta, phi, lam in (("theta", TINY_ANGLE, 0., 0.),
                                  ("phase", 0., TINY_ANGLE, 0.),
                                  ("two_phases", 0., TINY_ANGLE, TINY_ANGLE)):
        cases.append(repeat_case(f"tiny_u_{name}", 1,
                                 f"u3({theta:.17g},{phi:.17g},{lam:.17g}) q[0];\n",
                                 u_matrix(theta, phi, lam), "tiny_rotations"))

    # Each snapped angle is harmless in isolation at 1e-7, but a thousand
    # same-sign deviations are not. Test both sides of a period boundary.
    for axis, center in (("x", math.pi), ("y", math.pi),
                         ("z", math.pi / 2), ("z", math.pi / 4),
                         ("z", -math.pi / 2), ("z", -math.pi / 4),
                         ("z", math.pi), ("z", -math.pi)):
        angle = center + TINY_ANGLE
        cases.append(NumericalCase(f"special_r{axis}_{center:.6f}", 1,
                                  f"r{axis}({angle:.17g}) q[0];\n" * REPETITIONS,
                                  rotation(axis, angle * REPETITIONS), "special_angles"))
    for sign in (-1, 1):
        angle = 2 * math.pi + sign * TINY_ANGLE
        cases.append(NumericalCase(f"full_turn_offset_{sign}", 1,
                                  f"rx({angle:.17g}) q[0];\n" * REPETITIONS,
                                  rotation("x", angle * REPETITIONS), "special_angles"))
    hadamard = np.array([[1., 1.], [1., -1.]]) / math.sqrt(2)
    for name, theta, phi, lam in (("theta", math.pi / 2 + TINY_ANGLE, 0., math.pi),
                                  ("phi", math.pi / 2, TINY_ANGLE, math.pi),
                                  ("lambda", math.pi / 2, 0., math.pi + TINY_ANGLE)):
        cases.append(repeat_case(f"almost_hadamard_u_{name}", 1,
                                 f"u3({theta:.17g},{phi:.17g},{lam:.17g}) q[0];\nh q[0];\n",
                                 hadamard @ u_matrix(theta, phi, lam), "special_angles"))

    # Approximate inverse-pair removal must not silently erase the residual.
    for axis in "xyz":
        first, second = .37, -.37 + TINY_ANGLE
        cases.append(NumericalCase(f"almost_inverse_r{axis}", 1,
                                  (f"r{axis}({first:.17g}) q[0];\n"
                                   f"r{axis}({second:.17g}) q[0];\n") * REPETITIONS,
                                  rotation(axis, (first + second) * REPETITIONS),
                                  "inverse_pairs"))

    for axis in "xyz":
        first, second = .37, -.37 + TINY_ANGLE
        cases.append(NumericalCase(f"almost_inverse_cr{axis}", 2,
                                  (f"cr{axis}({first:.17g}) q[0],q[1];\n"
                                   f"cr{axis}({second:.17g}) q[0],q[1];\n") * REPETITIONS,
                                  controlled(rotation(axis, (first + second) * REPETITIONS)),
                                  "controlled_rotations"))
        cases.append(NumericalCase(f"controlled_full_turn_r{axis}", 2,
                                  f"cr{axis}({2 * math.pi:.17g}) q[0],q[1];\n",
                                  controlled(rotation(axis, 2 * math.pi)),
                                  "controlled_rotations"))
        cases.append(NumericalCase(f"tiny_cr{axis}", 2,
                                  f"cr{axis}({TINY_ANGLE:.17g}) q[0],q[1];\n" * REPETITIONS,
                                  controlled(rotation(axis, TINY_ANGLE * REPETITIONS)),
                                  "controlled_rotations"))
        for turns in (1, 2):
            angle = turns * 2 * math.pi + TINY_ANGLE
            count = REPETITIONS + 1  # Odd count retains the 2pi control phase.
            cases.append(NumericalCase(f"controlled_r{axis}_{turns}_turns_offset", 2,
                                      f"cr{axis}({angle:.17g}) q[0],q[1];\n" * count,
                                      controlled(rotation(axis, angle * count)),
                                      "controlled_rotations"))

    x, y = rotation("x", TINY_ANGLE), rotation("y", TINY_ANGLE)
    cases.append(repeat_case("noncommuting_xy", 1,
                             f"rx({TINY_ANGLE:.17g}) q[0];\nry({TINY_ANGLE:.17g}) q[0];\n",
                             y @ x, "interleavings"))
    cases.append(repeat_case("conjugated_tiny_rotation", 1,
                             f"rz(0.7) q[0];\nrx({TINY_ANGLE:.17g}) q[0];\nrz(-0.7) q[0];\n",
                             rotation("z", -.7) @ x @ rotation("z", .7), "interleavings"))
    cx = controlled(np.array([[0, 1], [1, 0]]))
    cases.append(repeat_case("entangled_tiny_rotation", 2,
                             f"cx q[0],q[1];\nry({TINY_ANGLE:.17g}) q[1];\ncx q[0],q[1];\n",
                             cx @ np.kron(y, np.eye(2)) @ cx, "interleavings"))

    # These gates are close to commuting with CX, but not exactly commuting.
    # Generic U parameters keep the small effect inside a nontrivial gate.
    control_u = u_matrix(TINY_ANGLE, .37, -.37)
    control_lift = np.kron(np.eye(2), control_u)
    cases.append(repeat_case("almost_diagonal_control_u", 2,
                             f"cx q[0],q[1];\nu3({TINY_ANGLE:.17g},0.37,-0.37) q[0];\ncx q[0],q[1];\n",
                             cx @ control_lift @ cx, "commutation"))
    target_u = u_matrix(.37, -math.pi / 2 + TINY_ANGLE, math.pi / 2)
    target_lift = np.kron(target_u, np.eye(2))
    cases.append(repeat_case("almost_x_axis_target_u", 2,
                             (f"cx q[0],q[1];\nu3(0.37,{-math.pi / 2 + TINY_ANGLE:.17g},"
                              f"{math.pi / 2:.17g}) q[1];\ncx q[0],q[1];\n"),
                             cx @ target_lift @ cx, "commutation"))
    # Echo away the large rotation so the noncommuting perturbation adds
    # coherently, rather than being hidden by its own repeated conjugation.
    for where, theta, phi, lam in (("control", TINY_ANGLE, .37, .21),
                                    ("target", 1.1, -math.pi / 2 + TINY_ANGLE,
                                     math.pi / 2)):
        q = 0 if where == "control" else 1
        matrix = u_matrix(theta, phi, lam)
        lift = np.kron(np.eye(2), matrix) if q == 0 else np.kron(matrix, np.eye(2))
        source = (f"cx q[0],q[1];\n"
                  f"u3({theta:.17g},{phi:.17g},{lam:.17g}) q[{q}];\n"
                  f"cx q[0],q[1];\n"
                  f"u3({-theta:.17g},{-lam:.17g},{-phi:.17g}) q[{q}];\n")
        cases.append(repeat_case(f"commutation_echo_{where}", 2, source,
                                 lift.conj().T @ cx @ lift @ cx, "commutation"))

    # A tiny rotation must survive across separate two-qubit blocks as well.
    unit = np.eye(8, dtype=complex)
    cx02 = unit.copy()
    for column in range(8):
        cx02[:, column] = unit[:, column ^ (4 if column & 1 else 0)]
    tiny2 = np.kron(y, np.eye(4))
    cases.append(repeat_case("routed_tiny_rotation", 3,
                             f"cx q[0],q[2];\nry({TINY_ANGLE:.17g}) q[2];\ncx q[0],q[2];\n",
                             cx02 @ tiny2 @ cx02, "interleavings", topology="line"))
    return cases


def phase_aligned_residual(expected, actual):
    """Maximum basis-column 2-norm with one phase for the complete operator."""
    overlap = np.vdot(expected, actual)
    phase = overlap / abs(overlap) if abs(overlap) else 1.
    return float(np.max(np.linalg.norm(actual - phase * expected, axis=0)))


def logical_output_matrix(case, result):
    width, gates = parse_qasm(result["out_qasm"])
    if width != case.width:
        raise AssertionError("these regressions require equal physical and logical widths")
    full = simulate(width, gates, np.eye(1 << width, dtype=complex))
    indices = []
    for key in ("initial_layout", "final_layout"):
        layout = result[key]
        if sorted(layout) != list(range(width)):
            raise AssertionError(f"invalid {key}: {layout}")
        indices.append([sum(((basis >> logical) & 1) << physical
                            for logical, physical in enumerate(layout))
                        for basis in range(1 << width)])
    initial, final = indices
    return full[np.ix_(final, initial)]


def compile_cases(binary, cases):
    """One bounded process per case avoids an unbounded hung worker test."""
    rows = []
    for case in cases:
        request = case.request()
        process = subprocess.run([str(binary)], input=json.dumps(request) + "\n",
                                 text=True, capture_output=True, timeout=30, check=True)
        lines = process.stdout.splitlines()
        if len(lines) != 2 or not json.loads(lines[0]).get("ready"):
            raise AssertionError(f"unexpected worker protocol: {process.stdout[:500]}")
        result = json.loads(lines[1])
        if result.get("status") != "compiled":
            raise AssertionError(f"{case.name}: {result}")
        output = logical_output_matrix(case, result)
        residual = phase_aligned_residual(case.expected, output)
        edges = [(a, b) for x, y in request["edges"] for a, b in ((x, y), (y, x))]
        validation = validate_output(case.qasm, result["out_qasm"], case.width, case.width,
                                     edges, result["initial_layout"], result["final_layout"])
        rows.append({"case": case.name, "category": case.category, "request": request,
                     "result": result, "analytic_max_column_residual": residual,
                     "validation": validation,
                     "passed": residual <= TOLERANCE and validation["passed"]})
    return rows


@unittest.skipIf(np is None, "optional NumPy dependency not installed")
class NumericalOracleTests(unittest.TestCase):
    def test_accumulated_rotation_is_detectably_nonidentity(self):
        expected = rotation("x", TINY_ANGLE * REPETITIONS)
        residual = phase_aligned_residual(expected, np.eye(2))
        self.assertGreater(residual, TOLERANCE)
        self.assertAlmostEqual(residual, 2 * math.sin(TINY_ANGLE * REPETITIONS / 4), places=18)

    def test_controlled_two_pi_is_not_global_phase(self):
        self.assertGreater(phase_aligned_residual(controlled(rotation("z", 2 * math.pi)),
                                                np.eye(4)), 1.)

    def test_closed_form_oracles_match_source_interpreter(self):
        for case in numerical_cases():
            with self.subTest(case=case.name):
                width, gates = parse_qasm(case.qasm)
                actual = simulate(width, gates, np.eye(1 << width))
                self.assertLess(phase_aligned_residual(case.expected, actual), 2e-11)


@unittest.skipIf(np is None, "optional NumPy dependency not installed")
@unittest.skipUnless(os.environ.get("QRUST_COMPARISON_BINARY"),
                     "set QRUST_COMPARISON_BINARY to a freshly built comparison_runner")
class PublicPipelineNumericalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        binary = Path(os.environ["QRUST_COMPARISON_BINARY"]).resolve()
        if not binary.is_file():
            raise AssertionError(f"configured comparison runner is missing: {binary}")
        cls.rows = compile_cases(binary, numerical_cases())

    def check_category(self, category):
        for row in self.rows:
            if row["category"] == category:
                with self.subTest(case=row["case"]):
                    self.assertTrue(row["validation"]["passed"], row["validation"])
                    self.assertLessEqual(row["analytic_max_column_residual"], TOLERANCE)

    def test_accumulated_tiny_rotations(self):
        self.check_category("tiny_rotations")

    def test_near_special_angles_and_period_boundaries(self):
        self.check_category("special_angles")

    def test_near_inverse_pairs(self):
        self.check_category("inverse_pairs")

    def test_controlled_rotations_preserve_relative_phase(self):
        self.check_category("controlled_rotations")

    def test_noncommuting_and_entangling_compositions(self):
        self.check_category("interleavings")

    def test_near_commutation_is_not_exact_commutation(self):
        self.check_category("commutation")


if __name__ == "__main__":
    unittest.main()
