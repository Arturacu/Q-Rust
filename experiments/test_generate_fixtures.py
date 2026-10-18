"""Algorithm-level checks independent of Q-Rust's simulator and decompositions.

Run with: python3 -m unittest discover -s experiments -p 'test_generate_fixtures.py'

The small statevector interpreter below uses analytic gate matrices and a
little-endian register: q[0] is the least significant bit of the basis index.
The Fourier oracle is the closed-form DFT, not another QFT gate construction.
Only Python's standard library is needed.
"""
from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
import re
import unittest

from generate_fixtures import FIXTURE_DIR, make_qft, make_wstate


_GATE = re.compile(r"([a-z][a-z0-9]*)(?:\(([^)]+)\))?\s+(.+);")
_QUBIT = re.compile(r"q\[(\d+)\]")


def assert_fixture_matches(actual: str, expected: str) -> None:
    """Require the same instructions, allowing only tiny libm angle rounding.

    Python delegates acos/sqrt to platform libm, so regenerated W-state angles
    need not have identical final decimal digits on macOS and Linux. Four ULPs
    allow a few final-bit differences in these composed calculations without a
    broad absolute tolerance: around a one-radian angle this is below 1e-15,
    and the allowance shrinks with the small QFT angles. This checks fixture
    synchronization only; the independent analytic state/matrix tests below
    still establish the intended algorithm-level semantics.
    """
    actual_lines, expected_lines = actual.splitlines(), expected.splitlines()
    if len(actual_lines) != len(expected_lines):
        raise AssertionError("fixture instruction count differs")
    for line_number, (left, right) in enumerate(zip(actual_lines, expected_lines), 1):
        parsed_left, parsed_right = _GATE.fullmatch(left), _GATE.fullmatch(right)
        if parsed_left is None or parsed_right is None:
            if left != right:
                raise AssertionError(f"fixture declaration/structure differs at line {line_number}")
            continue
        left_gate, left_params, left_operands = parsed_left.groups()
        right_gate, right_params, right_operands = parsed_right.groups()
        if (left_gate, left_operands) != (right_gate, right_operands):
            raise AssertionError(f"fixture gate/order/operands differ at line {line_number}")
        left_values = [float(p) for p in left_params.split(",")] if left_params is not None else []
        right_values = [float(p) for p in right_params.split(",")] if right_params is not None else []
        if len(left_values) != len(right_values):
            raise AssertionError(f"fixture parameter arity differs at line {line_number}")
        for parameter, (a, b) in enumerate(zip(left_values, right_values)):
            if not (math.isfinite(a) and math.isfinite(b)):
                raise AssertionError(f"non-finite fixture angle at line {line_number}")
            if abs(a - b) > 4 * max(math.ulp(a), math.ulp(b)):
                raise AssertionError(
                    f"fixture angle differs by more than four ULPs at line {line_number}, "
                    f"parameter {parameter}: {a!r} != {b!r}"
                )


def simulate(qasm: str, n: int, basis: int = 0) -> list[complex]:
    """Apply the fixture's QASM gates to one basis vector, without Q-Rust."""
    state = [0j] * (1 << n)
    state[basis] = 1.0 + 0j
    for source in qasm.splitlines():
        line = source.partition("//")[0].strip()
        if not line or line.startswith(("OPENQASM ", "include ", "qreg ")):
            continue
        parsed = _GATE.fullmatch(line)
        if parsed is None:
            raise ValueError(f"unsupported fixture syntax: {line}")
        gate, parameter, operands = parsed.groups()
        qubits = [int(q) for q in _QUBIT.findall(operands)]
        parameters = [float(value) for value in parameter.split(",")] if parameter is not None else []
        theta = parameters[0] if parameters else None
        if gate in {"h", "x", "ry", "rz", "u3"}:
            if len(qubits) != 1:
                raise ValueError(f"wrong arity: {line}")
            q = qubits[0]
            if gate == "h":
                a, b, c, d = (1 / math.sqrt(2),) * 3 + (-1 / math.sqrt(2),)
            elif gate == "x":
                a, b, c, d = 0, 1, 1, 0
            elif gate == "ry":
                a = d = math.cos(theta / 2)
                c = math.sin(theta / 2)
                b = -c
            elif gate == "rz":
                a, b, c, d = cmath.exp(-0.5j * theta), 0, 0, cmath.exp(0.5j * theta)
            else:
                theta, phi, lam = parameters
                a = math.cos(theta / 2)
                b = -cmath.exp(1j * lam) * math.sin(theta / 2)
                c = cmath.exp(1j * phi) * math.sin(theta / 2)
                d = cmath.exp(1j * (phi + lam)) * math.cos(theta / 2)
            for index in range(len(state)):
                if not (index >> q) & 1:
                    partner = index | (1 << q)
                    v0, v1 = state[index], state[partner]
                    state[index], state[partner] = a * v0 + b * v1, c * v0 + d * v1
        elif gate in {"cx", "crz", "swap"}:
            if len(qubits) != 2 or qubits[0] == qubits[1]:
                raise ValueError(f"invalid two-qubit gate: {line}")
            control, target = qubits
            if gate == "crz":
                for index in range(len(state)):
                    if (index >> control) & 1:
                        sign = 1 if (index >> target) & 1 else -1
                        state[index] *= cmath.exp(0.5j * theta * sign)
            elif gate == "cx":
                for index in range(len(state)):
                    if (index >> control) & 1 and not (index >> target) & 1:
                        partner = index ^ (1 << target)
                        state[index], state[partner] = state[partner], state[index]
            else:
                for index in range(len(state)):
                    if (index >> control) & 1 and not (index >> target) & 1:
                        partner = index ^ (1 << control) ^ (1 << target)
                        state[index], state[partner] = state[partner], state[index]
        else:
            raise ValueError(f"unsupported gate in independent fixture simulator: {gate}")
    return state


class AlgorithmFixtureTests(unittest.TestCase):
    def assert_fourier_transform(self, qasm: str, n: int) -> None:
        """Compare every matrix element using ONE global phase for all inputs."""
        dim = 1 << n
        phase = simulate(qasm, n)[0] * math.sqrt(dim)
        self.assertAlmostEqual(abs(phase), 1.0, delta=1e-11)
        for x in range(dim):
            state = simulate(qasm, n, basis=x)
            for y, actual in enumerate(state):
                expected = phase * cmath.exp(2j * math.pi * x * y / dim) / math.sqrt(dim)
                self.assertLess(
                    abs(actual - expected), 1e-11,
                    f"QFT-{n}, input {x}, output {y}: {actual} != {expected}",
                )

    def assert_w_state(self, qasm: str, n: int) -> None:
        state = simulate(qasm, n)
        for index, actual in enumerate(state):
            expected = 1 / math.sqrt(n) if index and index & (index - 1) == 0 else 0.0
            self.assertLess(abs(actual - expected), 1e-11, f"W-{n}, basis {index}")
        self.assertAlmostEqual(sum(abs(amplitude) ** 2 for amplitude in state), 1.0, delta=1e-11)

    def test_qft_matches_closed_form_fourier_matrix(self) -> None:
        for n in range(1, 7):
            with self.subTest(n=n):
                self.assert_fourier_transform(make_qft(n).qasm(), n)

    def test_w_state_has_equal_positive_single_excitation_amplitudes(self) -> None:
        for n in range(1, 7):
            with self.subTest(n=n):
                self.assert_w_state(make_wstate(n).qasm(), n)

    def test_stored_qft_fixtures_match_closed_form_fourier_matrix(self) -> None:
        for n in (3, 5):
            with self.subTest(n=n):
                self.assert_fourier_transform((Path(FIXTURE_DIR) / f"qft_{n}.qasm").read_text(), n)

    def test_rust_e2e_qft_fixtures_are_correct_fourier_transforms(self) -> None:
        fixture_dir = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
        for n in (3, 5):
            with self.subTest(n=n):
                qasm = (fixture_dir / f"qft_{n}.qasm").read_text()
                self.assert_fourier_transform(qasm, n)
                assert_fixture_matches(qasm, make_qft(n).qasm())

    def test_stored_w_state_fixtures_have_no_vacuum_or_multiple_excitations(self) -> None:
        for n in (5, 7, 10, 12, 14):
            with self.subTest(n=n):
                self.assert_w_state((Path(FIXTURE_DIR) / f"wstate_{n}.qasm").read_text(), n)

    def test_stored_qft_and_w_state_files_and_metadata_match_generator(self) -> None:
        for metadata in json.loads((Path(FIXTURE_DIR) / "manifest.json").read_text()):
            if metadata["family"] not in {"qft", "wstate"}:
                continue
            generator = make_qft if metadata["family"] == "qft" else make_wstate
            fixture = generator(metadata["nominal_n"])
            with self.subTest(fixture=fixture.name):
                assert_fixture_matches((Path(FIXTURE_DIR) / f"{fixture.name}.qasm").read_text(), fixture.qasm())
                for field in ("source_gates", "source_2q_gates", "source_cx_gates", "num_qubits", "notes"):
                    self.assertEqual(metadata[field], getattr(fixture, field), field)

    def test_qft_and_w_state_require_positive_register_width(self) -> None:
        for generator in (make_qft, make_wstate):
            for n in (0, -1):
                with self.subTest(generator=generator.__name__, n=n):
                    with self.assertRaises(ValueError):
                        generator(n)


class FixtureSynchronizationTests(unittest.TestCase):
    @staticmethod
    def fixture(angle: float) -> str:
        return f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\nry({angle:.17g}) q[1];\ncx q[0],q[1];\n'

    def test_accepts_a_few_ulp_of_platform_angle_rounding(self) -> None:
        for angle in (1.1071487177940904, -1.1071487177940904, 0.000003):
            for direction in (-math.inf, math.inf):
                changed = angle
                for steps in range(5):
                    with self.subTest(angle=angle, direction=direction, ulps=steps):
                        assert_fixture_matches(self.fixture(changed), self.fixture(angle))
                    changed = math.nextafter(changed, direction)

    def test_rejects_angle_changes_beyond_the_rounding_allowance(self) -> None:
        angle = 1.1071487177940904
        changed = angle
        for _ in range(8):
            changed = math.nextafter(changed, math.inf)
        for wrong in (changed, angle + 1e-10, angle + 0.01):
            with self.subTest(angle=wrong), self.assertRaises(AssertionError):
                assert_fixture_matches(self.fixture(wrong), self.fixture(angle))

    def test_rejects_gate_order_arity_operands_and_declaration_changes(self) -> None:
        original = self.fixture(1.0)
        changes = [
            original.replace("ry(1)", "rz(1)"),
            original.replace("ry(1)", "ry(1,0)"),
            original.replace("ry(1) q[1]", "ry(1) q[0]"),
            original.replace("cx q[0],q[1]", "cx q[1],q[0]"),
            original.replace("cx q[0],q[1]", "cx q[0],q[1],q[0]"),
            original.replace("ry(1) q[1];\ncx q[0],q[1];", "cx q[0],q[1];\nry(1) q[1];"),
            original.replace("cx q[0],q[1];\n", ""),
            original.replace("qreg q[2]", "qreg q[3]"),
        ]
        for changed in changes:
            with self.subTest(changed=changed), self.assertRaises(AssertionError):
                assert_fixture_matches(changed, original)

    def test_rejects_nonfinite_angles_even_when_both_strings_match(self) -> None:
        for angle in (math.nan, math.inf, -math.inf):
            with self.subTest(angle=angle), self.assertRaises(AssertionError):
                assert_fixture_matches(self.fixture(angle), self.fixture(angle))


if __name__ == "__main__":
    unittest.main()
