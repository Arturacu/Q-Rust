"""Shared numerical evidence and hardware metrics for compiler comparisons.

This deliberately does not import Q-Rust, Qiskit, Cirq, or TKET. The small QASM
interpreter uses analytic standard-gate matrices and little-endian statevectors.
It checks the *same output* that was timed, outside the timed region. Numerical
equivalence evidence is distinct from output quality (CX count and depth).

Up to six physical qubits, every logical basis input is checked. Through sixteen
physical qubits, three adversarial basis inputs, a uniform superposition, and
four seeded Haar inputs are checked. This finite sample is not a proof. Both
paths compare full physical isometries with one shared global phase; dirty
ancilla amplitudes are never projected out or normalized away.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from functools import lru_cache
import math
import re
import time

import numpy as np


DEFAULT_SEED = 0x5152555354
TOLERANCE = 1e-7
FULL_MAX_PHYSICAL = 6
SAMPLED_MAX_PHYSICAL = 16


@dataclass(frozen=True)
class Gate:
    name: str
    params: tuple[float, ...]
    qubits: tuple[int, ...]


_SIGNATURES = {
    "id": (0, 1), "x": (0, 1), "y": (0, 1), "z": (0, 1),
    "h": (0, 1), "s": (0, 1), "sdg": (0, 1), "t": (0, 1),
    "tdg": (0, 1), "sx": (0, 1), "sxdg": (0, 1),
    "rx": (1, 1), "ry": (1, 1), "rz": (1, 1), "p": (1, 1),
    "u1": (1, 1), "u2": (2, 1), "u3": (3, 1),
    "cx": (0, 2), "cy": (0, 2), "cz": (0, 2), "ch": (0, 2),
    "swap": (0, 2), "crx": (1, 2), "cry": (1, 2), "crz": (1, 2),
    "cp": (1, 2), "cu1": (1, 2), "rzz": (1, 2), "ccx": (0, 3),
}


def _angle(expression: str) -> float:
    """A deliberately bounded arithmetic parser, never Python eval."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id == "pi":
            return math.pi
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right
        raise ValueError(f"unsupported parameter expression: {expression}")
    if len(expression) > 200:
        raise ValueError("overlong parameter expression")
    value = visit(ast.parse(expression, mode="eval").body)
    if not math.isfinite(value):
        raise ValueError("non-finite gate parameter")
    return value


@lru_cache(maxsize=256)
def parse_qasm(text: str) -> tuple[int, tuple[Gate, ...]]:
    """Read the static OpenQASM 2 subset used by this benchmark, fail closed."""
    registers = {}
    width = 0
    gates = []
    text = re.sub(r"//[^\n]*", "", text)
    for statement in text.split(";"):
        statement = statement.strip()
        if not statement:
            continue
        if re.fullmatch(r"OPENQASM\s+2\.0", statement):
            continue
        if re.fullmatch(r'include\s+"qelib1\.inc"', statement):
            continue
        match = re.fullmatch(r"qreg\s+([a-zA-Z_]\w*)\[(\d+)\]", statement)
        if match:
            name, size = match.group(1), int(match.group(2))
            if name in registers or size < 1 or gates:
                raise ValueError("invalid, repeated, or late quantum register")
            registers[name] = (width, size)
            width += size
            continue
        # Unused classical storage is harmless; measurement is not supported.
        if re.fullmatch(r"creg\s+[a-zA-Z_]\w*\[\d+\]", statement):
            continue
        match = re.fullmatch(r"([a-zA-Z]\w*)(?:\((.*)\))?\s+(.+)", statement)
        if not match:
            raise ValueError(f"unsupported QASM statement: {statement[:120]}")
        name, parameters, arguments = match.groups()
        name = name.lower()
        if name == "u":
            name = "u3"  # OpenQASM's built-in U has the same analytic matrix.
        if name not in _SIGNATURES:
            raise ValueError(f"unsupported gate: {name}")
        params = tuple(_angle(x.strip()) for x in parameters.split(",")) if parameters else ()
        qubits = []
        for argument in arguments.split(","):
            operand = re.fullmatch(r"([a-zA-Z_]\w*)\[(\d+)\]", argument.strip())
            if operand is None or operand.group(1) not in registers:
                raise ValueError(f"invalid gate operand: {argument}")
            offset, size = registers[operand.group(1)]
            q = int(operand.group(2))
            if q >= size:
                raise ValueError("out-of-range gate operand")
            qubits.append(offset + q)
        if (len(params), len(qubits)) != _SIGNATURES[name] or len(set(qubits)) != len(qubits):
            raise ValueError(f"invalid parameter/operand arity: {name}")
        gates.append(Gate(name, params, tuple(qubits)))
    if not width:
        raise ValueError("no quantum register")
    return width, tuple(gates)


def _one_matrix(name: str, params: tuple[float, ...]) -> tuple[complex, ...]:
    if name == "id":
        return (1, 0, 0, 1)
    if name == "x":
        return (0, 1, 1, 0)
    if name == "y":
        return (0, -1j, 1j, 0)
    if name == "z":
        return (1, 0, 0, -1)
    if name == "h":
        h = math.sqrt(0.5)
        return (h, h, h, -h)
    phases = {"s": math.pi / 2, "sdg": -math.pi / 2,
              "t": math.pi / 4, "tdg": -math.pi / 4}
    if name in phases:
        return (1, 0, 0, np.exp(1j * phases[name]))
    if name in ("sx", "sxdg"):
        sign = 1 if name == "sx" else -1
        return ((1 + sign * 1j) / 2, (1 - sign * 1j) / 2,
                (1 - sign * 1j) / 2, (1 + sign * 1j) / 2)
    if name in ("p", "u1"):
        return (1, 0, 0, np.exp(1j * params[0]))
    if name == "rz":
        return (np.exp(-0.5j * params[0]), 0, 0, np.exp(0.5j * params[0]))
    if name in ("rx", "ry"):
        c, s = math.cos(params[0] / 2), math.sin(params[0] / 2)
        return (c, -1j * s, -1j * s, c) if name == "rx" else (c, -s, s, c)
    theta, phi, lam = (math.pi / 2, *params) if name == "u2" else params
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return (c, -np.exp(1j * lam) * s, np.exp(1j * phi) * s, np.exp(1j * (phi + lam)) * c)


@lru_cache(maxsize=256)
def _controlled_pairs(width: int, target: int, controls: tuple[int, ...]):
    indices = np.arange(1 << width, dtype=np.int64)
    mask = sum(1 << control for control in controls)
    zero = indices[((indices & mask) == mask) & ((indices & (1 << target)) == 0)]
    return zero, zero | (1 << target)


def simulate(width: int, gates: tuple[Gate, ...], inputs: np.ndarray) -> np.ndarray:
    """Apply analytic matrices to several columns at once without a dense operator."""
    state = np.array(inputs, dtype=np.complex128, copy=True)
    if state.ndim != 2 or state.shape[0] != 1 << width:
        raise ValueError("invalid statevector batch shape")
    for gate in gates:
        name, params, qubits = gate.name, gate.params, gate.qubits
        if name == "rzz":
            indices = np.arange(1 << width, dtype=np.int64)
            parity = ((indices >> qubits[0]) ^ (indices >> qubits[1])) & 1
            state *= np.exp(-0.5j * params[0] * (1 - 2 * parity))[:, None]
            continue
        if name == "swap":
            first, second = qubits
            indices = np.arange(1 << width, dtype=np.int64)
            left = indices[((indices >> first) & 1 == 1) & ((indices >> second) & 1 == 0)]
            right = left ^ (1 << first) ^ (1 << second)
            state[left], state[right] = state[right].copy(), state[left].copy()
            continue
        if len(qubits) > 1:
            controls, target = qubits[:-1], qubits[-1]
            one_name = "x" if name == "ccx" else name[1:]
            if name == "cu1":
                one_name = "u1"
            a, b, c, d = _one_matrix(one_name, params)
            zero, one = _controlled_pairs(width, target, controls)
            first, second = state[zero].copy(), state[one].copy()
            state[zero], state[one] = a * first + b * second, c * first + d * second
        else:
            a, b, c, d = _one_matrix(name, params)
            view = state.reshape(-1, 2, 1 << qubits[0], state.shape[1])
            first, second = view[:, 0].copy(), view[:, 1].copy()
            view[:, 0], view[:, 1] = a * first + b * second, c * first + d * second
    return state


def _checked_layout(n_logical: int, physical_width: int, layout):
    if layout is None:
        if physical_width != n_logical:
            raise ValueError("layout metadata required for added physical qubits")
        layout = list(range(n_logical))
    if len(layout) not in (n_logical, physical_width):
        raise ValueError("layout must describe every logical qubit or the full physical permutation")
    if any(type(q) is not int or q < 0 or q >= physical_width for q in layout):
        raise ValueError("layout contains a noninteger or out-of-range index")
    if len(set(layout)) != len(layout):
        raise ValueError("layout is not injective")
    return layout


def _layout_indices(n_logical: int, physical_width: int, layout) -> np.ndarray:
    layout = _checked_layout(n_logical, physical_width, layout)
    logical = np.arange(1 << n_logical, dtype=np.int64)
    physical = np.zeros_like(logical)
    for logical_q, physical_q in enumerate(layout[:n_logical]):
        physical |= ((logical >> logical_q) & 1) << physical_q
    return physical


def common_metrics(width: int, gates: tuple[Gate, ...]) -> dict:
    """Unit-cost ASAP DAG depth; CX-only depth removes all one-qubit vertices."""
    frontier = [0] * width
    two_frontier = [0] * width
    for gate in gates:
        level = 1 + max(frontier[q] for q in gate.qubits)
        for q in gate.qubits:
            frontier[q] = level
        if len(gate.qubits) > 1:
            two_level = 1 + max(two_frontier[q] for q in gate.qubits)
            for q in gate.qubits:
                two_frontier[q] = two_level
    return {"gate_count": len(gates), "cx_count": sum(g.name == "cx" for g in gates),
            "depth": max(frontier, default=0), "two_q_depth": max(two_frontier, default=0)}


def validate_output(source_qasm: str, output_qasm: str, n_logical: int,
                    physical_width: int, edges, initial_layout=None, final_layout=None,
                    seed: int = DEFAULT_SEED, samples: int = 4) -> dict:
    """Return auditable evidence; failures are records, never silently discarded.

    ``edges`` is a directed list of allowed CX pairs, or None for all-to-all.
    Both orientations must be listed for bidirectional devices. Layout entries
    map logical indices to physical indices at entry and exit respectively.
    """
    start = time.perf_counter_ns()
    result = {"passed": False, "evidence_class": "invalid", "tolerance": TOLERANCE,
              "logical_qubits": n_logical, "physical_qubits": physical_width,
              "seed": seed, "haar_samples_requested": samples,
              "checked_input_count": 0, "oracle": "independent_numpy_analytic_qasm",
              "basis_legal": False, "coupling_legal": False}
    try:
        if type(n_logical) is not int or type(physical_width) is not int or not 0 < n_logical <= physical_width:
            raise ValueError("invalid logical/physical widths")
        if type(samples) is not int or samples < 1:
            raise ValueError("at least one Haar sample is required")
        source_width, source = parse_qasm(source_qasm)
        output_width, output = parse_qasm(output_qasm)
        if source_width != n_logical or output_width != physical_width:
            raise ValueError("declared width does not match QASM register width")
        result["metrics"] = common_metrics(physical_width, output)
        if any(g.name not in ("u3", "cx") for g in output):
            raise ValueError("output contains gates outside the common U3/CX basis")
        result["basis_legal"] = True
        if edges is not None:
            allowed = set()
            for edge in edges:
                if len(edge) != 2 or any(type(q) is not int or q < 0 or q >= physical_width for q in edge) or edge[0] == edge[1]:
                    raise ValueError("invalid coupling edge")
                allowed.add(tuple(edge))
            if any(g.name == "cx" and g.qubits not in allowed for g in output):
                raise ValueError("output contains a CX outside the directed coupling graph")
        result["coupling_legal"] = True
        initial_layout = _checked_layout(n_logical, physical_width, initial_layout)
        final_layout = _checked_layout(n_logical, physical_width, final_layout)
        if physical_width > SAMPLED_MAX_PHYSICAL:
            result.update(evidence_class="unchecked_width", passed=None)
            return result
        initial_indices = _layout_indices(n_logical, physical_width, initial_layout)
        final_indices = _layout_indices(n_logical, physical_width, final_layout)
        dim = 1 << n_logical
        if physical_width <= FULL_MAX_PHYSICAL:
            inputs = np.eye(dim, dtype=np.complex128)
            result.update(evidence_class="full_isometry", haar_samples_checked=0,
                          input_scheme="every computational basis vector",
                          fidelity_interpretation="full isometry overlap; numerical correctness evidence, not output quality")
        else:
            basis = sorted({0, dim - 1, sum(1 << q for q in range(0, n_logical, 2))})
            inputs = np.zeros((dim, len(basis) + 1 + samples), dtype=np.complex128)
            for column, index in enumerate(basis):
                inputs[index, column] = 1
            inputs[:, len(basis)] = 1 / math.sqrt(dim)
            rng = np.random.default_rng(seed)
            random = rng.standard_normal((dim, samples)) + 1j * rng.standard_normal((dim, samples))
            random /= np.linalg.norm(random, axis=0)
            inputs[:, len(basis) + 1:] = random
            result.update(evidence_class="sampled_isometry", haar_samples_checked=samples,
                          input_scheme="zero, all-one, alternating-bit basis; uniform superposition; seeded Haar",
                          fidelity_interpretation="finite sampled overlap; not process fidelity, a proof, or output quality")
        count = inputs.shape[1]
        actual_inputs = np.zeros((1 << physical_width, count), dtype=np.complex128)
        actual_inputs[initial_indices] = inputs
        actual = simulate(physical_width, output, actual_inputs)
        expected = np.zeros_like(actual)
        expected[final_indices] = simulate(n_logical, source, inputs)
        overlap = np.vdot(expected, actual)
        phase = overlap / abs(overlap) if abs(overlap) else 1
        residuals = np.linalg.norm(actual - phase * expected, axis=0)
        leakage = np.sum(np.abs(actual) ** 2, axis=0) - np.sum(np.abs(actual[final_indices]) ** 2, axis=0)
        residual = float(np.max(residuals))
        result.update(
            passed=bool(math.isfinite(residual) and residual < TOLERANCE),
            checked_input_count=count,
            normalized_frobenius_residual=float(np.linalg.norm(residuals) / math.sqrt(count)),
            max_input_residual=residual,
            isometry_fidelity=float(abs(overlap / count) ** 2),
            max_ancilla_leakage_probability=float(max(0, np.max(leakage))),
            shared_global_phase_real=float(np.real(phase)),
            shared_global_phase_imag=float(np.imag(phase)),
        )
    except (ValueError, TypeError, SyntaxError, ZeroDivisionError, OverflowError) as error:
        result.update(error=str(error), evidence_class="invalid", passed=False)
    finally:
        result["validation_ns"] = time.perf_counter_ns() - start
    return result
