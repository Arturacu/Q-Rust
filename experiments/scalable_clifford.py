"""Independent Clifford-isometry checker for large, static benchmark circuits.

The comparison uses normalized Choi states: one reference wire is Bell-paired
with each logical input, and unused physical inputs start in |0>. Equality of
the final signed stabilizer groups establishes equality of the projected
Clifford isometries up to one global phase, including unused physical outputs.

Arbitrary 1q matrices are recognized against the 24 Clifford representatives.
The recorded sum of local phase-aligned Frobenius residuals bounds the ideal
real-arithmetic operator error by telescoping. This is numerical recognition,
not a certified floating-point error bound. Non-Clifford gates fail closed.
No compiler SDK or dense many-qubit simulation is used.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math

import numpy as np

from comparison_validation import parse_qasm, _one_matrix


LOCAL_TOLERANCE = 1e-10
TOTAL_TOLERANCE = 1e-7
H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
S = np.diag([1, 1j])


def phase_residual(actual, expected):
    overlap = np.vdot(expected, actual)
    if abs(overlap) < 1e-14:
        return math.inf
    return float(np.linalg.norm(actual - overlap / abs(overlap) * expected))


@lru_cache(maxsize=1)
def representatives():
    found = [(np.eye(2, dtype=complex), ())]
    cursor = 0
    while cursor < len(found):
        matrix, word = found[cursor]
        cursor += 1
        for symbol, gate in (("h", H), ("s", S)):
            candidate = gate @ matrix
            if all(phase_residual(candidate, other) > 1e-12 for other, _ in found):
                found.append((candidate, (*word, symbol)))
    if len(found) != 24:
        raise RuntimeError("Clifford representative construction failed")
    return tuple(found)


@lru_cache(maxsize=8192)
def recognize(name, params):
    matrix = np.asarray(_one_matrix(name, params), dtype=complex).reshape(2, 2)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Non-finite local matrix")
    choices = [(phase_residual(matrix, target), word) for target, word in representatives()]
    residual, word = min(choices)
    if residual > LOCAL_TOLERANCE:
        raise ValueError(f"Non-Clifford local gate: {name}{params}; residual={residual:.6g}")
    return word, residual


def multiply(a, b):
    """Multiply commuting Hermitian Paulis represented by (X bits,Z bits,sign)."""
    x1, z1, r1 = a
    x2, z2, r2 = b
    x, z = x1 ^ x2, z1 ^ z2
    phase = (2 * (r1 + r2) + (x1 & z1).bit_count() + (x2 & z2).bit_count()
             + 2 * (z1 & x2).bit_count() - (x & z).bit_count()) % 4
    if phase & 1:
        raise ValueError("Anticommuting stabilizer rows")
    return x, z, phase // 2


def canonical(rows, width):
    rows = list(rows)
    pivot_row = 0
    for column in range(2 * width):
        bit = 1 << column
        index = next((j for j in range(pivot_row, len(rows))
                      if (rows[j][0] | (rows[j][1] << width)) & bit), None)
        if index is None:
            continue
        rows[pivot_row], rows[index] = rows[index], rows[pivot_row]
        for j in range(len(rows)):
            if j != pivot_row and (rows[j][0] | (rows[j][1] << width)) & bit:
                rows[j] = multiply(rows[j], rows[pivot_row])
        pivot_row += 1
        if pivot_row == len(rows):
            break
    if pivot_row != width or len(rows) != width:
        raise ValueError("Stabilizer rank is not a pure state")
    return tuple(rows)


class State:
    def __init__(self, width):
        self.width = width
        self.rows = [(0, 1 << i, 0) for i in range(width)]

    def apply(self, name, *qubits):
        if name in {"h", "s"}:
            bit = 1 << qubits[0]
            result = []
            for x, z, r in self.rows:
                xb, zb = bool(x & bit), bool(z & bit)
                r ^= xb & zb
                if name == "h":
                    if xb != zb:
                        x ^= bit
                        z ^= bit
                elif xb:
                    z ^= bit
                result.append((x, z, int(r)))
            self.rows = result
        elif name == "cx":
            cb, tb = 1 << qubits[0], 1 << qubits[1]
            result = []
            for x, z, r in self.rows:
                xc, xt, zc, zt = bool(x & cb), bool(x & tb), bool(z & cb), bool(z & tb)
                r ^= xc & zt & (xt ^ zc ^ True)
                if xc:
                    x ^= tb
                if zt:
                    z ^= cb
                result.append((x, z, int(r)))
            self.rows = result
        elif name == "cz":
            self.apply("h", qubits[1])
            self.apply("cx", *qubits)
            self.apply("h", qubits[1])
        elif name == "cy":
            for _ in range(3):
                self.apply("s", qubits[1])
            self.apply("cx", *qubits)
            self.apply("s", qubits[1])
        elif name == "swap":
            self.apply("cx", *qubits)
            self.apply("cx", qubits[1], qubits[0])
            self.apply("cx", *qubits)
        else:
            raise ValueError("Unsupported symbolic Clifford: " + name)


def _map(mapping, logical, physical):
    if (len(mapping) != logical or any(type(q) is not int or not 0 <= q < physical for q in mapping)
            or len(set(mapping)) != logical):
        raise ValueError("Invalid logical-to-physical map")


def _choi(gates, logical, physical, bell_map, gate_map):
    state = State(logical + physical)
    for q, location in enumerate(bell_map):
        state.apply("h", q)
        state.apply("cx", q, logical + location)
    errors = []
    recognized = 0
    for gate in gates:
        qubits = tuple(logical + gate_map[q] for q in gate.qubits)
        if len(qubits) == 1:
            word, residual = recognize(gate.name, gate.params)
            for symbol in word:
                state.apply(symbol, *qubits)
            errors.append(residual)
            recognized += 1
        elif gate.name in {"cx", "cy", "cz", "swap"} and not gate.params:
            state.apply(gate.name, *qubits)
        else:
            raise ValueError("Unsupported/non-Clifford gate: " + gate.name)
    return canonical(state.rows, state.width), {
        "recognized_one_qubit_gates": recognized,
        "max_local_residual": max(errors, default=0),
        "local_residual_sum": math.fsum(errors),
    }


def validate_clifford(source_qasm, output_qasm, n_logical, physical_width,
                      initial_layout, final_layout, edges=None,
                      output_basis=("u3", "cx"), directed_edges=False):
    """Return fail-closed evidence for the emitted output and its two layouts."""
    result = {"passed": False, "evidence_class": "invalid_or_non_clifford",
              "all_logical_inputs": False, "n_logical": n_logical,
              "physical_width": physical_width, "local_tolerance": LOCAL_TOLERANCE,
              "total_tolerance": TOTAL_TOLERANCE,
              "recognition_is_numerical": True,
              "coupling_contract": "directed" if directed_edges else "bidirectional",
              "source_sha256": hashlib.sha256(source_qasm.encode()).hexdigest(),
              "output_sha256": hashlib.sha256(output_qasm.encode()).hexdigest()}
    try:
        if type(n_logical) is not int or type(physical_width) is not int or not 0 < n_logical <= physical_width:
            raise ValueError("Invalid widths")
        _map(initial_layout, n_logical, physical_width)
        _map(final_layout, n_logical, physical_width)
        source_width, source = parse_qasm(source_qasm)
        output_width, output = parse_qasm(output_qasm)
        if source_width != n_logical or output_width != physical_width:
            raise ValueError("QASM width mismatch")
        allowed = None
        if edges is not None:
            allowed = set()
            for edge in edges:
                if (len(edge) != 2 or any(type(q) is not int or not 0 <= q < physical_width for q in edge)
                        or edge[0] == edge[1]):
                    raise ValueError("Invalid coupling edge")
                allowed.add(tuple(edge) if directed_edges else tuple(sorted(edge)))
        for gate in output:
            if output_basis is not None and gate.name not in output_basis:
                raise ValueError("Output outside declared basis: " + gate.name)
            if len(gate.qubits) > 1 and allowed is not None:
                if directed_edges and gate.name in {"cx", "cy"}:
                    legal = gate.qubits in allowed
                else:
                    legal = (tuple(sorted(gate.qubits)) in allowed if not directed_edges
                             else gate.qubits in allowed or tuple(reversed(gate.qubits)) in allowed)
                if not legal:
                    raise ValueError("Output violates coupling graph")
        expected, input_residuals = _choi(source, n_logical, physical_width, final_layout, final_layout)
        actual, output_residuals = _choi(output, n_logical, physical_width, initial_layout, list(range(physical_width)))
        error_sum = input_residuals['local_residual_sum'] + output_residuals['local_residual_sum']
        result.update(input_recognition=input_residuals, output_recognition=output_residuals,
                      evidence_class="clifford_choi_isometry", all_logical_inputs=True,
                      reconstruction_error_sum=error_sum,
                      tableau_equal=actual == expected,
                      choi_qubits=n_logical + physical_width,
                      initial_layout=initial_layout, final_layout=final_layout,
                      expected_tableau_sha256=hashlib.sha256(json.dumps(expected).encode()).hexdigest(),
                      actual_tableau_sha256=hashlib.sha256(json.dumps(actual).encode()).hexdigest())
        result['passed'] = actual == expected and error_sum <= TOTAL_TOLERANCE
        if not result['passed']:
            result['error'] = "Clifford isometries differ" if actual != expected else "Recognition residual budget exceeded"
    except (ValueError, TypeError, IndexError, OverflowError, SyntaxError, ZeroDivisionError) as error:
        result['error'] = str(error)
    return result
