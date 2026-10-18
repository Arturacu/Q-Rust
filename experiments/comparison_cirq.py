"""Cirq adapter for the repeated, common-target compiler comparison.

This is a specified Cirq transformer pipeline, not a purported Cirq "O3" mode.
CZ and CX differ by conjugation with a Hadamard on the target, so full CZ
synthesis followed by that conjugation reaches the common {arbitrary 1q, CX}
target without changing the entangling-gate count. All synthesis, routing,
decomposition, and fusion are inside the compilation timer. Parsing and
one-for-one output serialization are outside it.

Primary API references:
https://quantumai.google/reference/python/cirq/CZTargetGateset
https://quantumai.google/reference/python/cirq/optimize_for_target_gateset
https://quantumai.google/reference/python/cirq/RouteCQC

The adapter accepts the study's unitary OpenQASM 2 subset with one q register.
The caller performs the SAME CRZ/RZZ normalization for every compiler.
"""
from __future__ import annotations

import math
import re
import time


SYNTHESIS_ATOL = 1e-10
MAX_SYNTHESIS_PASSES = 2


def _complete_initial_layout(initial_map, logical_qubits, physical_qubits):
    """Keep router placement and allocate dropped/inactive logical wires.

    Circuit optimizers can remove every operation on a wire. Its logical state
    still exists, so the reported layout must include it. Prefer its same-index
    physical wire where available and otherwise take the first unused wire.
    """
    completed = dict(initial_map)
    available = set(physical_qubits) - set(completed.values())
    for index, logical in enumerate(logical_qubits):
        if logical in completed:
            continue
        preferred = physical_qubits[index]
        physical = preferred if preferred in available else min(available)
        completed[logical] = physical
        available.remove(physical)
    return completed


def _to_cx(circuit, cirq):
    """Perform the CZ-to-CX lowering and resulting one-qubit fusion.

    CZTargetGateset can retain a full CZ represented by another gate type (for
    example the QASM importer's controlled U1). Recognize that case by its
    two-qubit unitary, up to ONE global phase, rather than by its Python class.
    This bounded canonicalization is part of compilation, never serialization.
    """
    cz_unitary = None

    def lower(op, _moment_index):
        nonlocal cz_unitary
        if len(op.qubits) != 2 or op.gate == cirq.CX:
            return op
        if isinstance(op.gate, cirq.CZPowGate):
            # Bound matrix error by |delta_angle|. Disable math.isclose's
            # default relative tolerance: it is larger than our absolute
            # synthesis tolerance and could silently round a fractional CZ.
            if not math.isclose(float(op.gate.exponent) % 2.0, 1.0,
                                rel_tol=0.0, abs_tol=SYNTHESIS_ATOL / math.pi):
                raise ValueError(f"Cirq emitted a fractional CZ: {op}")
        else:
            matrix = cirq.unitary(op, default=None)
            if cz_unitary is None:
                cz_unitary = cirq.unitary(cirq.CZ)
            if matrix is None or not cirq.linalg.allclose_up_to_global_phase(
                    matrix, cz_unitary, atol=SYNTHESIS_ATOL, rtol=0.0):
                raise ValueError(f"Cirq emitted a two-qubit operation outside full CZ/CX: {op}")
        control, target = op.qubits
        return [cirq.H(target), cirq.CX(control, target), cirq.H(target)]

    lowered = cirq.map_operations_and_unroll(circuit, lower)
    return cirq.merge_single_qubit_gates_to_phxz(lowered, atol=SYNTHESIS_ATOL)


def _serialize(circuit, physical_width, cirq):
    """Serialize one output gate per compiled gate, with no hidden lowering.

    A one-qubit unitary is expressed as u3 up to global phase; the study checks
    circuit equivalence up to a single global phase. Explicit global phases
    have no quantum gate cost and are omitted. Neither fusion nor any change
    to the two-qubit gate sequence is permitted in this function.
    """
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";',
             f"qreg q[{physical_width}];"]
    for op in circuit.all_operations():
        width = len(op.qubits)
        if width == 0:
            if not isinstance(op.gate, cirq.GlobalPhaseGate):
                raise ValueError(f"Unsupported zero-qubit operation: {op}")
        elif width == 1:
            # Euler-angle conversion is representational only: every
            # fused gate remains exactly one arbitrary single-qubit gate.
            lam, theta, phi = cirq.linalg.deconstruct_single_qubit_matrix_into_angles(cirq.unitary(op))
            if not all(math.isfinite(x) for x in (theta, phi, lam)):
                raise ValueError("Non-finite angle in Cirq output")
            lines.append(f"u3({theta:.17g},{phi:.17g},{lam:.17g}) "
                         f"q[{op.qubits[0].x}];")
        elif width == 2 and op.gate == cirq.CX:
            lines.append(f"cx q[{op.qubits[0].x}],q[{op.qubits[1].x}];")
        else:
            raise ValueError(f"Cirq output is not closed over arbitrary 1q + CX: {op}")
    return "\n".join(lines) + "\n"


def compile_one(qasm, n_logical, physical_width, edges, topology, seed):
    """Compile one trial and return the exact serialized result and layouts.

    ``seed`` is recorded for cross-tool trial pairing, but this pipeline exposes
    no random seed. Repeated Cirq trials measure timing variation, not independent
    randomized searches. The orchestrator fixes PYTHONHASHSEED across workers.
    """
    import cirq
    import networkx as nx
    from cirq.contrib.qasm_import import circuit_from_qasm

    if not 0 < n_logical <= physical_width:
        raise ValueError("Require 0 < logical width <= physical width")
    parse_start = time.perf_counter_ns()
    parsed = circuit_from_qasm(qasm)
    remap = {}
    for qubit in parsed.all_qubits():
        match = re.fullmatch(r"q_(\d+)", qubit.name)
        if match is None or int(match.group(1)) >= n_logical:
            raise ValueError(f"Unexpected input qubit: {qubit}")
        remap[qubit] = cirq.LineQubit(int(match.group(1)))
    parsed = parsed.transform_qubits(remap)
    if any(not cirq.has_unitary(op) for op in parsed.all_operations()):
        raise ValueError("Comparison adapter accepts only unitary circuits")
    parse_ns = time.perf_counter_ns() - parse_start

    compile_start = time.perf_counter_ns()
    circuit = parsed.copy()
    logical_qubits = cirq.LineQubit.range(n_logical)
    physical_qubits = cirq.LineQubit.range(physical_width)
    gateset = cirq.CZTargetGateset(
        atol=SYNTHESIS_ATOL,
        allow_partial_czs=False,
        preserve_moment_structure=False,
        reorder_operations=True,
    )
    circuit = cirq.optimize_for_target_gateset(
        circuit, gateset=gateset, max_num_passes=MAX_SYNTHESIS_PASSES,
        ignore_failures=False,
    )
    if topology == "all_to_all":
        initial_map = dict(zip(logical_qubits, physical_qubits))
        swap_map = {}
    else:
        graph = nx.Graph()
        graph.add_nodes_from(physical_qubits)
        for a, b in edges:
            if not (0 <= a < physical_width and 0 <= b < physical_width and a != b):
                raise ValueError(f"Invalid physical edge: {(a, b)}")
            graph.add_edge(physical_qubits[a], physical_qubits[b])
        if circuit.all_qubits():
            circuit, initial_map, swap_map = cirq.RouteCQC(graph).route_circuit(circuit)
        else:
            initial_map, swap_map = {}, {}
        initial_map = _complete_initial_layout(initial_map, logical_qubits, physical_qubits)
        # Synthesize SWAPs together with neighboring same-pair gates; this can
        # reduce routing overhead while preserving the router's final mapping.
        circuit = cirq.optimize_for_target_gateset(
            circuit, gateset=gateset, max_num_passes=MAX_SYNTHESIS_PASSES,
            ignore_failures=False,
        )
    circuit = _to_cx(circuit, cirq)
    compile_ns = time.perf_counter_ns() - compile_start

    initial_layout = [initial_map[qubit].x for qubit in logical_qubits]
    final_layout = [swap_map.get(initial_map[qubit], initial_map[qubit]).x
                    for qubit in logical_qubits]
    return {
        "compile_ns": compile_ns,
        "parse_ns": parse_ns,
        "out_qasm": _serialize(circuit, physical_width, cirq),
        "initial_layout": initial_layout,
        "final_layout": final_layout,
        "output_qubits": physical_width,
        "seed": seed,
        "seed_used": False,
        "pipeline": {
            "name": "Cirq CZ synthesis / RouteCQC / CZ synthesis / CX closure",
            "target": "arbitrary 1q + CX",
            "synthesis_atol": SYNTHESIS_ATOL,
            "max_synthesis_passes": MAX_SYNTHESIS_PASSES,
            "preserve_moment_structure": False,
            "reorder_operations": True,
            "routing": "omitted" if topology == "all_to_all" else "RouteCQC",
            "initial_mapper": "identity" if topology == "all_to_all" else "LineInitialMapper",
            "routing_lookahead_radius": 8,
            "randomized_search": False,
            "timing": "copy, pass/target construction, synthesis, mapping, routing, lowering, fusion",
        },
    }
