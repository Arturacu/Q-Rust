#!/usr/bin/env python3
"""Cross-tool baseline benchmark: Qiskit, Cirq, and tket vs. Q-Rust.

Transpiles the *same* fixtures (experiments/fixtures/) onto the *same* coupling
maps used by the Q-Rust runner, targeting the *same* logical level ({1-qubit
rotations, CX}), so CX count / depth / compile time are directly comparable
head-to-head with the qrust numbers in results/raw_benchmark_data.json.

Methodology notes:
  * Qiskit and tket receive the original OpenQASM fixture. Cirq receives a
    semantics-preserving normalized form in which `crz(theta)` and `rzz(theta)`
    are expanded into their exact, CX-optimal 2-CX decompositions because its
    OpenQASM importer does not accept those instructions. Q-Rust consumes the
    original form and lowers the same gates to two CX.
  * Multi-qubit gates (CCX) are decomposed to <=2-qubit gates *before* routing
    for every tool (qiskit/tket do this internally; for cirq we do it
    explicitly). This is the same discipline that the HighArityDecompositionPass
    enforces for Q-Rust.
  * Comparison is at each tool's best effort (qiskit optimization_level=3, tket
    FullPeepholeOptimise, cirq RouteCQC + single-qubit fusion), followed by
    rebasing to arbitrary one-qubit gates plus CX, i.e. the "O3" column.

Output: results/baselines.json  (one record per tool x fixture x topology)

Usage:
    python3 experiments/bench_baselines.py
"""
from __future__ import annotations

import json
import math
import os
import re
import statistics
import sys
import time
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")
RESULTS = os.path.join(HERE, "results")
# Gate-count, depth, and wall-time statistics use five independent timed
# repetitions for every tool. Qiskit varies the transpiler seed across trials;
# tket and Cirq are deterministic in circuit output, but wall-clock time still
# varies with the host, so their timing medians also require replication.
TRIALS_QISKIT = 5
TRIALS_TKET = 5
TRIALS_CIRQ = 5

# Coupling maps mirror src/backend.rs exactly (undirected unique edges).
QUITO = [(0, 1), (1, 2), (1, 3), (3, 4)]
NAIROBI = [(0, 1), (1, 2), (1, 3), (3, 5), (4, 5), (5, 6)]


def grid_edges(rows, cols):
    e = []
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            if c + 1 < cols:
                e.append((i, i + 1))
            if r + 1 < rows:
                e.append((i, i + cols))
    return e


def all_to_all_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def linear_edges(n):
    return [(i, i + 1) for i in range(n - 1)]


def ring_edges(n):
    return [(i, (i + 1) % n) for i in range(n)] if n > 2 else linear_edges(n)


def star_edges(n):
    return [(0, i) for i in range(1, n)]


def topologies_for(n):
    """Return {name: (num_qubits, undirected_edges)} that can host `n` qubits.

    Sized topologies (linear/ring/star/all_to_all) are built at the circuit's
    own width; the fixed IBM/grid devices are included only when they fit.
    """
    out = {}
    if n <= 5:
        out["ibm_quito"] = (5, QUITO)
    if n <= 7:
        out["ibm_nairobi"] = (7, NAIROBI)
    if n <= 16:
        out["grid_4x4"] = (16, grid_edges(4, 4))
    out["all_to_all"] = (n, all_to_all_edges(n))
    out["linear"] = (n, linear_edges(n))
    out["ring"] = (n, ring_edges(n))
    out["star"] = (n, star_edges(n))
    return out


# --------------------------------------------------------------------------- #
# QASM normalization: expand crz -> rz/cx/rz/cx (exact, 2-CX optimal).         #
# --------------------------------------------------------------------------- #
_CRZ = re.compile(r'crz\(([^)]+)\)\s*q\[(\d+)\]\s*,\s*q\[(\d+)\]\s*;')
_RZZ = re.compile(r'rzz\(([^)]+)\)\s*q\[(\d+)\]\s*,\s*q\[(\d+)\]\s*;')


def _eval_angle(expr: str) -> float:
    return float(eval(expr, {"pi": math.pi, "__builtins__": {}}, {}))


def normalize_qasm(text: str) -> str:
    """Expand gates the Cirq QASM importer does not recognise into their exact,
    CX-optimal decompositions. crz -> 2 CX; rzz -> 2 CX. (qiskit and tket
    receive the raw QASM and handle these natively.)"""
    def crz_repl(m):
        th = _eval_angle(m.group(1))
        c, t = m.group(2), m.group(3)
        return (f"rz({th / 2}) q[{t}];\ncx q[{c}],q[{t}];\n"
                f"rz({-th / 2}) q[{t}];\ncx q[{c}],q[{t}];")

    def rzz_repl(m):
        th = _eval_angle(m.group(1))
        a, b = m.group(2), m.group(3)
        return (f"cx q[{a}],q[{b}];\nrz({th}) q[{b}];\ncx q[{a}],q[{b}];")

    return _RZZ.sub(rzz_repl, _CRZ.sub(crz_repl, text))


def median_record(times_ms, gate_counts, cx_counts, depths):
    """Return per-metric medians and preserve samples for auditability."""
    return {
        "gate_count": int(statistics.median(gate_counts)),
        "cx_count": int(statistics.median(cx_counts)),
        "depth": int(statistics.median(depths)),
        "compile_ms_median": statistics.median(times_ms),
        "compile_ms_min": min(times_ms), "compile_ms_max": max(times_ms),
        "compile_ms_samples": times_ms,
        "gate_count_samples": gate_counts,
        "cx_count_samples": cx_counts,
        "depth_samples": depths,
    }


# --------------------------------------------------------------------------- #
# Qiskit                                                                       #
# --------------------------------------------------------------------------- #
def run_qiskit(qasm_raw, edges, opt=3):
    # LEGACY_CUSTOM_INSTRUCTIONS makes the full qelib1 set (crz, swap, ...)
    # available so we can load the original fixtures unmodified.
    from qiskit import qasm2, transpile
    from qiskit.transpiler import CouplingMap
    qc = qasm2.loads(qasm_raw, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    directed = [[a, b] for (a, b) in edges] + [[b, a] for (a, b) in edges]
    cm = CouplingMap(couplinglist=directed)
    times, gate_counts, cx_counts, depths = [], [], [], []
    for trial in range(TRIALS_QISKIT):
        t = time.perf_counter()
        out = transpile(qc, coupling_map=cm, basis_gates=["u", "cx"],
                        optimization_level=opt, seed_transpiler=42 + trial)
        times.append((time.perf_counter() - t) * 1000.0)
        ops = out.count_ops()
        gate_counts.append(len(out.data))
        cx_counts.append(ops.get("cx", 0))
        depths.append(out.depth())
    return median_record(times, gate_counts, cx_counts, depths)


# --------------------------------------------------------------------------- #
# tket                                                                         #
# --------------------------------------------------------------------------- #
def run_tket(qasm_raw, edges, opt=3):
    from pytket.qasm import circuit_from_qasm_str
    from pytket.architecture import Architecture
    from pytket.circuit import OpType
    from pytket.passes import (DefaultMappingPass, FullPeepholeOptimise,
                               DecomposeSwapsToCXs, SequencePass)
    try:
        from pytket.passes import AutoRebase
        rebase = AutoRebase({OpType.CX, OpType.TK1})
    except ImportError:
        from pytket.passes import auto_rebase_pass
        rebase = auto_rebase_pass({OpType.CX, OpType.TK1})

    arch = Architecture([tuple(e) for e in edges])
    times, gate_counts, cx_counts, depths = [], [], [], []
    for _ in range(TRIALS_TKET):
        c = circuit_from_qasm_str(qasm_raw)
        passes = []
        if opt >= 3:
            passes.append(FullPeepholeOptimise())
        passes += [DefaultMappingPass(arch), DecomposeSwapsToCXs(arch), rebase]
        t = time.perf_counter()
        SequencePass(passes).apply(c)
        times.append((time.perf_counter() - t) * 1000.0)
        gate_counts.append(c.n_gates)
        cx_counts.append(c.n_gates_of_type(OpType.CX))
        depths.append(c.depth())
    return median_record(times, gate_counts, cx_counts, depths)


# --------------------------------------------------------------------------- #
# Cirq                                                                         #
# --------------------------------------------------------------------------- #
def run_cirq(qasm_norm, edges, n):
    # Cirq 1.3 references the pre-NumPy-2.0 location of ComplexWarning.
    import numpy as np
    if not hasattr(np, "ComplexWarning"):
        np.ComplexWarning = np.exceptions.ComplexWarning
    import cirq
    import networkx as nx
    from cirq.contrib.qasm_import import circuit_from_qasm

    circ = circuit_from_qasm(qasm_norm)
    # Map NamedQubit('q_N') -> LineQubit(N).
    remap = {}
    for q in circ.all_qubits():
        m = re.match(r"q_(\d+)", q.name)
        remap[q] = cirq.LineQubit(int(m.group(1)))
    circ = circ.transform_qubits(remap)
    g = nx.Graph()
    g.add_nodes_from(cirq.LineQubit(i) for i in range(n))
    for a, b in edges:
        g.add_edge(cirq.LineQubit(a), cirq.LineQubit(b))

    def rebase_to_cx(circuit):
        """Lower Cirq's residual CZ gates to the common {1q, CX} level."""
        def lower_cz(op, _moment_index):
            if isinstance(op.gate, cirq.CZPowGate):
                exponent = float(op.gate.exponent)
                if abs((exponent % 2.0) - 1.0) > 1e-8:
                    raise ValueError(f"non-integral CZ exponent after normalization: {op}")
                control, target = op.qubits
                return [cirq.H(target), cirq.CNOT(control, target), cirq.H(target)]
            return op

        circuit = cirq.map_operations_and_unroll(circuit, lower_cz)
        circuit = cirq.merge_single_qubit_gates_to_phxz(circuit)
        for op in circuit.all_operations():
            if cirq.num_qubits(op) != 2:
                continue
            if not isinstance(op.gate, cirq.CNotPowGate):
                raise ValueError(f"Cirq output is not closed over CX: {op}")
            exponent = float(op.gate.exponent)
            if abs((exponent % 2.0) - 1.0) > 1e-8:
                raise ValueError(f"non-CX controlled-X exponent in output: {op}")
        return circuit

    times, gate_counts, cx_counts, depths = [], [], [], []
    for _ in range(TRIALS_CIRQ):
        t = time.perf_counter()
        # This is compilation work, so keep it inside the timed region just as
        # Qiskit and tket keep their high-arity lowering inside transpilation.
        lowered = cirq.Circuit(cirq.decompose(
            circ, keep=lambda op: cirq.num_qubits(op) <= 2))
        r = cirq.RouteCQC(g)(lowered)
        # Decompose routing SWAPs and merge single-qubit gates.
        r = cirq.Circuit(cirq.decompose(
            r, keep=lambda op: cirq.num_qubits(op) <= 2
            and not isinstance(op.gate, cirq.SwapPowGate)))
        r = rebase_to_cx(r)
        times.append((time.perf_counter() - t) * 1000.0)
        two_q = sum(cirq.num_qubits(op) == 2 for op in r.all_operations())
        gate_counts.append(sum(1 for _ in r.all_operations()))
        cx_counts.append(two_q)
        depths.append(len(r))

    return median_record(times, gate_counts, cx_counts, depths)


TOOLS = {"qiskit": run_qiskit, "tket": run_tket}  # cirq handled separately (diff signature)


def main():
    manifest = json.load(open(os.path.join(FIXTURE_DIR, "manifest.json")))
    os.makedirs(RESULTS, exist_ok=True)
    records = []
    total = sum(len(topologies_for(fx["num_qubits"])) for fx in manifest)
    idx = 0
    for fx in manifest:
        raw = open(os.path.join(HERE, fx["file"])).read()
        norm = normalize_qasm(raw)  # crz-expanded, for cirq's importer only
        n = fx["num_qubits"]
        for topo, (bq, edges) in topologies_for(n).items():
            idx += 1
            base = {"fixture": fx["name"], "family": fx["family"],
                    "nominal_n": fx["nominal_n"], "num_qubits": n,
                    "topology": topo, "opt_level": 3}
            for tool, fn in TOOLS.items():
                try:
                    rec = fn(raw, edges, 3)
                    records.append({**base, "tool": tool, "status": "ok", **rec})
                except Exception as e:
                    records.append({**base, "tool": tool, "status": "error",
                                    "error": f"{type(e).__name__}: {str(e)[:120]}"})
            try:
                rec = run_cirq(norm, edges, n)
                records.append({**base, "tool": "cirq", "status": "ok", **rec})
            except Exception as e:
                records.append({**base, "tool": "cirq", "status": "error",
                                "error": f"{type(e).__name__}: {str(e)[:120]}"})
            print(f"[{idx:>3}/{total}] {fx['name']:<12} {topo:<12} done", file=sys.stderr)

    out = os.path.join(RESULTS, "baselines.json")
    json.dump(records, open(out, "w"), indent=2)
    ok = sum(1 for r in records if r.get("status") == "ok")
    print(f"\n[baselines] {ok}/{len(records)} ok -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
