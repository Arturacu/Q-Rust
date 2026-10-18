"""Single-trial Qiskit/TKET adapters for the corrected-fixture comparison.

The caller supplies identical normalized QASM to every compiler. Parsing and
serialization are outside the compile timer. The timer includes copying the
parsed circuit, constructing target/pass objects, compilation and retrieving
the compiled circuit. Output is the actual measured circuit, in physical q[i]
coordinates, with logical input/output locations recorded separately.

References:
https://qiskit.qotlabs.org/docs/api/qiskit/qiskit.transpiler.TranspileLayout
https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html
https://docs.quantinuum.com/tket/user-guide/manual/manual_circuit.html
https://docs.quantinuum.com/tket/api-docs/passes.html
"""
from __future__ import annotations

import math
from time import perf_counter_ns


def _target(n_logical, physical_width, edges):
    if not 1 <= n_logical <= physical_width:
        raise ValueError("logical width must be positive and fit physical width")
    result = set()
    for a, b in edges:
        if not (0 <= a < physical_width and 0 <= b < physical_width) or a == b:
            raise ValueError(f"invalid target edge {(a, b)}")
        result.add(tuple(sorted((a, b))))
    return sorted(result)


def _layouts(initial, final, n_logical, physical_width):
    for label, mapping in (("initial", initial), ("final", final)):
        if (len(mapping) != n_logical or len(set(mapping)) != n_logical
                or any(not 0 <= i < physical_width for i in mapping)):
            raise ValueError(f"invalid {label} layout: {mapping}")


def _qasm_header(physical_width):
    return ['OPENQASM 2.0;', 'include "qelib1.inc";', f'qreg q[{physical_width}];']


def _number(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("non-finite gate parameter in compiled output")
    return format(number, ".17g")


def _qiskit(qasm, n_logical, physical_width, edges, topology, seed):
    from qiskit import qasm2, transpile
    from qiskit.transpiler import CouplingMap

    started = perf_counter_ns()
    circuit = qasm2.loads(qasm, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    parse_ns = perf_counter_ns() - started
    if circuit.num_qubits != n_logical:
        raise ValueError("declared logical width differs from parsed QASM")
    if circuit.num_clbits or any(x.operation.name in {"measure", "reset"} for x in circuit.data):
        raise ValueError("comparison supports static unitary circuits only")

    started = perf_counter_ns()
    source = circuit.copy()
    if topology == "all_to_all":
        coupling = None
    else:
        coupling = CouplingMap()
        for i in range(physical_width):
            coupling.add_physical_qubit(i)
        for a, b in edges:
            coupling.add_edge(a, b)
            coupling.add_edge(b, a)
    output = transpile(
        source, coupling_map=coupling, basis_gates=["u", "cx"],
        optimization_level=3, seed_transpiler=seed, approximation_degree=1.0,
        num_processes=1, qubits_initially_zero=False,
    )
    compile_ns = perf_counter_ns() - started

    if output.layout is None:
        initial = final = list(range(n_logical))
    else:
        initial = list(output.layout.initial_index_layout(filter_ancillas=True))
        final = list(output.layout.final_index_layout(filter_ancillas=True))
    _layouts(initial, final, n_logical, physical_width)
    if output.num_qubits > physical_width:
        raise ValueError("compiler exceeded physical width")
    lines = _qasm_header(physical_width)
    for instruction in output.data:
        op = instruction.operation
        qubits = [output.find_bit(q).index for q in instruction.qubits]
        if op.name == "barrier":
            continue
        if op.name == "cx" and len(qubits) == 2:
            lines.append(f"cx q[{qubits[0]}],q[{qubits[1]}];")
        elif op.name == "u" and len(qubits) == 1:
            lines.append(f"u3({','.join(_number(p) for p in op.params)}) q[{qubits[0]}];")
        else:
            raise ValueError(f"Qiskit output outside U3/CX basis: {op.name}")
    return {
        "parse_ns": parse_ns, "compile_ns": compile_ns,
        "out_qasm": "\n".join(lines) + "\n",
        "initial_layout": initial, "final_layout": final,
        "output_qubits": physical_width,
        "pipeline": {
            "compiler": "qiskit.transpile", "optimization_level": 3,
            "basis_gates": ["u", "cx"], "approximation_degree": 1.0,
            "seed_transpiler": seed, "num_processes": 1,
            "qubits_initially_zero": False,
            "routing": "none (unconstrained)" if coupling is None else "default preset, bidirectional coupling",
            "serialization": "measured U/CX instructions; global phase omitted",
        },
    }


def _tket(qasm, n_logical, physical_width, edges, topology, seed):
    from pytket.architecture import Architecture
    from pytket.circuit import OpType, Qubit
    from pytket.passes import (
        AutoRebase, DefaultMappingPass, DecomposeSwapsToCXs,
        FullPeepholeOptimise, SequencePass, SynthesiseTket,
    )
    from pytket.predicates import CompilationUnit
    from pytket.qasm import circuit_from_qasm_str

    started = perf_counter_ns()
    circuit = circuit_from_qasm_str(qasm)
    parse_ns = perf_counter_ns() - started
    originals = [Qubit("q", i) for i in range(n_logical)]
    if set(circuit.qubits) != set(originals):
        raise ValueError("expected one logical register q of the declared width")
    if circuit.bits or any(c.op.type in {OpType.Measure, OpType.Reset, OpType.Conditional}
                           for c in circuit.get_commands()):
        raise ValueError("comparison supports static unitary circuits only")

    started = perf_counter_ns()
    unit = CompilationUnit(circuit.copy())
    passes = [FullPeepholeOptimise(allow_swaps=True, target_2qb_gate=OpType.CX)]
    if topology != "all_to_all":
        architecture = Architecture(edges)
        # The TKET compilation manual recommends SynthesiseTket after routing:
        # its light optimizations remove routing redundancies while preserving
        # connectivity. FullPeepholeOptimise is intentionally kept before mapping.
        passes.extend([
            DefaultMappingPass(architecture), DecomposeSwapsToCXs(architecture),
            SynthesiseTket(),
        ])
    passes.append(AutoRebase({OpType.CX, OpType.U3}))
    SequencePass(passes).apply(unit)
    output = unit.circuit
    compile_ns = perf_counter_ns() - started

    def physical_index(qubit):
        if qubit.reg_name not in {"node", "q"} or len(qubit.index) != 1:
            raise ValueError(f"unresolved physical qubit {qubit}")
        value = int(qubit.index[0])
        if not 0 <= value < physical_width:
            raise ValueError(f"physical qubit outside target: {qubit}")
        return value

    initial = [physical_index(unit.initial_map[q]) for q in originals]
    permutation = output.implicit_qubit_permutation()
    inverse_permutation = {end: start for start, end in permutation.items()}
    # final_map names the circuit's output wires. get_commands() follows the
    # physical wires before the implicit output rewiring, so invert that
    # rewiring to locate each logical output in the emitted command sequence.
    final = [physical_index(inverse_permutation[unit.final_map[q]]) for q in originals]
    _layouts(initial, final, n_logical, physical_width)
    lines = _qasm_header(physical_width)
    # Commands are emitted directly: a QASM exporter may materialize implicit
    # wire permutations and thereby double-apply a permutation in final_layout.
    for command in output.get_commands():
        op = command.op
        qubits = [physical_index(q) for q in command.qubits]
        if op.type in {OpType.Barrier, OpType.Phase}:
            continue
        if op.type == OpType.CX and len(qubits) == 2:
            lines.append(f"cx q[{qubits[0]}],q[{qubits[1]}];")
        elif op.type == OpType.U3 and len(qubits) == 1:
            lines.append(f"u3({','.join(_number(float(p) * math.pi) for p in op.params)}) q[{qubits[0]}];")
        else:
            raise ValueError(f"TKET output outside U3/CX basis: {op.type}")
    return {
        "parse_ns": parse_ns, "compile_ns": compile_ns,
        "out_qasm": "\n".join(lines) + "\n",
        "initial_layout": initial, "final_layout": final,
        "output_qubits": physical_width,
        "pipeline": {
            "compiler": "pytket.SequencePass",
            "passes": ["FullPeepholeOptimise(allow_swaps=True,target_2qb_gate=CX)"]
                + ([] if topology == "all_to_all" else [
                    "DefaultMappingPass(GraphPlacement; undirected architecture)",
                    "DecomposeSwapsToCXs(architecture)",
                    "SynthesiseTket(connectivity-preserving postrouting cleanup)",
                ]) + ["AutoRebase({CX,U3})"],
            "seed": None, "seed_note": "chosen passes expose no stochastic seed",
            "implicit_permutation": {
                str(q): str(r) for q, r in permutation.items()
            },
            "serialization": "explicit measured commands; final layout is inverse(implicit_permutation) composed with final_map; global phase omitted",
        },
    }


def compile_one(tool, qasm, n_logical, physical_width, edges, topology, seed):
    """Compile one trial, without validation, retries or timed file I/O."""
    edges = _target(n_logical, physical_width, edges)
    if tool == "qiskit":
        return _qiskit(qasm, n_logical, physical_width, edges, topology, seed)
    if tool == "tket":
        return _tket(qasm, n_logical, physical_width, edges, topology, seed)
    raise ValueError(f"unsupported compiler: {tool}")
