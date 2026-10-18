#!/usr/bin/env python3
"""Secondary recipe ablations and stage probes; never alters primary records.

The five cases and six recipes are declared here before running. Each recipe
uses the original study's identical normalized input, common CX/1q output,
free layouts, and full physical-state oracle. Timings are warm compilation,
with parsing, serialization and validation outside the clock. This targeted,
post-hoc study explains mechanisms; it is not a new whole-corpus estimate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import time

CASE_IDS = ("qft_14__linear", "qft_8__linear", "ghz_8__all_to_all",
            "clifford_12__linear", "wstate_14__all_to_all")
RECIPES = ("cirq_original", "cirq_one_pass", "cirq_native_cx",
           "tket_original", "tket_two_qubit", "tket_light")
AVAILABLE_RECIPES = ("qrust", *RECIPES)
MASTER_SEED = 20261016


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def _native_cx(circuit, cirq):
    """Exact light path for the declared cases' 1q/CX/SWAP input subset.

    This intentionally does not claim to be a general Cirq transpiler. It
    preserves native CXs and lowers each input/routing SWAP to three CXs,
    then fuses one-qubit gates. No numerical two-qubit resynthesis occurs.
    """
    def lower(op, _):
        if len(op.qubits) <= 1 and cirq.has_unitary(op):
            return op
        if op.gate == cirq.CX:
            return op
        if op.gate == cirq.SWAP:
            a, b = op.qubits
            return [cirq.CX(a, b), cirq.CX(b, a), cirq.CX(a, b)]
        raise ValueError(f"native-CX diagnostic input outside 1q/CX/SWAP: {op}")
    return cirq.map_operations_and_unroll(circuit, lower)


def cirq_trial(case, recipe):
    import cirq
    import networkx as nx
    from cirq.contrib.qasm_import import circuit_from_qasm
    from comparison_cirq import _complete_initial_layout, _serialize, _to_cx, SYNTHESIS_ATOL

    parse_start = time.perf_counter_ns()
    parsed = circuit_from_qasm(case["qasm"])
    mapping = {}
    for qubit in parsed.all_qubits():
        match = re.fullmatch(r"q_(\d+)", qubit.name)
        if match is None or int(match.group(1)) >= case["n_logical"]:
            raise ValueError(f"Unexpected input wire {qubit}")
        mapping[qubit] = cirq.LineQubit(int(match.group(1)))
    parsed = parsed.transform_qubits(mapping)
    parse_ns = time.perf_counter_ns() - parse_start
    stages = {}
    start = mark = time.perf_counter_ns()
    circuit = parsed.copy()
    logical = cirq.LineQubit.range(case["n_logical"])
    physical = cirq.LineQubit.range(case["physical_width"])
    gateset = None
    if recipe != "cirq_native_cx":
        gateset = cirq.CZTargetGateset(atol=SYNTHESIS_ATOL, allow_partial_czs=False,
                                     preserve_moment_structure=False, reorder_operations=True)
    now = time.perf_counter_ns()
    stages["copy_and_configuration_ns"] = now - mark
    mark = now
    passes = 2 if recipe == "cirq_original" else 1
    if gateset is None:
        circuit = _native_cx(circuit, cirq)
    else:
        circuit = cirq.optimize_for_target_gateset(circuit, gateset=gateset,
                         max_num_passes=passes, ignore_failures=False)
    now = time.perf_counter_ns()
    stages["prerouting_synthesis_ns"] = now - mark
    mark = now
    if case["topology"] == "all_to_all":
        initial = dict(zip(logical, physical))
        swaps = {}
    else:
        graph = nx.Graph()
        graph.add_nodes_from(physical)
        graph.add_edges_from((physical[a], physical[b]) for a, b in case["edges"])
        if circuit.all_qubits():
            circuit, initial, swaps = cirq.RouteCQC(graph).route_circuit(circuit)
        else:
            initial, swaps = {}, {}
        initial = _complete_initial_layout(initial, logical, physical)
    now = time.perf_counter_ns()
    stages["graph_and_routing_ns"] = now - mark
    mark = now
    if case["topology"] != "all_to_all":
        if gateset is None:
            circuit = _native_cx(circuit, cirq)
        else:
            circuit = cirq.optimize_for_target_gateset(circuit, gateset=gateset,
                             max_num_passes=passes, ignore_failures=False)
    now = time.perf_counter_ns()
    stages["postrouting_synthesis_ns"] = now - mark
    mark = now
    if gateset is None:
        circuit = cirq.merge_single_qubit_gates_to_phxz(circuit, atol=SYNTHESIS_ATOL)
    else:
        circuit = _to_cx(circuit, cirq)
    finish = time.perf_counter_ns()
    stages["cx_closure_and_1q_fusion_ns"] = finish - mark
    return {"compile_ns": finish - start, "parse_ns": parse_ns, "stages_ns": stages,
            "initial_layout": [initial[q].x for q in logical],
            "final_layout": [swaps.get(initial[q], initial[q]).x for q in logical],
            "out_qasm": _serialize(circuit, case["physical_width"], cirq),
            "output_qubits": case["physical_width"]}


def tket_trial(case, recipe, stage_probe=False):
    from pytket.architecture import Architecture
    from pytket.circuit import OpType, Qubit
    from pytket.passes import (AutoRebase, DefaultMappingPass, DecomposeSwapsToCXs,
                              FullPeepholeOptimise, PeepholeOptimise2Q,
                              SequencePass, SynthesiseTket)
    from pytket.predicates import CompilationUnit
    from pytket.qasm import circuit_from_qasm_str
    from comparison_qiskit_tket import _layouts, _number, _qasm_header

    parse_start = time.perf_counter_ns()
    circuit = circuit_from_qasm_str(case["qasm"])
    originals = [Qubit("q", i) for i in range(case["n_logical"])]
    parse_ns = time.perf_counter_ns() - parse_start
    stages = {}
    start = time.perf_counter_ns()
    unit = CompilationUnit(circuit.copy())
    if recipe == "tket_original":
        passes = [("FullPeepholeOptimise", FullPeepholeOptimise(allow_swaps=True, target_2qb_gate=OpType.CX))]
    elif recipe == "tket_two_qubit":
        passes = [("PeepholeOptimise2Q", PeepholeOptimise2Q(allow_swaps=True))]
    elif recipe == "tket_light":
        passes = [("SynthesiseTket_pre", SynthesiseTket())]
    else:
        raise ValueError(f"Unknown recipe {recipe}")
    if case["topology"] != "all_to_all":
        architecture = Architecture(case["edges"])
        passes.extend([("DefaultMappingPass", DefaultMappingPass(architecture)),
                       ("DecomposeSwapsToCXs", DecomposeSwapsToCXs(architecture)),
                       ("SynthesiseTket_post", SynthesiseTket())])
    passes.append(("AutoRebase", AutoRebase({OpType.CX, OpType.U3})))
    sequence = SequencePass([p for _, p in passes])
    configured = time.perf_counter_ns()
    stages["copy_and_configuration_ns"] = configured - start
    if stage_probe:
        for name, compiler_pass in passes:
            before = time.perf_counter_ns()
            compiler_pass.apply(unit)
            stages[name + "_ns"] = time.perf_counter_ns() - before
    else:
        before = time.perf_counter_ns()
        sequence.apply(unit)
        stages["SequencePass_apply_ns"] = time.perf_counter_ns() - before
    applied = time.perf_counter_ns()
    output = unit.circuit
    finish = time.perf_counter_ns()
    stages["output_retrieval_ns"] = finish - applied

    def physical_index(qubit):
        if qubit.reg_name not in {"node", "q"} or len(qubit.index) != 1:
            raise ValueError(f"Unresolved physical wire {qubit}")
        value = int(qubit.index[0])
        if not 0 <= value < case["physical_width"]:
            raise ValueError(f"Physical wire outside target {qubit}")
        return value

    initial = [physical_index(unit.initial_map[q]) for q in originals]
    inverse = {end: begin for begin, end in output.implicit_qubit_permutation().items()}
    final = [physical_index(inverse[unit.final_map[q]]) for q in originals]
    _layouts(initial, final, case["n_logical"], case["physical_width"])
    lines = _qasm_header(case["physical_width"])
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
            raise ValueError(f"Output outside U3/CX basis: {op.type}")
    return {"compile_ns": finish - start, "parse_ns": parse_ns, "stages_ns": stages,
            "initial_layout": initial, "final_layout": final,
            "out_qasm": "\n".join(lines) + "\n", "output_qubits": case["physical_width"]}


def load_cases(protocol_path, all_cases=False):
    protocol = json.loads(protocol_path.read_text())
    indexed = {case["id"]: case for case in protocol["cases"]}
    cases = list(indexed.values()) if all_cases else [indexed[case_id] for case_id in CASE_IDS]
    for case in cases:
        if digest(case["qasm"]) != case["normalized_sha256"] or digest(case["source_qasm"]) != case["source_sha256"]:
            raise ValueError(f"Frozen fixture hash mismatch: {case['id']}")
    return cases


def compile_recipe(case, recipe, qrust_worker=None, seed=MASTER_SEED):
    """Dispatch a trial without moving IPC or serialization into native timers."""
    if recipe == "qrust":
        if qrust_worker is None:
            raise ValueError("qrust recipe requires a persistent native worker")
        request = {key: case[key] for key in
                   ("qasm", "n_logical", "physical_width", "edges", "topology")}
        request["seed"] = seed
        result = qrust_worker.call(request, 60)
        result.setdefault("stages_ns", {})
        return result
    if recipe not in RECIPES:
        raise ValueError(f"Unknown recipe {recipe}")
    return cirq_trial(case, recipe) if recipe.startswith("cirq_") else tket_trial(case, recipe)


def timing(args):
    import cirq
    import pytket
    cases = load_cases(args.protocol, args.all_cases)
    recipes = args.recipes or RECIPES
    # A new campaign must never overwrite an earlier protocol or raw records.
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = {
        "kind": ("post-hoc full-corpus recipe sensitivity; distinct from the unchanged-recipe primary run"
                 if args.all_cases else "post-hoc, targeted recipe ablations and stage diagnostics; excluded from primary estimates"),
        "case_ids": [c["id"] for c in cases], "recipes": list(recipes), "cases": cases,
        "selection": ("All original protocol cases; recipe sensitivity, selected after five-case diagnostics"
                      if args.all_cases else "Two QFT widths on a line plus simple GHZ, Clifford routing and W-state synthesis; chosen before this diagnostic run, after seeing primary results"),
        "versions": {"cirq": cirq.__version__, "pytket": pytket.__version__},
        "source_sha256": digest(Path(__file__).read_text()),
        "primary_protocol_sha256": hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
        "measured_repeats": args.trials, "warmups": 1, "stage_probe_repeats": args.probes,
        "master_seed": MASTER_SEED, "randomization": "shuffle cases once; within each case/trial shuffle recipes; all compilation sequential",
        "thread_environment": {k: os.environ.get(k) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS", "PYTHONHASHSEED"]},
        "timing_boundary": "copy, recipe/target construction, compile, result retrieval; excludes import, parse, output serialization and validation",
        "variants": {
            "qrust": "Canonical public O3 API through persistent comparison_runner; native timer and effective pipeline retained per output",
            "cirq_original": "Frozen CZ/RouteCQC/CZ/CX recipe, max 2 synthesis iterations",
            "cirq_one_pass": "Same recipe with max 1 synthesis iteration per phase",
            "cirq_native_cx": "Exact 1q/CX/SWAP-subset lowering, RouteCQC, exact 3-CX SWAP lowering, 1q fusion; no 2q resynthesis",
            "tket_original": "Frozen FullPeepholeOptimise + mapping + connectivity-preserving cleanup + U3/CX rebase",
            "tket_two_qubit": "Replace only initial FullPeepholeOptimise with PeepholeOptimise2Q",
            "tket_light": "Replace only initial FullPeepholeOptimise with SynthesiseTket",
        },
        "stage_probe_caveat": "TKET probes apply identical top-level passes individually to one CompilationUnit instead of SequencePass; separate from recipe timings and independently validated. Python dispatch/precondition overhead can differ.",
        "limitations": "Post-hoc recipe selection on the same local development corpus; no independent holdout. Recipe choices do not have equal search budgets. Does not estimate inherent language costs or each framework's best configuration. Native-CX Cirq recipe only supports the declared 1q/CX/SWAP subset.",
    }
    if "qrust" in recipes:
        if args.qrust_binary is None:
            raise ValueError("qrust recipe requires --qrust-binary")
        metadata["qrust_binary"] = {"path": str(args.qrust_binary.resolve()),
            "sha256": hashlib.sha256(args.qrust_binary.read_bytes()).hexdigest(),
            "worker": "comparison_study.Worker; persistent process; native compile_ns excludes IPC",
            "pipeline_evidence": "Each output retains pipeline, effective_search and seed_used from the native worker"}
    (args.output / "tradeoffs_protocol.json").write_text(json.dumps(metadata, indent=2) + "\n")
    rng = random.Random(MASTER_SEED)
    rng.shuffle(cases)
    with ExitStack() as resources, (args.output / "tradeoffs_timed.jsonl").open("w") as output:
        qrust_worker = None
        if "qrust" in recipes:
            from comparison_study import Worker
            qrust_worker = Worker("qrust", args.qrust_binary, args.output)
            resources.callback(qrust_worker.close)
        for case in cases:
            for trial in range(-1, args.trials):
                order = list(recipes)
                rng.shuffle(order)
                for recipe in order:
                    result = compile_recipe(case, recipe, qrust_worker, MASTER_SEED + max(trial, 0))
                    row = {"case_id": case["id"], "recipe": recipe, "trial": trial,
                           "kind": "warmup" if trial < 0 else "measured", **result}
                    output.write(json.dumps(row, allow_nan=False) + "\n")
                    output.flush()
            print(json.dumps({"timed_case": case["id"]}), flush=True)
        # Stage probes are a separately labeled, post-timing diagnostic.
        for case in cases if args.probes else []:
            for trial in range(-1, args.probes):
                order = [r for r in recipes if r.startswith("tket_")]
                rng.shuffle(order)
                for recipe in order:
                    row = {"case_id": case["id"], "recipe": recipe, "trial": trial,
                           "kind": "stage_warmup" if trial < 0 else "stage_probe",
                           **tket_trial(case, recipe, stage_probe=True)}
                    output.write(json.dumps(row, allow_nan=False) + "\n")
                    output.flush()
            print(json.dumps({"probed_case": case["id"]}), flush=True)


def validate_record_grid(metadata, rows):
    """Reject incomplete, duplicate or unexpected scheduled trial records."""
    case_ids = metadata["case_ids"]
    recipes = metadata["recipes"]
    if len(case_ids) != len(set(case_ids)) or len(recipes) != len(set(recipes)):
        raise ValueError("Protocol repeats a case or recipe")
    repeats, probes = metadata["measured_repeats"], metadata["stage_probe_repeats"]
    if (type(repeats) is not int or repeats < 1 or type(probes) is not int or probes < 0
            or metadata["warmups"] != 1):
        raise ValueError("Unsupported protocol repeat or warmup count")
    expected = set()
    for case_id in case_ids:
        for recipe in recipes:
            for trial in range(-1, repeats):
                expected.add((case_id, recipe, trial, "warmup" if trial == -1 else "measured"))
            if probes and recipe.startswith("tket_"):
                for trial in range(-1, probes):
                    expected.add((case_id, recipe, trial, "stage_warmup" if trial == -1 else "stage_probe"))
    actual = Counter()
    for row in rows:
        if type(row.get("trial")) is not int:
            raise ValueError("Trial index must be an integer")
        actual[(row["case_id"], row["recipe"], row["trial"], row["kind"])] += 1
    missing = expected - actual.keys()
    extra = actual.keys() - expected
    duplicate = {key for key, count in actual.items() if count != 1}
    if missing or extra or duplicate:
        raise ValueError(f"Incomplete or invalid trial grid: {len(missing)} missing, "
                         f"{len(extra)} unexpected, {len(duplicate)} duplicate")
    return {"expected_records": len(expected), "actual_records": len(rows),
            "missing": 0, "unexpected": 0, "duplicates": 0}


def validate(args):
    from comparison_validation import validate_output
    metadata = json.loads((args.output / "tradeoffs_protocol.json").read_text())
    cases = {c["id"]: c for c in metadata["cases"]}
    if set(cases) != set(metadata["case_ids"]) or len(cases) != len(metadata["cases"]):
        raise ValueError("Protocol case identifiers and fixture records disagree")
    for case in cases.values():
        if digest(case["qasm"]) != case["normalized_sha256"] or digest(case["source_qasm"]) != case["source_sha256"]:
            raise ValueError(f"Frozen fixture hash mismatch: {case['id']}")
    with (args.output / "tradeoffs_timed.jsonl").open() as source:
        scheduled = [json.loads(line) for line in source]
    grid = validate_record_grid(metadata, scheduled)
    cache = {}
    rows = []
    with (args.output / "tradeoffs_validated.jsonl").open("w") as dest:
        for row in scheduled:
            case = cases[row["case_id"]]
            if row.get("status", "compiled") != "compiled" or "out_qasm" not in row:
                row["validation"] = {"passed": False, "evidence_class": "compile_failure",
                                     "error": row.get("error", "Missing compiled output")}
                row["validation_cache_key"] = None
                dest.write(json.dumps(row, allow_nan=False) + "\n")
                rows.append(row)
                continue
            key = digest(json.dumps([row["case_id"], row["out_qasm"], row["initial_layout"], row["final_layout"]]))
            if key not in cache:
                edges = None if case["topology"] == "all_to_all" else [e for a, b in case["edges"] for e in [(a, b), (b, a)]]
                cache[key] = validate_output(case["source_qasm"], row["out_qasm"], case["n_logical"],
                       case["physical_width"], edges, row["initial_layout"], row["final_layout"])
            row["validation"] = cache[key]
            row["validation_cache_key"] = key
            dest.write(json.dumps(row, allow_nan=False) + "\n")
            rows.append(row)
    summaries = []
    for case_id in metadata["case_ids"]:
        for recipe in metadata["recipes"]:
            for kind in ("measured", "stage_probe"):
                group = [r for r in rows if r["case_id"] == case_id and r["recipe"] == recipe and r["kind"] == kind]
                if not group:
                    continue
                if any(r["validation"]["passed"] is not True for r in group):
                    summaries.append({"case_id": case_id, "recipe": recipe, "kind": kind,
                                      "trials": len(group), "all_passed": False,
                                      "failed_trials": sum(r["validation"]["passed"] is not True for r in group)})
                    continue
                stages = group[0]["stages_ns"]
                summaries.append({"case_id": case_id, "recipe": recipe, "kind": kind,
                    "trials": len(group), "median_compile_ns": statistics.median(r["compile_ns"] for r in group),
                    "median_cx": statistics.median(r["validation"]["metrics"]["cx_count"] for r in group),
                    "median_depth": statistics.median(r["validation"]["metrics"]["depth"] for r in group),
                    "stage_median_ns": {s: statistics.median(r["stages_ns"][s] for r in group) for s in stages},
                    "stage_median_fraction": {s: statistics.median(r["stages_ns"][s] / r["compile_ns"] for r in group) for s in stages},
                    "all_passed": all(r["validation"]["passed"] is True for r in group)})
    summary = {"rows": len(rows), "unique_artifacts": len(cache), "grid_integrity": grid,
               "all_passed": all(r["validation"]["passed"] is True for r in rows),
               "max_input_residual": max((r["validation"]["max_input_residual"] for r in rows
                                           if "max_input_residual" in r["validation"]), default=None),
               "cases": summaries}
    (args.output / "tradeoffs_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["all_passed"]:
        raise SystemExit("At least one diagnostic output failed independent validation")


def stage_counts(args):
    """Inspect exact Cirq intermediate circuits in separate untimed executions."""
    import cirq
    from unittest.mock import patch
    import comparison_cirq
    from comparison_validation import validate_output

    cases = load_cases(args.protocol, args.all_cases)
    records = []
    optimize = cirq.optimize_for_target_gateset
    route = cirq.RouteCQC.route_circuit
    lower = comparison_cirq._to_cx

    def counts(circuit):
        ops = list(circuit.all_operations())
        return {"operations": len(ops), "moments": len(circuit),
                "two_qubit": sum(len(op.qubits) == 2 for op in ops),
                "cx": sum(op.gate == cirq.CX for op in ops),
                "cz": sum(op.gate == cirq.CZ for op in ops),
                "swap": sum(op.gate == cirq.SWAP for op in ops),
                "gate_classes": dict(Counter(type(op.gate).__name__ for op in ops))}

    for case in cases:
        for recipe in (r for r in (args.recipes or RECIPES) if r.startswith("cirq_")):
            stages = []

            def observed_optimize(circuit, *positional, **options):
                before = counts(circuit)
                result = optimize(circuit, *positional, **options)
                stages.append({"stage": "CZ_synthesis", "max_passes": options["max_num_passes"],
                               "before": before, "after": counts(result)})
                return result

            def observed_route(router, circuit, *positional, **options):
                before = counts(circuit)
                result = route(router, circuit, *positional, **options)
                stages.append({"stage": "RouteCQC", "before": before, "after": counts(result[0])})
                return result

            def observed_lower(circuit, cirq_module):
                before = counts(circuit)
                result = lower(circuit, cirq_module)
                stages.append({"stage": "CZ_to_CX_and_1q_fusion", "before": before, "after": counts(result)})
                return result

            with patch.object(cirq, "optimize_for_target_gateset", observed_optimize), \
                 patch.object(cirq.RouteCQC, "route_circuit", observed_route), \
                 patch.object(comparison_cirq, "_to_cx", observed_lower):
                result = cirq_trial(case, recipe)
            # Timings from these instrumented executions are intentionally dropped.
            for key in ("compile_ns", "parse_ns", "stages_ns"):
                result.pop(key)
            edges = None if case["topology"] == "all_to_all" else [e for a, b in case["edges"] for e in [(a, b), (b, a)]]
            validation = validate_output(case["source_qasm"], result["out_qasm"],
                    case["n_logical"], case["physical_width"], edges,
                    result["initial_layout"], result["final_layout"])
            records.append({"case_id": case["id"], "recipe": recipe,
                            "stages": stages, "validation": validation, **result})
    result = {"kind": "separate untimed Cirq intermediate-gate-count diagnostic",
              "source_sha256": digest(Path(__file__).read_text()),
              "all_passed": all(r["validation"]["passed"] is True for r in records),
              "records": records}
    (args.output / "cirq_stage_counts.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"outputs": len(records), "all_passed": result["all_passed"]}))
    if not result["all_passed"]:
        raise SystemExit("At least one stage-inspection output failed independent validation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["time", "validate", "stage-counts"])
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=7)
    parser.add_argument("--probes", type=int, default=3)
    parser.add_argument("--all-cases", action="store_true",
                        help="Use every case in the original protocol for a separately labeled sensitivity run")
    parser.add_argument("--recipes", choices=AVAILABLE_RECIPES, nargs="+",
                        help="Select named recipes; omitted means all six")
    parser.add_argument("--qrust-binary", type=Path,
                        help="Persistent canonical Q-Rust worker, required when selecting recipe qrust")
    args = parser.parse_args()
    if args.trials < 1 or args.probes < 0:
        parser.error("trial count must be positive and probe count nonnegative")
    if args.recipes and len(args.recipes) != len(set(args.recipes)):
        parser.error("recipe names must be unique")
    if args.mode == "time" and args.recipes and "qrust" in args.recipes:
        if args.qrust_binary is None or not args.qrust_binary.is_file():
            parser.error("qrust recipe requires an existing --qrust-binary")
    if args.mode in {"time", "stage-counts"}:
        if args.protocol is None:
            parser.error("--protocol is required for timing or stage counts")
        if args.mode == "time":
            timing(args)
        else:
            stage_counts(args)
    else:
        validate(args)


if __name__ == "__main__":
    main()
