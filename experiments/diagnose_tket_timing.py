#!/usr/bin/env python3
"""Secondary timing-boundary diagnostic; never modifies primary study records.

Reconstructs the frozen TKET pass recipe, then splits its timed region into input
copy, pass/target construction, SequencePass.apply and compiled-circuit retrieval.
These post-study measurements are explanatory, not additional headline trials.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
from time import perf_counter_ns

from pytket import __version__ as tket_version
from pytket.architecture import Architecture
from pytket.circuit import OpType
from pytket.passes import (
    AutoRebase, DefaultMappingPass, DecomposeSwapsToCXs,
    FullPeepholeOptimise, SequencePass, SynthesiseTket,
)
from pytket.predicates import CompilationUnit, ConnectivityPredicate, GateSetPredicate
from pytket.qasm import circuit_from_qasm_str


def diagnostic_trial(circuit, case):
    started = perf_counter_ns()
    unit = CompilationUnit(circuit.copy())
    copied = perf_counter_ns()
    passes = [FullPeepholeOptimise(allow_swaps=True, target_2qb_gate=OpType.CX)]
    architecture = None
    if case["topology"] != "all_to_all":
        architecture = Architecture(case["edges"])
        passes.extend([DefaultMappingPass(architecture), DecomposeSwapsToCXs(architecture), SynthesiseTket()])
    passes.append(AutoRebase({OpType.CX, OpType.U3}))
    sequence = SequencePass(passes)
    configured = perf_counter_ns()
    sequence.apply(unit)
    applied = perf_counter_ns()
    output = unit.circuit
    retrieved = perf_counter_ns()
    # Predicates and structural diagnostics are outside all measured intervals.
    return {
        "copy_ns": copied - started,
        "configuration_ns": configured - copied,
        "apply_ns": applied - configured,
        "output_retrieval_ns": retrieved - applied,
        "total_ns": retrieved - started,
        "output_gates": output.n_gates,
        "output_cx": output.n_gates_of_type(OpType.CX),
        "basis_legal": GateSetPredicate({OpType.CX, OpType.U3}).verify(output),
        "connectivity_legal": architecture is None or ConnectivityPredicate(architecture).verify(output),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    protocol = json.loads((args.output_dir / "protocol.json").read_text())
    case_ids = ("ghz_8__all_to_all", "qft_8__linear")
    cases = {case["id"]: case for case in protocol["cases"]}
    result = {
        "kind": "secondary post-study timing diagnostic, excluded from primary estimates",
        "selection": "GHZ8 all-to-all (simple structured) and QFT8 linear (routed), selected for distinct workloads",
        "tket_version": tket_version,
        "warmups_per_case": 1, "measured_repeats_per_case": 3,
        "units": "nanoseconds from time.perf_counter_ns",
        "recipe": "identical pass sequence and options to frozen comparison_qiskit_tket.py; no pipeline tuning",
        "boundary": "copy + target/pass/SequencePass construction + SequencePass.apply + output retrieval; parsing/imports/export/predicates excluded",
        "limitations": "post-study diagnostic; host load may differ; stage clock calls add small overhead; not added to headline statistics",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": [],
    }
    for case_id in case_ids:
        case = cases[case_id]
        circuit = circuit_from_qasm_str(case["qasm"])
        warmup = diagnostic_trial(circuit, case)
        trials = [diagnostic_trial(circuit, case) for _ in range(3)]
        stages = ("copy_ns", "configuration_ns", "apply_ns", "output_retrieval_ns", "total_ns")
        medians = {stage: statistics.median(trial[stage] for trial in trials) for stage in stages}
        shares = {stage: statistics.median(trial[stage] / trial["total_ns"] for trial in trials) for stage in stages[:-1]}
        result["cases"].append({"case_id": case_id, "warmup": warmup, "trials": trials,
                                "medians_ns": medians, "median_per_trial_fraction": shares})
    path = args.output_dir / "tket_timing_diagnostic.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({case["case_id"]: {"medians_ns": case["medians_ns"], "median_per_trial_fraction": case["median_per_trial_fraction"]} for case in result["cases"]}, indent=2))


if __name__ == "__main__":
    main()
