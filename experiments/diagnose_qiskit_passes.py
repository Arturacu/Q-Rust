#!/usr/bin/env python3
"""Observe the unchanged Qiskit recipe with pass callbacks, then validate outputs.

Callback metric collection adds overhead, so total diagnostic latency is not a
headline runtime. Qiskit's own per-pass elapsed times and gate counts explain
where work is performed and where output costs change.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import statistics

from comparison_study import digest, THREAD_ENV


def main():
    import os
    os.environ.update(THREAD_ENV)
    import qiskit
    from comparison_qiskit_tket import compile_one
    from comparison_validation import validate_output
    from diagnose_compiler_tradeoffs import CASE_IDS

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protocol", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    protocol = json.loads(args.protocol.read_text())
    indexed = {c["id"]: c for c in protocol["cases"]}
    original_transpile = qiskit.transpile
    rows = []
    for case_id in CASE_IDS:
        case = indexed[case_id]
        assert digest(case["source_qasm"]) == case["source_sha256"]
        assert digest(case["qasm"]) == case["normalized_sha256"]
        for trial in range(-1, 7):
            passes = []

            def callback(**kwargs):
                compiler_pass, dag = kwargs["pass_"], kwargs["dag"]
                settings = {}
                for key in ("layout_trials", "swap_trials", "max_iterations", "seed", "trials", "heuristic"):
                    value = getattr(compiler_pass, key, None)
                    if isinstance(value, (str, int, float, bool)):
                        settings[key] = value
                passes.append({"name": compiler_pass.name(),
                    "elapsed_ns": round(kwargs["time"] * 1e9),
                    "operations": dict(dag.count_ops()), "depth": dag.depth(),
                    "settings": settings})

            def profiled(*positional, **keywords):
                return original_transpile(*positional, callback=callback, **keywords)

            qiskit.transpile = profiled
            try:
                output = compile_one("qiskit", **{key: case[key] for key in
                    ("qasm", "n_logical", "physical_width", "edges", "topology")},
                    seed=20261015+max(trial, 0))
            finally:
                qiskit.transpile = original_transpile
            rows.append({"case_id": case_id, "trial": trial, "passes": passes,
                         "output_sha256": digest(output["out_qasm"]), **output})
    # All checks run after compilation so they cannot perturb subsequent timing.
    cache = {}
    for row in rows:
        case = indexed[row["case_id"]]
        key = digest(json.dumps([case["source_sha256"], row["output_sha256"],
                                 row["initial_layout"], row["final_layout"]]))
        if key not in cache:
            edges = None if case["topology"] == "all_to_all" else [
                pair for a, b in case["edges"] for pair in ((a, b), (b, a))]
            cache[key] = validate_output(case["source_qasm"], row["out_qasm"],
                case["n_logical"], case["physical_width"], edges,
                row["initial_layout"], row["final_layout"])
        row["validation"] = cache[key]
    summaries = []
    for case_id in CASE_IDS:
        measured = [r for r in rows if r["case_id"] == case_id and r["trial"] >= 0]
        elapsed = defaultdict(list)
        for row in measured:
            by_name = defaultdict(int)
            for compiler_pass in row["passes"]:
                by_name[compiler_pass["name"]] += compiler_pass["elapsed_ns"]
            for name, value in by_name.items():
                elapsed[name].append(value)
        summaries.append({"case_id": case_id,
            "all_validated": all(r["validation"]["passed"] for r in measured),
            "median_cx": statistics.median(r["validation"]["metrics"]["cx_count"] for r in measured),
            "pass_median_elapsed_ns_when_executed": {
                name: statistics.median(values) for name, values in elapsed.items()}})
    result = {"kind": "secondary instrumented Qiskit pass diagnostic",
        "primary_protocol_sha256": digest(args.protocol.read_bytes()),
        "script_sha256": digest(Path(__file__).read_bytes()), "qiskit_version": qiskit.__version__,
        "cases": [indexed[i] for i in CASE_IDS], "warmups": 1, "trials": 7,
        "limitations": "Callbacks add gate-count/depth overhead; outer compile_ns must not be compared to primary timings. Individual pass times are Qiskit's own measurements, summed by name within each trial. Conditional pass medians include executions only. Five exploratory cases, not whole-corpus evidence.",
        "summaries": summaries, "rows": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(summaries, indent=2))
    if not all(row["validation"]["passed"] for row in rows):
        raise SystemExit("diagnostic output failed independent validation")


if __name__ == "__main__":
    main()
