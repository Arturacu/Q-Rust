#!/usr/bin/env python3
"""Compare routing changes against the preserved, fully supported external study.

All six recipes must retain the full predeclared cohort in both sessions.
No harness, corpus, recipe, schedule, numerical oracle or environment changes
are permitted. Sessions remain separate; unchanged compilers are drift controls.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import statistics

from analyze_comparison_study import METRICS
from analyze_external_comparison import PROFILES, digest
from analyze_external_followup import (
    VARIABLE_METADATA, before_after_recipes, bind_validated_rows, cells_for,
    changed_hashes, cohort_summary, common_cases, load_study, recipe_cells,
    trial_changes,
)

ALLOWED_COMPILER_CHANGES = {
    "src/transpiler/layout.rs", "src/transpiler/mod.rs",
    "src/transpiler/optimization.rs", "src/transpiler/routing.rs",
    "tests/cu1_routing_regressions.rs", "tests/logical_swap_regressions.rs",
}
FOCAL_FIXTURES = (
    "mqtbench__ae_n10", "mqtbench__qpeexact_n10", "mqtbench__qpeinexact_n10",
)
FOCAL_TOPOLOGIES = ("all_to_all", "linear")
DRIFT_CONTROLS = ("qiskit", "cirq", "tket", "cirq_one_pass", "tket_light")
REPORTED_METRICS = ("compile_ns", "request_wall_ns", *METRICS)


def assert_comparable(before, after):
    if before.keys() != after.keys():
        raise ValueError("Protocol fields differ")
    different = [key for key in before if key not in VARIABLE_METADATA and before[key] != after[key]]
    if different:
        raise ValueError("Frozen study settings differ: " + ", ".join(sorted(different)))
    harness = changed_hashes(before["source_sha256"], after["source_sha256"])
    compiler = changed_hashes(before["compiler_provenance"]["compiler_source_sha256"],
                              after["compiler_provenance"]["compiler_source_sha256"])
    if harness or set(compiler) - ALLOWED_COMPILER_CHANGES:
        raise ValueError(f"Unexpected source changes: harness={harness}, compiler={compiler}")
    return {
        "harness_changed": harness, "compiler_changed": compiler,
        "before_binary_sha256": before["qrust_binary_sha256"],
        "after_binary_sha256": after["qrust_binary_sha256"],
        "machine_notes": [before.get("machine_note"), after.get("machine_note")],
        "identical_corpus_schedules_seeds_recipes_environment_oracle": True,
    }


def focal_case_ids(protocol):
    expected = {(fixture, topology) for fixture in FOCAL_FIXTURES for topology in FOCAL_TOPOLOGIES}
    selected = [case for case in protocol["cases"] if case["fixture"] in FOCAL_FIXTURES]
    if len(selected) != 6 or {(c["fixture"], c["topology"]) for c in selected} != expected:
        raise ValueError("Expected exactly six predeclared AE/QPE-10 focal cases")
    if any(c["n_logical"] != 10 for c in selected):
        raise ValueError("Focal fixtures must have ten logical qubits")
    return {case["id"] for case in selected}


def full_cohort(protocol, sessions, expected_cases):
    if set(sessions) != {1, 2} or protocol["sessions"] != 2:
        raise ValueError("Exactly two separate timing sessions are required")
    planned = {case["id"] for case in protocol["cases"]}
    if len(protocol["cases"]) != expected_cases or len(planned) != expected_cases:
        raise ValueError(f"Expected {expected_cases} distinct predeclared cases")
    result = {}
    for session, records in sessions.items():
        bind_validated_rows(*records)
        if any(row["session"] != session for rows in records for row in rows):
            raise ValueError("Row belongs to a different session")
        result[session] = {}
        for profile in PROFILES:
            cells = cells_for(protocol, records, profile)
            eligible = common_cases(cells)
            if eligible != planned:
                raise ValueError(
                    f"Full cohort lost in session {session}/{profile}: "
                    f"{len(eligible)}/{expected_cases}; missing={sorted(planned - eligible)}"
                )
            result[session][profile] = cells
    return result, planned


def ratio_record(before, after):
    """Represent zero-cost transitions explicitly without infinities or epsilon."""
    return {
        "before": before, "after": after, "delta": after - before,
        "ratio": after / before if before else None,
        "ratio_class": ("both_zero" if after == 0 else "denominator_zero") if before == 0
                       else ("numerator_zero" if after == 0 else "positive"),
        "outcome": "improved" if after < before else "regressed" if after > before else "unchanged",
    }


def case_changes(before_cells, after_cells, before_rows, after_rows):
    left, right = recipe_cells(before_cells), recipe_cells(after_cells)
    if left.keys() != right.keys():
        raise ValueError("Before/after recipe-case grids differ")
    request_medians = []
    for rows in (before_rows, after_rows):
        grouped = defaultdict(list)
        for row in rows:
            grouped[row["case_id"], row["tool"]].append(row["request_wall_ns"])
        request_medians.append({key: statistics.median(values) for key, values in grouped.items()})
    detail = defaultdict(list)
    for key in sorted(left):
        old, new = left[key], right[key]
        if not old["eligible"] or not new["eligible"]:
            raise ValueError("Incomplete cell cannot enter the fixed-cohort comparison")
        metrics = {}
        for metric in REPORTED_METRICS:
            values = [side[key] for side in request_medians] if metric == "request_wall_ns" else [
                old["medians"][metric], new["medians"][metric]]
            metrics[metric] = ratio_record(*values)
        detail[key[1]].append({
            "case_id": key[0], "fixture": old["fixture"], "topology": old["topology"],
            "n_logical": old["n_logical"], "metrics": metrics,
        })
    result = {}
    for recipe, cells in detail.items():
        counts = {}
        for metric in REPORTED_METRICS:
            outcomes = Counter({"improved": 0, "unchanged": 0, "regressed": 0})
            outcomes.update(cell["metrics"][metric]["outcome"] for cell in cells)
            counts[metric] = dict(outcomes)
        result[recipe] = {
            "paired_cases": len(cells), "outcomes": counts, "cells": cells,
            "quality_adverse_cases": [cell for cell in cells if any(
                cell["metrics"][metric]["outcome"] == "regressed" for metric in METRICS)],
            "compile_adverse_cases": [cell for cell in cells
                                      if cell["metrics"]["compile_ns"]["outcome"] == "regressed"],
        }
    return result


def compare(before_protocol, after_protocol, before_sessions, after_sessions, draws=2000,
            expected_cases=100):
    if draws < 1:
        raise ValueError("Positive bootstrap draws required")
    integrity = assert_comparable(before_protocol, after_protocol)
    old_cells, planned = full_cohort(before_protocol, before_sessions, expected_cases)
    new_cells, new_planned = full_cohort(after_protocol, after_sessions, expected_cases)
    if planned != new_planned:
        raise ValueError("Before/after cohorts differ")
    focal = focal_case_ids(before_protocol)
    result = {
        "schema": 1, "analysis": "routing_followup", "integrity": integrity,
        "ratio_definition": "after / before; below one means shorter time or lower output cost",
        "timing_caveat": "Sessions remain separate. Unchanged Qiskit, Cirq and TKET recipes are drift controls. "
                         "Before/after timing alone does not establish causal attribution or remove host drift.",
        "quality_caveat": "CX count and depth are circuit-cost proxies, not execution fidelity or optimality. "
                          "A total-depth reduction does not imply a CX-depth reduction; both are retained.",
        "cohorts": {
            "original_common_case_ids": sorted(planned), "followup_common_case_ids": sorted(planned),
            "newly_common_case_ids": [], "lost_original_common_case_ids": [],
            "focal_case_ids": sorted(focal),
        },
        "sessions": {},
    }
    for session in (1, 2):
        before, after = before_sessions[session], after_sessions[session]
        recipes = before_after_recipes(before_protocol, old_cells[session], new_cells[session],
                                       before[0], after[0], draws, planned)
        focal_recipes = before_after_recipes(before_protocol, old_cells[session], new_cells[session],
                                             before[0], after[0], draws, focal)
        details = case_changes(old_cells[session], new_cells[session], before[0], after[0])
        profiles = {}
        for profile in PROFILES:
            old = cohort_summary(old_cells[session][profile], planned, f"s{session}/{profile}/old-fixed", draws)
            new = cohort_summary(new_cells[session][profile], planned, f"s{session}/{profile}/new-fixed", draws)
            profiles[profile] = {
                "original_on_original_common": old, "followup_on_original_common": new,
                "followup_on_all_common": new,
                "followup_on_newly_common": cohort_summary(new_cells[session][profile], set(), "empty", draws),
                "original_on_focal": cohort_summary(old_cells[session][profile], focal,
                                                    f"s{session}/{profile}/old-focal", draws),
                "followup_on_focal": cohort_summary(new_cells[session][profile], focal,
                                                    f"s{session}/{profile}/new-focal", draws),
            }
        result["sessions"][str(session)] = {
            "recipes": recipes, "recipes_on_original_common": recipes,
            "recipes_on_focal": focal_recipes, "profiles": profiles,
            "trial_changes": trial_changes(before[1], after[1]), "case_changes": details,
            "focal_cases": {case_id: {recipe: next(cell for cell in value["cells"] if cell["case_id"] == case_id)
                                      for recipe, value in details.items()} for case_id in sorted(focal)},
            "unchanged_recipe_drift_controls": list(DRIFT_CONTROLS),
        }
    return result


def markdown(result):
    count = len(result["cohorts"]["original_common_case_ids"])
    lines = ["# External routing follow-up comparison", "",
             f"The same {count} cases passed every trial for all six recipes in both studies and both sessions.", "",
             "Ratios are after/before geometric means of complete case medians. "
             "Zero-cost pairs remain in outcome counts and outside positive-ratio geometric means.", "",
             "| Cohort | Session | Recipe | Cases | Compile | Request | CX | Depth | CX depth |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for session, info in result["sessions"].items():
        for field, label in (("recipes", "Full fixed"), ("recipes_on_focal", "AE/QPE-10")):
            for recipe, stats in info[field].items():
                values = [stats["ratios"][metric]["geomean"] for metric in
                          ("compile_ns", "request_wall_ns", "cx_count", "depth", "two_q_depth")]
                shown = ["—" if value is None else f"{value:.5g}" for value in values]
                lines.append(f"| {label} | {session} | {recipe} | {stats['paired_complete_cases']} | "
                             + " | ".join(shown) + " |")
    lines += ["", result["timing_caveat"], "", result["quality_caveat"], "",
              "The JSON retains all case-level medians, adverse cases, trial changes, fixture-cluster intervals, "
              "and primary/light baseline comparisons. Primary/light profiles reuse Q-Rust/Qiskit observations.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("followup", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    args = parser.parse_args()
    old, old_sha, old_sessions = load_study(args.original)
    new, new_sha, new_sessions = load_study(args.followup)
    result = compare(old, new, old_sessions, new_sessions, args.bootstrap_draws)
    result.update(original_protocol_sha256=old_sha, followup_protocol_sha256=new_sha,
                  analysis_source_sha256=digest(Path(__file__)),
                  analysis_dependency_sha256={name: digest(Path(__file__).with_name(name)) for name in (
                      "analyze_external_followup.py", "analyze_external_comparison.py",
                      "analyze_comparison_study.py", "external_comparison.py", "external_corpus.py")})
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "routing_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (args.output / "ROUTING_COMPARISON.md").write_text(markdown(result))
    print(json.dumps({key: len(value) for key, value in result["cohorts"].items()}))


if __name__ == "__main__":
    main()
