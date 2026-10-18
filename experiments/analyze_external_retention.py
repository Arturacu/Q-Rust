#!/usr/bin/env python3
"""Compare candidate-retention changes against the preserved external routing study.

The full cohort and all six recipes stay fixed. Prespecified recovery/retention
ceilings are evaluated on every Q-Rust trial, separately from numerical validity.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from analyze_comparison_study import METRICS
from analyze_external_comparison import PROFILES, digest
from analyze_external_followup import (
    VARIABLE_METADATA, before_after_recipes, changed_hashes, cohort_summary,
    keyed, load_study, trial_changes,
)
from analyze_external_routing import case_changes, full_cohort

ALLOWED_COMPILER_CHANGES = {
    "src/transpiler/candidate_selection.rs", "src/transpiler/mod.rs",
    "src/transpiler/numerical_policy.rs", "src/transpiler/routing.rs",
    "src/transpiler/layout.rs", "tests/routing_candidate_regressions.rs",
}
FOCAL_WIDTHS = {
    "mqtbench__ae_n10": 10, "mqtbench__qpeexact_n10": 10,
    "mqtbench__qpeinexact_n10": 10, "qasmbench__basis_trotter_n4": 4,
    "qasmbench__pea_n5": 5, "qasmbench__qft_n4": 4,
}
FOCAL_TOPOLOGIES = ("all_to_all", "linear")
DRIFT_CONTROLS = ("qiskit", "cirq", "tket", "cirq_one_pass", "tket_light")
QUALITY_TARGET_METRICS = ("cx_count", "depth", "two_q_depth")
PRESPECIFIED_PLAN_SHA256 = "29799fcc11fd5b4b74cc095b0a8e9b00209bbf76220aa1ecb9e9d8511fd27416"
PRESERVATION_SHA256 = "61818c58b3bcb043dd427d0bba7173ff1edca6c00d9aef71dc3fb1c03410d38c"
IDENTITY_FIELDS = (
    "status", "out_qasm", "initial_layout", "final_layout", "output_qubits",
    "artifact_sha256",
)


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
    if before["tools"] != ["qrust", *DRIFT_CONTROLS] or before["trials"] != 7:
        raise ValueError("Expected the six frozen recipes and seven trials")
    return {
        "harness_changed": harness, "compiler_changed": compiler,
        "before_binary_sha256": before["qrust_binary_sha256"],
        "after_binary_sha256": after["qrust_binary_sha256"],
        "before_compiler_provenance": before["compiler_provenance"],
        "after_compiler_provenance": after["compiler_provenance"],
        "frozen_harness_sha256": before["source_sha256"],
        "corpus_manifest_sha256": before["corpus_manifest_sha256"],
        "machine_notes": [before.get("machine_note"), after.get("machine_note")],
        "identical_corpus_schedules_seeds_recipes_environment_oracle": True,
    }


def focal_case_ids(protocol):
    selected = [case for case in protocol["cases"] if case["fixture"] in FOCAL_WIDTHS]
    expected = {(fixture, topology) for fixture in FOCAL_WIDTHS for topology in FOCAL_TOPOLOGIES}
    if len(selected) != 12 or {(c["fixture"], c["topology"]) for c in selected} != expected:
        raise ValueError("Expected exactly twelve predeclared retention/recovery focal cases")
    if any(c["n_logical"] != FOCAL_WIDTHS[c["fixture"]] for c in selected):
        raise ValueError("Focal fixture logical widths differ")
    return {case["id"] for case in selected}


def assert_plan(plan, protocol, protocol_sha256, before_sessions):
    for key, value in (("baseline_protocol_sha256", protocol_sha256),
                       ("baseline_binary_sha256", protocol["qrust_binary_sha256"]),
                       ("corpus_manifest_sha256", protocol["corpus_manifest_sha256"])):
        if plan.get(key) != value:
            raise ValueError("Prespecified plan differs from baseline: " + key)
    focal = focal_case_ids(protocol)
    cases = {case["id"]: case for case in protocol["cases"]}
    declared = plan["focal_cases"]
    if len(declared) != 12 or {case["id"] for case in declared} != focal:
        raise ValueError("Plan focal case grid differs")
    for case in declared:
        if any(case[field] != cases[case["id"]][field] for field in (
                "fixture", "n_logical", "topology", "source_sha256", "normalized_sha256")):
            raise ValueError("Plan focal source/width/topology differs")
    ceilings = plan["predeclared_quality_ceilings"]
    if len(ceilings) != 12 or {item["case_id"] for item in ceilings} != focal:
        raise ValueError("Plan quality ceiling grid differs")
    for item in ceilings:
        case = cases[item["case_id"]]
        recover = case["topology"] == "linear" and case["fixture"].startswith("qasmbench__")
        if (item["purpose"], item["reference_study"]) != (
                ("recover_previous_quality", "compatibility") if recover else
                ("retain_current_quality", "routing")):
            raise ValueError("Plan quality ceiling reference differs")
        metrics = item["reference_metrics"]
        if set(metrics) != set(QUALITY_TARGET_METRICS) or any(
                type(value) is not int or value < 0 for value in metrics.values()):
            raise ValueError("Plan ceilings require three nonnegative integer metrics")
        for _, checked in before_sessions.values():
            for row in checked:
                if row["tool"] != "qrust" or row["case_id"] != item["case_id"]:
                    continue
                reference = row["validation"]["metrics"]
                if not recover and any(reference[m] != metrics[m] for m in metrics):
                    raise ValueError("Retention ceiling differs from measured routing baseline")
                if recover and (any(metrics[m] > reference[m] for m in metrics)
                                or all(metrics[m] == reference[m] for m in metrics)):
                    raise ValueError("Recovery ceiling does not recover a baseline regression")
    return {item["case_id"]: item for item in ceilings}


def validate_trial_grid(protocol, sessions):
    expected = {(case["id"], tool, trial) for case in protocol["cases"]
                for tool in protocol["tools"] for trial in range(protocol["trials"])}
    for session, pair in sessions.items():
        for rows in pair:
            if set(keyed(rows)) != expected:
                raise ValueError(f"Missing or extra trial cells in session {session}")
            for row in rows:
                if any(field not in row for field in IDENTITY_FIELDS):
                    raise ValueError("Missing emitted artifact identity field")


def identity_changes(left, right):
    fields = [field for field in IDENTITY_FIELDS if left[field] != right[field]]
    if left.get("output_sha256") != right.get("output_sha256"):
        fields.append("output_sha256")
    if left["validation"]["metrics"] != right["validation"]["metrics"]:
        fields.append("validated_metrics")
    return fields


def unchanged_baseline_identity(before_sessions, after_sessions):
    counts = Counter()
    for session in (1, 2):
        left, right = keyed(before_sessions[session][1]), keyed(after_sessions[session][1])
        for key in sorted(left):
            old, new = left[key], right[key]
            if (old["seed"], old["session"]) != (new["seed"], new["session"]):
                raise ValueError("Before/after seed or session pairing differs")
            if key[1] in DRIFT_CONTROLS:
                changed = identity_changes(old, new)
                if changed:
                    raise ValueError(f"Unchanged compiler baseline changed: session={session}, {key}, {changed}")
                counts[key[1]] += 1
    return {"paired_trials": sum(counts.values()), "identical_trials": sum(counts.values()),
            "by_recipe": dict(sorted(counts.items())), "all_identical": True,
            "checked_fields": [*IDENTITY_FIELDS, "output_sha256", "validated_metrics"],
            "timings_excluded_from_identity": True}


def session_identity(sessions):
    left, right = keyed(sessions[1][1]), keyed(sessions[2][1])
    if left.keys() != right.keys():
        raise ValueError("Session grids differ")
    counts, changed = Counter(), []
    for key in sorted(left):
        if left[key]["seed"] != right[key]["seed"]:
            raise ValueError("Cross-session same-seed pairing differs")
        fields = identity_changes(left[key], right[key])
        counts[key[1]] += 1
        if fields:
            changed.append({"case_id": key[0], "tool": key[1], "trial": key[2],
                            "seed": left[key]["seed"], "changed_fields": fields})
    return {"paired_trials": sum(counts.values()), "identical_trials": sum(counts.values()) - len(changed),
            "by_recipe": dict(sorted(counts.items())), "all_identical": not changed,
            "changed_trials": changed}


def quality_acceptance(ceilings, before_sessions, after_sessions):
    focal_results, all_regressions = [], []
    for session in (1, 2):
        before, after = keyed(before_sessions[session][1]), keyed(after_sessions[session][1])
        for key in sorted(after):
            if key[1] != "qrust":
                continue
            old = before[key]["validation"]["metrics"]
            new = after[key]["validation"]["metrics"]
            exceeded = {m: {"before": old[m], "after": new[m]} for m in QUALITY_TARGET_METRICS
                        if new[m] > old[m]}
            if exceeded:
                all_regressions.append({"session": session, "case_id": key[0], "trial": key[2],
                                        "seed": after[key]["seed"], "regressed_metrics": exceeded})
        for case_id, ceiling in sorted(ceilings.items()):
            rows = [row for key, row in after.items() if key[0] == case_id and key[1] == "qrust"]
            failed = [{"trial": row["trial"], "seed": row["seed"],
                       "metrics": row["validation"]["metrics"],
                       "exceeded_metrics": [m for m in QUALITY_TARGET_METRICS
                                            if row["validation"]["metrics"][m] > ceiling["reference_metrics"][m]]}
                      for row in rows if any(row["validation"]["metrics"][m] > ceiling["reference_metrics"][m]
                                             for m in QUALITY_TARGET_METRICS)]
            focal_results.append({**ceiling, "session": session, "trials_checked": len(rows),
                                  "observed_ranges": {m: {
                                      "min": min(row["validation"]["metrics"][m] for row in rows),
                                      "max": max(row["validation"]["metrics"][m] for row in rows)}
                                      for m in QUALITY_TARGET_METRICS},
                                  "passed": not failed, "failed_trials": failed})
    return {
        "metrics": list(QUALITY_TARGET_METRICS),
        "timing_has_no_predeclared_ceiling": True,
        "numerical_correctness_is_a_separate_required_gate": True,
        "focal_case_count": len(ceilings), "focal_session_cells": len(focal_results),
        "focal_trials_checked": sum(item["trials_checked"] for item in focal_results),
        "focal_targets": focal_results, "all_focal_targets_passed": all(item["passed"] for item in focal_results),
        "full_cohort_trials_checked": sum(row["tool"] == "qrust" for pair in after_sessions.values() for row in pair[1]),
        "full_cohort_regressed_trials": all_regressions,
        "full_cohort_nonworsening_passed": not all_regressions,
        "gate_count_reported_separately_from_selection_objective": True,
    }


def compare(before_protocol, after_protocol, before_sessions, after_sessions, *, plan,
            baseline_protocol_sha256, draws=2000, expected_cases=100):
    if draws < 1:
        raise ValueError("Positive bootstrap draws required")
    integrity = assert_comparable(before_protocol, after_protocol)
    for protocol, sessions in ((before_protocol, before_sessions), (after_protocol, after_sessions)):
        validate_trial_grid(protocol, sessions)
    old_cells, planned = full_cohort(before_protocol, before_sessions, expected_cases)
    new_cells, new_planned = full_cohort(after_protocol, after_sessions, expected_cases)
    if planned != new_planned:
        raise ValueError("Before/after cohorts differ")
    ceilings = assert_plan(plan, before_protocol, baseline_protocol_sha256, before_sessions)
    focal = set(ceilings)
    result = {
        "schema": 1, "analysis": "routing_retention_followup", "integrity": integrity,
        "original_protocol_sha256": baseline_protocol_sha256,
        "ratio_definition": "after / before; below one means shorter time or lower output cost",
        "timing_caveat": "Sessions remain separate. Unchanged Qiskit, Cirq and TKET recipes are drift controls. "
                         "Before/after timing alone does not establish causal attribution or remove host drift.",
        "quality_caveat": "CX count and depth are circuit-cost proxies, not execution fidelity or optimality. "
                          "This corpus informed the fixes and is not an untouched generalization holdout.",
        "cohorts": {
            "original_common_case_ids": sorted(planned), "followup_common_case_ids": sorted(planned),
            "newly_common_case_ids": [], "lost_original_common_case_ids": [], "focal_case_ids": sorted(focal),
        },
        "non_qrust_identity": unchanged_baseline_identity(before_sessions, after_sessions),
        "same_seed_session_identity": {"original": session_identity(before_sessions),
                                       "followup": session_identity(after_sessions)},
        "quality_acceptance": quality_acceptance(ceilings, before_sessions, after_sessions),
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


def check_preservation(repo, receipt_path):
    if digest(receipt_path) != PRESERVATION_SHA256:
        raise ValueError("Preservation receipt differs from the prespecified receipt")
    receipt = json.loads(receipt_path.read_text())
    if len(receipt["sha256"]) != 43:
        raise ValueError("Expected all 43 preserved artifacts")
    for relative, checksum in receipt["sha256"].items():
        path = repo / relative
        if not path.is_file() or digest(path) != checksum:
            raise ValueError("Previously retained artifact changed: " + relative)
    return {"files_checked": 43, "all_unchanged": True, "receipt_sha256": digest(receipt_path),
            "source_commit": receipt["source_commit"]}


def markdown(result):
    lines = ["# External candidate-retention comparison", "",
             f"All {len(result['cohorts']['original_common_case_ids'])} cases retained all six recipes in both sessions.", "",
             "| Cohort | Session | Recipe | Cases | Compile | Request | CX | Depth | CX depth |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for session, info in result["sessions"].items():
        for field, label in (("recipes", "Full fixed"), ("recipes_on_focal", "Retention/recovery")):
            for recipe, stats in info[field].items():
                values = [stats["ratios"][metric]["geomean"] for metric in
                          ("compile_ns", "request_wall_ns", "cx_count", "depth", "two_q_depth")]
                shown = ["—" if value is None else f"{value:.5g}" for value in values]
                lines.append(f"| {label} | {session} | {recipe} | {stats['paired_complete_cases']} | "
                             + " | ".join(shown) + " |")
    acceptance = result["quality_acceptance"]
    lines += ["", result["ratio_definition"] + ".", "", result["timing_caveat"], "", result["quality_caveat"], "",
              f"All twelve predeclared quality targets passed: {acceptance['all_focal_targets_passed']}. "
              f"Full-cohort trial-level nonworsening passed: {acceptance['full_cohort_nonworsening_passed']}.", "",
              "The JSON retains every case, adverse result, trial-level target check, artifact identity check, "
              "fixture-cluster interval, and primary/light comparison. Timing has no predeclared acceptance ceiling.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("followup", type=Path)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    parser.add_argument("--preservation", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    if digest(args.plan) != PRESPECIFIED_PLAN_SHA256:
        raise ValueError("Investigation plan differs from its prespecified SHA256")
    plan = json.loads(args.plan.read_text())
    preservation = check_preservation(args.repo, args.preservation) if args.preservation else None
    old, old_sha, old_sessions = load_study(args.original)
    new, new_sha, new_sessions = load_study(args.followup)
    result = compare(old, new, old_sessions, new_sessions, plan=plan,
                     baseline_protocol_sha256=old_sha, draws=args.bootstrap_draws)
    result.update(followup_protocol_sha256=new_sha, investigation_plan_sha256=digest(args.plan),
                  original_preservation=preservation, analysis_source_sha256=digest(Path(__file__)),
                  analysis_dependency_sha256={name: digest(Path(__file__).with_name(name)) for name in (
                      "analyze_external_routing.py", "analyze_external_followup.py", "analyze_external_comparison.py",
                      "analyze_comparison_study.py", "external_comparison.py", "external_corpus.py")})
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "retention_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (args.output / "RETENTION_COMPARISON.md").write_text(markdown(result))
    print(json.dumps({key: len(value) for key, value in result["cohorts"].items()}))


if __name__ == "__main__":
    main()
