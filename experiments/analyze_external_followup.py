#!/usr/bin/env python3
"""Compare preserved external studies, separating old-common and newly supported cases."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import statistics

from analyze_comparison_study import BASELINES, METRICS, TOOLS, build_cells, geometric_ratio_summary, summarize_pairs
from analyze_external_comparison import PROFILES, digest, project_rows, read_rows, verify_session_files

ALLOWED_HARNESS_CHANGES = {"experiments/comparison_cirq.py"}
ALLOWED_COMPILER_CHANGES = {
    "src/ir/gates.rs", "src/ir/gate_def.rs", "src/transpiler/mod.rs",
    "tests/cu1_regressions.rs",
}
VARIABLE_METADATA = {
    "created_utc", "compiler_commit", "compiler_provenance", "qrust_binary_sha256",
    "source_sha256", "machine_note",
}


def changed_hashes(before, after):
    return sorted(key for key in before.keys() | after.keys() if before.get(key) != after.get(key))


def assert_comparable(before, after):
    """Reject changes to the corpus, schedule, seeds, recipes, environment or oracle."""
    if before.keys() != after.keys():
        raise ValueError("Protocol fields differ")
    different = [key for key in before if key not in VARIABLE_METADATA and before[key] != after[key]]
    if different:
        raise ValueError("Frozen study settings differ: " + ", ".join(sorted(different)))
    harness = changed_hashes(before["source_sha256"], after["source_sha256"])
    compiler = changed_hashes(before["compiler_provenance"]["compiler_source_sha256"],
                              after["compiler_provenance"]["compiler_source_sha256"])
    if set(harness) - ALLOWED_HARNESS_CHANGES or set(compiler) - ALLOWED_COMPILER_CHANGES:
        raise ValueError(f"Unexpected source changes: harness={harness}, compiler={compiler}")
    return {"harness_changed": harness, "compiler_changed": compiler,
            "before_binary_sha256": before["qrust_binary_sha256"],
            "after_binary_sha256": after["qrust_binary_sha256"],
            "machine_notes": [before.get("machine_note"), after.get("machine_note")],
            "identical_corpus_schedules_seeds_recipes_environment_oracle": True}


def bind_validated_rows(timed, validated):
    """Validation may add evidence, but may not change any measured field."""
    if len(timed) != len(validated):
        raise ValueError("Timed/validated row counts differ")
    for raw, checked in zip(timed, validated):
        if any(key not in checked or checked[key] != value for key, value in raw.items()):
            raise ValueError("Validation changed a measured field")
        if set(checked) - set(raw) - {"validation", "validation_cache_key"}:
            raise ValueError("Unexpected validation metadata")


def load_study(path):
    from external_comparison import check_rows, load_protocol
    if (path / "active.lock").exists():
        raise ValueError("Study is still running")
    protocol, checksum = load_protocol(path)
    if protocol["sessions"] != 2:
        raise ValueError("Expected two separately retained sessions")
    receipt = json.loads((path / "validation_summary.json").read_text())
    if receipt["protocol_sha256"] != checksum:
        raise ValueError("Validation receipt belongs to another protocol")
    if receipt["oracle_sha256"] != protocol["source_sha256"]["experiments/comparison_validation.py"]:
        raise ValueError("Validation oracle differs from frozen protocol")
    root = Path(__file__).resolve().parent.parent
    for name in ("experiments/external_comparison.py", "experiments/external_corpus.py",
                 "experiments/analyze_external_comparison.py", "experiments/analyze_comparison_study.py"):
        if digest(root / name) != protocol["source_sha256"][name]:
            raise ValueError("Analysis dependency differs from frozen source: " + name)
    result = {}
    for session in (1, 2):
        directory = path / f"session_{session:02d}"
        matches = [row for row in receipt["sessions"] if row["session"] == directory.name]
        if len(matches) != 1:
            raise ValueError("Missing or duplicate validation receipt")
        verify_session_files(directory, checksum, session, matches[0])
        timed = read_rows(directory / "timed.jsonl")
        validated = read_rows(directory / "validated.jsonl")
        check_rows(protocol, session, timed)
        check_rows(protocol, session, validated)
        bind_validated_rows(timed, validated)
        result[session] = (timed, validated)
    return protocol, checksum, result


def cells_for(protocol, records, profile):
    timed, validated = records
    return build_cells(protocol, project_rows(timed, profile), project_rows(validated, profile),
                       expected_trials=protocol["trials"])


def common_cases(cells):
    eligible = defaultdict(set)
    for cell in cells:
        if cell["eligible"]:
            eligible[cell["case_id"]].add(cell["tool"])
    return {case for case, tools in eligible.items() if tools == set(TOOLS)}


def keyed(rows):
    result = {}
    for row in rows:
        key = (row["case_id"], row["tool"], row["trial"])
        if key in result:
            raise ValueError("Duplicate trial identity")
        result[key] = row
    return result


def trial_changes(before, after):
    a, b = keyed(before), keyed(after)
    if a.keys() != b.keys():
        raise ValueError("Before/after trial grids differ")
    counts, changes = defaultdict(Counter), []
    for key, left in a.items():
        right = b[key]
        if left["seed"] != right["seed"] or left["session"] != right["session"]:
            raise ValueError("Before/after seed or session pairing differs")
        fields = []
        if left["status"] != right["status"]:
            fields.append("compile_status")
        lv, rv = left.get("validation") or {}, right.get("validation") or {}
        if lv.get("passed") != rv.get("passed"):
            fields.append("validation_status")
        if left["status"] == right["status"] == "compiled":
            if lv.get("metrics") != rv.get("metrics"):
                fields.append("quality")
            if left["artifact_sha256"] != right["artifact_sha256"]:
                fields.append("artifact")
        counts[key[1]]["paired_trials"] += 1
        for field in fields:
            counts[key[1]][field + "_changed"] += 1
        if not fields:
            counts[key[1]]["unchanged_outcome"] += 1
        else:
            changes.append({"case_id": key[0], "tool": key[1], "trial": key[2], "seed": left["seed"],
                            "changed": fields, "before_status": left["status"], "after_status": right["status"],
                            "before_metrics": lv.get("metrics"), "after_metrics": rv.get("metrics")})
    return {"by_recipe": {tool: dict(value) for tool, value in counts.items()}, "changed_trials": changes}


def recipe_cells(profile_cells):
    result = {}
    for profile, cells in profile_cells.items():
        reverse = {value: key for key, value in PROFILES[profile].items()}
        for cell in cells:
            result[cell["case_id"], reverse[cell["tool"]]] = cell
    return result


def before_after_recipes(protocol, before_cells, after_cells, before_rows, after_rows, draws, case_ids=None):
    a, b = recipe_cells(before_cells), recipe_cells(after_cells)
    if case_ids is not None:
        a = {key: cell for key, cell in a.items() if key[0] in case_ids}
        b = {key: cell for key, cell in b.items() if key[0] in case_ids}
    grouped = []
    for rows in (before_rows, after_rows):
        index = defaultdict(list)
        for row in rows:
            index[row["case_id"], row["tool"]].append(row)
        grouped.append(index)
    results = {}
    for tool in protocol["tools"]:
        pairs = [(key, cell, b[key]) for key, cell in a.items()
                 if key[1] == tool and cell["eligible"] and b[key]["eligible"]]
        stats = {"planned_cases": len(protocol["cases"]) if case_ids is None else len(case_ids),
                 "before_eligible_cases": sum(c["eligible"] for k, c in a.items() if k[1] == tool),
                 "after_eligible_cases": sum(c["eligible"] for k, c in b.items() if k[1] == tool),
                 "paired_complete_cases": len(pairs), "ratios": {}}
        for metric in ("compile_ns", "request_wall_ns", *METRICS):
            observations = []
            for key, left, right in pairs:
                if metric == "request_wall_ns":
                    values = [[row.get(metric) for row in index[key]] for index in grouped]
                    if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0
                           for side in values for v in side):
                        raise ValueError("Invalid complete request timing")
                    old, new = map(statistics.median, values)
                else:
                    old, new = left["medians"][metric], right["medians"][metric]
                observations.append({"fixture": left["fixture"], "numerator": new, "denominator": old})
            stats["ratios"][metric] = geometric_ratio_summary(observations, label=f"before-after/{tool}/{metric}", draws=draws)
        results[tool] = stats
    return results


def cohort_summary(cells, case_ids, label, draws):
    selected = [cell for cell in cells if cell["case_id"] in case_ids]
    return {"planned_cases": len(case_ids), "complete_all_four_cases": len(common_cases(selected)),
            "pairwise": {tool: summarize_pairs(selected, tool, label=label, draws=draws) for tool in BASELINES}}


def compare(before_protocol, after_protocol, before_sessions, after_sessions, draws=2000):
    integrity = assert_comparable(before_protocol, after_protocol)
    # This cohort is determined solely from the original study, before using
    # any follow-up outcomes. Both original sessions must be complete/valid.
    old_cells = {s: {p: cells_for(before_protocol, rows, p) for p in PROFILES}
                 for s, rows in before_sessions.items()}
    old_common = set.intersection(*(common_cases(value["primary"]) for value in old_cells.values()))
    new_cells = {s: {p: cells_for(after_protocol, rows, p) for p in PROFILES}
                 for s, rows in after_sessions.items()}
    new_common = set.intersection(*(common_cases(value["primary"]) for value in new_cells.values()))
    newly_common = new_common - old_common
    result = {"schema": 1, "integrity": integrity,
              "ratio_definition": "after / before; below one means shorter time or lower output cost",
              "timing_caveat": "Separate sessions and fixed corpora; unchanged Qiskit/TKET recipes are drift controls. No causal timing attribution or correction is made.",
              "cohorts": {"original_common_case_ids": sorted(old_common),
                          "followup_common_case_ids": sorted(new_common),
                          "newly_common_case_ids": sorted(newly_common),
                          "lost_original_common_case_ids": sorted(old_common - new_common)},
              "sessions": {}}
    for session in before_sessions:
        before, after = before_sessions[session], after_sessions[session]
        changes = trial_changes(before[1], after[1])
        recipes = before_after_recipes(before_protocol, old_cells[session], new_cells[session], before[0], after[0], draws)
        fixed_recipes = before_after_recipes(before_protocol, old_cells[session], new_cells[session],
                                            before[0], after[0], draws, case_ids=old_common)
        profiles = {}
        for profile in PROFILES:
            profiles[profile] = {
                "original_on_original_common": cohort_summary(old_cells[session][profile], old_common, f"s{session}/{profile}/old-fixed", draws),
                "followup_on_original_common": cohort_summary(new_cells[session][profile], old_common, f"s{session}/{profile}/new-fixed", draws),
                "followup_on_all_common": cohort_summary(new_cells[session][profile], new_common, f"s{session}/{profile}/new-all", draws),
                "followup_on_newly_common": cohort_summary(new_cells[session][profile], newly_common, f"s{session}/{profile}/newly-common", draws),
            }
        result["sessions"][str(session)] = {"recipes": recipes, "recipes_on_original_common": fixed_recipes,
                                             "profiles": profiles, "trial_changes": changes,
                                             "unchanged_recipe_drift_controls": ["qiskit", "tket", "tket_light"]}
    return result


def markdown(result):
    cohorts = result["cohorts"]
    lines = ["# External benchmark compatibility follow-up comparison", "",
             f"Original common cohort: {len(cohorts['original_common_case_ids'])} cases. "
             f"Follow-up common cohort: {len(cohorts['followup_common_case_ids'])}; "
             f"newly common: {len(cohorts['newly_common_case_ids'])}.", "",
             "Cohorts require all four primary recipes to pass every trial in both sessions. "
             "The original cohort uses only original results. Profiles reuse Q-Rust/Qiskit observations.", "",
             "Ratios below are follow-up/original geometric means of paired complete-case medians. "
             "Full request time includes parsing, compilation, serialization and IPC.", "",
             "Each recipe's supported cohort can differ; the fixed original cohort holds the case set constant for all recipes, including drift controls.", "",
             "| Cohort | Session | Recipe | Original / follow-up complete cases | Paired cases | Compile ratio | Request ratio | CX ratio | Depth ratio |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for session, info in result["sessions"].items():
        for field, label in (("recipes", "Recipe-supported"), ("recipes_on_original_common", "Fixed original")):
            for tool, stats in info[field].items():
                values = [stats["ratios"][m]["geomean"] for m in ("compile_ns", "request_wall_ns", "cx_count", "depth")]
                formatted = ["—" if v is None else f"{v:.4g}" for v in values]
                lines.append(f"| {label} | {session} | {tool} | {stats['before_eligible_cases']} / {stats['after_eligible_cases']} | {stats['paired_complete_cases']} | " + " | ".join(formatted) + " |")
    lines += ["", result["timing_caveat"], "",
              "The JSON retains confidence intervals, zero-cost cases, all metric and artifact changes, "
              "and separate primary/light baseline comparisons for the fixed original cohort, full follow-up cohort and newly common cohort.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("followup", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    parser.add_argument("--expected-original-common", type=int, default=82)
    args = parser.parse_args()
    if args.bootstrap_draws < 1:
        raise ValueError("Positive bootstrap draws required")
    old, old_sha, old_sessions = load_study(args.original)
    new, new_sha, new_sessions = load_study(args.followup)
    result = compare(old, new, old_sessions, new_sessions, args.bootstrap_draws)
    if len(result["cohorts"]["original_common_case_ids"]) != args.expected_original_common:
        raise ValueError("Original common cohort differs from the predeclared size")
    result.update(original_protocol_sha256=old_sha, followup_protocol_sha256=new_sha,
                  analysis_source_sha256=digest(Path(__file__)))
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "followup_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (args.output / "FOLLOWUP_COMPARISON.md").write_text(markdown(result))
    print(json.dumps({key: len(value) for key, value in result["cohorts"].items()}))


if __name__ == "__main__":
    main()
