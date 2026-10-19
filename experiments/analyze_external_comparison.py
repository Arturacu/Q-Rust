#!/usr/bin/env python3
"""Analyze interleaved external recipes without pooling timing sessions."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from analyze_comparison_study import (
    BASELINES, METRICS, TOOLS, analyze, geometric_ratio_summary,
    narrative, summarize_pairs, write_csv,
)

PROFILES = {
    "primary": {"qrust": "qrust", "qiskit": "qiskit", "cirq": "cirq", "tket": "tket"},
    "light": {"qrust": "qrust", "qiskit": "qiskit", "cirq_one_pass": "cirq", "tket_light": "tket"},
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def verify_session_files(directory, protocol_sha, session, validation_receipt):
    info = json.loads((directory / "session.json").read_text())
    if not info["completed"] or info["session"] != session or info["protocol_sha256"] != protocol_sha:
        raise ValueError("Session completion or identity differs from the frozen protocol")
    for name in ("timed", "warmup"):
        if digest(directory / f"{name}.jsonl") != info[f"{name}_sha256"]:
            raise ValueError(f"Session {name} records changed after timing")
    if digest(directory / "validated.jsonl") != validation_receipt["validated_sha256"]:
        raise ValueError("Validated rows changed after validation")


def project_rows(rows, profile):
    """Keep exact shared observations; only relabel recipe identities for analysis."""
    mapping = PROFILES[profile]
    result = []
    seen = set()
    for row in rows:
        if row["tool"] not in mapping:
            continue
        projected = {**row, "raw_tool": row["tool"], "raw_trial_id": row["id"],
                     "tool": mapping[row["tool"]]}
        projected["id"] = f"{row['case_id']}__{projected['tool']}__{row['trial']}"
        if projected["id"] in seen:
            raise ValueError("Duplicate trial in recipe profile")
        seen.add(projected["id"])
        result.append(projected)
    return result


def suite_of(case):
    suite = case.get("suite") or case.get("source_suite") or case.get("community_source")
    aliases = {"QASMBench": "qasmbench", "MQT Bench": "mqt", "mqtbench": "mqt"}
    suite = aliases.get(suite, suite)
    if suite is None:
        prefix = case["fixture"].split("_")[0].lower()
        prefix = aliases.get(prefix, prefix)
        if prefix not in {"qasmbench", "mqt"}:
            raise ValueError(f"Missing suite metadata: {case['id']}")
        suite = prefix
    return suite


def width_bin(width):
    return "2-6" if width <= 6 else "7-10" if width <= 10 else "11-14"


def analyze_profile(protocol, timed, validated, profile, draws=2000):
    view = {**protocol, "tools": list(TOOLS)}
    summary, cells = analyze(view, project_rows(timed, profile), project_rows(validated, profile),
                             expected_trials=protocol["trials"], draws=draws)
    cases = {case["id"]: case for case in protocol["cases"]}
    for cell in cells:
        cell["suite"] = suite_of(cases[cell["case_id"]])
        cell["width_bin"] = width_bin(cell["n_logical"])
    summary["recipe_profile"] = profile
    summary["recipe_mapping"] = PROFILES[profile]
    summary["shared_observations"] = "Q-Rust and Qiskit rows are identical observations in both profile views, not additional trials."
    for field in ("suite", "width_bin"):
        summary["pairwise"]["by_" + field] = {
            group: {tool: summarize_pairs([cell for cell in cells if cell[field] == group], tool,
                        label=f"{profile}/{field}/{group}", draws=draws) for tool in BASELINES}
            for group in sorted({cell[field] for cell in cells})}
    return summary, cells


def recipe_sensitivity(primary, light, draws=2000):
    """Compare baseline recipes on the same session and paired successful cells."""
    a = {(c["case_id"], c["tool"]): c for c in primary if c["eligible"]}
    b = {(c["case_id"], c["tool"]): c for c in light if c["eligible"]}
    results = {}
    for tool in ("cirq", "tket"):
        pairs = [(cell, b[key]) for key, cell in a.items() if key[1] == tool and key in b]
        stats = {"paired_cells": len(pairs)}
        for metric in ("compile_ns", *METRICS):
            timing = metric == "compile_ns"
            observations = [{"fixture": x["fixture"],
                "numerator": (x if timing else y)["medians"][metric],
                "denominator": (y if timing else x)["medians"][metric]} for x, y in pairs]
            stats[metric] = geometric_ratio_summary(observations,
                label=f"recipe/{tool}/{metric}", draws=draws)
            stats[metric]["definition"] = "primary / light runtime" if timing else "light / primary output cost"
        results[tool] = stats
    return results


def session_agreement(protocol, first_rows, second_rows, profile_cells, draws=2000):
    def keyed(rows):
        result = {}
        for row in rows:
            key = (row["case_id"], row["tool"], row["trial"])
            if key in result:
                raise ValueError("Duplicate raw session trial")
            result[key] = row
        return result
    left, right = keyed(first_rows), keyed(second_rows)
    if left.keys() != right.keys():
        raise ValueError("Session grids differ")
    outcome = Counter()
    per_tool = {tool: Counter() for tool in protocol["tools"]}
    changed = []
    for key, a in left.items():
        b = right[key]
        if a["status"] != b["status"]:
            state = "compile_status_changed"
        elif a["status"] != "compiled":
            state = "same_nonsuccess_status"
        elif a["validation"].get("passed") != b["validation"].get("passed"):
            state = "validation_status_changed"
        elif a["validation"].get("metrics") != b["validation"].get("metrics"):
            state = "quality_changed"
        elif a["artifact_sha256"] != b["artifact_sha256"]:
            state = "same_quality_different_artifact"
        else:
            state = "identical_artifact"
        outcome[state] += 1
        per_tool[key[1]][state] += 1
        if state in {"compile_status_changed", "validation_status_changed", "quality_changed"}:
            changed.append({"case_id": key[0], "tool": key[1], "trial": key[2], "difference": state})
    timing = {}
    for profile, (first, second) in profile_cells.items():
        one = {(c["case_id"], c["tool"]): c for c in first if c["eligible"]}
        two = {(c["case_id"], c["tool"]): c for c in second if c["eligible"]}
        timing[profile] = {tool: geometric_ratio_summary([
            {"fixture": c["fixture"], "numerator": two[key]["medians"]["compile_ns"],
             "denominator": c["medians"]["compile_ns"]}
            for key, c in one.items() if key[1] == tool and key in two],
            label=f"sessions/{profile}/{tool}", draws=draws) for tool in TOOLS}
    return {"same_seed_comparisons": len(left), "counts": dict(outcome),
            "by_recipe": {key: dict(value) for key, value in per_tool.items()},
            "changed_trials": changed, "session_2_over_1_compile_time": timing,
            "interpretation": "Sessions are retained separately; ratios describe session variation, not an intervention or independent workload coverage."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--bootstrap-draws", type=int, default=2000)
    args = parser.parse_args()
    from external_comparison import check_rows, load_protocol
    protocol, protocol_sha = load_protocol(args.study)
    if protocol["sessions"] != 2 or args.bootstrap_draws < 1:
        raise ValueError("This analysis requires two complete sessions and positive resampling draws")
    if protocol["recipe_profiles"] != {name: list(mapping) for name, mapping in PROFILES.items()}:
        raise ValueError("Recipe profiles differ from the frozen analysis contract")
    root = Path(__file__).resolve().parent.parent
    for name in ("experiments/analyze_external_comparison.py", "experiments/analyze_comparison_study.py"):
        if digest(root / name) != protocol["source_sha256"][name]:
            raise ValueError(f"Analysis source differs from frozen protocol: {name}")
    receipt = json.loads((args.study / "validation_summary.json").read_text())
    if receipt["protocol_sha256"] != protocol_sha:
        raise ValueError("Validation receipt belongs to another protocol")
    outdir = args.study / "analysis"
    outdir.mkdir(exist_ok=False)
    result = {"protocol_sha256": protocol_sha, "profiles_share_qrust_qiskit_observations": True,
              "sessions": {}, "source_sha256": digest(Path(__file__))}
    session_rows, profile_cells = [], {profile: [] for profile in PROFILES}
    for session in (1, 2):
        directory = args.study / f"session_{session:02d}"
        matching = next(x for x in receipt["sessions"] if x["session"] == directory.name)
        verify_session_files(directory, protocol_sha, session, matching)
        timed, validated = read_rows(directory / "timed.jsonl"), read_rows(directory / "validated.jsonl")
        check_rows(protocol, session, timed)
        check_rows(protocol, session, validated)
        session_rows.append(validated)
        result["sessions"][str(session)] = {}
        for profile in PROFILES:
            summary, cells = analyze_profile(protocol, timed, validated, profile, args.bootstrap_draws)
            profile_cells[profile].append(cells)
            result["sessions"][str(session)][profile] = summary
            write_csv(outdir / f"session_{session}_{profile}_cells.csv", cells)
            (outdir / f"session_{session}_{profile}.md").write_text(narrative(summary))
        result["sessions"][str(session)]["recipe_sensitivity"] = recipe_sensitivity(
            profile_cells["primary"][-1], profile_cells["light"][-1], args.bootstrap_draws)
    result["session_agreement"] = session_agreement(protocol, *session_rows, profile_cells, args.bootstrap_draws)
    (outdir / "external_summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"sessions": 2, "profiles": list(PROFILES),
                      "session_agreement": result["session_agreement"]["counts"]}))


if __name__ == "__main__":
    main()
