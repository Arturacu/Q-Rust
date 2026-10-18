#!/usr/bin/env python3
"""Analyze paired binary revisions or routing ablations after exact-output validation.

Complete prescribed cells only; failures remain visible. Timing ratios are
reference/candidate and quality ratios candidate/reference. Confidence intervals
resample fixtures with their topology cells together, never trial repetitions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from analyze_comparison_study import (
    BOOTSTRAP_DRAWS, BOOTSTRAP_SEED, family_balanced_ratio_summary,
    geometric_ratio_summary, read_jsonl,
)

METRICS = ("compile_ns", "parse_compile_ns", "cx_count", "gate_count", "depth", "two_q_depth")


def numeric(value, positive=False):
    return type(value) in (int, float) and math.isfinite(value) and (value > 0 if positive else value >= 0)


def build_cells(protocol, timed, validated):
    cases = {c["id"]: c for c in protocol["cases"]}
    tools = protocol["tools"]
    trials = protocol["trials"]
    if len(cases) != len(protocol["cases"]) or len(set(tools)) != len(tools):
        raise ValueError("duplicate case or variant identifier")
    if type(trials) is not int or trials < 1:
        raise ValueError("prescribed trial count must be positive")
    raw = {r["id"]: r for r in timed}
    checked = {r["id"]: r for r in validated}
    if len(raw) != len(timed) or len(checked) != len(validated):
        raise ValueError("duplicate timed or validated record identifier")
    if set(checked) - set(raw):
        raise ValueError("validation contains unmeasured records")
    groups = defaultdict(list)
    for record_id, row in raw.items():
        if row["case_id"] not in cases or row["tool"] not in tools:
            raise ValueError("unknown case or variant")
        result = checked.get(record_id)
        if result is not None and any(k not in result or result[k] != v for k, v in row.items()):
            raise ValueError(f"validation changed measured record: {record_id}")
        reasons = []
        validation = result.get("validation", {}) if result else {}
        if row.get("status") != "compiled":
            reasons.append("compilation failed")
        if validation.get("passed") is not True:
            reasons.append("independent validation missing or failed")
        if not numeric(row.get("compile_ns"), positive=True) or not numeric(row.get("parse_ns")):
            reasons.append("invalid timing")
        metrics = validation.get("metrics", {})
        if any(not numeric(metrics.get(m)) for m in METRICS if not m.endswith("_ns")):
            reasons.append("invalid output metrics")
        if row.get("status") == "compiled":
            output = row.get("out_qasm")
            if not isinstance(output, str) or hashlib.sha256(output.encode()).hexdigest() != row.get("output_sha256"):
                raise ValueError(f"output hash mismatch: {record_id}")
        groups[row["case_id"], row["tool"]].append((row, metrics, reasons))
    cells = []
    for case_id, case in cases.items():
        for tool in tools:
            rows = groups[case_id, tool]
            reasons = []
            if len(rows) != trials or {r[0].get("trial") for r in rows} != set(range(trials)):
                reasons.append("trial identities do not exactly match prescribed repetitions")
            failures = [{"id": r["id"], "status": r.get("status"), "reasons": why} for r, _, why in rows if why]
            if failures:
                reasons.append("one or more trials failed")
            cell = {"case_id": case_id, "fixture": case["fixture"], "family": case["family"],
                    "topology": case["topology"], "n_logical": case["n_logical"], "tool": tool,
                    "eligible": not reasons, "expected_trials": trials, "received_trials": len(rows),
                    "exclusion_reasons": reasons, "failed_trials": failures,
                    "status_counts": dict(Counter(r.get("status") for r, _, _ in rows)),
                    "medians": None, "ranges": None}
            if not reasons:
                values = {
                    "compile_ns": [r["compile_ns"] for r, _, _ in rows],
                    "parse_compile_ns": [r["parse_ns"] + r["compile_ns"] for r, _, _ in rows],
                    **{m: [metrics[m] for _, metrics, _ in rows] for m in METRICS if not m.endswith("_ns")},
                }
                cell["medians"] = {m: statistics.median(v) for m, v in values.items()}
                cell["ranges"] = {m: {"min": min(v), "max": max(v)} for m, v in values.items()}
            cells.append(cell)
    return cells


def summarize_pairs(cells, reference, candidate, *, label, draws=BOOTSTRAP_DRAWS):
    index = {(c["case_id"], c["tool"]): c for c in cells if c["eligible"]}
    pairs = [(c, index[case, reference]) for (case, tool), c in index.items()
             if tool == candidate and (case, reference) in index]
    result = {"paired_cases": len(pairs), "fixtures": len({c["fixture"] for c, _ in pairs}), "per_case": []}
    for metric in METRICS:
        timing = metric.endswith("_ns")
        observations = [{"fixture": c["fixture"], "family": c["family"],
                         "numerator": (r if timing else c)["medians"][metric],
                         "denominator": (c if timing else r)["medians"][metric]} for c, r in pairs]
        key = f"{label}/{candidate}/{metric}"
        stats = geometric_ratio_summary(observations, label=key, draws=draws)
        stats["family_balanced"] = family_balanced_ratio_summary(observations, label=key, draws=draws)
        stats["candidate_win_tie_loss"] = dict(Counter({"win": 0, "tie": 0, "loss": 0}) + Counter(
            "win" if c["medians"][metric] < r["medians"][metric] else
            "loss" if c["medians"][metric] > r["medians"][metric] else "tie" for c, r in pairs))
        # Keep zero-count outcomes explicit instead of Counter dropping them.
        stats["candidate_win_tie_loss"] = {k: stats["candidate_win_tie_loss"].get(k, 0) for k in ("win", "tie", "loss")}
        stats["ratio_definition"] = "reference / candidate; >1 means candidate faster" if timing else "candidate / reference; <1 means candidate fewer/shallower"
        result[metric] = stats
    for c, r in pairs:
        result["per_case"].append({"case_id": c["case_id"], "fixture": c["fixture"],
            "family": c["family"], "topology": c["topology"],
            "reference_medians": r["medians"], "candidate_medians": c["medians"],
            "ratios": {m: (r["medians"][m] / c["medians"][m] if m.endswith("_ns") else c["medians"][m] / r["medians"][m])
                       if (c if m.endswith("_ns") else r)["medians"][m] else None for m in METRICS}})
    return result


def analyze(protocol, timed, validated, reference, draws=BOOTSTRAP_DRAWS):
    if reference not in protocol["tools"]:
        raise ValueError("reference variant is absent from protocol")
    cells = build_cells(protocol, timed, validated)
    comparisons = {}
    for candidate in protocol["tools"]:
        if candidate == reference:
            continue
        comparisons[candidate] = {"aggregate": summarize_pairs(cells, reference, candidate, label="aggregate", draws=draws)}
        for grouping in ("topology", "family"):
            comparisons[candidate][f"by_{grouping}"] = {
                group: summarize_pairs([c for c in cells if c[grouping] == group], reference, candidate,
                                       label=f"{grouping}/{group}", draws=draws)
                for group in sorted({c[grouping] for c in cells})}
    return {
        "schema": 1, "reference": reference,
        "method": {"cell_estimator": "median of all prescribed complete validated trials",
                   "aggregate_estimator": "equal-cell geometric mean of positive ratios",
                   "bootstrap": {"unit": "fixture with related topology cells together", "draws": draws,
                                 "seed": BOOTSTRAP_SEED, "interval": "95% percentile; none for fewer than two fixtures"},
                   "zero_policy": "zero numerator/denominator counted separately; excluded from log means",
                   "limitations": "Fixed development corpus. Intervals do not capture machine/session uncertainty or resample trials. Profile-instrumented variants are diagnostic."},
        "coverage": {"planned_cases": len(protocol["cases"]), "variants": protocol["tools"],
                     "expected_trials": len(protocol["cases"])*len(protocol["tools"])*protocol["trials"],
                     "timed_trials": len(timed), "validation_rows": len(validated),
                     "eligible_cells": {t: sum(c["eligible"] and c["tool"] == t for c in cells) for t in protocol["tools"]},
                     "excluded_cell_count": sum(not c["eligible"] for c in cells)},
        "comparisons": comparisons, "cells": cells,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--draws", type=int, default=BOOTSTRAP_DRAWS)
    args = parser.parse_args()
    protocol = json.loads((args.directory / "protocol.json").read_text())
    result = analyze(protocol, read_jsonl(args.directory / "timed.jsonl"),
                     read_jsonl(args.directory / "validated.jsonl"), args.reference, args.draws)
    result["provenance"] = {name: hashlib.sha256((args.directory/name).read_bytes()).hexdigest()
                            for name in ("protocol.json", "timed.jsonl", "validated.jsonl")}
    (args.directory / "summary.json").write_text(json.dumps(result, indent=2)+"\n")
    with (args.directory / "cell_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case_id", "fixture", "family", "topology", "n_logical", "tool", "eligible", *METRICS])
        writer.writeheader()
        for cell in result["cells"]:
            writer.writerow({key: cell[key] for key in writer.fieldnames if key not in METRICS} | (cell["medians"] or {}))
    print(json.dumps({"coverage": result["coverage"], "aggregate": {
        c: {m: {k: v[k] for k in ("geomean", "ci95", "candidate_win_tie_loss")} for m, v in data["aggregate"].items() if m in METRICS}
        for c, data in result["comparisons"].items()}}, indent=2))


if __name__ == "__main__":
    main()
