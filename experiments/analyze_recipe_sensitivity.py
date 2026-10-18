#!/usr/bin/env python3
"""Aggregate a separately labeled, interleaved recipe sensitivity campaign.

Q-Rust speed ratios are alternative recipe / Q-Rust; output ratios are Q-Rust /
alternative. Never substitutes alternative rows into the primary comparison.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from analyze_comparison_study import BOOTSTRAP_DRAWS, read_jsonl, summarize_pairs
from analyze_routing_revisions import build_cells


def convert_record(row):
    result = dict(row)
    result.update(id=f"{row['case_id']}__{row['recipe']}__{row['trial']}", tool=row["recipe"])
    result.setdefault("status", "compiled")
    if "out_qasm" in result:
        result.setdefault("output_sha256", hashlib.sha256(result["out_qasm"].encode()).hexdigest())
    return result


def sensitivity_quality(cells, tools):
    grouped = defaultdict(dict)
    for c in cells:
        if c["eligible"]:
            grouped[c["case_id"]][c["tool"]] = c
    complete = {key: row for key, row in grouped.items() if set(row) == set(tools)}
    counts = {}
    for metric in ("cx_count", "depth", "two_q_depth"):
        bands = Counter({"at_best": 0, "within_5_percent": 0, "within_10_percent": 0,
                         "best_zero": 0, "all_zero": 0})
        for row in complete.values():
            qrust = row["qrust"]["medians"][metric]
            best = min(c["medians"][metric] for t, c in row.items() if t != "qrust")
            for label, factor in (("at_best", 1), ("within_5_percent", 1.05), ("within_10_percent", 1.1)):
                bands[label] += qrust <= best * factor
            bands["best_zero"] += best == 0
            bands["all_zero"] += best == 0 and qrust == 0
        counts[metric] = dict(bands)
    return {"complete_cases": len(complete),
            "baseline_recipes": [t for t in tools if t != "qrust"],
            "definition": "Q-Rust versus best of the selected alternative recipes ONLY; excludes Qiskit and must not be called all-compiler near-best coverage",
            "near_best_counts": counts}


def analyze(metadata, timed, validated, draws=BOOTSTRAP_DRAWS):
    tools = metadata["recipes"]
    if "qrust" not in tools or len(tools) < 2:
        raise ValueError("sensitivity aggregation requires Q-Rust and at least one alternative recipe")
    protocol = {"cases": metadata["cases"], "tools": tools, "trials": metadata["measured_repeats"]}
    measured = [convert_record(r) for r in timed if r["kind"] == "measured"]
    checked = [convert_record(r) for r in validated if r["kind"] == "measured"]
    cells = build_cells(protocol, measured, checked)
    baselines = [t for t in tools if t != "qrust"]
    pairwise = {"aggregate": {t: summarize_pairs(cells, t, label="sensitivity/aggregate", draws=draws) for t in baselines}}
    for grouping in ("topology", "family"):
        pairwise[f"by_{grouping}"] = {group: {
            t: summarize_pairs([c for c in cells if c[grouping] == group], t,
                               label=f"sensitivity/{grouping}/{group}", draws=draws) for t in baselines}
            for group in sorted({c[grouping] for c in cells})}
    artifacts = {r["validation_cache_key"]: r["validation"] for r in validated
                 if r.get("validation", {}).get("passed") is True and "validation_cache_key" in r}
    return {"schema": 1, "kind": "secondary exploratory full-corpus recipe sensitivity",
            "method": {"estimator": "complete-cell medians; paired geometric mean; fixture-cluster percentile intervals",
                       "draws": draws, "speed_ratio": "alternative recipe / Q-Rust; >1 means Q-Rust faster",
                       "quality_ratio": "Q-Rust / alternative recipe; <1 means fewer gates/shallower",
                       "limitations": "Recipes selected after diagnostic cases. Same development corpus, no holdout. Intervals do not capture machine/session variation. Primary comparator recipes remain separately reported."},
            "coverage": {"cases": len(protocol["cases"]), "fixtures": len({c["fixture"] for c in protocol["cases"]}),
                         "recipes": tools, "measured_rows": len(measured), "all_rows_including_warmups_probes": len(timed),
                         "expected_measured_rows": len(protocol["cases"])*len(tools)*protocol["trials"],
                         "eligible_cells": {t: sum(c["eligible"] and c["tool"] == t for c in cells) for t in tools},
                         "excluded_cells": sum(not c["eligible"] for c in cells),
                         "measured_evidence_counts": dict(Counter(r.get("validation", {}).get("evidence_class", "missing") for r in checked)),
                         "unique_semantic_artifacts_including_warmups": len(artifacts),
                         "unique_artifact_evidence_counts": dict(Counter(v["evidence_class"] for v in artifacts.values()))},
            "pairwise": pairwise, "selected_recipe_quality": sensitivity_quality(cells, tools), "cells": cells}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--draws", type=int, default=BOOTSTRAP_DRAWS)
    args = parser.parse_args()
    metadata = json.loads((args.directory/"tradeoffs_protocol.json").read_text())
    result = analyze(metadata, read_jsonl(args.directory/"tradeoffs_timed.jsonl"),
                     read_jsonl(args.directory/"tradeoffs_validated.jsonl"), args.draws)
    result["provenance"] = {name: hashlib.sha256((args.directory/name).read_bytes()).hexdigest()
                            for name in ("tradeoffs_protocol.json", "tradeoffs_timed.jsonl", "tradeoffs_validated.jsonl")}
    (args.directory/"sensitivity_aggregate.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"coverage": result["coverage"], "aggregate": result["pairwise"]["aggregate"],
                     "selected_recipe_quality": result["selected_recipe_quality"]}, indent=2))


if __name__ == "__main__":
    main()
