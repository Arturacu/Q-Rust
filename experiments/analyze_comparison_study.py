#!/usr/bin/env python3
"""Analyze repeated, paired, independently validated compiler measurements.

Only complete cells with the prescribed number of successful, validated trials
are eligible. Confidence intervals resample fixtures, retaining correlated
topologies within each fixture. They describe variability within this fixed
corpus, not a claim that the fixtures are a random sample of all circuits.
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

import numpy as np


TOOLS = ("qrust", "qiskit", "cirq", "tket")
BASELINES = TOOLS[1:]
TOOL_LABELS = {"qrust": "Q-Rust", "qiskit": "Qiskit", "cirq": "Cirq", "tket": "TKET"}
METRICS = ("cx_count", "gate_count", "depth", "two_q_depth")
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20261015


def read_jsonl(path):
    if not path.exists():
        return []
    rows = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.strip():
                try:
                    rows.append(json.loads(line))
                except ValueError as error:
                    raise ValueError(f"{path}:{line_number}: invalid JSON") from error
    return rows


def _numeric(value, *, positive=False):
    return (type(value) in (int, float) and math.isfinite(value)
            and (value > 0 if positive else value >= 0))


def _percentile(values, percentile):
    return float(np.percentile(values, percentile)) if values else None


def _brief_failure(row, reasons):
    validation = row.get("validation") or {}
    return {
        "id": row.get("id"), "trial": row.get("trial"),
        "status": row.get("status"), "reasons": reasons,
        "error": row.get("error"),
        "validation": {key: value for key, value in validation.items() if key != "metrics"},
    }


def build_cells(protocol, timed, validated, expected_trials=7):
    """Retain incomplete/failed cells; never replace them with successful subsets."""
    cases = protocol["cases"]
    if len({case["id"] for case in cases}) != len(cases):
        raise ValueError("duplicate case id in protocol")
    case_ids = {case["id"] for case in cases}
    timed_by_id = {}
    for row in timed:
        if row.get("id") in timed_by_id:
            raise ValueError(f"duplicate timed trial id: {row.get('id')}")
        timed_by_id[row.get("id")] = row
    validated_by_id = {}
    for row in validated:
        if row.get("id") in validated_by_id:
            raise ValueError(f"duplicate validation trial id: {row.get('id')}")
        validated_by_id[row.get("id")] = row
    if set(validated_by_id) - set(timed_by_id):
        raise ValueError("validated trials are absent from timed.jsonl")
    grouped = defaultdict(list)
    for trial_id, timed_row in timed_by_id.items():
        row = dict(timed_row)
        checked = validated_by_id.get(trial_id)
        if checked:
            # Validation must refer to the exact timed output, not a recompile
            # or a re-timed row. Reject conflicting measurement identities.
            for key in ("case_id", "tool", "trial", "compile_ns", "parse_ns", "out_qasm", "output_sha256",
                        "initial_layout", "final_layout", "output_qubits"):
                if key in checked and key in timed_row and checked[key] != timed_row[key]:
                    raise ValueError(f"validation changed timed field {key} for {trial_id}")
            row["validation"] = checked.get("validation")
        if row.get("case_id") not in case_ids or row.get("tool") not in TOOLS:
            raise ValueError(f"trial has unknown case/tool: {trial_id}")
        grouped[row["case_id"], row["tool"]].append(row)

    cells = []
    for case in cases:
        for tool in TOOLS:
            rows = grouped[case["id"], tool]
            failures = []
            passed = []
            duplicate_trials = len({row.get("trial") for row in rows}) != len(rows)
            for row in rows:
                reasons = []
                validation = row.get("validation") or {}
                if row.get("status") != "compiled":
                    reasons.append("compilation did not succeed")
                if validation.get("passed") is not True:
                    reasons.append("independent validation absent or failed")
                if not _numeric(row.get("compile_ns"), positive=True):
                    reasons.append("missing or invalid compile_ns")
                if not _numeric(row.get("parse_ns")):
                    reasons.append("missing or invalid parse_ns")
                metrics = validation.get("metrics") or {}
                if any(not _numeric(metrics.get(key)) for key in METRICS):
                    reasons.append("missing or invalid common metrics")
                if reasons:
                    failures.append(_brief_failure(row, reasons))
                else:
                    passed.append(row)
            eligible = (len(rows) == expected_trials and len(passed) == expected_trials
                        and not duplicate_trials)
            reasons = []
            if len(rows) != expected_trials:
                reasons.append(f"received {len(rows)} of exactly {expected_trials} expected trials")
            if failures:
                reasons.append(f"{len(failures)} trials failed eligibility checks")
            if duplicate_trials:
                reasons.append("duplicate trial number")
            cell = {
                "case_id": case["id"], "fixture": case["fixture"],
                "family": case["family"], "topology": case["topology"],
                "n_logical": case["n_logical"], "physical_width": case["physical_width"],
                "tool": tool, "eligible": eligible,
                "expected_trials": expected_trials, "received_trials": len(rows),
                "passed_trials": len(passed), "exclusion_reasons": reasons,
                "status_counts": dict(Counter(str(row.get("status")) for row in rows)),
                "evidence_counts": dict(Counter((row.get("validation") or {}).get("evidence_class", "missing") for row in rows)),
                "failed_trials": failures,
                "medians": None, "ranges": None,
            }
            if eligible:
                values = {
                    "compile_ns": [row["compile_ns"] for row in rows],
                    "parse_ns": [row["parse_ns"] for row in rows],
                    # Median of each trial's sum, never sum/subtraction of
                    # separately aggregated timing statistics.
                    "parse_compile_ns": [row["parse_ns"] + row["compile_ns"] for row in rows],
                    **{key: [row["validation"]["metrics"][key] for row in rows] for key in METRICS},
                }
                cell["medians"] = {key: statistics.median(v) for key, v in values.items()}
                cell["ranges"] = {key: {"min": min(v), "q25": _percentile(v, 25),
                                       "q75": _percentile(v, 75), "max": max(v)}
                                  for key, v in values.items()}
            cells.append(cell)
    return cells


def geometric_ratio_summary(observations, *, label, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    """observations contain numerator, denominator and fixture for each cell."""
    positive = []
    zeros = Counter({"both_zero": 0, "numerator_zero": 0, "denominator_zero": 0})
    for item in observations:
        numerator, denominator = item["numerator"], item["denominator"]
        if numerator == 0 and denominator == 0:
            zeros["both_zero"] += 1
        elif numerator == 0:
            zeros["numerator_zero"] += 1
        elif denominator == 0:
            zeros["denominator_zero"] += 1
        else:
            positive.append((item["fixture"], math.log(numerator / denominator)))
    result = {
        "paired_cells": len(observations), "positive_ratio_cells": len(positive),
        "zero_cases": dict(zeros), "geomean": None, "ci95": None,
        "median_ratio": None, "p10_ratio": None, "p90_ratio": None,
        "bootstrap_fixture_clusters": len({fixture for fixture, _ in positive}),
    }
    if not positive:
        return result
    ratios = [math.exp(log_ratio) for _, log_ratio in positive]
    result.update({"geomean": math.exp(statistics.mean(log_ratio for _, log_ratio in positive)),
                   "median_ratio": statistics.median(ratios),
                   "p10_ratio": _percentile(ratios, 10), "p90_ratio": _percentile(ratios, 90)})
    groups = defaultdict(list)
    for fixture, log_ratio in positive:
        groups[fixture].append(log_ratio)
    if len(groups) < 2 or draws < 1:
        return result
    names = sorted(groups)
    totals = np.array([sum(groups[name]) for name in names])
    counts = np.array([len(groups[name]) for name in names])
    label_seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")
    rng = np.random.default_rng(seed ^ label_seed)
    indices = rng.integers(0, len(names), size=(draws, len(names)))
    bootstraps = np.exp(totals[indices].sum(axis=1) / counts[indices].sum(axis=1))
    result["ci95"] = [float(np.percentile(bootstraps, 2.5)), float(np.percentile(bootstraps, 97.5))]
    return result


def _wins(pairs, metric):
    counts = Counter({"win": 0, "tie": 0, "loss": 0})
    for qrust, baseline in pairs:
        q, b = qrust["medians"][metric], baseline["medians"][metric]
        counts["win" if q < b else "loss" if q > b else "tie"] += 1
    return dict(counts)


def family_balanced_ratio_summary(observations, *, label, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    """Equal family weights after averaging cell log ratios within each family.

    The secondary interval keeps the represented families fixed and resamples
    fixture clusters independently within each family. A family with only zero
    ratio cases cannot contribute a finite log-ratio estimate and is explicit.
    """
    grouped = defaultdict(lambda: defaultdict(list))
    all_families = {item["family"] for item in observations}
    for item in observations:
        if item["numerator"] > 0 and item["denominator"] > 0:
            grouped[item["family"]][item["fixture"]].append(math.log(item["numerator"] / item["denominator"]))
    means = {family: statistics.mean(value for values in fixtures.values() for value in values)
             for family, fixtures in grouped.items()}
    result = {
        "geomean": math.exp(statistics.mean(means.values())) if means else None,
        "families": len(means), "family_geomeans": {key: math.exp(value) for key, value in sorted(means.items())},
        "families_without_positive_ratios": sorted(all_families - set(means)), "ci95": None,
        "interval_method": "fixed families, fixture-cluster resampling within each family",
    }
    if not means or draws < 1 or not any(len(fixtures) > 1 for fixtures in grouped.values()):
        return result
    label_seed = int.from_bytes(hashlib.sha256((label + "/family-balanced").encode()).digest()[:8], "big")
    rng = np.random.default_rng(seed ^ label_seed)
    bootstrapped_log_means = []
    for family in sorted(grouped):
        fixtures = grouped[family]
        names = sorted(fixtures)
        totals = np.array([sum(fixtures[name]) for name in names])
        counts = np.array([len(fixtures[name]) for name in names])
        indices = rng.integers(0, len(names), size=(draws, len(names)))
        bootstrapped_log_means.append(totals[indices].sum(axis=1) / counts[indices].sum(axis=1))
    values = np.exp(np.mean(bootstrapped_log_means, axis=0))
    result["ci95"] = [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]
    return result


def summarize_pairs(cells, baseline, *, label="aggregate", draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    index = {(cell["case_id"], cell["tool"]): cell for cell in cells if cell["eligible"]}
    pairs = [(q, index[case, baseline]) for (case, tool), q in index.items()
             if tool == "qrust" and (case, baseline) in index]
    result = {"paired_cells": len(pairs), "fixtures": len({q["fixture"] for q, _ in pairs})}
    for metric in ("compile_ns", "parse_compile_ns", *METRICS):
        speed = metric.endswith("_ns")
        observations = [{"fixture": q["fixture"], "family": q["family"],
                         "numerator": (b if speed else q)["medians"][metric],
                         "denominator": (q if speed else b)["medians"][metric]}
                        for q, b in pairs]
        result[metric] = geometric_ratio_summary(observations, label=f"{label}/{baseline}/{metric}", draws=draws, seed=seed)
        result[metric]["family_balanced"] = family_balanced_ratio_summary(observations, label=f"{label}/{baseline}/{metric}", draws=draws, seed=seed)
        result[metric]["qrust_win_tie_loss"] = _wins(pairs, metric)
        result[metric]["ratio_definition"] = "baseline / Q-Rust; >1 means Q-Rust faster" if speed else "Q-Rust / baseline; <1 means Q-Rust fewer/shallower"
    return result


def _dominates(a, b, metrics):
    return all(a[key] <= b[key] for key in metrics) and any(a[key] < b[key] for key in metrics)


def complete_case_quality(cells):
    grouped = defaultdict(dict)
    for cell in cells:
        if cell["eligible"]:
            grouped[cell["case_id"]][cell["tool"]] = cell
    entries = []
    frontier_counts = {name: Counter({tool: 0 for tool in TOOLS})
                       for name in ("speed_cx", "speed_depth", "speed_cx_depth")}
    for case_id, tools in grouped.items():
        if set(tools) != set(TOOLS):
            continue
        qrust = tools["qrust"]
        quality = {}
        for metric in ("cx_count", "depth", "two_q_depth"):
            values = {tool: cell["medians"][metric] for tool, cell in tools.items()}
            best = min(values.values())
            best_baseline = min(values[tool] for tool in BASELINES)
            q = values["qrust"]
            quality[metric] = {
                "qrust": q, "best_observed_all_four": best,
                "best_observed_baseline": best_baseline,
                "best_all_four_tools": [tool for tool in TOOLS if values[tool] == best],
                "best_baseline_tools": [tool for tool in BASELINES if values[tool] == best_baseline],
                "ratio_to_best_baseline": q / best_baseline if best_baseline > 0 else None,
                "ratio_to_best_all_four": q / best if best > 0 else None,
                "best_zero": best_baseline == 0, "all_zero": all(value == 0 for value in values.values()),
                "within_5_percent": q <= best_baseline * 1.05,
                "within_10_percent": q <= best_baseline * 1.10,
                "at_best": q <= best_baseline,
                "at_all_four_best": q == best,
                "within_5_percent_all_four": q <= best * 1.05,
                "within_10_percent_all_four": q <= best * 1.10,
            }
        frontiers = {}
        for name, metrics in (("speed_cx", ("compile_ns", "cx_count")),
                              ("speed_depth", ("compile_ns", "depth")),
                              ("speed_cx_depth", ("compile_ns", "cx_count", "depth"))):
            frontiers[name] = [tool for tool in TOOLS if not any(
                _dominates(tools[other]["medians"], tools[tool]["medians"], metrics)
                for other in TOOLS if other != tool)]
            frontier_counts[name].update(frontiers[name])
        entries.append({"case_id": case_id, "fixture": qrust["fixture"],
                        "family": qrust["family"], "topology": qrust["topology"],
                        "quality": quality, "pareto_frontiers": frontiers})
    totals = {}
    for metric in ("cx_count", "depth", "two_q_depth"):
        totals[metric] = {key: sum(entry["quality"][metric][key] for entry in entries)
                         for key in ("at_best", "within_5_percent", "within_10_percent", "best_zero", "all_zero")}
        totals[metric]["percentage_of_cells"] = {
            key: 100 * totals[metric][key] / len(entries) if entries else None
            for key in ("at_best", "within_5_percent", "within_10_percent")
        }
    return {"complete_all_four_cells": len(entries), "near_best_counts": totals,
            "pareto_counts": {key: dict(value) for key, value in frontier_counts.items()},
            "cells": entries}


def worst_tails(cells, limit=10):
    index = {(cell["case_id"], cell["tool"]): cell for cell in cells if cell["eligible"]}
    result = {}
    for baseline in BASELINES:
        rows = []
        for (case, tool), q in index.items():
            if tool != "qrust" or (case, baseline) not in index:
                continue
            b = index[case, baseline]
            rows.append({"case_id": case, "fixture": q["fixture"], "family": q["family"], "topology": q["topology"],
                         "qrust": q["medians"], "baseline": b["medians"],
                         "compile_speed_ratio": b["medians"]["compile_ns"] / q["medians"]["compile_ns"]})
        result[baseline] = {"slowest_relative_cases": sorted(rows, key=lambda row: row["compile_speed_ratio"])[:limit]}
        for metric in ("cx_count", "depth", "two_q_depth"):
            def ratio(row):
                q, b = row["qrust"][metric], row["baseline"][metric]
                return q / b if b else (math.inf if q else 1.0)
            tail = []
            for row in sorted(rows, key=lambda row: (ratio(row), row["qrust"][metric] - row["baseline"][metric]), reverse=True)[:limit]:
                r = ratio(row)
                tail.append({**row, "quality_ratio": r if math.isfinite(r) else None,
                             "ratio_kind": "positive numerator / zero baseline" if math.isinf(r) else "both zero" if row["baseline"][metric] == 0 else "finite",
                             "absolute_excess": row["qrust"][metric] - row["baseline"][metric]})
            result[baseline][f"worst_{metric}"] = tail
    return result


def analyze(protocol, timed, validated, *, expected_trials=7, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    cells = build_cells(protocol, timed, validated, expected_trials)
    summary = {
        "schema_version": 1,
        "method": {
            "expected_trials_per_cell_tool": expected_trials,
            "cell_estimator": "median of all prescribed trials; only complete validated cells eligible",
            "aggregate_estimator": "equal-cell-weight geometric mean of positive paired ratios",
            "secondary_family_balanced_estimator": "mean of within-family mean log ratios; equal weight for each family with positive ratios; interval resamples fixture clusters within fixed families",
            "bootstrap": {"unit": "fixture; all its topology cells remain together", "draws": draws,
                          "seed": seed, "interval": "percentile 95%; no interval with fewer than two fixtures"},
            "zero_policy": "zero ratios and zero denominators excluded from log/geomean summaries and counted separately; all pairs retained for win/tie/loss",
            "timing": "warm API compilation; parse+compile is secondary; process startup, imports, serialization and validation excluded",
            "near_best": "primary bands compare Q-Rust against the best of Qiskit/Cirq/TKET cell medians; all-four-complete cells only; at_best means at or better than best baseline; separate all-four ratios/bands also retained",
            "scope": "descriptive fixed-corpus comparison; no causal language-speed attribution or optimality claim",
        },
        "coverage": {
            "planned_cases": len(protocol["cases"]), "planned_fixtures": len({case["fixture"] for case in protocol["cases"]}),
            "expected_total_trials": len(protocol["cases"]) * len(TOOLS) * expected_trials,
            "timed_trials": len(timed), "validation_rows": len(validated),
            "eligible_cells": dict(Counter(cell["tool"] for cell in cells if cell["eligible"])),
            "excluded_cell_count": sum(not cell["eligible"] for cell in cells),
            "validation_evidence_counts": dict(Counter((row.get("validation") or {}).get("evidence_class", "missing") for row in validated)),
        },
        "pairwise": {"aggregate": {tool: summarize_pairs(cells, tool, draws=draws, seed=seed) for tool in BASELINES}},
        "complete_case_quality": complete_case_quality(cells),
        "worst_tails": worst_tails(cells),
        "excluded_cells": [cell for cell in cells if not cell["eligible"]],
    }
    for field, destination in (("topology", "by_topology"), ("family", "by_family")):
        summary["pairwise"][destination] = {
            group: {tool: summarize_pairs([cell for cell in cells if cell[field] == group], tool,
                                           label=f"{field}/{group}", draws=draws, seed=seed) for tool in BASELINES}
            for group in sorted({cell[field] for cell in cells})}
    return summary, cells


def _ratio_text(stat):
    if stat["geomean"] is None:
        return "undefined (no positive pairs)"
    value = f"{stat['geomean']:.3g}"
    if stat["ci95"]:
        value += f" [{stat['ci95'][0]:.3g}, {stat['ci95'][1]:.3g}]"
    return value


def _wlt(stat):
    values = stat["qrust_win_tie_loss"]
    return f"{values['win']}/{values['tie']}/{values['loss']}"


def narrative(summary):
    coverage = summary["coverage"]
    lines = ["# Corrected-fixture compiler comparison", "",
             f"The prescribed corpus contains {coverage['planned_fixtures']} fixtures and {coverage['planned_cases']} fixture–topology cells. "
             f"The study retained {coverage['timed_trials']} timed trials and {coverage['validation_rows']} validation records; "
             f"{coverage['excluded_cell_count']} compiler–cell combinations were ineligible for paired analysis.", "",
             "Only cells with every prescribed trial successful and independently validated contribute to comparisons. "
             "Results below use per-cell medians. Brackets are 95% fixture-cluster bootstrap intervals; "
             "they preserve the relationship between multiple topologies of the same fixture.", "",
             "## Compilation speed", "",
             "Ratios are baseline time divided by Q-Rust time; a value above 1 favors Q-Rust. "
             "These are warm API compilation measurements. Parsing is reported separately; startup, imports, serialization and validation are excluded.", "",
             "| Baseline | Paired cells | Compile speed ratio [95% CI] | Parse + compile ratio [95% CI] | Q-Rust compile win/tie/loss |",
             "|---|---:|---:|---:|---:|"]
    aggregate = summary["pairwise"]["aggregate"]
    for tool in BASELINES:
        stat = aggregate[tool]
        lines.append(f"| {tool} | {stat['paired_cells']} | {_ratio_text(stat['compile_ns'])} | {_ratio_text(stat['parse_compile_ns'])} | {_wlt(stat['compile_ns'])} |")
    lines += ["", "## Output quality", "",
              "Ratios are Q-Rust divided by the baseline; a value below 1 favors Q-Rust. CX count and depth are structural cost measures, not execution fidelity. "
              "Zero/zero and one-sided-zero comparisons remain in win/tie/loss counts but are excluded from geometric means.", "",
              "| Baseline | CX ratio [95% CI] | CX win/tie/loss | Depth ratio [95% CI] | Depth win/tie/loss | CX zero cases (both / Q-Rust / baseline) |",
              "|---|---:|---:|---:|---:|---:|"]
    for tool in BASELINES:
        stat = aggregate[tool]
        zero = stat["cx_count"]["zero_cases"]
        lines.append(f"| {tool} | {_ratio_text(stat['cx_count'])} | {_wlt(stat['cx_count'])} | {_ratio_text(stat['depth'])} | {_wlt(stat['depth'])} | {zero['both_zero']}/{zero['numerator_zero']}/{zero['denominator_zero']} |")
    lines += ["", "## Secondary family-balanced comparison", "",
              "The primary tables weight each eligible cell equally. This secondary summary first averages log ratios within each family, "
              "then weights the represented families equally. Its intervals resample fixture clusters within each fixed family. "
              "Families with no positive ratios for a metric are excluded from that metric and recorded in the JSON.", "",
              "| Baseline | Compile speed ratio [95% CI] | CX ratio [95% CI] | Depth ratio [95% CI] | Families (speed / CX / depth) |",
              "|---|---:|---:|---:|---:|"]
    for tool in BASELINES:
        stats = [aggregate[tool][metric]["family_balanced"] for metric in ("compile_ns", "cx_count", "depth")]
        lines.append(f"| {tool} | {' | '.join(_ratio_text(stat) for stat in stats)} | {' / '.join(str(stat['families']) for stat in stats)} |")
    complete = summary["complete_case_quality"]
    count = complete["complete_all_four_cells"]
    lines += ["", f"Among the {count} cells eligible for all four tools:", "",
              "| Metric | Q-Rust at or better than best baseline | Within 5% of best baseline | Within 10% of best baseline | Best baseline is zero |",
              "|---|---:|---:|---:|---:|"]
    for metric in ("cx_count", "depth", "two_q_depth"):
        stat = complete["near_best_counts"][metric]
        values = [f"{stat[key]}/{count}" + (f" ({stat['percentage_of_cells'][key]:.1f}%)" if count else "")
                  for key in ("at_best", "within_5_percent", "within_10_percent")]
        lines.append(f"| {metric} | {' | '.join(values)} | {stat['best_zero']}/{count} |")
    lines += ["", "The observed best baseline is the minimum median among Qiskit, Cirq and TKET, not a proven optimum. "
              "When the best baseline count is zero, only a zero-count Q-Rust result is within 5% or 10%. "
              "The JSON additionally records ratios and bands against the best of all four tools.", "",
              "## Topology and family breakdowns", "",
              "| Group | Baseline | Cells | Compile speed ratio [95% CI] | CX ratio [95% CI] | Depth ratio [95% CI] |",
              "|---|---|---:|---:|---:|---:|"]
    for grouping in ("by_topology", "by_family"):
        for group, tools in summary["pairwise"][grouping].items():
            for tool in BASELINES:
                stat = tools[tool]
                lines.append(f"| {group} | {tool} | {stat['paired_cells']} | {_ratio_text(stat['compile_ns'])} | {_ratio_text(stat['cx_count'])} | {_ratio_text(stat['depth'])} |")
    lines += ["", "## Tradeoffs and failure tails", "",
              "A pipeline is on the Pareto frontier when no other measured pipeline is at least as good on every listed metric and strictly better on one. "
              "Frontiers use per-cell medians and count tied pipelines separately.", "",
              "| Objectives | Q-Rust | Qiskit | Cirq | TKET |", "|---|---:|---:|---:|---:|"]
    for name, counts in complete["pareto_counts"].items():
        lines.append(f"| {name} | " + " | ".join(f"{counts[tool]}/{count}" for tool in TOOLS) + " |")
    lines += ["", "| Baseline | Worst CX excess case | Q-Rust CX / baseline CX | Slowest relative case | Baseline / Q-Rust time |",
              "|---|---|---:|---|---:|"]
    for tool in BASELINES:
        tails = summary["worst_tails"][tool]
        if not tails["worst_cx_count"]:
            continue
        quality = tails["worst_cx_count"][0]
        speed = tails["slowest_relative_cases"][0]
        lines.append(f"| {tool} | {quality['case_id']} | {quality['qrust']['cx_count']} / {quality['baseline']['cx_count']} | {speed['case_id']} | {speed['compile_speed_ratio']:.3g} |")
    lines += ["", "## What the evidence supports", "",
              "The measurements support comparisons of the specified compiler recipes on these corrected, static circuits and target graphs, under the recorded single-machine conditions. "
              "Any speed claim should quote the paired ratio and name this scope. Any compactness claim should name CX count or depth and include the loss cases. "
              "The per-family and per-topology tables can identify useful subsets, but a subset selected after viewing these results is exploratory and needs a separate confirmation corpus.", "",
              "The experiment does not isolate the effect of Rust: the tools differ in algorithms, search effort, representations and API overhead, and competitors have native implementations. "
              "Seven trials equalize repetitions, not optimization work or runtime budgets. These recipes are documented choices, not a search over every available competitor pipeline. "
              "The fixed corpus and limited widths do not establish general scalability, hardware fidelity, or closeness to a mathematical optimum. "
              "Numerical validation is evidence for the measured outputs; sampled checks are not proofs of equivalence.", "",
              "## Eligibility and reproducibility", "",
              f"Expected trial records: {coverage['expected_total_trials']}. Eligible cells by tool: "
              + ", ".join(f"{tool} {coverage['eligible_cells'].get(tool, 0)}" for tool in TOOLS) + ".", "",
              "Validation evidence classes: `" + json.dumps(coverage["validation_evidence_counts"], sort_keys=True) + "`.", "",
              "Every excluded cell and its individual trial failures are listed in `summary.json`; missing trials are explicit. "
              "`cell_summary.csv` contains all compiler–cell combinations, including exclusions. "
              "The JSON also retains zero cases, subgroup intervals, Pareto membership and the ten worst tails per baseline. "
              "Aggregate geometric means give equal weight to eligible cells, so fixtures with more target topologies contribute more cells. "
              "Bootstrap intervals resample fixture clusters and condition on positive-ratio cells; they do not correct corpus-selection bias or multiple comparisons.", ""]
    if summary["excluded_cells"]:
        lines += ["| Excluded cell | Tool | Passed / received / expected | Reasons |", "|---|---|---:|---|"]
        for cell in summary["excluded_cells"]:
            lines.append(f"| {cell['case_id']} | {cell['tool']} | {cell['passed_trials']} / {cell['received_trials']} / {cell['expected_trials']} | {'; '.join(cell['exclusion_reasons'])} |")
        lines.append("")
    return "\n".join(lines)


def write_csv(path, cells):
    timings = ("compile_ns", "parse_ns", "parse_compile_ns")
    fields = ["case_id", "fixture", "family", "topology", "n_logical", "physical_width", "tool", "eligible",
              "expected_trials", "received_trials", "passed_trials", "exclusion_reasons", "status_counts", "evidence_counts"]
    fields += [f"median_{key}" for key in (*timings, *METRICS)]
    fields += [f"{key}_{quantile}" for key in timings for quantile in ("min", "q25", "q75", "max")]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cell in cells:
            row = {key: cell.get(key) for key in fields if key in cell}
            for key in ("exclusion_reasons", "status_counts", "evidence_counts"):
                row[key] = json.dumps(row[key], sort_keys=True)
            if cell["eligible"]:
                row.update({f"median_{key}": value for key, value in cell["medians"].items()})
                row.update({f"{key}_{quantile}": cell["ranges"][key][quantile] for key in timings for quantile in ("min", "q25", "q75", "max")})
            writer.writerow(row)


def write_plots(output_dir, summary, cells):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 130, "savefig.dpi": 180})
    colors = {"qiskit": "#2463A6", "cirq": "#C47620", "tket": "#36846B"}
    index = {(cell["case_id"], cell["tool"]): cell for cell in cells if cell["eligible"]}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1))
    for ax, metric, title in zip(axes, ("compile_ns", "cx_count", "depth"),
                                 ("Warm compile speed", "CX count ratio", "Circuit depth ratio")):
        for position, tool in enumerate(BASELINES):
            stat = summary["pairwise"]["aggregate"][tool][metric]
            if stat["geomean"] is None:
                continue
            value = stat["geomean"]
            interval = stat["ci95"]
            ax.plot(position, value, "o", color=colors[tool], markersize=7)
            if interval:
                ax.vlines(position, interval[0], interval[1], color=colors[tool], linewidth=2)
                ax.hlines(interval, position - 0.07, position + 0.07, color=colors[tool])
            ax.annotate(f"{value:.3g}", (position, value), xytext=(7, 0), textcoords="offset points", va="center")
        ax.axhline(1, color="#666666", linestyle="--", linewidth=1)
        ax.set_yscale("log")
        ax.set_xticks(range(3), [TOOL_LABELS[tool] for tool in BASELINES])
        ax.set_xlim(-0.45, 2.65)
        ax.set_title(title)
        ax.set_ylabel("Baseline / Q-Rust (>1 favors Q-Rust)" if metric == "compile_ns" else "Q-Rust / baseline (<1 favors Q-Rust)")
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Paired geometric means with 95% fixture-cluster intervals", y=1.03)
    fig.tight_layout()
    fig.savefig(output_dir / "aggregate_ratios.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, metric, title in zip(axes, ("cx_count", "depth"), ("CX count", "Circuit depth")):
        omitted = 0
        for tool in BASELINES:
            points = []
            for (case, name), q in index.items():
                if name != "qrust" or (case, tool) not in index:
                    continue
                b = index[case, tool]
                if not q["medians"][metric] or not b["medians"][metric]:
                    omitted += 1
                    continue
                points.append((b["medians"]["compile_ns"] / q["medians"]["compile_ns"],
                               q["medians"][metric] / b["medians"][metric]))
            if points:
                x, y = zip(*points)
                ax.scatter(x, y, s=16, alpha=.55, color=colors[tool], label=TOOL_LABELS[tool])
        ax.axhline(1, color="#777777", linewidth=1, linestyle="--")
        ax.axvline(1, color="#777777", linewidth=1, linestyle="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Compile speed ratio: baseline / Q-Rust")
        ax.set_ylabel(f"{title} ratio: Q-Rust / baseline")
        ax.set_title(f"{title}; {omitted} zero comparisons omitted")
        ax.grid(alpha=.15)
    axes[0].legend(fontsize=9)
    fig.suptitle("Each point is a paired fixture–topology cell; lower right favors Q-Rust", y=1.02)
    fig.tight_layout()
    fig.savefig(output_dir / "speed_quality_tradeoffs.png", bbox_inches="tight")
    plt.close(fig)

    groups = summary["pairwise"]["by_topology"]
    matrix = np.array([[groups[group][tool]["compile_ns"]["geomean"] or np.nan for tool in BASELINES] for group in groups])
    fig, ax = plt.subplots(figsize=(6.8, max(3.0, len(groups) * .75)))
    finite = matrix[np.isfinite(matrix)]
    bound = max(1.0, float(np.max(np.abs(np.log2(finite))))) if len(finite) else 1.0
    graphic = ax.imshow(np.log2(matrix), cmap="RdBu", vmin=-bound, vmax=bound, aspect="auto")
    ax.set_xticks(range(3), [TOOL_LABELS[tool] for tool in BASELINES])
    ax.set_yticks(range(len(groups)), list(groups))
    for row in range(len(groups)):
        for column in range(3):
            value = matrix[row, column]
            ax.text(column, row, "—" if not math.isfinite(value) else f"{value:.3g}×", ha="center", va="center",
                    color="white" if math.isfinite(value) and abs(math.log2(value)) > bound * .55 else "black")
    ax.set_title("Warm compile speed by target topology\nBaseline / Q-Rust; >1 favors Q-Rust")
    fig.colorbar(graphic, ax=ax, label="log₂ speed ratio", shrink=.8)
    fig.tight_layout()
    fig.savefig(output_dir / "topology_speed.png", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--expected-trials", type=int)
    parser.add_argument("--bootstrap-draws", type=int, default=BOOTSTRAP_DRAWS)
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    protocol = json.loads((args.output_dir / "protocol.json").read_text())
    expected = protocol.get("trials", 7) if args.expected_trials is None else args.expected_trials
    if not isinstance(expected, int) or expected < 1 or args.bootstrap_draws < 1:
        parser.error("trial count and bootstrap draw count must be positive integers")
    summary, cells = analyze(protocol, read_jsonl(args.output_dir / "timed.jsonl"),
                             read_jsonl(args.output_dir / "validated.jsonl"),
                             expected_trials=expected, draws=args.bootstrap_draws, seed=args.seed)
    summary["provenance"] = {
        name + "_sha256": hashlib.sha256((args.output_dir / name).read_bytes()).hexdigest()
        for name in ("protocol.json", "timed.jsonl", "validated.jsonl")
    }
    summary["provenance"]["analysis_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv(args.output_dir / "cell_summary.csv", cells)
    (args.output_dir / "results.md").write_text(narrative(summary))
    if not args.no_plots:
        write_plots(args.output_dir, summary, cells)
    print(json.dumps(summary["coverage"], indent=2))


if __name__ == "__main__":
    main()
