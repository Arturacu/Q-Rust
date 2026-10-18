#!/usr/bin/env python3
"""Present candidate-retention follow-up results without modifying either preserved study.

The renderer copies existing statistics; it does not reselect cases or refit
intervals. The before/after analyzer defines the common comparison cohort.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path

SPECS = (
    ("Qiskit O3", "primary", "qiskit"),
    ("Cirq (max 2)", "primary", "cirq"),
    ("Cirq (max 1)", "light", "cirq"),
    ("TKET (full)", "primary", "tket"),
    ("TKET (light)", "light", "tket"),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trim_family_members(value):
    if isinstance(value, dict):
        value.pop("family_geomeans", None)
        for child in value.values():
            trim_family_members(child)
    elif isinstance(value, list):
        for child in value:
            trim_family_members(child)


def compact_summary(full, comparison, protocol, validation):
    result = copy.deepcopy(full)
    trim_family_members(result)
    for session in result["sessions"].values():
        for profile in ("primary", "light"):
            value = session[profile]
            value["pairwise"].pop("by_family", None)
            value["complete_case_quality"].pop("cells", None)
            for cell in value["excluded_cells"]:
                trials = cell.pop("failed_trials", [])
                cell["distinct_errors"] = sorted({t.get("error") or "" for t in trials})
            for tails in value["worst_tails"].values():
                for metric in tails:
                    tails[metric] = tails[metric][:3]
    changes = result["session_agreement"].pop("changed_trials", [])
    result["session_agreement"]["changed_trial_count"] = len(changes)
    result["session_agreement"]["first_changed_trials"] = changes[:10]
    followup = copy.deepcopy(comparison)
    trim_family_members(followup)
    for session in followup["sessions"].values():
        changes = session["trial_changes"].pop("changed_trials", [])
        session["trial_changes"]["changed_trial_count"] = len(changes)
        session["trial_changes"]["first_changed_trials"] = changes[:10]
    # Keep the public view compact; the archive retains every paired cell.
    for session in followup["sessions"].values():
        session.pop("recipes_on_original_common", None)  # same fixed cohort as recipes
        for recipe, changes in session.get("case_changes", {}).items():
            changes.pop("cells", None)
            if recipe != "qrust":
                changes.pop("compile_adverse_cases", None)
        for profile in session.get("profiles", {}).values():
            profile.pop("followup_on_original_common", None)  # same as followup_on_all_common
            profile.pop("followup_on_newly_common", None)  # empty in this unchanged cohort
    result["before_after"] = followup
    result["validation_receipt"] = validation
    result["study_metadata"] = {key: protocol[key] for key in (
        "compiler_commit", "compiler_provenance", "qrust_binary_sha256", "corpus_manifest_sha256",
        "versions", "python", "platform", "machine", "logical_cpu_count", "machine_note",
        "thread_environment", "recipes", "recipe_profiles", "trials", "sessions", "master_seed",
        "timeout_seconds", "timing_boundary", "source_sha256",
    )}
    result["presentation_note"] = (
        "Compact presentation of preserved analyses. Individual family ratios, complete-case "
        "quality records and individual failed-trial records are omitted; worst-tail lists retain "
        "three entries and changed-trial lists retain ten. Counts and aggregate statistics remain. "
        "The evidence archive retains complete records. This view also omits duplicated fixed-cohort summaries, "
        "full case-change cells and unchanged-SDK timing-tail lists. Public CSV copies use LF line endings. "
        "Primary/light views reuse the same Q-Rust and Qiskit observations. Two sessions are "
        "timing repetitions, not additional independent circuits. The figure uses the follow-up "
        "primary all-four common cohort in both sessions; any profile-specific missing pairs "
        "remain explicit in the comparison statistics. Original results are separate artifacts."
    )
    return result


def span(values):
    values = sorted(set(values))
    return str(values[0]) if len(values) == 1 else f"{values[0]}–{values[-1]}"


def figure_notes(full, comparison, protocol, validation):
    cohort = comparison["cohorts"]["followup_common_case_ids"]
    fixtures = {case["fixture"] for case in protocol["cases"] if case["id"] in set(cohort)}
    planned = len(protocol["cases"])
    complete = []
    for session in ("1", "2"):
        primary = full["sessions"][session]["primary"]["coverage"]["eligible_cells"]
        light = full["sessions"][session]["light"]["coverage"]["eligible_cells"]
        complete.append(sum(primary.values()) + light.get("cirq", 0) + light.get("tket", 0))
    evidence = validation["unique_evidence_classes"]
    counts = Counter()
    for session in validation["sessions"]:
        counts.update(session["counts"])
    passed = sum(value for name, value in counts.items() if name.endswith("_passed"))
    failed = sum(value for name, value in counts.items() if name.endswith("_failed"))
    other = sum(counts.values()) - passed - failed
    positives = {metric: [comparison["sessions"][s]["profiles"][profile]["followup_on_all_common"]
                         ["pairwise"][tool][metric]["positive_ratio_cells"]
                         for s in ("1", "2") for _, profile, tool in SPECS]
                 for metric in ("compile_ns", "cx_count", "depth")}
    return [
        f"Cohort: {len(cohort)} primary-common cases / {len(fixtures)} fixtures. Complete recipe–case cells: "
        f"S1 {complete[0]}/{planned * len(protocol['tools'])}; S2 {complete[1]}/{planned * len(protocol['tools'])}.",
        "Geometric means of case medians. Filled circles + bars: session 1 and 95% fixture-cluster intervals; "
        "open diamonds: session 2 point estimates.",
        f"Positive ratios: compile {span(positives['compile_ns'])}; CX {span(positives['cx_count'])}; "
        f"depth {span(positives['depth'])}. Zero-cost pairs remain in counts, outside geometric means.",
        f"Numerical checks: {passed}/{passed + failed} emitted trial outputs passed; {other} non-emitting outcomes. "
        f"Unique artifacts: {validation['unique_artifacts']} = {evidence.get('full_isometry', 0)} full + "
        f"{evidence.get('sampled_isometry', 0)} sampled isometries.",
        "Intervals describe this fixed corpus. Sampled checks are not proofs; CX count and depth are not execution fidelity.",
    ]


def render_focal(comparison, output):
    """Keep the twelve prespecified focal cells explicit; this view shows the line."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    cases = (
        ("Trotter-4", "qasmbench__basis_trotter_n4__linear"),
        ("QFT-4", "qasmbench__qft_n4__linear"),
        ("PEA-5", "qasmbench__pea_n5__linear"),
        ("AE-10", "mqtbench__ae_n10__linear"),
        ("QPE-10 exact", "mqtbench__qpeexact_n10__linear"),
        ("QPE-10 inexact", "mqtbench__qpeinexact_n10__linear"),
    )
    metrics = (
        ("compile_ns", "Warm compilation (ms)", 1e6),
        ("cx_count", "CX count", 1),
        ("depth", "Total depth", 1),
        ("two_q_depth", "CX depth", 1),
    )
    series = (
        ("Q-Rust before", "qrust", "before", "#6c737b", -.16),
        ("Q-Rust after", "qrust", "after", "#137d80", 0),
        ("Qiskit current", "qiskit", "after", "#9e4f8d", .16),
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.2))
    fig.subplots_adjust(left=.15, right=.97, top=.83, bottom=.20, wspace=.37, hspace=.45)
    for axis, (metric, title, scale) in zip(axes.flat, metrics):
        for row, (_, case_id) in enumerate(cases):
            for label, recipe, when, color, offset in series:
                for session, jitter in (("1", -.025), ("2", .025)):
                    value = comparison["sessions"][session]["focal_cases"][case_id][recipe]["metrics"][metric][when] / scale
                    axis.plot(value, row + offset + jitter, color=color,
                              marker="o" if session == "1" else "s", linestyle="none",
                              markerfacecolor=color if session == "1" else "white", markersize=5)
        axis.set(title=title, yticks=range(len(cases)), yticklabels=[label for label, _ in cases],
                 ylim=(len(cases)-.5, -.5), xlim=(0, None))
        axis.grid(axis="x", alpha=.17)
    handles = [Line2D([], [], color=color, marker="o", linestyle="none", label=label)
               for label, _, _, color, _ in series]
    handles += [Line2D([], [], color="#333333", marker="s", markerfacecolor="white",
                       linestyle="none", label="Session 2 (session 1: filled circles)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .10), ncol=2, frameon=False)
    fig.suptitle("Regression recovery and AE/QPE retention on linear graphs", fontsize=17, weight="bold", y=.965)
    fig.text(.5, .91, "Preserved before / after Q-Rust results and the current Qiskit O3 comparison",
             ha="center", fontsize=11)
    fig.text(.035, .025,
             "Each point is one case's median of seven trials; sessions are separate. Lower values are better.\n"
             "Times are from separate historical/current runs, with unchanged baselines as drift controls. "
             "CX and depth are structural costs.\n"
             "This focal view is not the full corpus; comparison.png reports all 100 cases and all five baseline recipes.",
             fontsize=8.4, va="bottom", linespacing=1.5)
    fig.savefig(output / "focal_linear.png", dpi=180)
    fig.savefig(output / "focal_linear.svg", metadata={"Date": None})
    plt.close(fig)
    svg = output / "focal_linear.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")


def render(study, comparison_path, output):
    # Import the numerical/graphics runtime only when explicitly rendering,
    # after every timed session has finished.
    if (study / "active.lock").exists():
        raise ValueError("Cannot render while the study is active")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    full_path = study / "analysis/external_summary.json"
    full = json.loads(full_path.read_text())
    comparison = json.loads(comparison_path.read_text())
    if comparison.get("analysis") != "routing_retention_followup":
        raise ValueError("Expected the candidate-retention follow-up analysis")
    protocol = json.loads((study / "protocol.json").read_text())
    validation = json.loads((study / "validation_summary.json").read_text())
    checksum = digest(study / "protocol.json")
    if not (checksum == full["protocol_sha256"] == comparison["followup_protocol_sha256"]
            == validation["protocol_sha256"] == (study / "protocol.sha256").read_text().strip()):
        raise ValueError("Presentation inputs refer to different frozen protocols")
    if not comparison["cohorts"]["followup_common_case_ids"]:
        raise ValueError("No common follow-up cohort to plot")
    compact = compact_summary(full, comparison, protocol, validation)
    compact["presentation_source_sha256"] = digest(Path(__file__))
    compact["presentation_input_sha256"] = {
        "external_summary.json": digest(full_path),
        "retention_comparison.json": digest(comparison_path),
        "validation_summary.json": digest(study / "validation_summary.json"),
    }
    # An exclusive destination prevents accidentally replacing original figures
    # or a previous presentation. Render revisions into a new directory.
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(json.dumps(compact, indent=2, allow_nan=False) + "\n")
    for session in (1, 2):
        for profile in ("primary", "light"):
            name = f"session_{session}_{profile}_cells.csv"
            (output / name).write_text((study / "analysis" / name).read_text(), newline="\n")
    (output / "validation_summary.json").write_bytes((study / "validation_summary.json").read_bytes())

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "qrust-external-retention"})
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.8))
    fig.subplots_adjust(left=.12, right=.98, top=.84, bottom=.31, wspace=.40)
    labels = [label for label, _, _ in SPECS]
    colors = {"compile_ns": "#137d80", "cx_count": "#196caa", "depth": "#9e4f8d"}
    for axis, metrics, title, xlabel in (
        (axes[0], ("compile_ns",), "Warm compilation", "Baseline / Q-Rust time  (higher favors Q-Rust)"),
        (axes[1], ("cx_count", "depth"), "Output CX count and total depth", "Q-Rust / baseline cost  (lower favors Q-Rust)"),
    ):
        for metric in metrics:
            metric_offset = 0 if len(metrics) == 1 else (-.14 if metric == "cx_count" else .14)
            for session, offset in (("1", -.045), ("2", .045)):
                for row, (_, profile, tool) in enumerate(SPECS):
                    stats = comparison["sessions"][session]["profiles"][profile]["followup_on_all_common"]["pairwise"][tool][metric]
                    value, ci = stats["geomean"], stats["ci95"]
                    if value is None:
                        continue
                    y = row + metric_offset + offset
                    if session == "1" and ci is not None:
                        axis.hlines(y, ci[0], ci[1], color=colors[metric], linewidth=1.7, alpha=.8)
                        axis.vlines(ci, y - .035, y + .035, color=colors[metric], linewidth=1)
                    axis.plot(value, y, marker="o" if session == "1" else "D", linestyle="none",
                              markersize=5, color=colors[metric],
                              markerfacecolor=colors[metric] if session == "1" else "white", markeredgewidth=1.1)
        axis.axvline(1, color="#888888", linestyle="--", linewidth=1)
        axis.set(yticks=range(len(labels)), yticklabels=labels, title=title, xlabel=xlabel,
                 ylim=(len(labels) - .5, -.5), xscale="log")
        axis.grid(axis="x", which="both", alpha=.17)
    axes[0].legend(handles=[
        Line2D([], [], color=colors["compile_ns"], marker="o", label="Session 1 + 95% interval"),
        Line2D([], [], color=colors["compile_ns"], marker="D", markerfacecolor="white", linestyle="none", label="Session 2"),
    ], loc="upper left", bbox_to_anchor=(0, -.17), frameon=False, fontsize=9)
    axes[1].legend(handles=[
        Line2D([], [], color=colors["cx_count"], marker="o", label="CX count"),
        Line2D([], [], color=colors["depth"], marker="o", label="Total depth"),
    ], loc="upper left", bbox_to_anchor=(0, -.17), frameon=False, ncol=2, fontsize=9)
    count = len(comparison["cohorts"]["followup_common_case_ids"])
    fig.suptitle("External compiler candidate-retention follow-up", fontsize=17, weight="bold", y=.965)
    fig.text(.5, .91, f"QASMBench + MQT Bench · {count} common cases · fixed recipes on one host", ha="center", fontsize=11)
    notes = figure_notes(full, comparison, protocol, validation)
    fig.text(.035, .025, "\n".join(notes), va="bottom", fontsize=8.2, linespacing=1.6)
    fig.savefig(output / "comparison.png", dpi=180)
    fig.savefig(output / "comparison.svg", metadata={"Date": None})
    plt.close(fig)
    svg = output / "comparison.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    (output / "FIGURE_NOTES.md").write_text("# External candidate-retention follow-up figure\n\n" + "\n\n".join(notes) + "\n")
    render_focal(comparison, output)
    print(json.dumps({"output": str(output), "common_cases": count,
                      "compact_summary_bytes": (output / "summary.json").stat().st_size}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("comparison", type=Path, help="retention_comparison.json from the separate before/after analyzer")
    parser.add_argument("output", type=Path, help="new destination; existing directories are rejected")
    args = parser.parse_args()
    render(args.study, args.comparison, args.output)
