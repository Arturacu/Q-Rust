#!/usr/bin/env python3
"""Render compact public tables and a figure from frozen external-study results.

This presentation step does not select cases or recompute the frozen statistics.
The evidence archive retains the complete original analysis and raw records.
"""
import argparse
import copy
import json
from pathlib import Path
import shutil


def render(study, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    full = json.loads((study / "analysis/external_summary.json").read_text())
    compact = copy.deepcopy(full)
    def trim_family_members(value):
        if isinstance(value, dict):
            value.pop("family_geomeans", None)
            for child in value.values():
                trim_family_members(child)
        elif isinstance(value, list):
            for child in value:
                trim_family_members(child)
    trim_family_members(compact)
    for session in compact["sessions"].values():
        for profile in ("primary", "light"):
            value = session[profile]
            value["pairwise"].pop("by_family")
            value["complete_case_quality"].pop("cells")
            for cell in value["excluded_cells"]:
                trials = cell.pop("failed_trials")
                cell["distinct_errors"] = sorted({t.get("error", "") for t in trials})
            for tails in value["worst_tails"].values():
                for metric in tails:
                    tails[metric] = tails[metric][:3]
    compact["presentation_note"] = (
        "Compact copy of the frozen analysis: individual family ratios, complete-case "
        "quality records and individual failed trials are omitted; worst-tail lists "
        "show their first three entries. Public CSV copies use LF line endings. "
        "All byte-exact original records remain in evidence.tar.gz. "
        "All 82 paired cells are conditional on Q-Rust accepting the input; 18 of "
        "100 selected cells failed. Session/profile repetitions are not independent circuits.")
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(json.dumps(compact, indent=2) + "\n")
    for source in sorted((study / "analysis").glob("*_cells.csv")):
        (output / source.name).write_text(source.read_text())
    shutil.copyfile(study / "validation_summary.json", output / "validation_summary.json")

    labels = ["Qiskit O3", "Cirq · max 2", "Cirq · max 1", "TKET · full", "TKET · light"]
    specs = [("primary", "qiskit"), ("primary", "cirq"), ("light", "cirq"),
             ("primary", "tket"), ("light", "tket")]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), layout="constrained")
    ax = axes[0, 0]
    coverage = [82, 100, 97, 96, 100, 100]
    names = ["Q-Rust O3", *labels]
    y = np.arange(len(names))
    ax.barh(y, coverage, color=["#137d80", *["#70889d"] * 5], height=.65)
    ax.barh(y, 100 - np.array(coverage), left=coverage, color="#e5e8ec", height=.65)
    for i, count in enumerate(coverage):
        ax.text(count - 2, i, f"{count}/100", va="center", ha="right", color="white", weight="bold")
    ax.set(yticks=y, yticklabels=names, xlim=(0, 104), xlabel="Cases with all seven trials validated",
           title="Coverage · same in both sessions")
    ax.invert_yaxis()

    for ax, metric, title, xlabel in [
        (axes[0, 1], "compile_ns", "Warm compilation · 82 paired cases", "Baseline / Q-Rust time  (higher favors Q-Rust)"),
        (axes[1, 0], "cx_count", "CX count · 80 positive-ratio cases", "Q-Rust / baseline CX  (lower favors Q-Rust)"),
        (axes[1, 1], "depth", "Total depth · 82 paired cases", "Q-Rust / baseline depth  (lower favors Q-Rust)"),
    ]:
        for session, color, offset in [("1", "#137d80", -.10), ("2", "#c87529", .10)]:
            if metric != "compile_ns" and session == "2":
                continue  # Every successful same-seed output matched exactly.
            stats = [full["sessions"][session][profile]["pairwise"]["aggregate"][tool][metric]
                     for profile, tool in specs]
            values = np.array([s["geomean"] for s in stats])
            errors = np.array([[s["geomean"] - s["ci95"][0], s["ci95"][1] - s["geomean"]]
                               for s in stats]).T
            positions = np.arange(len(labels)) + (offset if metric == "compile_ns" else 0)
            ax.errorbar(values, positions, xerr=errors, fmt="o", color=color, capsize=3,
                        label=f"Session {session}", markersize=5)
        ax.axvline(1, color="#999999", linestyle="--", linewidth=1)
        ax.set(yticks=np.arange(len(labels)), yticklabels=labels, title=title, xlabel=xlabel)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=.2)
        if metric == "compile_ns":
            ax.set_xscale("log")
            ax.legend(loc="upper right", frameon=False)
        else:
            ax.set_xlim(.45, 1.23)
    fig.suptitle("External compiler comparison · QASMBench + MQT Bench", fontsize=16, weight="bold")
    fig.supxlabel("Geometric means of per-case median ratios; bars: 95% fixture-cluster resampling intervals.\n"
                  "Quality is identical across sessions. Two zero-CX cases remain in counts, outside CX geometric means.\n"
                  "41 of 50 fixtures contribute to paired results; input/adapter failures are not repaired or replaced.",
                  fontsize=9)
    fig.savefig(output / "comparison.png", dpi=180)
    fig.savefig(output / "comparison.svg")
    svg = output / "comparison.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.study, args.output)
