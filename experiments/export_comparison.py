#!/usr/bin/env python3
"""Cross-tool comparison deliverables: Q-Rust vs Qiskit vs Cirq vs tket.

Merges Q-Rust's O3 results (results/raw_benchmark_data.json) with the baseline
tool results (results/baselines.json) and emits:

  results/comparison_tables.tex   head-to-head LaTeX tables (CX, compile time)
  results/fig_cmp_cx_all2all.png   CX count by tool, ideal connectivity
  results/fig_cmp_cx_hardware.png  CX count by tool, IBM Nairobi (with caveat)
  results/fig_cmp_speed.png        compile-time comparison (log scale)

Honest methodology caveats baked into the captions:
  * Comparison is at each tool's best effort (qrust O3, qiskit optimization_level=3,
    tket FullPeepholeOptimise, cirq RouteCQC).
  * CX counts are CX-equivalent: every tool's SWAPs are decomposed to 3 CX.
  * On hardware topologies, cirq's RouteCQC performs no smart initial placement,
    so its routed CX is inflated relative to the Sabre-based placers in
    qrust/qiskit/tket — the all-to-all table is the cleanest apples-to-apples
    decomposition/optimization comparison.
"""
from __future__ import annotations

import json
import os
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

TOOLS = ["qrust", "qiskit", "tket", "cirq"]
OTHERS = ["qiskit", "tket", "cirq"]
TOOL_LABEL = {"qrust": "Q-Rust", "qiskit": "Qiskit", "tket": "tket", "cirq": "Cirq"}
TOOL_COLOR = {"qrust": "#2563eb", "qiskit": "#dc2626", "tket": "#16a34a", "cirq": "#9333ea"}
FAMILY_ORDER = ["ghz", "qft", "adder", "clifford", "qaoa", "bv", "ising", "wstate", "graphstate"]
FAMILY_LABEL = {"ghz": "GHZ", "qft": "QFT", "adder": "Adder", "clifford": "Clifford",
                "qaoa": "QAOA", "bv": "BV", "ising": "Ising", "wstate": "W-state",
                "graphstate": "GraphState"}
TOPO_ORDER = ["all_to_all", "linear", "ring", "star", "grid_4x4", "ibm_nairobi", "ibm_quito"]
TOPO_LABEL = {"all_to_all": "All2All", "linear": "Linear", "ring": "Ring", "star": "Star",
              "grid_4x4": "Grid4x4", "ibm_nairobi": "Nairobi", "ibm_quito": "Quito"}


def load():
    # Prefer the expanded O3 comparison sweep if present.
    qfile = "cmp_raw_benchmark_data.json"
    if not os.path.exists(os.path.join(RESULTS, qfile)):
        qfile = "raw_benchmark_data.json"
    qrust = [r for r in json.load(open(os.path.join(RESULTS, qfile)))
             if r.get("status") == "ok" and r["opt_level"] == 3 and not r["ablation"]]
    base = [r for r in json.load(open(os.path.join(RESULTS, "baselines.json")))
            if r.get("status") == "ok"]
    # Index: (fixture, topology, tool) -> {cx, depth, compile_ms}
    data = {}
    for r in qrust:
        # Use transpile-only time (exclude QASM parsing) to match the baselines,
        # which parse before starting their timer.
        t = r["timing_us"]
        transpile_us = t["total"] - t["parse"]
        data[(r["fixture"], r["topology"], "qrust")] = {
            "cx": r["post"]["cx"], "depth": r["post"]["depth"],
            "compile_ms": transpile_us / 1000.0,
            "cx_min": r["post"].get("cx_min"), "cx_max": r["post"].get("cx_max"),
            "n": r["num_qubits_logical"], "family": r["family"],
        }
    for r in base:
        data[(r["fixture"], r["topology"], r["tool"])] = {
            "cx": r["cx_count"], "depth": r["depth"],
            "compile_ms": r["compile_ms_median"],
            "n": r["num_qubits"], "family": r["family"],
        }
    return data


def fixtures_on(data, topo):
    seen = {}
    for (fx, tp, tool), v in data.items():
        if tp == topo:
            seen.setdefault(fx, (v["family"], v["n"]))
    return sorted(seen.items(), key=lambda kv: (FAMILY_ORDER.index(kv[1][0]), kv[1][1]))


# --------------------------------------------------------------------------- #
# Tables                                                                       #
# --------------------------------------------------------------------------- #
def tabular(caption, label, header, rows, colspec=None):
    colspec = colspec or ("l" + "r" * (len(header) - 1))
    out = ["\\begin{table}[htbp]", "  \\centering", f"  \\caption{{{caption}}}",
           f"  \\label{{{label}}}", f"  \\begin{{tabular}}{{{colspec}}}", "    \\hline",
           "    " + " & ".join(header) + " \\\\", "    \\hline"]
    out += ["    " + " & ".join(str(c) for c in r) + " \\\\" for r in rows]
    out += ["    \\hline", "  \\end{tabular}", "\\end{table}", ""]
    return "\n".join(out)


def _bold_min(vals):
    """Format a list of (tool,value); bold the minimum."""
    present = [(t, v) for t, v in vals if v is not None]
    if not present:
        return {t: "--" for t, _ in vals}
    best = min(v for _, v in present)
    return {t: (f"\\textbf{{{v:g}}}" if v == best else f"{v:g}") if v is not None else "--"
            for t, v in vals}


def t_cx(data, topo, label_suffix):
    rows = []
    for fx, (fam, n) in fixtures_on(data, topo):
        vals = [(t, data.get((fx, topo, t), {}).get("cx")) for t in TOOLS]
        fmt = _bold_min(vals)
        rows.append([fx.replace("_", "\\_"), n] + [fmt[t] for t in TOOLS])
    return tabular(
        f"CX-count comparison at best effort, {label_suffix}. Lowest per row in "
        f"bold. CX-equivalent (SWAP = 3 CX).",
        f"tab:cmp_cx_{topo}",
        ["Circuit", "$n$"] + [TOOL_LABEL[t] for t in TOOLS], rows)


def t_speed(data, topo):
    rows = []
    speedups = []
    for fx, (fam, n) in fixtures_on(data, topo):
        vals = [(t, data.get((fx, topo, t), {}).get("compile_ms")) for t in TOOLS]
        cells = []
        for t, v in vals:
            cells.append("--" if v is None else (f"\\textbf{{{v:.3f}}}" if v == min(
                x for _, x in vals if x is not None) else f"{v:.3f}"))
        qr = data.get((fx, topo, "qrust"), {}).get("compile_ms")
        qk = data.get((fx, topo, "qiskit"), {}).get("compile_ms")
        if qr and qk:
            speedups.append(qk / qr)
        rows.append([fx.replace("_", "\\_"), n] + cells)
    cap = (f"Compile-time comparison (ms, median) at best effort, all-to-all. "
           f"Q-Rust is native code; the others are Python. ")
    if speedups:
        cap += f"Median Q-Rust speedup vs.\\ Qiskit: {statistics.median(speedups):.0f}$\\times$."
    return tabular(cap, "tab:cmp_speed",
                   ["Circuit", "$n$"] + [TOOL_LABEL[t] for t in TOOLS], rows)


def t_geomean(data):
    """Summary: geometric-mean CX ratio of each tool vs Q-Rust, all-to-all."""
    rows = []
    for t in TOOLS:
        ratios = []
        for fx, (fam, n) in fixtures_on(data, "all_to_all"):
            qr = data.get((fx, "all_to_all", "qrust"), {}).get("cx")
            tv = data.get((fx, "all_to_all", t), {}).get("cx")
            if qr and tv:
                ratios.append(tv / qr)
        if ratios:
            gm = statistics.geometric_mean(ratios)
            rows.append([TOOL_LABEL[t], f"{gm:.2f}",
                         "fewer CX than Q-Rust" if gm < 1 else "more CX than Q-Rust"])
    return tabular(
        "Geometric-mean CX ratio vs.\\ Q-Rust on all-to-all at O3 ($<1$ means the "
        "tool emits fewer CX than Q-Rust on average).",
        "tab:cmp_geomean", ["Tool", "GM(CX/CX$_{qrust}$)", "Interpretation"], rows,
        colspec="lrl")


# --------------------------------------------------------------------------- #
# Figures                                                                      #
# --------------------------------------------------------------------------- #
def fig_cx(data, topo, fname, title):
    fxs = fixtures_on(data, topo)
    labels = [fx.replace("_", "") for fx, _ in fxs]
    x = range(len(fxs))
    width = 0.2
    fig, ax = plt.subplots(figsize=(max(8, len(fxs) * 0.55), 4.6))
    for i, t in enumerate(TOOLS):
        ys = [data.get((fx, topo, t), {}).get("cx", 0) for fx, _ in fxs]
        ax.bar([xi + (i - 1.5) * width for xi in x], ys, width,
               label=TOOL_LABEL[t], color=TOOL_COLOR[t])
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("CX count (CX-equivalent)")
    ax.set_title(title)
    ax.legend(frameon=False, ncol=4, fontsize=9)
    fig.tight_layout()
    p = os.path.join(RESULTS, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_speed(data):
    fxs = fixtures_on(data, "all_to_all")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for t in TOOLS:
        pts = [(v[1], data.get((fx, "all_to_all", t), {}).get("compile_ms"))
               for fx, v in fxs]
        pts = [(n, ms) for (n, ms) in pts if ms is not None]
        pts.sort()
        if pts:
            xs = [n for n, _ in pts]
            ys = [ms for _, ms in pts]
            ax.scatter(xs, ys, label=TOOL_LABEL[t], color=TOOL_COLOR[t], s=36, alpha=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("Number of qubits $n$")
    ax.set_ylabel("Compile time (ms, log)")
    ax.set_title("Compile-time comparison at O3 (all-to-all)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(RESULTS, "fig_cmp_speed.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def _cx_ratio_grid(data, against="best"):
    """family x topology -> geomean(qrust_cx / reference_cx) over sizes.

    `against` is either "best" (min over Qiskit/tket/Cirq) or a specific tool
    name. <1 means Q-Rust emits fewer CX than the reference.
    Returns (grid dict, families_present, topos_present).
    """
    by_cell = defaultdict(list)  # (family, topo) -> [ratio,...]
    fams, topos = set(), set()
    keyset = {(fx, tp) for (fx, tp, t) in data}
    for (fx, tp) in keyset:
        qr = data.get((fx, tp, "qrust"), {}).get("cx")
        if not qr or qr <= 0:
            continue
        if against == "best":
            ref_vals = [data.get((fx, tp, t), {}).get("cx") for t in OTHERS]
            ref_vals = [c for c in ref_vals if c and c > 0]
            ref = min(ref_vals) if ref_vals else None
        else:
            ref = data.get((fx, tp, against), {}).get("cx")
            ref = ref if ref and ref > 0 else None
        if ref is None:
            continue
        fam = data[(fx, tp, "qrust")]["family"]
        by_cell[(fam, tp)].append(qr / ref)
        fams.add(fam)
        topos.add(tp)
    grid = {k: statistics.geometric_mean(v) for k, v in by_cell.items()}
    fams = [f for f in FAMILY_ORDER if f in fams]
    topos = [t for t in TOPO_ORDER if t in topos]
    return grid, fams, topos


def t_winloss(data):
    grid, fams, topos = _cx_ratio_grid(data)
    rows = []
    wins = ties = losses = 0
    for fam in fams:
        cells = []
        for tp in topos:
            r = grid.get((fam, tp))
            if r is None:
                cells.append("--")
                continue
            if r < 0.98:
                wins += 1
                cells.append(f"\\textbf{{{r:.2f}}}")
            elif r <= 1.02:
                ties += 1
                cells.append(f"{r:.2f}")
            else:
                losses += 1
                cells.append(f"{r:.2f}")
        rows.append([FAMILY_LABEL[fam]] + cells)
    cap = (f"Q-Rust CX-count ratio vs.\\ the best of Qiskit/tket/Cirq, by circuit "
           f"family $\\times$ topology at O3 (geometric mean over sizes). "
           f"Values $<1$ (bold) = Q-Rust wins; $\\approx 1$ = tie; $>1$ = trails. "
           f"Across {wins + ties + losses} combinations: {wins} wins, {ties} ties, "
           f"{losses} losses.")
    return tabular(cap, "tab:cmp_winloss",
                   ["Family"] + [TOPO_LABEL[t] for t in topos], rows)


def t_speed_winloss(data):
    """Median/max compile-time speedup of Q-Rust vs. each baseline (all
    topologies/sizes at O3). Reported per-tool because Cirq is fast only by
    skipping heavy optimization, so a 'vs fastest' figure would understate the
    gap against the real optimizing transpilers (Qiskit, tket)."""
    rows = []
    for tool in OTHERS:
        sp = []
        for (fx, tp, t), v in data.items():
            if t != "qrust":
                continue
            o = data.get((fx, tp, tool), {}).get("compile_ms")
            if v["compile_ms"] and v["compile_ms"] > 0 and o:
                sp.append(o / v["compile_ms"])
        if sp:
            rows.append([TOOL_LABEL[tool], f"{statistics.median(sp):.0f}$\\times$",
                         f"{max(sp):.0f}$\\times$", len(sp)])
    return tabular(
        "Compile-time speedup of Q-Rust vs.\\ each baseline (median and max over "
        "all topologies/sizes at O3). Q-Rust is native code; baselines are Python.",
        "tab:cmp_speed_winloss",
        ["vs.\\ Tool", "Median", "Max", "$n$"], rows, colspec="lrrr")


def _draw_heatmap(ax, data, against, title, show_cbar=True, fig=None):
    """Render one CX-ratio heatmap onto `ax`. Returns the image handle."""
    import numpy as np
    import matplotlib.colors as mcolors
    grid, fams, topos = _cx_ratio_grid(data, against)
    M = np.full((len(fams), len(topos)), np.nan)
    for i, fam in enumerate(fams):
        for j, tp in enumerate(topos):
            if (fam, tp) in grid:
                M[i, j] = grid[(fam, tp)]
    norm = mcolors.TwoSlopeNorm(vmin=0.6, vcenter=1.0, vmax=1.6)
    im = ax.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(len(topos)))
    ax.set_xticklabels([TOPO_LABEL[t] for t in topos], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([FAMILY_LABEL[f] for f in fams], fontsize=8)
    wins = 0
    total = 0
    for i in range(len(fams)):
        for j in range(len(topos)):
            if M[i, j] == M[i, j]:  # not NaN
                total += 1
                wins += M[i, j] < 0.98
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                        fontsize=7, color="black")
    ax.set_title(f"{title}  ({wins}/{total} wins)", fontsize=10)
    if show_cbar and fig is not None:
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label="CX$_{qrust}$ / CX$_{ref}$")
    return im


def fig_winloss(data, against, fname, title):
    grid, fams, topos = _cx_ratio_grid(data, against)
    fig, ax = plt.subplots(figsize=(1.1 * len(topos) + 2.5, 0.55 * len(fams) + 2))
    _draw_heatmap(ax, data, against, title, show_cbar=True, fig=fig)
    fig.tight_layout()
    p = os.path.join(RESULTS, fname)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def fig_winloss_panel(data):
    """2x2 panel: vs best, vs Qiskit, vs tket, vs Cirq."""
    specs = [("best", "vs. best baseline"), ("qiskit", "vs. Qiskit"),
             ("tket", "vs. tket"), ("cirq", "vs. Cirq")]
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    im = None
    for ax, (against, title) in zip(axes.flat, specs):
        im = _draw_heatmap(ax, data, against, title, show_cbar=False)
    fig.suptitle("Q-Rust CX-count ratio by family $\\times$ topology "
                 "(blue $<1$ = Q-Rust fewer CX)", fontsize=13)
    fig.colorbar(im, ax=axes, fraction=0.025, pad=0.04, label="CX$_{qrust}$ / CX$_{ref}$")
    p = os.path.join(RESULTS, "fig_cmp_winloss_panel.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    data = load()
    blocks = [
        "% Auto-generated by export_comparison.py — Q-Rust vs Qiskit/tket/Cirq.",
        "", t_winloss(data), t_speed_winloss(data),
        t_cx(data, "all_to_all", "all-to-all (ideal connectivity)"),
        t_geomean(data),
        t_cx(data, "ibm_nairobi", "IBM Nairobi (cirq lacks placement --- see text)"),
        t_speed(data, "all_to_all"),
    ]
    tex = os.path.join(RESULTS, "comparison_tables.tex")
    open(tex, "w").write("\n".join(blocks))

    figs = [
        fig_winloss(data, "best", "fig_cmp_winloss_heatmap.png",
                    "Q-Rust CX vs. best baseline (<1 = Q-Rust wins)"),
        fig_winloss(data, "qiskit", "fig_cmp_winloss_qiskit.png", "Q-Rust CX vs. Qiskit"),
        fig_winloss(data, "tket", "fig_cmp_winloss_tket.png", "Q-Rust CX vs. tket"),
        fig_winloss(data, "cirq", "fig_cmp_winloss_cirq.png", "Q-Rust CX vs. Cirq"),
        fig_winloss_panel(data),
        fig_cx(data, "all_to_all", "fig_cmp_cx_all2all.png",
               "CX count by transpiler at O3 (all-to-all, ideal)"),
        fig_cx(data, "ibm_nairobi", "fig_cmp_cx_hardware.png",
               "CX count by transpiler at O3 (IBM Nairobi; cirq has no placement)"),
        fig_speed(data),
    ]
    print(f"comparison tables -> {tex}")
    for f in figs:
        print(f"figure            -> {f}")


if __name__ == "__main__":
    main()
