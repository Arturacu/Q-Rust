#!/usr/bin/env python3
"""Generate booktabs \\input fragments and publication-quality figures.

Emits one self-contained `\\begin{tabular}...\\end{tabular}` per table (drop
into the draft's existing table floats via \\input) and the three figures the
draft references. Table/figure names match the \\label names used by the
accompanying LaTeX report.

Data sources:
  results/raw_benchmark_data.json      characterisation suite, 4 topologies, O0-O3
  results/char3_raw_benchmark_data.json  same suite on linear/ring/star (for tab:verify)
  results/abl_raw_benchmark_data.json    QFT ablation, median of 5 (for tab:ablation)
  results/cmp_raw_benchmark_data.json    comparison suite, median of 5, O3, 7 topologies
  results/baselines.json                 Qiskit/tket/Cirq on the same cells

Outputs land in results/ as tab_*.tex and fig_*.png.
"""
from __future__ import annotations

import json
import os
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Serif body font + Computer Modern math, so matplotlib
# figures blend with the LaTeX text. Vector PDF is the primary output (crisp in
# print); a PNG twin is kept for quick preview.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "CMU Serif", "Times New Roman"],
    "mathtext.fontset": "cm",
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "figure.dpi": 150, "savefig.bbox": "tight",
})

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")


def save_fig(fig, stem):
    """Write both a vector PDF (for LaTeX inclusion) and a PNG (for preview)."""
    fig.savefig(os.path.join(R, stem + ".pdf"))
    fig.savefig(os.path.join(R, stem + ".png"), dpi=150)
    plt.close(fig)
    print("  wrote", stem + ".pdf (+ .png)")

FAM_ORDER = ["ghz", "qft", "adder", "clifford", "qaoa", "bv", "ising", "wstate", "graphstate"]
FAM_LABEL = {"ghz": "GHZ", "qft": "QFT", "adder": "Adder", "clifford": "Rand.\\ Clifford",
             "qaoa": "QAOA", "bv": "BV", "ising": "Ising", "wstate": "W-state",
             "graphstate": "Graph state"}
TOPO_ORDER = ["all_to_all", "linear", "ring", "star", "grid_4x4", "ibm_nairobi", "ibm_quito"]
TOPO_LABEL = {"all_to_all": "All-to-all", "linear": "Linear", "ring": "Ring", "star": "Star",
              "grid_4x4": "Grid $4{\\times}4$", "ibm_nairobi": "Nairobi", "ibm_quito": "Quito"}
TOOLS = ["qrust", "qiskit", "tket", "cirq"]
OTHERS = ["qiskit", "tket", "cirq"]
TIE_BAND = 0.04  # Ratios in [0.96, 1.04] are ties.
TOOL_LABEL = {"qrust": "Q-Rust", "qiskit": "Qiskit", "tket": "t$\\vert$ket$\\rangle$", "cirq": "Cirq"}


def load(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else []


def ok(recs, o3_only=False):
    out = [r for r in recs if r.get("status") == "ok"]
    if o3_only:
        out = [r for r in out if r["opt_level"] == 3 and not r.get("ablation")]
    return out


def tabular(colspec, header, rows):
    L = ["\\begin{tabular}{%s}" % colspec, "  \\toprule",
         "  " + " & ".join(header) + " \\\\", "  \\midrule"]
    L += ["  " + " & ".join(str(c) for c in r) + " \\\\" for r in rows]
    L += ["  \\bottomrule", "\\end{tabular}"]
    return "\n".join(L) + "\n"


def write(name, text):
    open(os.path.join(R, name), "w").write(text)
    print("  wrote", name)


# =========================================================================== #
# Comparison data index                                                       #
# =========================================================================== #
def comparison_index():
    cmp = ok(load("cmp_raw_benchmark_data.json"), o3_only=True)
    base = ok(load("baselines.json"))
    d = {}
    meta = {}
    for r in cmp:
        k = (r["fixture"], r["topology"])
        d[(k, "qrust")] = {"cx": r["post"]["cx"], "depth": r["post"]["depth"],
                           "gates": r["post"]["gates"],
                           "ms": (r["timing_us"]["total"] - r["timing_us"]["parse"]) / 1000.0,
                           "cxmax": r["post"].get("cx_max", r["post"]["cx"])}
        meta[k] = {"family": r["family"], "n": r["num_qubits_logical"]}
    for r in base:
        k = (r["fixture"], r["topology"])
        d[(k, r["tool"])] = {"cx": r["cx_count"], "depth": r["depth"], "gates": r["gate_count"],
                             "ms": r["compile_ms_median"]}
    return d, meta


# =========================================================================== #
# CHARACTERISATION TABLES                                                      #
# =========================================================================== #
def t_suite():
    """Family-summary of the fixture suite (too many fixtures for one row each)."""
    man = json.load(open(os.path.join(HERE, "fixtures", "manifest.json")))
    byf = defaultdict(list)
    for f in man:
        byf[f["family"]].append(f)
    rows = []
    for fam in FAM_ORDER:
        fs = byf.get(fam, [])
        if not fs:
            continue
        ws = sorted(f["num_qubits"] for f in fs)
        g = [f["source_gates"] for f in fs]
        cx = [f["source_2q_gates"] for f in fs]
        wr = f"{ws[0]}--{ws[-1]}" if len(set(ws)) > 1 else str(ws[0])
        rows.append([FAM_LABEL[fam], len(fs), wr,
                     f"{min(g)}--{max(g)}" if min(g) != max(g) else str(g[0]),
                     f"{min(cx)}--{max(cx)}" if min(cx) != max(cx) else str(cx[0])])
    n_core = sum(1 for f in man if f["num_qubits"] <= 18)
    rows.append(["\\textbf{Total}", f"\\textbf{{{len(man)}}}",
                 f"{n_core} core, {len(man) - n_core} scaling", "", ""])
    write("tab_suite.tex",
          tabular("lrrrr",
                  ["Family", "\\#", "$n$", "Gates", "2q gates"], rows))


def _verify_recs():
    """Characterisation-suite verification records over all 7 topologies.
    Prefers the consolidated verify7 run; falls back to the 4-topo + char3 pair."""
    v7 = load("verify7_raw_benchmark_data.json")
    recs = ok(v7) if v7 else ok(load("raw_benchmark_data.json")) + ok(load("char3_raw_benchmark_data.json"))
    return [r for r in recs if not r.get("ablation")]


def t_verify():
    recs = _verify_recs()
    rows = []
    tot = defaultdict(int)
    for tp in TOPO_ORDER:
        rs = [r for r in recs if r["topology"] == tp]
        if not rs:
            continue
        c = defaultdict(int)
        fids = []
        for r in rs:
            v = r["verification"]
            c[v["verdict_class"]] += 1
            if v["verified"] and v["fidelity"] is not None:
                fids.append(v["fidelity"])
        for k in c:
            tot[k] += c[k]
        fmin = f"{min(fids):.6f}" if fids else "--"
        rows.append([TOPO_LABEL[tp], c["Exact"], c["Statistical"],
                     c["Unverifiable"], c["NotEquivalent"], fmin])
    rows.append(["\\textbf{Total}", f"\\textbf{{{tot['Exact']}}}",
                 f"\\textbf{{{tot['Statistical']}}}", f"\\textbf{{{tot['Unverifiable']}}}",
                 f"\\textbf{{{tot['NotEquivalent']}}}", "--"])
    write("tab_verify.tex",
          tabular("lrrrrr",
                  ["Topology", "Exact", "Statistical", "Unverifiable", "Not eq.", "$F_{\\min}$"],
                  rows))


def _char_o3(topo):
    return [r for r in ok(load("raw_benchmark_data.json"), o3_only=True) if r["topology"] == topo]


def t_cxred(topo, name):
    recs = _char_o3(topo)
    rows = []
    for fam in ["ghz", "qft", "adder", "clifford"]:
        rs = [r for r in recs if r["family"] == fam]
        if not rs:
            continue
        pre = sum(r["pre"]["cx"] for r in rs)
        post = sum(r["post"]["cx"] for r in rs)
        factor = "n/a" if pre == 0 else f"{post / pre:.2f}$\\times$"
        rows.append([FAM_LABEL[fam], pre, post, factor])
    write(name, tabular("lrrr", ["Family", "CX$_{\\text{pre}}$", "CX$_{\\text{post}}$",
                                 "Post/Pre"], rows))


def t_timing():
    recs = _char_o3("all_to_all")
    rows = []
    for fam in ["ghz", "qft", "adder", "clifford"]:
        rs = [r for r in recs if r["family"] == fam]
        if not rs:
            continue
        m = lambda k: statistics.mean(r["timing_us"][k] for r in rs)
        rows.append([FAM_LABEL[fam], f"{m('parse'):.0f}", f"{m('logic_opt'):.0f}",
                     f"{m('routing'):.0f}", f"{m('basis_decomp'):.0f}",
                     f"{m('cleanup'):.0f}", f"{m('total'):.0f}"])
    write("tab_timing.tex",
          tabular("lrrrrrr",
                  ["Family", "Parse", "Logic opt.", "Routing", "Decomp.", "Cleanup", "Total"],
                  rows))


def t_ram():
    recs = _char_o3("all_to_all")
    rows = []
    for fam in ["ghz", "qft", "adder", "clifford"]:
        rs = sorted([r for r in recs if r["family"] == fam], key=lambda r: r["num_qubits_logical"])
        mibs = [r["peak_rss_bytes"] / 2**20 for r in rs]
        if not mibs:
            continue
        rows.append([FAM_LABEL[fam], f"{rs[0]['num_qubits_logical']}--{rs[-1]['num_qubits_logical']}",
                     f"{min(mibs):.1f}", f"{max(mibs):.1f}"])
    write("tab_ram.tex", tabular("lrrr", ["Family", "$n$ range", "Min (MiB)", "Max (MiB)"], rows))


def t_ablation():
    """Controlled ablation: pre_cleanup vs post on the SAME routed circuit
    (median of 5), so the delta isolates the cleanup effect from routing noise."""
    recs = ok(load("abl_raw_benchmark_data.json"))
    cells = {}
    for r in recs:
        if r["family"] != "qft" or r["opt_level"] != 3 or r.get("ablation") \
           or "pre_cleanup" not in r:
            continue
        cells[(r["fixture"], r["topology"])] = r
    if not cells:
        write("tab_ablation.tex", "% abl_raw_benchmark_data.json (with pre_cleanup) not found\n")
        return
    rows = []
    for tp in ["all_to_all", "grid_4x4", "ibm_nairobi", "linear", "ring"]:
        for key in sorted((k for k in cells if k[1] == tp),
                          key=lambda k: cells[k]["num_qubits_logical"]):
            r = cells[key]
            off_cx, on_cx = r["pre_cleanup"]["cx"], r["post"]["cx"]
            off_d, on_d = r["pre_cleanup"]["depth"], r["post"]["depth"]
            dcx, dd = off_cx - on_cx, off_d - on_d  # positive = cleanup removed
            rows.append([key[0].replace("_", "\\_"), TOPO_LABEL[tp],
                         int(on_cx), int(off_cx), f"{int(dcx):+d}",
                         int(on_d), int(off_d), f"{int(dd):+d}"])
    write("tab_ablation.tex",
          tabular("llrrrrrr",
                  ["Circuit", "Topology", "CX$_{\\text{on}}$", "CX$_{\\text{off}}$", "$\\Delta$CX",
                   "D$_{\\text{on}}$", "D$_{\\text{off}}$", "$\\Delta$D"], rows))


# =========================================================================== #
# COMPARISON TABLES                                                           #
# =========================================================================== #
def _ratio_grid(d, meta, metric="cx", against="best"):
    cell = defaultdict(list)
    fams, topos = set(), set()
    for k in {kk for (kk, t) in d}:
        q = d.get((k, "qrust"), {}).get(metric)
        if not q or q <= 0:
            continue
        if against == "best":
            ov = [d.get((k, t), {}).get(metric) for t in OTHERS]
            ov = [v for v in ov if v and v > 0]
            ref = min(ov) if ov else None
        else:
            ref = d.get((k, against), {}).get(metric)
            ref = ref if ref and ref > 0 else None
        if not ref:
            continue
        cell[(meta[k]["family"], k[1])].append(q / ref)
        fams.add(meta[k]["family"])
        topos.add(k[1])
    grid = {kk: statistics.geometric_mean(v) for kk, v in cell.items()}
    return grid, [f for f in FAM_ORDER if f in fams], [t for t in TOPO_ORDER if t in topos]


def t_cmp_winloss(d, meta):
    grid, fams, topos = _ratio_grid(d, meta)
    rows = []
    for fam in fams:
        cells = []
        for tp in topos:
            r = grid.get((fam, tp))
            if r is None:
                cells.append("--")
            elif r < 1.0 - TIE_BAND:
                cells.append(f"\\textbf{{{r:.2f}}}")
            else:
                cells.append(f"{r:.2f}")
        rows.append([FAM_LABEL[fam]] + cells)
    write("tab_cmp_winloss.tex",
          tabular("l" + "r" * len(topos), ["Family"] + [TOPO_LABEL[t] for t in topos], rows))


def t_cmp_geomean(d, meta):
    rows = []
    for tool in OTHERS:
        rs = [d[(k, tool)]["cx"] / d[(k, "qrust")]["cx"]
              for k in {kk for (kk, t) in d} if k[1] == "all_to_all"
              and d.get((k, tool)) and d.get((k, "qrust"))]
        if rs:
            rows.append([TOOL_LABEL[tool], f"{statistics.geometric_mean(rs):.2f}"])
    write("tab_cmp_geomean.tex",
          tabular("lr", ["Tool", "GM(CX$_{\\text{tool}}$/CX$_{\\text{Q-Rust}}$)"], rows))


def _representative(d, meta, topo, per_family=3):
    """Up to `per_family` fixtures per family on `topo`, spanning the widths
    (smallest, a middle, largest) so per-circuit tables stay readable."""
    byf = defaultdict(list)
    for fx in {k[0] for (k, t) in d if k[1] == topo}:
        byf[meta[(fx, topo)]["family"]].append(fx)
    out = []
    for fam in FAM_ORDER:
        fs = sorted(byf.get(fam, []), key=lambda fx: meta[(fx, topo)]["n"])
        if not fs:
            continue
        if len(fs) <= per_family:
            pick = fs
        else:
            pick = [fs[0], fs[len(fs) // 2], fs[-1]]
        out.extend(pick)
    return out


def t_cmp_cx(d, meta, topo, name):
    rows = []
    fxs = _representative(d, meta, topo)
    for fx in fxs:
        k = (fx, topo)
        vals = {t: d.get((k, t), {}).get("cx") for t in TOOLS}
        present = [v for v in vals.values() if v is not None]
        best = min(present) if present else None
        cells = []
        for t in TOOLS:
            v = vals[t]
            if v is None:
                cells.append("--")
            elif v == best:
                cells.append(f"\\textbf{{{v}}}")
            else:
                cells.append(str(v))
        label = fx.replace("_", "\\_") + ("$^{\\dagger}$" if fx == "ghz_25" else "")
        rows.append([label, meta[k]["n"]] + cells)
    write(name, tabular("lr" + "r" * 4, ["Circuit", "$n$"] + [TOOL_LABEL[t] for t in TOOLS], rows))


def t_cmp_speed_winloss(d, meta):
    rows = []
    keys = {kk for (kk, t) in d}
    for tool in OTHERS:
        sp = []
        for k in keys:
            q = d.get((k, "qrust"), {}).get("ms")
            o = d.get((k, tool), {}).get("ms")
            if q and q > 0 and o:
                sp.append(o / q)
        if sp:
            rows.append([TOOL_LABEL[tool], f"{statistics.median(sp):.0f}$\\times$",
                         f"{max(sp):.0f}$\\times$"])
    write("tab_cmp_speed_winloss.tex",
          tabular("lrr", ["Baseline", "Median speedup", "Max speedup"], rows))


def t_cmp_speed(d, meta, topo="all_to_all"):
    rows = []
    fxs = _representative(d, meta, topo)
    for fx in fxs:
        k = (fx, topo)
        vals = {t: d.get((k, t), {}).get("ms") for t in TOOLS}
        present = [v for v in vals.values() if v is not None]
        best = min(present) if present else None
        cells = []
        for t in TOOLS:
            v = vals[t]
            if v is None:
                cells.append("--")
            elif v == best:
                cells.append(f"\\textbf{{{v:.3f}}}")
            else:
                cells.append(f"{v:.3f}")
        label = fx.replace("_", "\\_") + ("$^{\\dagger}$" if fx == "ghz_25" else "")
        rows.append([label, meta[k]["n"]] + cells)
    write("tab_cmp_speed.tex",
          tabular("lr" + "r" * 4, ["Circuit", "$n$"] + [TOOL_LABEL[t] + " (ms)" for t in TOOLS], rows))


def t_cmp_depth(d, meta):
    """Depth and total-gate-count ratios vs best baseline, by family (fills the gap)."""
    dg, fams, _ = _ratio_grid(d, meta, metric="depth")
    gg, _, _ = _ratio_grid(d, meta, metric="gates")
    # aggregate to family (geomean over topologies)
    def fam_geo(grid):
        by = defaultdict(list)
        for (fam, tp), v in grid.items():
            by[fam].append(v)
        return {f: statistics.geometric_mean(v) for f, v in by.items()}
    dgm, ggm = fam_geo(dg), fam_geo(gg)
    rows = []
    for fam in FAM_ORDER:
        if fam in dgm:
            rows.append([FAM_LABEL[fam], f"{dgm[fam]:.2f}", f"{ggm.get(fam, float('nan')):.2f}"])
    write("tab_cmp_depth.tex",
          tabular("lrr", ["Family", "Depth ratio", "Total-gate ratio"], rows))
    # also print win/tie/loss for the prose (XX/YY)
    wins = sum(1 for v in dg.values() if v < 1.0 - TIE_BAND)
    ties = sum(1 for v in dg.values() if 1.0 - TIE_BAND <= v <= 1.0 + TIE_BAND)
    print(f"    [depth] vs-best: {wins} wins, {ties} ties, {len(dg) - wins - ties} losses "
          f"of {len(dg)} combos")


# =========================================================================== #
# FIGURES                                                                      #
# =========================================================================== #
import numpy as np  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402

TOOL_COLOR = {"qrust": "#5b2c9e", "qiskit": "#d62728", "tket": "#2ca02c", "cirq": "#7f7f7f"}
MARK = {"qrust": "o", "qiskit": "s", "tket": "^", "cirq": "D"}


def _tlab(tool):
    return TOOL_LABEL[tool].replace("$\\vert$", "|").replace("$\\rangle$", ">")


def _topo_txt(t):
    return TOPO_LABEL[t].replace("$", "").replace("\\times", "x").replace("{", "").replace("}", "")


def _fam_txt(f):
    return FAM_LABEL[f].replace("\\ ", " ").replace("Rand. ", "")


def _draw_ratio_heatmap(ax, d, meta, tool, big=False):
    """Draw the Q-Rust/baseline CX-ratio heatmap onto ax; return (im, wins, n)."""
    grid, fams, topos = _ratio_grid(d, meta, against=tool)
    M = np.full((len(fams), len(topos)), np.nan)
    for i, fam in enumerate(fams):
        for j, tp in enumerate(topos):
            if (fam, tp) in grid:
                M[i, j] = grid[(fam, tp)]
    norm = mcolors.TwoSlopeNorm(vmin=0.85, vcenter=1.0, vmax=1.5)
    im = ax.imshow(np.clip(M, 0.85, 1.5), cmap="RdBu_r", norm=norm, aspect="auto")
    fs = 9 if big else 6.5
    ax.set_xticks(range(len(topos)))
    ax.set_xticklabels([_topo_txt(t) for t in topos], rotation=40, ha="right", fontsize=fs + 1)
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([_fam_txt(f) for f in fams], fontsize=fs + 1)
    wins = n = 0
    for i in range(len(fams)):
        for j in range(len(topos)):
            if M[i, j] == M[i, j]:
                n += 1
                tie = 1.0 - TIE_BAND <= M[i, j] <= 1.0 + TIE_BAND
                wins += M[i, j] < 1.0 - TIE_BAND
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=fs,
                        color="0.5" if tie else "black")
            else:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, hatch="///",
                                           fill=False, edgecolor="0.6", lw=0))
    return im, wins, n


def fig_cmp_panel(d, meta):
    """Combined 1x3 overview (kept as an optional compact figure)."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    im = None
    for ax, tool in zip(axes, OTHERS):
        im, wins, n = _draw_ratio_heatmap(ax, d, meta, tool)
        ax.set_title(f"vs. {_tlab(tool)}  ({wins}/{n} Q-Rust wins)", fontsize=11)
    fig.suptitle("Q-Rust CX-count ratio by family (rows) $\\times$ topology (columns): "
                 "blue $<1$ = Q-Rust fewer CX", fontsize=12)
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="CX ratio")
    save_fig(fig, "fig_cmp_panel_3")


def fig_cmp_vs(d, meta, tool):
    """Standalone one-to-one heatmap: Q-Rust vs a single baseline (legible)."""
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    im, wins, n = _draw_ratio_heatmap(ax, d, meta, tool, big=True)
    ax.set_title(f"Q-Rust CX count relative to {_tlab(tool)} "
                 f"({wins} of {n} combinations in Q-Rust's favour)", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, label="CX ratio (Q-Rust / baseline)")
    fig.tight_layout()
    save_fig(fig, f"fig_cmp_vs_{tool}")


def fig_cmp_tradeoff(d, meta):
    """All-tools summary over the common successful benchmark instances.

    Quality is the geometric mean of cell-level tool/Q-Rust CX ratios; time is
    the median of each cell's median-over-sizes compile time. Restricting the
    input to instances completed by every tool prevents missing cells from
    changing the workload, while cell-level aggregation prevents families with
    more fixture widths from receiving disproportionate weight.
    """
    common_keys = [
        k for k in meta
        if all(
            d.get((k, tool), {}).get("cx", 0) > 0
            and d.get((k, tool), {}).get("ms", 0) > 0
            for tool in TOOLS
        )
    ]
    if not common_keys:
        raise ValueError("no benchmark instances succeeded for all tools")

    common_cells = {(meta[k]["family"], k[1]) for k in common_keys}
    pts = {}
    for tool in TOOLS:
        quality_by_cell = defaultdict(list)
        times_by_cell = defaultdict(list)
        for k in common_keys:
            cell = (meta[k]["family"], k[1])
            quality_by_cell[cell].append(d[(k, tool)]["cx"] / d[(k, "qrust")]["cx"])
            times_by_cell[cell].append(d[(k, tool)]["ms"])

        cell_quality = [statistics.geometric_mean(values)
                        for values in quality_by_cell.values()]
        cell_medians = [statistics.median(values) for values in times_by_cell.values()]
        pts[tool] = (statistics.geometric_mean(cell_quality),
                     statistics.median(cell_medians))
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for tool, (x, y) in pts.items():
        ax.scatter(x, y, s=150, color=TOOL_COLOR[tool], marker=MARK[tool],
                   zorder=3, edgecolor="black", linewidth=0.6)
        ax.annotate(_tlab(tool), (x, y), textcoords="offset points", xytext=(9, 5), fontsize=11)
    ax.axvline(1.0, color="0.6", ls="--", lw=1, label="Q-Rust reference")
    ax.set_yscale("log")
    ax.set_xlabel("Geometric-mean CX ratio to Q-Rust  (left = fewer two-qubit gates)")
    ax.set_ylabel("Median cell-level compile time (ms, log)  (down = faster)")
    ax.set_title("Quality--speed trade-off on common successful cells "
                 f"({len(common_keys)} instances; {len(common_cells)} cells)")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    save_fig(fig, "fig_cmp_tradeoff")


def fig_cmp_record(d, meta):
    """Overall head-to-head record: Q-Rust win/tie/loss per baseline."""
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    labels, W, T, L = [], [], [], []
    for tool in OTHERS:
        grid, _, _ = _ratio_grid(d, meta, against=tool)
        labels.append(_tlab(tool))
        W.append(sum(v < 1.0 - TIE_BAND for v in grid.values()))
        T.append(sum(1.0 - TIE_BAND <= v <= 1.0 + TIE_BAND for v in grid.values()))
        L.append(sum(v > 1.0 + TIE_BAND for v in grid.values()))
    ax.barh(labels, W, color="#2ca02c", label="Q-Rust fewer CX")
    ax.barh(labels, T, left=W, color="#b0b0b0", label="tie (within $\\pm$4%)")
    ax.barh(labels, L, left=[a + b for a, b in zip(W, T)], color="#d62728", label="Q-Rust more CX")
    for i, (w, t, l) in enumerate(zip(W, T, L)):
        ax.text(w / 2, i, str(w), ha="center", va="center", fontsize=9, color="white")
        ax.text(w + t / 2, i, str(t), ha="center", va="center", fontsize=9)
        ax.text(w + t + l / 2, i, str(l), ha="center", va="center", fontsize=9, color="white")
    ax.set_xlabel("Family $\\times$ topology combinations")
    ax.set_title("Q-Rust head-to-head CX record vs. each baseline")
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.invert_yaxis()
    fig.tight_layout()
    save_fig(fig, "fig_cmp_record")


def fig_cmp_depth_bar(d, meta):
    """Table tab:cmp_depth as a grouped bar: depth and total-gate ratio/family."""
    dg, _, _ = _ratio_grid(d, meta, metric="depth")
    gg, _, _ = _ratio_grid(d, meta, metric="gates")

    def fam_geo(grid):
        by = defaultdict(list)
        for (fam, tp), v in grid.items():
            by[fam].append(v)
        return {f: statistics.geometric_mean(v) for f, v in by.items()}
    dgm, ggm = fam_geo(dg), fam_geo(gg)
    fams = [f for f in FAM_ORDER if f in dgm]
    y = np.arange(len(fams))
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.barh(y - 0.2, [dgm[f] for f in fams], height=0.38, color="#9b6dd1", label="Depth")
    ax.barh(y + 0.2, [ggm[f] for f in fams], height=0.38, color="#f28e2b", label="Total gates")
    ax.axvline(1.0, color="0.5", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([_fam_txt(f) for f in fams])
    ax.set_xlabel("Ratio to best baseline (1.0 = matches best; $<1$ = fewer)")
    ax.set_title("Depth and total gate count vs. the best baseline (O3)")
    ax.legend(frameon=False)
    ax.invert_yaxis()
    fig.tight_layout()
    save_fig(fig, "fig_cmp_depth")


def fig_cmp_scaling(d, meta):
    """Two panels: QFT (superlinear gate growth) and GHZ (linear), so the
    Qiskit 5->7 step is shown to be family-independent."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, fam, title in zip(axes, ("qft", "ghz"), ("QFT", "GHZ")):
        for tool in TOOLS:
            pts = sorted((meta[k]["n"], d[(k, tool)]["ms"])
                         for k in {kk for (kk, t) in d}
                         if k[1] == "all_to_all" and meta[k]["family"] == fam and d.get((k, tool)))
            if pts:
                xs, ys = zip(*pts)
                ax.plot(xs, ys, MARK[tool] + "-", color=TOOL_COLOR[tool],
                        label=_tlab(tool), lw=1.8, ms=6)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Number of qubits $n$")
        ax.set_title(f"{title} (all-to-all, O3)")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Compile time (ms, log)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_fig(fig, "fig_cmp_scaling")


def fig_stage_breakdown():
    recs = _char_o3("all_to_all")
    stages = ["parse", "logic_opt", "routing", "basis_decomp", "cleanup"]
    labels = ["Parse", "Logic opt.", "Routing", "Decomp.", "Cleanup"]
    colors = ["#c9b6e4", "#9b6dd1", "#f28e2b", "#b07aa1", "#d4a6c8"]
    fams = ["ghz", "qft", "adder", "clifford"]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    left = np.zeros(len(fams))
    for st, lab, col in zip(stages, labels, colors):
        vals = []
        for fam in fams:
            rs = [r for r in recs if r["family"] == fam]
            tot = statistics.mean(r["timing_us"]["total"] for r in rs)
            vals.append(100 * statistics.mean(r["timing_us"][st] for r in rs) / tot)
        ax.barh([FAM_LABEL[f].replace("\\ ", " ") for f in fams], vals, left=left,
                label=lab, color=col)
        if st == "routing":
            for i, v in enumerate(vals):
                ax.text(left[i] + v / 2, i, f"{v:.0f}%", ha="center", va="center", fontsize=8)
        left += np.array(vals)
    ax.set_xlabel("Share of transpile time (%)")
    ax.set_xlim(0, 100)
    ax.set_title("Per-stage compile-time breakdown (all-to-all, O3)")
    ax.legend(ncol=5, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    fig.tight_layout()
    save_fig(fig, "fig_stage_breakdown")


def fig_ablation():
    """Cleanup CX-reduction (off - on) vs circuit size, one line per topology;
    visualizes that the benefit grows with size and topology sparsity."""
    recs = ok(load("abl_raw_benchmark_data.json"))
    cells = {}
    for r in recs:
        if r["family"] != "qft" or r["opt_level"] != 3 or r.get("ablation") \
           or "pre_cleanup" not in r:
            continue
        cells[(r["fixture"], r["topology"])] = r
    if not cells:
        print("  (skip fig_ablation: no abl data)")
        return
    fig, ax = plt.subplots(figsize=(7, 4.4))
    styles = {"all_to_all": ("#9b6dd1", "o"), "grid_4x4": ("#f28e2b", "s"),
              "ibm_nairobi": ("#2ca02c", "^"), "linear": ("#d62728", "D"),
              "ring": ("#1f77b4", "v")}
    for tp, (col, mk) in styles.items():
        pts = sorted((cells[k]["num_qubits_logical"],
                      cells[k]["pre_cleanup"]["cx"] - cells[k]["post"]["cx"])
                     for k in cells if k[1] == tp)
        if pts:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, mk + "-", color=col, lw=1.8, ms=6,
                    label=TOPO_LABEL[tp].replace("$", "").replace("\\times", "x")
                    .replace("{", "").replace("}", ""))
    ax.set_xlabel("Number of qubits $n$")
    ax.set_ylabel("CX gates removed by cleanup")
    ax.set_title("Post-routing cleanup benefit on QFT (O3, median of 5)")
    ax.legend(frameon=False, fontsize=9, title="Topology")
    fig.tight_layout()
    save_fig(fig, "fig_ablation")


def fig_verify():
    """Stacked bar of verdict distribution per topology (characterisation suite)."""
    recs = _verify_recs()
    classes = ["Exact", "Statistical", "Unverifiable"]
    cols = {"Exact": "#2ca02c", "Statistical": "#f0c000", "Unverifiable": "#b0b0b0"}
    topos = [t for t in TOPO_ORDER if any(r["topology"] == t for r in recs)]
    counts = {c: [] for c in classes}
    for tp in topos:
        rs = [r for r in recs if r["topology"] == tp]
        for c in classes:
            counts[c].append(sum(1 for r in rs if r["verification"]["verdict_class"] == c))
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    left = np.zeros(len(topos))
    labels = [TOPO_LABEL[t].replace("$", "").replace("\\times", "x").replace("{", "").replace("}", "")
              for t in topos]
    for c in classes:
        ax.barh(labels, counts[c], left=left, color=cols[c], label=c, edgecolor="white")
        left += np.array(counts[c])
    ax.set_xlabel("Runs (21 circuits $\\times$ 4 optimisation levels, minus width skips)")
    ax.set_title("Verification verdicts by topology")
    ax.legend(frameon=False, ncol=3, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.15))
    ax.invert_yaxis()
    fig.tight_layout()
    save_fig(fig, "fig_verify")


def main():
    print("Characterisation tables:")
    t_suite()
    t_verify()
    t_cxred("all_to_all", "tab_cxred_all_to_all.tex")
    t_cxred("ibm_nairobi", "tab_cxred_ibm_nairobi.tex")
    t_timing()
    t_ram()
    t_ablation()
    print("Comparison tables:")
    d, meta = comparison_index()
    t_cmp_winloss(d, meta)
    t_cmp_geomean(d, meta)
    t_cmp_cx(d, meta, "all_to_all", "tab_cmp_cx_all_to_all.tex")
    t_cmp_cx(d, meta, "ibm_nairobi", "tab_cmp_cx_ibm_nairobi.tex")
    t_cmp_speed_winloss(d, meta)
    t_cmp_speed(d, meta)
    t_cmp_depth(d, meta)
    print("Figures:")
    fig_cmp_panel(d, meta)                 # combined overview (optional)
    for tool in OTHERS:                    # standalone one-to-one heatmaps
        fig_cmp_vs(d, meta, tool)
    fig_cmp_tradeoff(d, meta)              # overall quality-speed trade-off
    fig_cmp_record(d, meta)                # overall win/tie/loss record
    fig_cmp_depth_bar(d, meta)             # tab:cmp_depth as a graph
    fig_cmp_scaling(d, meta)
    fig_stage_breakdown()
    fig_ablation()
    fig_verify()


if __name__ == "__main__":
    main()
