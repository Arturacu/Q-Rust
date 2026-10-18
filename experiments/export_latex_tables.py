#!/usr/bin/env python3
"""Turn raw benchmark data into LaTeX tables and figures.

Reads experiments/results/raw_benchmark_data.json (written by bench_runner.py)
and emits:

  results/tables.tex   pre-formatted \\begin{tabular} blocks, one per results
                       section (suite overview, CX reduction, stage-timing
                       breakdown, peak RAM, verification summary, cleanup ablation)
  results/fig_*.png    publication-style figures (matplotlib)

The script is defensive: crashed/skipped cells (e.g. the adder-on-hardware
Toffoli-routing defect) are excluded from quantitative tables but reported in a
dedicated "failures" table so downstream reports can account for them honestly.

Usage:
    python3 experiments/export_latex_tables.py
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
RAW = os.path.join(RESULTS, "raw_benchmark_data.json")

FAMILY_ORDER = ["ghz", "qft", "adder", "clifford"]
FAMILY_LABEL = {"ghz": "GHZ", "qft": "QFT", "adder": "Adder", "clifford": "Rand.\\ Clifford"}
TOPO_LABEL = {
    "ibm_quito": "IBM Quito (5q)", "ibm_nairobi": "IBM Nairobi (7q)",
    "all_to_all": "All-to-all", "grid_4x4": "Grid $4\\times4$ (16q)",
}
TOPO_ORDER = ["ibm_quito", "ibm_nairobi", "grid_4x4", "all_to_all"]


def load() -> list[dict]:
    if not os.path.exists(RAW):
        raise SystemExit(f"missing {RAW}; run bench_runner.py first")
    return json.load(open(RAW))


def ok(recs):
    return [r for r in recs if r.get("status") == "ok"]


# --------------------------------------------------------------------------- #
# LaTeX helpers                                                                #
# --------------------------------------------------------------------------- #
def tabular(caption, label, header, rows, colspec=None):
    ncol = len(header)
    colspec = colspec or ("l" + "r" * (ncol - 1))
    lines = [
        "\\begin{table}[htbp]", "  \\centering",
        f"  \\caption{{{caption}}}", f"  \\label{{{label}}}",
        f"  \\begin{{tabular}}{{{colspec}}}", "    \\hline",
        "    " + " & ".join(header) + " \\\\", "    \\hline",
    ]
    for row in rows:
        lines.append("    " + " & ".join(str(c) for c in row) + " \\\\")
    lines += ["    \\hline", "  \\end{tabular}", "\\end{table}", ""]
    return "\n".join(lines)


def fmt_pct(v):
    return "--" if v is None or v == "" else f"{float(v):.1f}"


def fmt_f(v, p=4):
    return "--" if v is None or v == "" else f"{float(v):.{p}f}"


# --------------------------------------------------------------------------- #
# Tables                                                                       #
# --------------------------------------------------------------------------- #
def t_suite(recs):
    """Table: fixture suite overview (one row per distinct fixture)."""
    seen = {}
    for r in ok(recs):
        k = r["fixture"]
        if k not in seen:
            seen[k] = r
    rows = []
    def sortkey(r):
        return (FAMILY_ORDER.index(r["family"]), r["num_qubits_logical"])
    for r in sorted(seen.values(), key=sortkey):
        rows.append([
            r["fixture"].replace("_", "\\_"), FAMILY_LABEL[r["family"]],
            r["num_qubits_logical"], r["pre"]["gates"], r["pre"]["two_q"],
            r["pre"]["cx"], r["pre"]["depth"],
        ])
    return tabular(
        "Benchmark fixture suite: source-level metrics (pre-transpilation).",
        "tab:suite",
        ["Circuit", "Family", "$n$", "Gates", "2q", "CX", "Depth"], rows)


def t_cx_reduction(recs):
    """Table: source vs post-transpilation CX, signed reduction %, AND the
    post/pre inflation factor (which makes the decomposition cost explicit:
    for {U,CX}, ccx->6 CX and crz->2 CX, so the ratio is typically >1).

    Computed at O3 on a representative hardware topology (IBM Nairobi) and on
    the ideal all-to-all, using ratio-of-means (sum post / sum pre) rather than
    mean-of-ratios, so widely-different fixture sizes aggregate honestly.
    """
    blocks = []
    for topo in ("ibm_nairobi", "all_to_all"):
        rows = []
        for fam in FAMILY_ORDER:
            rs = [r for r in ok(recs) if r["family"] == fam
                  and r["topology"] == topo and r["opt_level"] == 3 and not r["ablation"]]
            if not rs:
                continue
            sum_pre = sum(r["pre"]["cx"] for r in rs)
            sum_post = sum(r["post"]["cx"] for r in rs)
            if sum_pre == 0:
                red = "--"        # no source CX (QFT): reduction undefined
                ratio = "$\\infty$"  # CX created from scratch by decomposition
            else:
                red = f"{(sum_pre - sum_post) / sum_pre * 100:+.0f}\\%"
                ratio = f"{sum_post / sum_pre:.2f}$\\times$"
            rows.append([FAMILY_LABEL[fam], sum_pre, sum_post, red, ratio])
        blocks.append(tabular(
            f"Source vs.\\ post-transpilation CX count (summed over fixtures), "
            f"signed reduction, and inflation factor (post/pre) at O3 --- "
            f"{TOPO_LABEL[topo]}. Negative reduction / ratio $>1$ reflects basis "
            f"decomposition into \\{{U,CX\\}} (ccx$\\to$6, crz$\\to$2 CX).",
            f"tab:cxred_{topo}",
            ["Family", "CX$_{pre}$", "CX$_{post}$", "Reduction", "Post/Pre"], rows))
    return "\n".join(blocks)


def t_timing(recs):
    """Table: mean per-stage timing breakdown (us) at O3, all-to-all."""
    rows = []
    for fam in FAMILY_ORDER:
        rs = [r for r in ok(recs)
              if r["family"] == fam and r["topology"] == "all_to_all"
              and r["opt_level"] == 3 and not r["ablation"]]
        if not rs:
            continue
        def m(stage):
            return statistics.mean(r["timing_us"][stage] for r in rs)
        rows.append([
            FAMILY_LABEL[fam], f"{m('parse'):.0f}", f"{m('logic_opt'):.0f}",
            f"{m('routing'):.0f}", f"{m('basis_decomp'):.0f}",
            f"{m('cleanup'):.0f}", f"{m('total'):.0f}",
        ])
    return tabular(
        "Per-stage compilation time (mean, microseconds) at O3 on the "
        "all-to-all topology.",
        "tab:timing",
        ["Family", "Parse", "LogicOpt", "Routing", "Decomp", "Cleanup", "Total"], rows)


def t_ram(recs):
    """Table: peak RSS (MiB) by family and width at O3, all-to-all."""
    rows = []
    for r in sorted([r for r in ok(recs)
                     if r["topology"] == "all_to_all" and r["opt_level"] == 3
                     and not r["ablation"]],
                    key=lambda r: (FAMILY_ORDER.index(r["family"]), r["num_qubits_logical"])):
        rows.append([
            r["fixture"].replace("_", "\\_"), r["num_qubits_logical"],
            f"{r['peak_rss_bytes'] / (1024 * 1024):.1f}",
        ])
    return tabular(
        "Peak resident-set size (getrusage) during transpilation, O3, all-to-all.",
        "tab:ram", ["Circuit", "$n$", "Peak RSS (MiB)"], rows)


def t_verification(recs):
    """Table: verification verdict distribution by topology + fidelity stats."""
    rows = []
    for topo in TOPO_ORDER:
        rs = [r for r in ok(recs) if r["topology"] == topo and not r["ablation"]]
        if not rs:
            continue
        cnt = defaultdict(int)
        fids = []
        for r in rs:
            cnt[r["verification"]["verdict_class"]] += 1
            f = r["verification"]["fidelity"]
            if f is not None:
                fids.append(f)
        minf = f"{min(fids):.6f}" if fids else "--"
        rows.append([
            TOPO_LABEL[topo], cnt["Exact"], cnt["Statistical"],
            cnt["Unverifiable"], cnt["NotEquivalent"], minf,
        ])
    return tabular(
        "Verification verdicts by topology (layout-aware harness). $F_{\\min}$ "
        "is the minimum recorded fidelity: process fidelity for exact checks "
        "and sampled output-state fidelity for statistical checks.",
        "tab:verify",
        ["Topology", "Exact", "Stat.", "Unver.", "NotEq.", "$F_{\\min}$"], rows)


def t_ablation(recs):
    """Cleanup-ablation table: QFT O3 cleanup ON vs OFF (Stage-5 bypass) deltas."""
    on = {}
    off = {}
    for r in ok(recs):
        if r["family"] != "qft" or r["opt_level"] != 3:
            continue
        key = (r["fixture"], r["topology"])
        if r["ablation"]:
            off[key] = r
        elif r["cleanup_enabled"]:
            on[key] = r
    rows = []
    for key in sorted(on.keys() & off.keys(), key=lambda k: (k[1], on[k]["num_qubits_logical"])):
        a, b = on[key], off[key]
        d_cx = b["post"]["cx"] - a["post"]["cx"]
        d_depth = b["post"]["depth"] - a["post"]["depth"]
        rows.append([
            key[0].replace("_", "\\_"), TOPO_LABEL[key[1]].split(" (")[0],
            a["post"]["cx"], b["post"]["cx"], f"{d_cx:+d}",
            a["post"]["depth"], b["post"]["depth"], f"{d_depth:+d}",
        ])
    return tabular(
        "Ablation: effect of Stage-5 post-routing cleanup on QFT at "
        "O3. `ON' = full pipeline; `OFF' = cleanup bypassed. $\\Delta$ = OFF $-$ ON "
        "(positive means cleanup removed that many).",
        "tab:ablation",
        ["Circuit", "Topo", "CX$_{on}$", "CX$_{off}$", "$\\Delta$CX",
         "D$_{on}$", "D$_{off}$", "$\\Delta$D"], rows)


def t_failures(recs):
    """Table: runs that crashed/were excluded, for honest accounting."""
    fails = [r for r in recs if r.get("status") != "ok"]
    if not fails:
        return "% (no failed runs)\n"
    agg = defaultdict(lambda: [0, ""])
    for r in fails:
        key = (r.get("fixture", "?"), r.get("topology", "?"))
        agg[key][0] += 1
        agg[key][1] = r.get("error", "").splitlines()[0][:60] if r.get("error") else r.get("status")
    rows = []
    for (fx, topo), (n, msg) in sorted(agg.items()):
        rows.append([fx.replace("_", "\\_"), topo, n, "\\texttt{" + msg.replace("_", "\\_") + "}"])
    return tabular(
        "Excluded runs (transpiler defects surfaced by the suite). The adder "
        "fixtures crash when routed on limited-connectivity devices: Q-Rust "
        "routes the native CCX before decomposing it, corrupting its operand "
        "list (\\texttt{gate\\_def.rs:245}).",
        "tab:failures",
        ["Circuit", "Topology", "\\#runs", "First error"], rows,
        colspec="llrl")


# --------------------------------------------------------------------------- #
# Figures                                                                      #
# --------------------------------------------------------------------------- #
plt.rcParams.update({
    "figure.dpi": 150, "font.size": 11, "axes.grid": True,
    "grid.alpha": 0.3, "axes.axisbelow": True,
})
PALETTE = {"ghz": "#2563eb", "qft": "#dc2626", "adder": "#16a34a", "clifford": "#9333ea"}


def fig_cx_reduction(recs):
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for fam in FAMILY_ORDER:
        xs, ys = [], []
        for o in range(4):
            rs = [r for r in ok(recs) if r["family"] == fam
                  and r["topology"] == "all_to_all" and r["opt_level"] == o
                  and not r["ablation"] and r["cx_reduction_pct"] is not None]
            if rs:
                xs.append(o)
                ys.append(statistics.mean(r["cx_reduction_pct"] for r in rs))
        if xs:
            ax.plot(xs, ys, "o-", lw=2, ms=7, color=PALETTE[fam], label=FAMILY_LABEL[fam].replace("\\", ""))
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"O{i}" for i in range(4)])
    ax.set_xlabel("Optimization level")
    ax.set_ylabel("Mean CX reduction (%)")
    ax.set_title("CX-gate reduction vs. optimization level (all-to-all)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(RESULTS, "fig_cx_reduction.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_stage_timing(recs):
    stages = ["parse", "logic_opt", "routing", "basis_decomp", "cleanup"]
    colors = ["#94a3b8", "#2563eb", "#dc2626", "#16a34a", "#f59e0b"]
    fams = [f for f in FAMILY_ORDER
            if any(r["family"] == f and r["topology"] == "all_to_all" and r["opt_level"] == 3
                   for r in ok(recs))]
    fig, ax = plt.subplots(figsize=(7, 4.3))
    bottom = [0] * len(fams)
    for st, col in zip(stages, colors):
        vals = []
        for fam in fams:
            rs = [r for r in ok(recs) if r["family"] == fam
                  and r["topology"] == "all_to_all" and r["opt_level"] == 3 and not r["ablation"]]
            vals.append(statistics.mean(r["timing_us"][st] for r in rs) if rs else 0)
        ax.bar([FAMILY_LABEL[f].replace("\\", "") for f in fams], vals, bottom=bottom,
               label=st, color=col)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_ylabel("Time (microseconds)")
    ax.set_title("Per-stage compile time at O3 (all-to-all)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    p = os.path.join(RESULTS, "fig_stage_timing.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_scalability(recs):
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for fam in FAMILY_ORDER:
        pts = sorted([(r["num_qubits_logical"], r["timing_us"]["total"])
                      for r in ok(recs) if r["family"] == fam
                      and r["topology"] == "all_to_all" and r["opt_level"] == 3 and not r["ablation"]])
        if pts:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, "o-", lw=2, ms=6, color=PALETTE[fam], label=FAMILY_LABEL[fam].replace("\\", ""))
    ax.set_yscale("log")
    ax.set_xlabel("Number of qubits $n$")
    ax.set_ylabel("Total compile time (microseconds, log)")
    ax.set_title("Transpilation-time scalability at O3 (all-to-all)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(RESULTS, "fig_scalability.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_routing_overhead(recs):
    """SWAP/CX inflation from routing: post CX by topology for GHZ across widths."""
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for topo in TOPO_ORDER:
        pts = sorted([(r["num_qubits_logical"], r["post"]["cx"])
                      for r in ok(recs) if r["family"] == "ghz"
                      and r["topology"] == topo and r["opt_level"] == 3 and not r["ablation"]])
        if pts:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, "o-", lw=2, ms=6, label=TOPO_LABEL[topo].split(" (")[0])
    ax.set_xlabel("Number of qubits $n$")
    ax.set_ylabel("Post-transpilation CX count")
    ax.set_title("Routing overhead: GHZ CX count by topology (O3)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p = os.path.join(RESULTS, "fig_routing_overhead.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def main():
    recs = load()
    os.makedirs(RESULTS, exist_ok=True)

    blocks = [
        "% Auto-generated by export_latex_tables.py — do not edit by hand.",
        "% Requires \\usepackage{booktabs} is NOT assumed; plain \\hline used.",
        "", t_suite(recs), t_cx_reduction(recs), t_timing(recs), t_ram(recs),
        t_verification(recs), t_ablation(recs), t_failures(recs),
    ]
    tex_path = os.path.join(RESULTS, "tables.tex")
    with open(tex_path, "w") as fh:
        fh.write("\n".join(blocks))

    figs = [fig_cx_reduction(recs), fig_stage_timing(recs),
            fig_scalability(recs), fig_routing_overhead(recs)]

    print(f"LaTeX tables -> {tex_path}")
    for f in figs:
        print(f"figure       -> {f}")
    nok = len(ok(recs))
    print(f"\nsummarized {nok}/{len(recs)} successful runs "
          f"({len(recs) - nok} excluded; see tab:failures)")


if __name__ == "__main__":
    main()
