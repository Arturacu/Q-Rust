#!/usr/bin/env python3
"""Sweep orchestrator for the Q-Rust benchmark suite.

Drives the compiled Rust `runner` binary across the full experiment matrix and
collects one structured record per cell:

  Circuits     all fixtures in experiments/fixtures/manifest.json
  Topologies   ibm_quito (5q), ibm_nairobi (7q), all_to_all (ideal), grid_4x4 (16q)
  Opt levels   O0, O1, O2, O3
  Target basis {U, CX}   (the common cross-tool comparison level)

Plus the cleanup ablation: every QFT fixture is additionally run at O3 with
Stage-5 (post-routing cleanup) bypassed, so the cleanup delta can be measured
against the standard O3 (cleanup-enabled) run.

Each cell captures: pre/post gate-CX-depth counts, CX reduction %, total and
per-stage µs timings, peak RSS (getrusage from the runner; psutil cross-check),
and the layout-aware verification verdict + fidelity.

Usage:
    python3 experiments/bench_runner.py --dry-run      # ghz_3, all opts, quito
    python3 experiments/bench_runner.py                # full sweep -> results/
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXTURE_DIR = os.path.join(HERE, "fixtures")
RESULTS_DIR = os.path.join(HERE, "results")
DEFAULT_RUNNER = os.path.join(ROOT, "target", "release", "runner")

TOPOLOGIES = ["ibm_quito", "ibm_nairobi", "all_to_all", "grid_4x4"]
TOPOLOGY_QUBITS = {
    "ibm_quito": 5, "ibm_nairobi": 7, "grid_4x4": 16,
    "all_to_all": None, "linear": None, "ring": None, "star": None,
}
OPT_LEVELS = [0, 1, 2, 3]

try:
    import psutil  # optional RAM cross-check
except ImportError:
    psutil = None


def load_manifest() -> list[dict]:
    path = os.path.join(FIXTURE_DIR, "manifest.json")
    if not os.path.exists(path):
        sys.exit(f"manifest not found: {path}\nRun generate_fixtures.py first.")
    with open(path) as fh:
        return json.load(fh)


def fits(topology: str, num_qubits: int) -> bool:
    cap = TOPOLOGY_QUBITS[topology]
    return cap is None or num_qubits <= cap


def _run_once(runner: str, fx: dict, topology: str, opt: int,
              cleanup: bool, ablation: bool, no_verify: bool = False) -> dict:
    """Invoke the runner once for a single matrix cell; return the parsed record.

    The runner's getrusage peak RSS is authoritative (kernel-tracked, no
    sampling races). When psutil is available we additionally poll the child's
    RSS as an independent cross-check, stored under `peak_rss_psutil_bytes`.
    """
    cmd = [
        runner,
        "--input", os.path.join(HERE, fx["file"]),
        "--fixture", fx["name"],
        "--family", fx["family"],
        "--nominal-n", str(fx["nominal_n"]),
        "--topology", topology,
        "--opt", str(opt),
    ]
    if not cleanup:
        cmd.append("--no-cleanup")
    if ablation:
        cmd.append("--ablation")
    if no_verify:
        cmd.append("--no-verify")

    psutil_peak = None
    if psutil is not None:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            p = psutil.Process(proc.pid)
            peak = 0
            while proc.poll() is None:
                try:
                    rss = p.memory_info().rss
                    peak = max(peak, rss)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
                time.sleep(0.002)
            psutil_peak = peak or None
        except psutil.NoSuchProcess:
            pass
        stdout, stderr = proc.communicate()
        rc = proc.returncode
    else:
        cp = subprocess.run(cmd, capture_output=True, text=True)
        stdout, stderr, rc = cp.stdout, cp.stderr, cp.returncode

    if rc != 0 and not stdout.strip():
        return {
            "fixture": fx["name"], "family": fx["family"], "nominal_n": fx["nominal_n"],
            "topology": topology, "opt_level": opt, "cleanup_enabled": cleanup,
            "ablation": ablation, "status": "runner_failed",
            "error": stderr.strip()[:500],
        }
    try:
        rec = json.loads(stdout.strip().splitlines()[-1])
    except (json.JSONDecodeError, IndexError):
        return {
            "fixture": fx["name"], "topology": topology, "opt_level": opt,
            "status": "bad_output", "error": (stdout + stderr).strip()[:500],
        }
    if psutil_peak is not None:
        rec["peak_rss_psutil_bytes"] = psutil_peak
    return rec


def run_one(runner, fx, topology, opt, cleanup, ablation, reps=1, no_verify=False):
    """Run a cell `reps` times and return a record with median metrics.

    Q-Rust's SABRE layout/routing is randomized (no fixed seed), so a single
    measurement of CX/depth/time is noisy (~±2-4%). For a defensible comparison
    we repeat each cell and report the median, plus min/max for the headline CX
    count and total time so variance is visible.
    """
    runs = [_run_once(runner, fx, topology, opt, cleanup, ablation, no_verify) for _ in range(reps)]
    ok = [r for r in runs if r.get("status") == "ok"]
    if not ok:
        return runs[0]
    base = dict(ok[0])  # keep verification, layout, ids from a representative run
    if reps > 1:
        med = statistics.median
        cxs = [r["post"]["cx"] for r in ok]
        for key in ("cx", "depth", "gates", "two_q", "swaps"):
            base["post"][key] = med([r["post"][key] for r in ok])
        if "pre_cleanup" in base:
            for key in ("cx", "depth", "gates"):
                base["pre_cleanup"][key] = med([r["pre_cleanup"][key] for r in ok])
        for key in ("parse", "logic_opt", "routing", "basis_decomp", "cleanup", "total"):
            base["timing_us"][key] = med([r["timing_us"][key] for r in ok])
        base["post"]["cx_min"] = min(cxs)
        base["post"]["cx_max"] = max(cxs)
        base["reps"] = len(ok)
    return base


def build_matrix(manifest, dry_run, topologies, opt_levels, ablation=True):
    """Return a list of (fixture, topology, opt, cleanup, ablation) cells."""
    cells = []
    if dry_run:
        fx = next(f for f in manifest if f["name"] == "ghz_3")
        for opt in opt_levels:
            cells.append((fx, "ibm_quito", opt, True, False))
        return cells

    for fx in manifest:
        nq = fx["num_qubits"]
        for topology in topologies:
            if not fits(topology, nq):
                continue
            for opt in opt_levels:
                cells.append((fx, topology, opt, True, False))
        # Cleanup ablation: QFT at O3, cleanup bypassed (paired with the
        # standard O3 cleanup-enabled run above).
        if ablation and fx["family"] == "qft" and 3 in opt_levels:
            for topology in topologies:
                if fits(topology, nq):
                    cells.append((fx, topology, 3, False, True))
    return cells


def main() -> None:
    ap = argparse.ArgumentParser(description="Q-Rust benchmark orchestrator")
    ap.add_argument("--dry-run", action="store_true",
                    help="ghz_3 only, all opt levels, ibm_quito; pretty-print payloads")
    ap.add_argument("--runner", default=DEFAULT_RUNNER, help="path to the release runner binary")
    ap.add_argument("--output-dir", default=RESULTS_DIR)
    ap.add_argument("--opts", default="0,1,2,3",
                    help="comma-separated optimization levels to sweep")
    ap.add_argument("--topologies", default=",".join(TOPOLOGIES),
                    help="comma-separated topology names")
    ap.add_argument("--out-prefix", default="",
                    help="filename prefix for outputs (avoids clobbering prior runs)")
    ap.add_argument("--reps", type=int, default=1,
                    help="repetitions per cell; report median (qrust routing is randomized)")
    ap.add_argument("--families", default="",
                    help="comma-separated family filter (default: all)")
    ap.add_argument("--no-verify", action="store_true",
                    help="skip equivalence verification (faster; for CX/depth-only sweeps)")
    args = ap.parse_args()
    opt_levels = [int(x) for x in args.opts.split(",")]
    topologies = [t.strip() for t in args.topologies.split(",")]
    fam_filter = {f.strip() for f in args.families.split(",") if f.strip()}

    if not os.path.exists(args.runner):
        sys.exit(f"runner not found: {args.runner}\nBuild it: cargo build --release -p qrust-bench")

    manifest = load_manifest()
    if fam_filter:
        manifest = [f for f in manifest if f["family"] in fam_filter]
    cells = build_matrix(manifest, args.dry_run, topologies, opt_levels)
    print(f"[bench] {len(cells)} runs queued "
          f"(psutil cross-check: {'on' if psutil else 'off'})", file=sys.stderr)

    records = []
    skipped = 0
    t0 = time.time()
    for idx, (fx, topology, opt, cleanup, ablation) in enumerate(cells, 1):
        rec = run_one(args.runner, fx, topology, opt, cleanup, ablation,
                      reps=args.reps, no_verify=args.no_verify)
        records.append(rec)
        status = rec.get("status", "?")
        if status == "skipped_too_wide":
            skipped += 1
        tag = " [ablation]" if ablation else ""
        verdict = rec.get("verification", {}).get("verdict_class", "-")
        print(f"[{idx:>3}/{len(cells)}] {fx['name']:<12} {topology:<12} O{opt}{tag} "
              f"-> {status:<16} {verdict}", file=sys.stderr)
        if args.dry_run:
            print(json.dumps(rec, indent=2))

    if args.dry_run:
        print(f"\n[bench] dry run complete ({len(records)} payloads above).", file=sys.stderr)
        return

    os.makedirs(args.output_dir, exist_ok=True)
    raw_path = os.path.join(args.output_dir, f"{args.out_prefix}raw_benchmark_data.json")
    with open(raw_path, "w") as fh:
        json.dump(records, fh, indent=2)

    csv_path = os.path.join(args.output_dir, f"{args.out_prefix}summary.csv")
    write_csv(records, csv_path)

    dt = time.time() - t0
    ok = sum(1 for r in records if r.get("status") == "ok")
    print(f"\n[bench] done in {dt:.1f}s: {ok} ok, {skipped} skipped, "
          f"{len(records) - ok - skipped} other", file=sys.stderr)
    print(f"[bench] raw  -> {raw_path}", file=sys.stderr)
    print(f"[bench] csv  -> {csv_path}", file=sys.stderr)


CSV_COLUMNS = [
    "fixture", "family", "nominal_n", "num_qubits_logical", "topology", "opt_level",
    "cleanup_enabled", "ablation", "status",
    "pre_gates", "pre_cx", "pre_two_q", "pre_depth",
    "post_gates", "post_cx", "post_two_q", "post_swaps", "post_depth", "post_qubits",
    "cx_reduction_pct", "two_q_reduction_pct", "depth_reduction_pct",
    "parse_us", "logic_opt_us", "routing_us", "basis_decomp_us", "cleanup_us", "total_us",
    "peak_rss_bytes", "peak_rss_psutil_bytes",
    "verdict_class", "fidelity", "verified", "verify_method",
]


def write_csv(records: list[dict], path: str) -> None:
    import csv

    def g(d, *keys, default=""):
        for k in keys:
            if d is None:
                return default
            d = d.get(k)
        return default if d is None else d

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for r in records:
            pre, post = r.get("pre", {}), r.get("post", {})
            tim, ver = r.get("timing_us", {}), r.get("verification", {})
            w.writerow([
                g(r, "fixture"), g(r, "family"), g(r, "nominal_n"), g(r, "num_qubits_logical"),
                g(r, "topology"), g(r, "opt_level"), g(r, "cleanup_enabled"), g(r, "ablation"),
                g(r, "status"),
                g(pre, "gates"), g(pre, "cx"), g(pre, "two_q"), g(pre, "depth"),
                g(post, "gates"), g(post, "cx"), g(post, "two_q"), g(post, "swaps"),
                g(post, "depth"), g(post, "num_qubits"),
                g(r, "cx_reduction_pct"), g(r, "two_q_reduction_pct"), g(r, "depth_reduction_pct"),
                g(tim, "parse"), g(tim, "logic_opt"), g(tim, "routing"),
                g(tim, "basis_decomp"), g(tim, "cleanup"), g(tim, "total"),
                g(r, "peak_rss_bytes"), g(r, "peak_rss_psutil_bytes"),
                g(ver, "verdict_class"), g(ver, "fidelity"), g(ver, "verified"), g(ver, "method"),
            ])


if __name__ == "__main__":
    main()
