#!/usr/bin/env python3
"""Interleave saved Q-Rust binaries or request-level ablations on the fixed corpus.

The main four-compiler study remains unchanged. This separate experiment isolates
Q-Rust revision/setting effects. Timing and correctness validation are separate.
Each variant is a JSON object with a binary path and optional request overrides.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

import comparison_study as study

SEED = 20261016
REQUEST_FIELDS = ("qasm", "n_logical", "physical_width", "edges", "topology")


def load_variants(path):
    variants = json.loads(path.read_text())
    if not isinstance(variants, dict) or len(variants) < 2:
        raise ValueError("at least two named variants are required")
    for label, variant in variants.items():
        if not label or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789_-" for c in label):
            raise ValueError("variant labels must be lowercase path-safe identifiers")
        if not isinstance(variant, dict):
            raise ValueError("each variant must describe a binary and optional request controls")
        binary = Path(variant["binary"]).resolve(strict=True)
        overrides = variant.get("request", {})
        if not isinstance(overrides, dict) or set(overrides) & {*REQUEST_FIELDS, "seed"}:
            raise ValueError("variants must not alter circuit, target, layout contract or paired trial seed")
        variant.update(binary=str(binary), binary_sha256=study.digest(binary.read_bytes()),
                       request=overrides)
    return variants


def timing(args):
    variants = load_variants(args.variants)
    cases = study.build_cases()
    if args.case:
        selected = set(args.case)
        known = {c["id"] for c in cases}
        if selected - known:
            raise ValueError(f"unknown cases: {sorted(selected-known)}")
        cases = [c for c in cases if c["id"] in selected]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    protocol = {
        "schema": 1, "experiment": "interleaved_qrust_revision_ablation",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "master_seed": SEED, "trials": args.trials, "warmups_per_cell_tool": 1,
        "tools": list(variants), "variants": variants, "cases": cases,
        "timeout_seconds": args.timeout, "thread_environment": study.THREAD_ENV,
        "python": sys.version, "platform": platform.platform(),
        "machine_note": args.machine_note,
        "harness_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=study.HERE, text=True).strip(),
        "script_sha256": study.digest(Path(__file__).read_bytes()),
        "timing_boundary": "Same warm public-API boundary as comparison_study.py; variants interleaved sequentially; validation runs separately.",
        "limitations": "Exploratory changes motivated by this internal corpus; no external holdout. Repetitions are not independent circuits. Stage instrumentation, if requested, has overhead and is diagnostic only.",
    }
    (out / "protocol.json").write_text(json.dumps(protocol, indent=2)+"\n")
    rng = random.Random(SEED)
    rng.shuffle(cases)
    workers = {}
    started = time.monotonic()
    try:
        for label, variant in variants.items():
            logdir = out / label
            logdir.mkdir()
            workers[label] = study.Worker("qrust", Path(variant["binary"]), logdir)
        with (out/"timed.jsonl").open("w") as measured, (out/"warmup.jsonl").open("w") as warmups:
            for ci, case in enumerate(cases):
                for trial in range(-1, args.trials):
                    order = list(variants)
                    rng.shuffle(order)
                    for position, label in enumerate(order):
                        request = {key: case[key] for key in REQUEST_FIELDS}
                        request.update(seed=SEED+max(trial, 0))
                        request.update(variants[label]["request"])
                        worker = workers[label]
                        if trial >= 0 and worker.process is None:
                            recovery = worker.call(request, args.timeout)
                            warmups.write(json.dumps({"case_id": case["id"], "tool": label,
                                "reason": "worker_restart", **recovery})+"\n")
                            warmups.flush()
                        result = worker.call(request, args.timeout)
                        row = {"id": f"{case['id']}__{label}__{trial}", "case_id": case["id"],
                               "tool": label, "trial": trial, "seed": request["seed"], "case_order": ci,
                               "tool_order": position, **result}
                        if "out_qasm" in row:
                            row["output_sha256"] = study.digest(row["out_qasm"])
                        stream = warmups if trial == -1 else measured
                        stream.write(json.dumps(row, separators=(",", ":"))+"\n")
                        stream.flush()
                print(json.dumps({"completed_cases": ci+1, "total_cases": len(cases),
                    "elapsed_seconds": round(time.monotonic()-started, 1)}), flush=True)
    finally:
        for worker in workers.values():
            worker.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    timer = commands.add_parser("time")
    timer.add_argument("output", type=Path)
    timer.add_argument("--variants", type=Path, required=True)
    timer.add_argument("--case", action="append")
    timer.add_argument("--trials", type=int, default=7)
    timer.add_argument("--timeout", type=float, default=60)
    timer.add_argument("--machine-note", required=True)
    validator = commands.add_parser("validate")
    validator.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "time":
        if args.trials < 1:
            parser.error("trials must be positive")
        timing(args)
    else:
        study.validate(args)


if __name__ == "__main__":
    main()
