#!/usr/bin/env python3
"""Revalidate a frozen comparison corpus after a numerical correctness repair.

This is not a timing campaign. Saved baseline outputs are checked again without
recompilation. Q-Rust is recompiled once per frozen case with the original request
and its output is checked with the unchanged independent oracle. No previous
evidence is edited, and output directories must be new.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from comparison_study import Worker, THREAD_ENV, build_cases
from comparison_validation import validate_output
import comparison_validation

HERE = Path(__file__).resolve().parent
METRICS = ("cx_count", "depth", "two_q_depth", "gate_count")


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def source_provenance(root):
    """Hash actual compiler inputs, including files not tracked by Git.

    Source archives are supported; an enclosing unrelated Git checkout must not
    be mistaken for the compiler repository. Git adds context but is not needed
    to identify source contents. Build/cache directories are excluded.
    """
    root = root.resolve()
    hashes = {}
    for directory, names, files in os.walk(root):
        names[:] = sorted(name for name in names if name not in {
            ".git", "target", ".venv", "venv", "__pycache__"})
        for name in sorted(files):
            if name.endswith(".rs") or name in {"Cargo.toml", "Cargo.lock"}:
                path = Path(directory) / name
                hashes[str(path.relative_to(root))] = digest(path.read_bytes())
    result = {"compiler_source_sha256": hashes,
              "compiler_source_manifest_sha256": digest(json.dumps(hashes, sort_keys=True)),
              "git_available": False, "git_head": None, "git_diff_sha256": None}
    try:
        git_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=root,
                                           stderr=subprocess.DEVNULL, text=True).strip()
        if Path(git_root).resolve() != root:
            result["git_unavailable_reason"] = "Source export is not the root of its enclosing Git repository"
            return result
        head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                       stderr=subprocess.DEVNULL, text=True).strip()
        diff = subprocess.check_output(["git", "diff", "HEAD"], cwd=root, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError) as error:
        result["git_unavailable_reason"] = type(error).__name__
    else:
        result.update(git_available=True, git_head=head, git_diff_sha256=digest(diff))
    return result


def artifact_key(case, row):
    """Include both layouts: equal text does not imply equal logical semantics."""
    return digest(json.dumps([case["source_sha256"], case["physical_width"], case["edges"],
                              row["output_sha256"], row["initial_layout"], row["final_layout"]]))


def check_frozen_inputs(protocol, rows, oracle_sha):
    expected_oracle = protocol["source_sha256"]["experiments/comparison_validation.py"]
    if oracle_sha != expected_oracle:
        raise ValueError("Independent oracle differs from frozen study; refusing tolerance/semantic drift")
    cases = {case["id"]: case for case in protocol["cases"]}
    if len(cases) != len(protocol["cases"]):
        raise ValueError("Duplicate case IDs")
    for case in cases.values():
        if digest(case["source_qasm"]) != case["source_sha256"] or digest(case["qasm"]) != case["normalized_sha256"]:
            raise ValueError(f"Fixture checksum mismatch: {case['id']}")
    expected = {(case_id, tool, trial) for case_id in cases for tool in protocol["tools"]
                for trial in range(protocol["trials"])}
    seen = set()
    for row in rows:
        cell = (row["case_id"], row["tool"], row["trial"])
        if cell in seen or cell not in expected:
            raise ValueError(f"Duplicate/unexpected trial: {cell}")
        seen.add(cell)
        if row["status"] == "compiled" and digest(row["out_qasm"]) != row["output_sha256"]:
            raise ValueError(f"Output checksum mismatch: {row['id']}")
    if seen != expected:
        raise ValueError(f"Missing {len(expected - seen)} frozen trial(s)")
    return cases


def validate_row(case, row, cache):
    if row["status"] != "compiled":
        return {"passed": False, "evidence_class": "not_compiled", "error": row.get("error", row["status"])}, None
    key = artifact_key(case, row)
    if key not in cache:
        edges = None if case["topology"] == "all_to_all" else [
            pair for a, b in case["edges"] for pair in ((a, b), (b, a))]
        cache[key] = validate_output(case["source_qasm"], row["out_qasm"], case["n_logical"],
                                     case["physical_width"], edges, row["initial_layout"], row["final_layout"])
    return cache[key], key


def summarize(rows, cache):
    groups = {}
    for tool in sorted({row["tool"] for row in rows}):
        selected = [r for r in rows if r["tool"] == tool]
        groups[tool] = {
            "rows": len(selected), "passed": sum(r["validation"]["passed"] is True for r in selected),
            "failed": sum(r["validation"]["passed"] is False for r in selected),
            "unchecked": sum(r["validation"]["passed"] is None for r in selected),
            "unique_artifacts": len({r["validation_cache_key"] for r in selected if r["validation_cache_key"]}),
            "max_residual": max((r["validation"].get("max_input_residual", 0) for r in selected), default=0),
        }
    return {"by_tool": groups, "unique_artifacts": len(cache),
            "unique_evidence_classes": dict(Counter(v["evidence_class"] for v in cache.values())),
            "failed_or_unchecked_ids": [r["id"] for r in rows if r["validation"]["passed"] is not True],
            "all_passed": all(r["validation"]["passed"] is True for r in rows)}


def provenance(previous, protocol):
    oracle_path = Path(comparison_validation.__file__)
    return {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "oracle_sha256": digest(oracle_path.read_bytes()),
        "prior_files_sha256": {name: digest((previous / name).read_bytes())
                               for name in ("protocol.json", "timed.jsonl", "validated.jsonl")},
        "frozen_compiler_commit": protocol["compiler_commit"],
        "frozen_binary_sha256": protocol["qrust_binary_sha256"],
        "python": sys.version, "platform": platform.platform(),
        "thread_environment": THREAD_ENV,
        "numpy_version": __import__("numpy").__version__,
        "scope": "Correctness and output quality only. No speed measurement or claim.",
        "validation_contract": protocol["validation"],
        "limitations": "Full numerical isometry through six physical qubits; otherwise the unchanged finite eight-input sample. Not a proof. Same development corpus, not an external holdout.",
    }


def load_previous(previous):
    protocol = json.loads((previous / "protocol.json").read_text())
    rows = read_rows(previous / "timed.jsonl")
    meta = provenance(previous, protocol)
    cases = check_frozen_inputs(protocol, rows, meta["oracle_sha256"])
    # The live builder uses tuples for edges; the frozen JSON uses arrays.
    current = {c["id"]: c for c in json.loads(json.dumps(build_cases()))}
    if current != cases:
        raise ValueError("Current corrected fixtures/normalization differ from frozen case corpus")
    old_validated = {r["id"]: r for r in read_rows(previous / "validated.jsonl")}
    if len(old_validated) != len(rows) or set(old_validated) != {r["id"] for r in rows}:
        raise ValueError("Saved timed/validated row grid differs")
    for row in rows:
        prior = old_validated[row["id"]]
        if any(prior.get(k) != row.get(k) for k in row):
            raise ValueError(f"Saved validated row differs from timed artifact: {row['id']}")
    return protocol, rows, cases, meta, old_validated


def revalidate_saved(args):
    protocol, rows, cases, meta, old_validated = load_previous(args.previous)
    args.output.mkdir(parents=True, exist_ok=False)
    dump(args.output / "protocol.json", protocol)
    dump(args.output / "provenance.json", meta)
    cache, checked, changes = {}, [], []
    with (args.output / "revalidated.jsonl").open("w") as stream:
        for i, row in enumerate(rows):
            validation, key = validate_row(cases[row["case_id"]], row, cache)
            record = {"id": row["id"], "case_id": row["case_id"], "tool": row["tool"],
                      "trial": row["trial"], "output_sha256": row.get("output_sha256"),
                      "initial_layout": row.get("initial_layout"), "final_layout": row.get("final_layout"),
                      "validation": validation, "validation_cache_key": key}
            checked.append(record)
            previous = old_validated[row["id"]]["validation"]
            if any(previous.get(k) != validation.get(k) for k in ("passed", "evidence_class", "metrics")):
                changes.append({"id": row["id"], "prior_validation": previous, "fresh_validation": validation})
            stream.write(json.dumps(record, separators=(",", ":")) + "\n")
            if (i + 1) % 100 == 0:
                stream.flush()
                print(json.dumps({"saved_rows_checked": i + 1, "unique_artifacts": len(cache)}), flush=True)
    summary = {**summarize(checked, cache), "changed_verdict_or_metrics": changes, "provenance": meta}
    dump(args.output / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if summary["all_passed"] and not changes else 1


def recompile_corrected(args):
    protocol, rows, cases, meta, _ = load_previous(args.previous)
    old_qrust = defaultdict(list)
    for row in rows:
        if row["tool"] == "qrust":
            old_qrust[row["case_id"]].append(row)
    for case_id, previous in old_qrust.items():
        if any(r["status"] != "compiled" for r in previous) or len({artifact_key(cases[case_id], r) for r in previous}) != 1:
            raise ValueError(f"Frozen Q-Rust output not uniquely determined: {case_id}")
    args.output.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    meta.update(corrected_binary_sha256=digest(binary.read_bytes()),
                worker_adapter_sha256=digest((HERE / "src/bin/comparison_runner.rs").read_bytes()),
                **source_provenance(HERE.parent),
                requests_per_case=1, timing_telemetry="Renamed diagnostic fields; not a timing trial and excluded from comparisons")
    dump(args.output / "protocol.json", protocol)
    dump(args.output / "provenance.json", meta)
    worker = Worker("qrust", binary, args.output)
    results = []
    try:
        with (args.output / "compiled.jsonl").open("w") as stream:
            for i, case in enumerate(sorted(cases.values(), key=lambda c: c["id"])):
                request = {key: case[key] for key in ("qasm", "n_logical", "physical_width", "edges", "topology")}
                request["seed"] = protocol["master_seed"]
                result = worker.call(request, args.timeout)
                for field in ("compile_ns", "parse_ns"):
                    if field in result:
                        result["diagnostic_" + field] = result.pop(field)
                row = {"id": case["id"] + "__qrust_corrected", "case_id": case["id"], "tool": "qrust", **result}
                if "out_qasm" in row:
                    row["output_sha256"] = digest(row["out_qasm"])
                results.append(row)
                stream.write(json.dumps(row, separators=(",", ":")) + "\n")
                stream.flush()
                print(json.dumps({"cases_compiled": i + 1, "total_cases": len(cases), "case": case["id"], "status": row["status"]}), flush=True)
    finally:
        worker.close()
    cache, checked, changes = {}, [], []
    with (args.output / "validated.jsonl").open("w") as stream:
        for i, row in enumerate(results):
            case = cases[row["case_id"]]
            row["validation"], row["validation_cache_key"] = validate_row(case, row, cache)
            previous = old_qrust[case["id"]][0]
            old_validation, _ = validate_row(case, previous, {})
            old_metrics = old_validation.get("metrics", {})
            new_metrics = row["validation"].get("metrics", {})
            change = {"case_id": case["id"], "family": case["family"], "topology": case["topology"],
                      "before_passed": old_validation["passed"], "after_passed": row["validation"]["passed"],
                      "before": old_metrics, "after": new_metrics,
                      "delta": {k: new_metrics[k] - old_metrics[k] for k in METRICS if k in old_metrics and k in new_metrics},
                      "same_output_text": row.get("output_sha256") == previous["output_sha256"],
                      "same_layouts": all(row.get(k) == previous[k] for k in ("initial_layout", "final_layout"))}
            changes.append(change)
            checked.append(row)
            stream.write(json.dumps(row, separators=(",", ":")) + "\n")
            stream.flush()
            print(json.dumps({"cases_validated": i + 1, "case": case["id"], "passed": row["validation"]["passed"], "delta": change["delta"]}), flush=True)
    meta["corrected_binary_sha256_after"] = digest(binary.read_bytes())
    if meta["corrected_binary_sha256_after"] != meta["corrected_binary_sha256"]:
        raise RuntimeError("Worker binary changed during campaign")
    summary = {**summarize(checked, cache), "provenance": meta, "case_comparisons": changes,
               "quality_counts": {k: dict(Counter("improved" if c["delta"][k] < 0 else "worse" if c["delta"][k] > 0 else "equal"
                                                    for c in changes if k in c["delta"])) for k in METRICS},
               "same_output_text_cases": sum(c["same_output_text"] for c in changes),
               "same_layout_cases": sum(c["same_layouts"] for c in changes)}
    dump(args.output / "summary.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k != "case_comparisons"}, indent=2), flush=True)
    return 0 if summary["all_passed"] and all(c["before_passed"] is True for c in changes) else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("saved", "corrected"):
        command = sub.add_parser(name)
        command.add_argument("--previous", type=Path, required=True, help="Frozen directory containing protocol/timed/validated")
        command.add_argument("--output", type=Path, required=True, help="New evidence directory; must not exist")
        if name == "corrected":
            command.add_argument("--binary", type=Path, required=True)
            command.add_argument("--timeout", type=float, default=60)
    args = parser.parse_args()
    return revalidate_saved(args) if args.command == "saved" else recompile_corrected(args)


if __name__ == "__main__":
    raise SystemExit(main())
