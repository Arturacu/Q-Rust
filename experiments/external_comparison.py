#!/usr/bin/env python3
"""Frozen external-corpus compiler comparison, with timing and checking separated.

``prepare`` admits and freezes a previously selected community corpus; ``time``
runs one fresh persistent-worker session; ``validate`` checks retained outputs
only after every prescribed session has completed. Existing compiler adapters
and the independent numerical oracle are reused without changing their recipes.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

from comparison_study import Worker, TOOLS as PRIMARY_TOOLS, THREAD_ENV

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PACKAGES = ("qiskit", "cirq-core", "pytket", "numpy", "scipy", "networkx", "rustworkx", "sympy", "ply", "lark")
TOOLS = (*PRIMARY_TOOLS, "cirq_one_pass", "tket_light")
PROFILES = {"primary": list(PRIMARY_TOOLS),
            "light": ["qrust", "qiskit", "cirq_one_pass", "tket_light"]}
CONTRACT_FILES = (
    "experiments/external_comparison.py", "experiments/comparison_study.py",
    "experiments/comparison_qiskit_tket.py", "experiments/comparison_cirq.py",
    "experiments/comparison_validation.py", "experiments/bench_baselines.py",
    "experiments/analyze_comparison_study.py",
    "experiments/diagnose_compiler_tradeoffs.py",
    "experiments/external_corpus.py", "experiments/analyze_external_comparison.py",
    "experiments/requirements-comparison.txt",
)


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


class ProfileWorker(Worker):
    """Change dispatch for two prespecified light recipes, never their passes."""
    def __init__(self, tool, binary, outdir):
        self.recipe = tool
        canonical_tool = {"cirq_one_pass": "cirq", "tket_light": "tket"}.get(tool, tool)
        try:
            super().__init__(canonical_tool, binary, outdir)
        except Exception:
            self.close()
            raise

    def close(self):
        super().close()
        # Popen itself can fail after the log was opened but before process
        # assignment. Such a startup failure still must release the handle.
        if hasattr(self, "log") and not self.log.closed:
            self.log.close()

    def start(self):
        if self.recipe in PRIMARY_TOOLS:
            try:
                return super().start()
            except Exception:
                self.close()
                raise
        env = dict(os.environ, **THREAD_ENV)
        env.setdefault("MPLCONFIGDIR", str(self.outdir / "matplotlib_cache"))
        self.log = (self.outdir / f"worker_{self.recipe}.log").open("a")
        try:
            self.process = subprocess.Popen(
                [sys.executable, str(Path(__file__).resolve()), "worker", self.tool],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log,
                text=True, bufsize=1, env=env)
            if not self.receive(120).get("ready"):
                raise RuntimeError(f"{self.tool} light worker did not become ready")
        except Exception:
            self.close()
            raise


def light_worker(tool):
    from diagnose_compiler_tradeoffs import cirq_trial, tket_trial
    if tool == "cirq":
        import cirq  # noqa: F401: load SDK outside the first worker request
    else:
        import pytket  # noqa: F401
    recipe = "cirq_one_pass" if tool == "cirq" else "tket_light"
    print(json.dumps({"ready": True, "tool": tool, "recipe": recipe}), flush=True)
    for line in sys.stdin:
        try:
            request = json.loads(line)
            result = (cirq_trial(request, recipe) if tool == "cirq"
                      else tket_trial(request, recipe, stage_probe=False))
            result.update(status="compiled", seed_used=False,
                          pipeline={"recipe": recipe, "source": "diagnose_compiler_tradeoffs.py",
                                    "stage_probe": False, "instrumentation_in_compile_ns": True})
        except Exception as error:
            result = {"status": "error", "error": f"{type(error).__name__}: {error}"}
        print(json.dumps(result), flush=True)


def read_rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def source_provenance(root=ROOT):
    """Hash build inputs without importing the simulator or compiler SDKs."""
    root = root.resolve()
    paths = []
    for directory, names, files in os.walk(root):
        names[:] = sorted(name for name in names if name not in {
            ".git", "target", ".venv", "venv", "__pycache__"})
        paths.extend(Path(directory) / name for name in sorted(files)
                     if name.endswith(".rs") or name in {"Cargo.toml", "Cargo.lock"})
    hashes = {str(path.relative_to(root)): digest(path.read_bytes()) for path in paths}
    result = {"compiler_source_sha256": hashes,
              "compiler_source_manifest_sha256": digest(canonical(hashes)),
              "git_head": None, "git_diff_sha256": None}
    try:
        git_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"],
                                           cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
        if Path(git_root).resolve() == root:
            result["git_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"],
                cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
            result["git_diff_sha256"] = digest(subprocess.check_output(["git", "diff", "HEAD"],
                cwd=root, stderr=subprocess.DEVNULL))
    except (OSError, subprocess.CalledProcessError):
        pass
    return result


def artifact_hash(case, row):
    """Bind the entire physical output and both logical wire maps."""
    return digest(canonical({"source": case["source_sha256"],
        "normalized": case["normalized_sha256"], "physical_width": case["physical_width"],
        "topology": case["topology"], "edges": case["edges"],
        "output": row["output_sha256"], "initial_layout": row["initial_layout"],
        "final_layout": row["final_layout"], "output_qubits": row["output_qubits"]}))


def build_cases(manifest):
    """Add a target grid without selecting on any compiler outcome."""
    fixtures = manifest["cases"]
    if not fixtures or len({case["id"] for case in fixtures}) != len(fixtures):
        raise ValueError("Corpus must contain unique nonempty fixture IDs")
    cases = []
    for fixture in fixtures:
        width = fixture["n_logical"]
        if type(width) is not int or not 1 <= width <= 16:
            raise ValueError("Corpus width must fit the unchanged oracle (1..16)")
        for field, hashed in (("source_qasm", "source_sha256"), ("qasm", "normalized_sha256")):
            if digest(fixture[field]) != fixture[hashed]:
                raise ValueError(f"Corpus {field} checksum mismatch: {fixture['id']}")
        evidence = fixture.get("normalization_evidence", {})
        mode = "full_unitary_numerical" if width <= 6 else "sampled_8_inputs_numerical"
        if (evidence.get("passed") is not True or evidence.get("mode") != mode
                or evidence.get("input_count") != (2 ** width if width <= 6 else 8)
                or evidence.get("tolerance") != 1e-7
                or any(type(evidence.get("residuals", {}).get(key)) not in (int, float)
                       or not 0 <= evidence["residuals"][key] <= 1e-7
                       for key in ("source", "normalized"))):
            raise ValueError(f"Missing or invalid pre-timing normalization evidence: {fixture['id']}")
        for topology in ("all_to_all", "linear"):
            cases.append({**fixture, "id": f"{fixture['id']}__{topology}",
                          "fixture": fixture.get("fixture", fixture["id"]),
                          "family": fixture["family"], "topology": topology,
                          "physical_width": width,
                          # The native worker consumes explicit connectivity
                          # even when Python adapters interpret the topology
                          # name as unconstrained. An empty list disconnects it.
                          "edges": ([[a, b] for a in range(width) for b in range(a + 1, width)]
                                    if topology == "all_to_all" else [[i, i + 1] for i in range(width - 1)])})
    return sorted(cases, key=lambda case: case["id"])


def schedule(cases, trials, master_seed, session):
    """Same trial seeds across sessions, independently shuffled request order."""
    rng = random.Random(master_seed + 1_000_003 * session)
    ordered = list(cases)
    rng.shuffle(ordered)
    rows = []
    for case_order, case in enumerate(ordered):
        trials_in_order = list(range(trials))
        rng.shuffle(trials_in_order)
        for trial in [-1, *trials_in_order]:
            tools = list(TOOLS)
            rng.shuffle(tools)
            for tool_order, tool in enumerate(tools):
                rows.append({"case_id": case["id"], "tool": tool, "trial": trial,
                             "seed": master_seed + max(trial, 0),
                             "case_order": case_order, "tool_order": tool_order})
    return rows


def prepare(args):
    if args.trials < 1 or args.sessions < 2 or not math.isfinite(args.timeout) or args.timeout <= 0:
        raise ValueError("Require positive trials/timeout and at least two sessions")
    manifest_bytes = args.corpus.read_bytes()
    manifest = json.loads(manifest_bytes)
    cases = build_cases(manifest)
    sources = {name: digest((ROOT / name).read_bytes()) for name in CONTRACT_FILES}
    protocol = {
        "schema": 1, "created_utc": utc_now(), "study": "external_community_corpus",
        "recipe_profiles": PROFILES,
        "recipes": {"qrust": "public O3 comparison_runner", "qiskit": "transpile optimization_level=3",
                    "cirq": "CZ/RouteCQC/CZ/CX, max two iterations",
                    "tket": "FullPeepholeOptimise/DefaultMappingPass/cleanup/U3,CX",
                    "cirq_one_pass": "same Cirq structure, max one synthesis iteration per phase",
                    "tket_light": "SynthesiseTket/DefaultMappingPass/cleanup/U3,CX"},
        "profile_caveat": "All six recipes are interleaved within the same sessions. Primary and light profile analyses reuse the SAME Q-Rust and Qiskit observations. The light recipes reuse existing diagnostic functions, including tiny stage-clock/dictionary bookkeeping inside compile_ns, with no timing subtraction. Neither profile claims universal SDK best performance or equal search effort.",
        "master_seed": args.seed, "trials": args.trials, "sessions": args.sessions,
        "warmups_per_cell_tool": 1, "timeout_seconds": args.timeout, "tools": list(TOOLS),
        "cases": cases, "corpus_manifest_sha256": digest(manifest_bytes),
        "corpus_exclusions": manifest.get("exclusions", []),
        "corpus_provenance": manifest.get("provenance", {}),
        "compiler_provenance": source_provenance(),
        "compiler_commit": args.compiler_commit,
        "qrust_binary_sha256": digest(args.binary.read_bytes()),
        "source_sha256": sources,
        "versions": {name: importlib.metadata.version(name) for name in PACKAGES},
        "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
        "logical_cpu_count": os.cpu_count(), "thread_environment": THREAD_ENV,
        "worker_executable": str(Path(sys.executable).resolve()), "machine_note": args.machine_note,
        "timing_boundary": "Warm in-process compile_ns: cloning, configuration, passes and compiled-circuit retrieval. Parsing is parse_ns. Imports/startup, serialization, IPC and checking are outside compile_ns. request_wall_ns covers the complete worker request/response. noncompile_nonparse_request_ns is wall minus parse minus compile: it combines serialization, layouts, diagnostics, IPC and interpreter overhead, and is NOT isolated serialization time.",
        "seed_policy": "Trial seed = master_seed + trial, reused in both sessions. Qiskit consumes the seed; Q-Rust/Cirq/TKET report deterministic recipes. Seven trials describe seven Qiskit seed outcomes, with two session repetitions per seed; other tools measure timing repetition. Distinct seeds are not independent circuits.",
        "target": "Arbitrary one-qubit U3 and CX; bidirectional ideal all-to-all and line targets; equal physical/logical width; free initial/final logical permutations; arbitrary logical inputs.",
        "normalization": "Community extraction and exact gate-definition expansion precede freezing; all tools receive identical normalized QASM. No baseline compiler optimizes another tool's input. Original unitary source is the independent correctness reference. Raw upstream source and transformation lineage remain in the corpus manifest.",
        "validation": "Unchanged independent NumPy analytic gate oracle, AFTER ALL timing sessions. Full numerical isometry through six physical qubits; three basis, uniform and four seeded Haar inputs otherwise through sixteen qubits. Shared global phase, full physical amplitudes, residual tolerance 1e-7. Sampled evidence is not a proof; CX and depth are not execution fidelity.",
        "analysis": "Analyze each session separately with the same complete-cell criteria and paired ratios. Retain all selected fixtures, warmup outcomes, failures/timeouts and ineligible cells. No replacement by successful subsets. Session repetition is not additional independent circuit coverage. Fixed recipes have unequal optimization/search effort and are not each SDK's universal best setting.",
        "runtime_limit": "Single host with thread environment caps. No CPU pinning/frequency control. No simulation, builds or other benchmark sessions may run concurrently. Runtime source and binary hashes checked before each session.",
        "schedules": {str(session): schedule(cases, args.trials, args.seed, session)
                      for session in range(1, args.sessions + 1)},
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "corpus_manifest.json").write_bytes(manifest_bytes)
    dump(args.output / "protocol.json", protocol)
    (args.output / "protocol.sha256").write_text(digest((args.output / "protocol.json").read_bytes()) + "\n")
    print(json.dumps({"cases": len(cases), "fixtures": len(manifest["cases"]),
        "sessions": args.sessions, "timed_rows_per_session": len(cases) * len(TOOLS) * args.trials,
        "protocol_sha256": (args.output / "protocol.sha256").read_text().strip()}))


def load_protocol(study):
    protocol_bytes = (study / "protocol.json").read_bytes()
    protocol_sha = digest(protocol_bytes)
    if protocol_sha != (study / "protocol.sha256").read_text().strip():
        raise ValueError("Frozen protocol checksum mismatch")
    protocol = json.loads(protocol_bytes)
    if digest((study / "corpus_manifest.json").read_bytes()) != protocol["corpus_manifest_sha256"]:
        raise ValueError("Frozen corpus checksum mismatch")
    manifest = json.loads((study / "corpus_manifest.json").read_bytes())
    if build_cases(manifest) != protocol["cases"]:
        raise ValueError("Frozen cases differ from the frozen corpus")
    for session in range(1, protocol["sessions"] + 1):
        expected = schedule(protocol["cases"], protocol["trials"], protocol["master_seed"], session)
        if expected != protocol["schedules"][str(session)]:
            raise ValueError("Frozen schedule differs from the prespecified seed policy")
    return protocol, protocol_sha


def check_runtime(protocol, binary):
    if digest(binary.read_bytes()) != protocol["qrust_binary_sha256"]:
        raise ValueError("Compiler binary differs from frozen protocol")
    for name, checksum in protocol["source_sha256"].items():
        if digest((ROOT / name).read_bytes()) != checksum:
            raise ValueError(f"Comparison source differs from frozen protocol: {name}")
    if source_provenance()["compiler_source_sha256"] != protocol["compiler_provenance"]["compiler_source_sha256"]:
        raise ValueError("Compiler source differs from frozen protocol")
    versions = {name: importlib.metadata.version(name) for name in PACKAGES}
    if versions != protocol["versions"] or sys.version != protocol["python"]:
        raise ValueError("Python or dependency versions differ from frozen protocol")
    if THREAD_ENV != protocol["thread_environment"]:
        raise ValueError("Worker thread environment differs from frozen protocol")


def timed_request(worker, request, timeout):
    started = time.perf_counter_ns()
    result = worker.call(request, timeout)
    wall = time.perf_counter_ns() - started
    row = {**result, "request_wall_ns": wall}
    if row.get("status") == "compiled":
        components = [row.get(key) for key in ("compile_ns", "parse_ns")]
        if any(type(value) is not int or value < 0 for value in components) or components[0] == 0:
            row.update(status="error", error="Worker returned invalid timing components")
        elif sum(components) > wall:
            row.update(status="error", error="Worker component timing exceeds complete request time")
        else:
            row["noncompile_nonparse_request_ns"] = wall - sum(components)
    return row


def timing(args):
    protocol, protocol_sha = load_protocol(args.study)
    if not 1 <= args.session <= protocol["sessions"]:
        raise ValueError("Session is outside frozen protocol")
    check_runtime(protocol, args.binary)
    # A lock covers the whole study, preventing overlapping sessions/checking.
    lock = args.study / "active.lock"
    descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(descriptor)
    workers = {}
    outdir = args.study / f"session_{args.session:02d}"
    try:
        outdir.mkdir(exist_ok=False)
        (outdir / "protocol.json").write_bytes((args.study / "protocol.json").read_bytes())
        session_info = {"session": args.session, "protocol_sha256": protocol_sha,
                        "started_utc": utc_now(), "finished_utc": None, "completed": False}
        dump(outdir / "session.json", session_info)
        startup_errors = {}
        for tool in TOOLS:
            try:
                workers[tool] = ProfileWorker(tool, args.binary, outdir)
            except Exception as error:
                workers[tool] = None
                startup_errors[tool] = f"{type(error).__name__}: {error}"
        dump(outdir / "startup.json", {"errors": startup_errors})
        cases = {case["id"]: case for case in protocol["cases"]}
        count, completed_cases, started = 0, set(), time.monotonic()
        with (outdir / "timed.jsonl").open("x") as measured, (outdir / "warmup.jsonl").open("x") as warmups:
            for entry in protocol["schedules"][str(args.session)]:
                case, tool = cases[entry["case_id"]], entry["tool"]
                request = {key: case[key] for key in ("qasm", "n_logical", "physical_width", "edges", "topology")}
                request["seed"] = entry["seed"]
                worker = workers[tool]
                if worker is None:
                    result = {"status": "error", "error": startup_errors[tool],
                              "phase": "worker_startup", "request_wall_ns": None}
                else:
                    replacement_failed = False
                    if entry["trial"] >= 0 and worker.process is None:
                        warmed = timed_request(worker, request, protocol["timeout_seconds"])
                        warmups.write(canonical({**entry, "session": args.session,
                            "reason": "worker_restart", **warmed}) + "\n")
                        warmups.flush()
                        if warmed.get("status") != "compiled":
                            result = {**warmed, "phase": "replacement_warmup",
                                      "error": "Replacement warmup failed; measured request not attempted: " + warmed.get("error", "unknown failure")}
                            replacement_failed = True
                    if not replacement_failed:
                        result = timed_request(worker, request, protocol["timeout_seconds"])
                row = {**result, **entry, "session": args.session,
                       "id": f"{case['id']}__{tool}__{entry['trial']}"}
                if row.get("status") == "compiled":
                    try:
                        row["output_sha256"] = digest(row["out_qasm"])
                        row["artifact_sha256"] = artifact_hash(case, row)
                    except (KeyError, TypeError) as error:
                        row.update(status="error", error=f"Incomplete worker output: {error}")
                stream = warmups if entry["trial"] == -1 else measured
                stream.write(canonical(row) + "\n")
                stream.flush()
                if entry["trial"] >= 0:
                    count += 1
                    if count % (protocol["trials"] * len(TOOLS)) == 0:
                        completed_cases.add(case["id"])
                        print(json.dumps({"session": args.session, "completed_cases": len(completed_cases),
                            "total_cases": len(cases), "last_case": case["id"],
                            "elapsed_seconds": round(time.monotonic() - started, 1)}), flush=True)
        session_info.update(finished_utc=utc_now(), completed=True, timed_rows=count,
            timed_sha256=digest((outdir / "timed.jsonl").read_bytes()),
            warmup_sha256=digest((outdir / "warmup.jsonl").read_bytes()))
        dump(outdir / "session.json", session_info)
    finally:
        for worker in workers.values():
            if worker is not None:
                worker.close()
        lock.unlink()


def check_rows(protocol, session, rows):
    expected = {(case["id"], tool, trial) for case in protocol["cases"]
                for tool in TOOLS for trial in range(protocol["trials"])}
    cases = {case["id"]: case for case in protocol["cases"]}
    scheduled = [entry for entry in protocol["schedules"][str(session)] if entry["trial"] >= 0]
    seen = set()
    for index, row in enumerate(rows):
        key = (row["case_id"], row["tool"], row["trial"])
        if key in seen or key not in expected or row["session"] != session:
            raise ValueError("Duplicate, unexpected, or wrong-session timed trial")
        if row["seed"] != protocol["master_seed"] + row["trial"]:
            raise ValueError("Recorded seed differs from frozen trial seed")
        if row["id"] != f"{row['case_id']}__{row['tool']}__{row['trial']}":
            raise ValueError("Timed row ID differs from its trial identity")
        if any(row[key] != scheduled[index][key] for key in scheduled[index]):
            raise ValueError("Timed request order differs from frozen schedule")
        seen.add(key)
        if row["status"] == "compiled":
            if digest(row["out_qasm"]) != row["output_sha256"]:
                raise ValueError("Timed output checksum mismatch")
            if artifact_hash(cases[row["case_id"]], row) != row["artifact_sha256"]:
                raise ValueError("Timed artifact/layout checksum mismatch")
    if seen != expected:
        raise ValueError(f"Missing {len(expected - seen)} prescribed timed trials")


def validate(args):
    protocol, protocol_sha = load_protocol(args.study)
    lock = args.study / "active.lock"
    descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(descriptor)
    try:
        sessions = []
        # Check every session BEFORE importing/running the numerical simulator.
        for session in range(1, protocol["sessions"] + 1):
            directory = args.study / f"session_{session:02d}"
            info = json.loads((directory / "session.json").read_text())
            if not info["completed"] or info["protocol_sha256"] != protocol_sha:
                raise ValueError("All prescribed timing sessions must finish before validation")
            for name in ("timed", "warmup"):
                if digest((directory / f"{name}.jsonl").read_bytes()) != info[f"{name}_sha256"]:
                    raise ValueError(f"Session {name} file checksum mismatch")
            rows = read_rows(directory / "timed.jsonl")
            check_rows(protocol, session, rows)
            if (directory / "validated.jsonl").exists():
                raise FileExistsError(directory / "validated.jsonl")
            sessions.append((directory, rows))
        oracle = HERE / "comparison_validation.py"
        if digest(oracle.read_bytes()) != protocol["source_sha256"]["experiments/comparison_validation.py"]:
            raise ValueError("Independent oracle differs from frozen protocol")
        from comparison_validation import validate_output
        cases = {case["id"]: case for case in protocol["cases"]}
        cache, summary = {}, []
        for directory, rows in sessions:
            counts = Counter()
            with (directory / "validated.jsonl").open("x") as output:
                for index, original in enumerate(rows):
                    row = dict(original)
                    if row["status"] == "compiled":
                        case = cases[row["case_id"]]
                        key = row["artifact_sha256"]
                        if key not in cache:
                            edges = None if case["topology"] == "all_to_all" else [
                                pair for a, b in case["edges"] for pair in ((a, b), (b, a))]
                            cache[key] = validate_output(case["source_qasm"], row["out_qasm"],
                                case["n_logical"], case["physical_width"], edges,
                                row["initial_layout"], row["final_layout"])
                        row["validation"] = cache[key]
                        row["validation_cache_key"] = key
                        outcome = "passed" if row["validation"]["passed"] is True else "failed"
                    else:
                        outcome = row["status"]
                    counts[f"{row['tool']}_{outcome}"] += 1
                    output.write(canonical(row) + "\n")
                    if (index + 1) % 100 == 0:
                        output.flush()
                        print(json.dumps({"session": directory.name, "validated_rows": index + 1,
                                          "unique_artifacts": len(cache)}), flush=True)
            summary.append({"session": directory.name, "counts": dict(counts),
                            "validated_sha256": digest((directory / "validated.jsonl").read_bytes())})
        dump(args.study / "validation_summary.json", {"protocol_sha256": protocol_sha,
            "oracle_sha256": digest(oracle.read_bytes()), "sessions": summary,
            "unique_artifacts": len(cache),
            "unique_evidence_classes": dict(Counter(value["evidence_class"] for value in cache.values())),
            "max_residual": max((value.get("max_input_residual", 0) for value in cache.values()), default=0)})
    finally:
        lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("corpus", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--compiler-commit", required=True)
    p.add_argument("--trials", type=int, default=7)
    p.add_argument("--sessions", type=int, default=2)
    p.add_argument("--timeout", type=float, default=60.0)
    p.add_argument("--seed", type=int, default=20261018)
    p.add_argument("--machine-note", default="")
    t = commands.add_parser("time")
    t.add_argument("study", type=Path)
    t.add_argument("--binary", required=True, type=Path)
    t.add_argument("--session", required=True, type=int)
    v = commands.add_parser("validate")
    v.add_argument("study", type=Path)
    w = commands.add_parser("worker")
    w.add_argument("tool", choices=("cirq", "tket"))
    args = parser.parse_args()
    if args.command == "worker":
        light_worker(args.tool)
    else:
        {"prepare": prepare, "time": timing, "validate": validate}[args.command](args)


if __name__ == "__main__":
    main()
