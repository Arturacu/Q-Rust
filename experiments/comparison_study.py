#!/usr/bin/env python3
"""Repeated, matched compilation study with independently checked artifacts.

Run timing and validation as separate commands. No numerical simulation runs
during timing. Workers persist, requests are sequential, and every outcome is
retained. See protocol.json in each result directory for the frozen contract.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import select
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
TOOLS = ("qrust", "qiskit", "cirq", "tket")
MASTER_SEED = 20261015
THREAD_ENV = {name: "1" for name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RAYON_NUM_THREADS",
    "QISKIT_NUM_PROCS")}
THREAD_ENV.update(QISKIT_PARALLEL="FALSE", PYTHONHASHSEED="20261015")


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def worker(tool):
    if tool == "cirq":
        from comparison_cirq import compile_one
    else:
        from comparison_qiskit_tket import compile_one
    print(json.dumps({"ready": True, "tool": tool}), flush=True)
    for line in sys.stdin:
        try:
            request = json.loads(line)
            args = {key: request[key] for key in
                    ("qasm", "n_logical", "physical_width", "edges", "topology", "seed")}
            result = compile_one(**args) if tool == "cirq" else compile_one(tool, **args)
            result.setdefault("status", "compiled")
            result.setdefault("seed_used", tool == "qiskit")
        except Exception as error:
            result = {"status": "error", "error": f"{type(error).__name__}: {error}"}
        print(json.dumps(result), flush=True)


class Worker:
    def __init__(self, tool, binary, outdir):
        self.tool, self.binary, self.outdir = tool, binary, outdir
        self.process = None
        self.start()

    def start(self):
        env = dict(os.environ, **THREAD_ENV)
        env.setdefault("MPLCONFIGDIR", str(self.outdir / "matplotlib_cache"))
        command = [str(self.binary)] if self.tool == "qrust" else [
            sys.executable, str(Path(__file__).resolve()), "worker", self.tool]
        self.log = (self.outdir / f"worker_{self.tool}.log").open("a")
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=self.log, text=True, bufsize=1, env=env)
        ready = self.receive(120)
        if not ready.get("ready"):
            raise RuntimeError(f"{self.tool} failed to start: {ready}")

    def receive(self, timeout):
        if not select.select([self.process.stdout], [], [], timeout)[0]:
            raise TimeoutError(f"{self.tool} exceeded {timeout} seconds")
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError(f"{self.tool} worker exited ({self.process.poll()})")
        return json.loads(line)

    def call(self, request, timeout):
        try:
            if self.process is None:
                self.start()
            self.process.stdin.write(json.dumps(request) + "\n")
            self.process.stdin.flush()
            return self.receive(timeout)
        except Exception as error:
            self.close()
            return {"status": "timeout" if isinstance(error, TimeoutError) else "error",
                    "error": f"{type(error).__name__}: {error}", "worker_restart_required": True}

    def close(self):
        if self.process:
            self.process.terminate()
            try:
                self.process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
            self.process = None
            self.log.close()


def build_cases(pilot=False):
    from bench_baselines import normalize_qasm, topologies_for
    fixtures = json.loads((HERE / "fixtures/manifest.json").read_text())
    cases = []
    pilot_names = {"qft_5", "wstate_5", "adder_6", "clifford_7", "ghz_3", "qft_14"}
    for fixture in fixtures:
        n = fixture["num_qubits"]
        if n > 14 or (pilot and fixture["name"] not in pilot_names):
            continue
        original = (HERE / fixture["file"]).read_text()
        normalized = normalize_qasm(original)
        for topology, (width, edges) in topologies_for(n).items():
            if topology not in {"all_to_all", "linear", "ibm_nairobi"}:
                continue
            cases.append({"id": f"{fixture['name']}__{topology}", "fixture": fixture["name"],
                          "family": fixture["family"], "topology": topology,
                          "n_logical": n, "physical_width": width, "edges": edges,
                          "source_qasm": original, "qasm": normalized,
                          "source_sha256": digest(original), "normalized_sha256": digest(normalized)})
    return sorted(cases, key=lambda c: c["id"])


def timing(args):
    outdir = args.output.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    cases = build_cases(args.pilot)
    sources = list(HERE.glob("comparison_*.py")) + [HERE / "src/bin/comparison_runner.rs",
               HERE / "bench_baselines.py", HERE.parent / "Cargo.lock"]
    versions = {name: importlib.metadata.version(name) for name in
                ("qiskit", "cirq-core", "pytket", "numpy", "scipy")}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    protocol = {
        "schema": 1, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "master_seed": MASTER_SEED, "pilot": args.pilot, "trials": args.trials,
        "warmups_per_cell_tool": 1, "timeout_seconds": args.timeout, "tools": TOOLS,
        "cases": cases, "compiler_commit": commit, "versions": versions,
        "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
        "logical_cpu_count": os.cpu_count(), "thread_environment": THREAD_ENV,
        "machine_note": args.machine_note,
        "worker_executable": str(Path(sys.executable).resolve()),
        "qrust_binary_sha256": digest(args.binary.read_bytes()),
        "source_sha256": {str(p.relative_to(HERE.parent)): digest(p.read_bytes()) for p in sources},
        "timing_boundary": "Warm in-process API: input clone, configuration and compiler passes included; imports, source parsing, serialization, IPC and validation excluded. parse_ns separately retained.",
        "target": "Arbitrary one-qubit U3 and CX, bidirectional ideal connectivity. Free input/output permutation. Clean zero ancillas, arbitrary logical input.",
        "normalization": "Every tool receives identical exact CRZ and RZZ expansion into RZ/CX. CCX remains for timed lowering. Validation compares against the original unnormalized fixture.",
        "validation": "Independent NumPy analytic gate simulator. Full isometry through six physical qubits; otherwise three basis, uniform and four seeded Haar inputs. One shared global phase; full physical amplitudes; residual threshold 1e-7. Not proof or hardware fidelity.",
        "sampling_limit": "51 in-house fixtures <=14 logical qubits, nine families; no external holdout. Seven repetitions are not seven independent circuits. Preset recipes have unequal search effort.",
        "runtime_limit": "Single host, single thread environment caps, no CPU affinity/frequency control. Q-Rust internal routing RNG is not externally seedable. Qiskit uses recorded seeds; Cirq/TKET recipes are deterministic.",
    }
    (outdir / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rng = random.Random(MASTER_SEED)
    rng.shuffle(cases)
    workers = {}
    started = time.monotonic()
    try:
        for tool in TOOLS:
            workers[tool] = Worker(tool, args.binary, outdir)
        with (outdir / "timed.jsonl").open("w") as records, (outdir / "warmup.jsonl").open("w") as warmups:
            for ci, case in enumerate(cases):
                for trial in range(-1, args.trials):
                    order = list(TOOLS)
                    rng.shuffle(order)
                    for position, tool in enumerate(order):
                        seed = MASTER_SEED + max(trial, 0)
                        request = {key: case[key] for key in
                                   ("qasm", "n_logical", "physical_width", "edges", "topology")}
                        request["seed"] = seed
                        if trial >= 0 and workers[tool].process is None:
                            replacement = workers[tool].call(request, args.timeout)
                            warmups.write(json.dumps({"case_id": case["id"], "tool": tool,
                                                      "reason": "worker_restart", **replacement}) + "\n")
                            warmups.flush()
                        result = workers[tool].call(request, args.timeout)
                        row = {"id": f"{case['id']}__{tool}__{trial}", "case_id": case["id"],
                               "tool": tool, "trial": trial, "seed": seed,
                               "case_order": ci, "tool_order": position, **result}
                        if "out_qasm" in row:
                            row["output_sha256"] = digest(row["out_qasm"])
                        stream = warmups if trial == -1 else records
                        stream.write(json.dumps(row, separators=(",", ":")) + "\n")
                        stream.flush()
                print(json.dumps({"completed_cases": ci + 1, "total_cases": len(cases),
                                  "last_case": case["id"], "elapsed_seconds": round(time.monotonic()-started, 1)}), flush=True)
    finally:
        for process in workers.values():
            process.close()


def validate(args):
    from comparison_validation import validate_output
    protocol = json.loads((args.output / "protocol.json").read_text())
    cases = {c["id"]: c for c in protocol["cases"]}
    for case in cases.values():
        if digest(case["source_qasm"]) != case["source_sha256"] or digest(case["qasm"]) != case["normalized_sha256"]:
            raise ValueError("fixture checksum mismatch")
    cache, counts = {}, {}
    destination = args.output / "validated.jsonl"
    if destination.exists():
        raise FileExistsError(destination)
    with destination.open("w") as output, (args.output / "timed.jsonl").open() as timed:
        for index, line in enumerate(timed):
            row = json.loads(line)
            case = cases[row["case_id"]]
            if row["status"] == "compiled":
                if digest(row["out_qasm"]) != row["output_sha256"]:
                    raise ValueError("output checksum mismatch")
                key = digest(json.dumps([case["source_sha256"], case["physical_width"], case["edges"],
                                        row["output_sha256"], row["initial_layout"], row["final_layout"]]))
                if key not in cache:
                    edges = None if case["topology"] == "all_to_all" else [
                        pair for a, b in case["edges"] for pair in ((a, b), (b, a))]
                    cache[key] = validate_output(case["source_qasm"], row["out_qasm"],
                                                 case["n_logical"], case["physical_width"], edges,
                                                 row["initial_layout"], row["final_layout"])
                row["validation"] = cache[key]
                row["validation_cache_key"] = key
                status = "passed" if row["validation"]["passed"] is True else "failed"
            else:
                status = row["status"]
            counts[(row["tool"], status)] = counts.get((row["tool"], status), 0) + 1
            output.write(json.dumps(row, separators=(",", ":")) + "\n")
            if (index + 1) % 100 == 0:
                output.flush()
                print(f"validated {index+1}; unique artifacts {len(cache)}", flush=True)
    print(json.dumps({"counts": {f"{k[0]}_{k[1]}": v for k, v in counts.items()},
                      "unique_artifacts": len(cache)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    w = commands.add_parser("worker")
    w.add_argument("tool", choices=("qiskit", "cirq", "tket"))
    t = commands.add_parser("time")
    t.add_argument("output", type=Path)
    t.add_argument("--binary", type=Path, required=True)
    t.add_argument("--trials", type=int, default=7)
    t.add_argument("--timeout", type=float, default=60)
    t.add_argument("--pilot", action="store_true")
    t.add_argument("--machine-note", default="CPU model and memory not separately recorded")
    v = commands.add_parser("validate")
    v.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "worker":
        worker(args.tool)
    elif args.command == "time":
        if args.trials < 1:
            parser.error("trials must be positive")
        timing(args)
    else:
        validate(args)


if __name__ == "__main__":
    main()
