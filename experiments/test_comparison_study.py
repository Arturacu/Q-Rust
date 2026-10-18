"""Lightweight protocol and artifact-integrity tests; no compilers are timed."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import comparison_study as study
from bench_baselines import normalize_qasm

try:
    import numpy as np
    from comparison_validation import parse_qasm, simulate
except ImportError:
    np = None


def text(width, operations=""):
    return f'OPENQASM 2.0; include "qelib1.inc"; qreg q[{width}]; {operations}'


def case():
    qasm = text(1)
    return {"id": "empty__all_to_all", "fixture": "empty", "family": "test",
            "topology": "all_to_all", "n_logical": 1, "physical_width": 1,
            "edges": [], "source_qasm": qasm, "qasm": qasm,
            "source_sha256": study.digest(qasm), "normalized_sha256": study.digest(qasm)}


def record(tool="qrust", output=None):
    output = text(1) if output is None else output
    return {"id": f"empty__all_to_all__{tool}__0", "case_id": "empty__all_to_all",
            "tool": tool, "trial": 0, "status": "compiled", "out_qasm": output,
            "output_sha256": study.digest(output), "initial_layout": [0], "final_layout": [0]}


class StudyProtocolTests(unittest.TestCase):
    def test_prespecified_main_scope(self):
        cases = study.build_cases()
        self.assertEqual(len(cases), 123)
        self.assertEqual(len({row["fixture"] for row in cases}), 51)
        self.assertEqual(len({row["family"] for row in cases}), 9)
        self.assertEqual(len({row["id"] for row in cases}), len(cases))
        self.assertEqual(sum(c["topology"] == "all_to_all" for c in cases), 51)
        self.assertEqual(sum(c["topology"] == "linear" for c in cases), 51)
        self.assertEqual(sum(c["topology"] == "ibm_nairobi" for c in cases), 21)
        self.assertTrue(all(0 < c["n_logical"] <= c["physical_width"] <= 14 for c in cases))

    def test_every_tool_receives_same_normalized_content_with_original_retained(self):
        for row in study.build_cases():
            with self.subTest(case=row["id"]):
                self.assertEqual(row["source_sha256"], study.digest(row["source_qasm"]))
                self.assertEqual(row["normalized_sha256"], study.digest(row["qasm"]))
                self.assertEqual(row["qasm"], normalize_qasm(row["source_qasm"]))
                self.assertNotIn("crz(", row["qasm"])
                self.assertNotIn("rzz(", row["qasm"])

    @unittest.skipIf(np is None, "optional NumPy checker dependency not installed")
    def test_normalization_exact_for_adversarial_angles_and_both_wire_orders(self):
        for gate in ("crz", "rzz"):
            for angle in ("0", "0.37", "-0.81", "pi", "2*pi", "-2*pi", "4*pi"):
                for first, second in ((0, 1), (1, 0)):
                    with self.subTest(gate=gate, angle=angle, order=(first, second)):
                        original = text(2, f"{gate}({angle}) q[{first}],q[{second}];")
                        normalized = normalize_qasm(original)
                        n, source = parse_qasm(original)
                        _, expanded = parse_qasm(normalized)
                        before = simulate(n, source, np.eye(4))
                        after = simulate(n, expanded, np.eye(4))
                        self.assertLess(float(np.max(np.abs(before - after))), 1e-12)

    def _timing_args(self, directory):
        directory = Path(directory)
        binary = directory / "dummy_binary"
        binary.write_bytes(b"test binary")
        return SimpleNamespace(output=directory / "results", binary=binary, trials=2,
                               timeout=1.0, pilot=True, machine_note="test host")

    def test_workers_already_started_are_closed_if_later_start_fails(self):
        started = []

        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                if tool == "qiskit":
                    raise RuntimeError("test startup failure")
                self.closed = False
                started.append(self)

            def close(self):
                self.closed = True

        with tempfile.TemporaryDirectory() as temporary, \
                patch.object(study, "build_cases", return_value=[case()]), \
                patch.object(study, "Worker", FakeWorker), \
                patch.object(study.importlib.metadata, "version", return_value="test"), \
                patch.object(study.platform, "platform", return_value="test-platform"), \
                patch.object(study.subprocess, "check_output", return_value="test-commit"):
            with self.assertRaisesRegex(RuntimeError, "test startup failure"):
                study.timing(self._timing_args(temporary))
        self.assertEqual(len(started), 1)
        self.assertTrue(started[0].closed)

    def test_failed_trial_retained_and_restarted_worker_warmed_before_next_trial(self):
        workers = {}

        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                self.tool, self.calls, self.closed = tool, 0, False
                self.process = object()
                workers[tool] = self

            def call(self, request, timeout):
                self.calls += 1
                self.process = object()
                if self.tool == "qrust" and self.calls == 2:
                    self.process = None
                    return {"status": "timeout", "worker_restart_required": True, "error": "test timeout"}
                return {"status": "compiled", "compile_ns": 10, "parse_ns": 1,
                        "out_qasm": request["qasm"], "initial_layout": [0], "final_layout": [0]}

            def close(self):
                self.closed = True

        with tempfile.TemporaryDirectory() as temporary, \
                patch.object(study, "build_cases", return_value=[case()]), \
                patch.object(study, "Worker", FakeWorker), \
                patch.object(study.importlib.metadata, "version", return_value="test"), \
                patch.object(study.platform, "platform", return_value="test-platform"), \
                patch.object(study.subprocess, "check_output", return_value="test-commit"), \
                redirect_stdout(io.StringIO()):
            args = self._timing_args(temporary)
            study.timing(args)
            measured = [json.loads(line) for line in (args.output / "timed.jsonl").read_text().splitlines()]
            warmed = [json.loads(line) for line in (args.output / "warmup.jsonl").read_text().splitlines()]
        self.assertEqual(len(measured), 8)
        self.assertEqual(sum(row["status"] == "timeout" for row in measured), 1)
        self.assertEqual(workers["qrust"].calls, 4)  # warmup, failure, replacement warmup, trial
        self.assertEqual(len(warmed), 5)
        self.assertTrue(all(worker.closed for worker in workers.values()))


@unittest.skipIf(np is None, "optional NumPy checker dependency not installed")
class ArtifactIntegrityTests(unittest.TestCase):
    def validate(self, rows, cases=None):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "protocol.json").write_text(json.dumps({"cases": cases or [case()]}))
            (output / "timed.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
            with redirect_stdout(io.StringIO()):
                study.validate(SimpleNamespace(output=output))
            return [json.loads(line) for line in (output / "validated.jsonl").read_text().splitlines()]

    def test_identical_artifact_cache_may_cross_compiler_identity(self):
        with patch("comparison_validation.validate_output", return_value={"passed": True}) as oracle:
            result = self.validate([record(), record("qiskit")])
        self.assertEqual(oracle.call_count, 1)
        self.assertEqual(result[0]["validation_cache_key"], result[1]["validation_cache_key"])
        self.assertTrue(all(row["validation"]["passed"] for row in result))

    def test_changed_output_with_stale_hash_is_rejected_before_cache_reuse(self):
        changed = record("qiskit", text(1, "u3(0,0,pi) q[0];"))
        changed["output_sha256"] = record()["output_sha256"]
        with patch("comparison_validation.validate_output", return_value={"passed": True}) as oracle:
            with self.assertRaises(ValueError):
                self.validate([record(), changed])
        self.assertEqual(oracle.call_count, 1)

    def test_changed_original_source_with_stale_hash_is_rejected(self):
        changed_case = case()
        changed_case["source_qasm"] = text(1, "h q[0];")
        with patch("comparison_validation.validate_output", return_value={"passed": True}) as oracle:
            with self.assertRaises(ValueError):
                self.validate([record()], [changed_case])
        oracle.assert_not_called()

    def test_failed_compilation_is_retained_without_validation(self):
        failure = {"case_id": case()["id"], "tool": "qrust", "status": "timeout", "error": "timeout"}
        with patch("comparison_validation.validate_output") as oracle:
            result = self.validate([failure])
        self.assertEqual(result, [failure])
        oracle.assert_not_called()


if __name__ == "__main__":
    unittest.main()
