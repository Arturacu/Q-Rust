"""Frozen external protocol and evidence integrity; no compiler benchmarks."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import external_comparison as study


def fixture(name="external_example"):
    text = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[2]; h q[0]; cx q[0],q[1];'
    return {"id": name, "fixture": name, "family": "test", "community_source": "test-suite",
            "n_logical": 2, "source_qasm": text, "qasm": text,
            "source_sha256": study.digest(text), "normalized_sha256": study.digest(text),
            "normalization_evidence": {"passed": True, "mode": "full_unitary_numerical", "input_count": 4,
                                       "tolerance": 1e-7, "residuals": {"source": 0.0, "normalized": 0.0}}}


def manifest():
    return {"schema": 1, "cases": [fixture()], "exclusions": [{"id": "excluded", "reason": "input cap"}],
            "provenance": {"source_revision": "frozen-test"}}


class ProtocolTests(unittest.TestCase):
    def test_target_grid_preserves_same_input_and_source_for_every_recipe(self):
        cases = study.build_cases(manifest())
        self.assertEqual([case["topology"] for case in cases], ["all_to_all", "linear"])
        self.assertEqual(cases[0]["edges"], [[0, 1]])
        self.assertEqual(cases[1]["edges"], [[0, 1]])
        self.assertTrue(all(case["qasm"] == fixture()["qasm"] for case in cases))
        self.assertTrue(all(case["n_logical"] == case["physical_width"] == 2 for case in cases))

    def test_unconstrained_topology_explicitly_connects_every_native_worker_pair(self):
        value = manifest()
        source = value["cases"][0]
        source["n_logical"] = 4
        source["normalization_evidence"]["input_count"] = 16
        case = next(c for c in study.build_cases(value) if c["topology"] == "all_to_all")
        self.assertEqual(case["edges"], [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])

    def test_duplicate_fixture_ids_and_unsupported_width_fail_closed(self):
        invalid = manifest()
        invalid["cases"] *= 2
        with self.assertRaisesRegex(ValueError, "unique"):
            study.build_cases(invalid)
        for width in (0, 17, True):
            invalid = manifest()
            invalid["cases"][0]["n_logical"] = width
            with self.assertRaisesRegex(ValueError, "width"):
                study.build_cases(invalid)

    def test_source_or_normalized_hash_mismatch_rejected(self):
        for key in ("source_qasm", "qasm"):
            invalid = manifest()
            invalid["cases"][0][key] += " x q[1];"
            with self.assertRaisesRegex(ValueError, "checksum"):
                study.build_cases(invalid)

    def test_sessions_have_same_seed_set_but_independent_order(self):
        cases = study.build_cases(manifest())
        first = study.schedule(cases, 7, 123, 1)
        second = study.schedule(cases, 7, 123, 2)
        self.assertEqual(first, study.schedule(cases, 7, 123, 1))
        self.assertNotEqual(first, second)
        def identities(rows):
            return sorted((row["case_id"], row["tool"], row["trial"], row["seed"]) for row in rows)
        self.assertEqual(identities(first), identities(second))
        self.assertEqual(len(first), 2 * 8 * 6)
        self.assertEqual(len({(r["case_id"], r["tool"], r["trial"]) for r in first}), len(first))
        for case in cases:
            selected = [row for row in first if row["case_id"] == case["id"]]
            self.assertEqual([row["trial"] for row in selected[:6]], [-1] * 6)

    def test_both_profiles_share_exact_qrust_qiskit_recipe_ids(self):
        self.assertEqual(study.PROFILES["primary"], ["qrust", "qiskit", "cirq", "tket"])
        self.assertEqual(study.PROFILES["light"], ["qrust", "qiskit", "cirq_one_pass", "tket_light"])

    def test_normalization_evidence_required_and_must_cover_declared_inputs(self):
        for change in ({"passed": False}, {"input_count": 1}, {"mode": "unknown"}, {"tolerance": 0.1},
                       {"residuals": {"source": float("nan"), "normalized": 0}}):
            value = manifest()
            value["cases"][0]["normalization_evidence"].update(change)
            with self.assertRaisesRegex(ValueError, "normalization evidence"):
                study.build_cases(value)
        value = manifest()
        del value["cases"][0]["normalization_evidence"]
        with self.assertRaisesRegex(ValueError, "normalization evidence"):
            study.build_cases(value)


class FrozenStudyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.binary = self.directory / "binary"
        self.binary.write_bytes(b"frozen worker")
        self.corpus = self.directory / "manifest.json"
        self.corpus.write_text(json.dumps(manifest()))
        self.output = self.directory / "study"
        args = SimpleNamespace(corpus=self.corpus, output=self.output, binary=self.binary,
            trials=2, sessions=2, timeout=1.0, seed=123, compiler_commit="test-commit", machine_note="test")
        with patch.object(study, "CONTRACT_FILES", ()), \
                patch.object(study, "source_provenance", return_value={"compiler_source_sha256": {}}), \
                patch.object(study.importlib.metadata, "version", return_value="test-version"), \
                redirect_stdout(io.StringIO()):
            study.prepare(args)
        self.protocol, self.protocol_sha = study.load_protocol(self.output)

    def timing(self, session, worker):
        with patch.object(study, "check_runtime"), patch.object(study, "ProfileWorker", worker), \
                redirect_stdout(io.StringIO()):
            study.timing(SimpleNamespace(study=self.output, binary=self.binary, session=session))

    def test_freeze_retains_excluded_inputs_and_source_manifest(self):
        self.assertEqual(self.protocol["corpus_exclusions"], manifest()["exclusions"])
        self.assertEqual(json.loads((self.output / "corpus_manifest.json").read_text()), manifest())

    def test_protocol_file_tampering_rejected(self):
        path = self.output / "protocol.json"
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "protocol checksum"):
            study.load_protocol(self.output)

    def test_corpus_file_tampering_rejected(self):
        path = self.output / "corpus_manifest.json"
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "corpus checksum"):
            study.load_protocol(self.output)

    def test_schedule_tampering_rejected_even_with_recomputed_protocol_hash(self):
        self.protocol["schedules"]["1"][0]["seed"] += 1
        study.dump(self.output / "protocol.json", self.protocol)
        (self.output / "protocol.sha256").write_text(study.digest((self.output / "protocol.json").read_bytes()))
        with self.assertRaisesRegex(ValueError, "schedule"):
            study.load_protocol(self.output)

    def test_changed_binary_rejected_before_worker_start(self):
        self.binary.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "binary"):
            study.check_runtime(self.protocol, self.binary)

    def test_timeout_retained_and_restarted_worker_gets_an_extra_warmup(self):
        workers = {}
        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                self.tool, self.calls, self.closed, self.process = tool, 0, False, object()
                workers[tool] = self

            def call(self, request, timeout):
                self.calls += 1
                self.process = object()
                if self.tool == "qrust" and self.calls == 2:
                    self.process = None
                    return {"status": "timeout", "error": "preserved timeout"}
                return {"status": "compiled", "compile_ns": 1, "parse_ns": 0,
                        "out_qasm": request["qasm"], "initial_layout": [0, 1],
                        "final_layout": [0, 1], "output_qubits": 2}

            def close(self):
                self.closed = True
        self.timing(1, FakeWorker)
        rows = study.read_rows(self.output / "session_01/timed.jsonl")
        warmups = study.read_rows(self.output / "session_01/warmup.jsonl")
        self.assertEqual(len(rows), 2 * 2 * 6)
        self.assertEqual(sum(row["status"] == "timeout" for row in rows), 1)
        self.assertEqual(len(warmups), 2 * 6 + 1)
        self.assertTrue(all(worker.closed for worker in workers.values()))
        self.assertFalse((self.output / "active.lock").exists())
        study.check_rows(self.protocol, 1, rows)
        changed = json.loads(json.dumps(rows))
        row = next(row for row in changed if row["status"] == "compiled")
        row["final_layout"] = [1, 0]
        with self.assertRaisesRegex(ValueError, "artifact/layout"):
            study.check_rows(self.protocol, 1, changed)
        changed = json.loads(json.dumps(rows))
        changed[0]["trial"] = 7
        with self.assertRaisesRegex(ValueError, "unexpected"):
            study.check_rows(self.protocol, 1, changed)
        with self.assertRaisesRegex(ValueError, "Missing"):
            study.check_rows(self.protocol, 1, rows[:-1])

    def test_startup_failures_retained_for_every_selected_trial(self):
        class FailedWorker:
            def __init__(self, tool, binary, outdir):
                raise RuntimeError("startup failure")
        self.timing(1, FailedWorker)
        rows = study.read_rows(self.output / "session_01/timed.jsonl")
        self.assertEqual(len(rows), 24)
        self.assertTrue(all(row["phase"] == "worker_startup" and row["status"] == "error" for row in rows))
        study.check_rows(self.protocol, 1, rows)

    def test_failed_replacement_warmup_cannot_become_a_cold_measured_success(self):
        workers = {}
        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                self.tool, self.calls, self.process = tool, 0, object()
                workers[tool] = self

            def call(self, request, timeout):
                self.calls += 1
                self.process = object()
                if self.tool == "qrust" and self.calls in {2, 3}:
                    self.process = None
                    return {"status": "timeout", "error": "test timeout"}
                return {"status": "compiled", "compile_ns": 1, "parse_ns": 0,
                        "out_qasm": request["qasm"], "initial_layout": [0, 1],
                        "final_layout": [0, 1], "output_qubits": 2}

            def close(self):
                pass
        self.timing(1, FakeWorker)
        rows = study.read_rows(self.output / "session_01/timed.jsonl")
        qrust = [row for row in rows if row["tool"] == "qrust"]
        self.assertEqual([row["status"] for row in qrust], ["timeout", "timeout", "compiled", "compiled"])
        self.assertEqual(qrust[1]["phase"], "replacement_warmup")
        self.assertEqual(workers["qrust"].calls, 6)

    def test_primary_startup_exception_releases_the_process(self):
        closed = []
        def failed_start(worker):
            raise RuntimeError("readiness timeout")
        def close(worker):
            closed.append(worker)
        with patch.object(study.Worker, "start", failed_start), patch.object(study.ProfileWorker, "close", close):
            with self.assertRaisesRegex(RuntimeError, "readiness"):
                study.ProfileWorker("qrust", self.binary, self.directory)
        self.assertGreaterEqual(len(closed), 1)

    def test_active_session_lock_blocks_timing_and_validation(self):
        (self.output / "active.lock").touch()
        with patch.object(study, "check_runtime"), patch.object(study, "ProfileWorker") as worker:
            with self.assertRaises(FileExistsError):
                study.timing(SimpleNamespace(study=self.output, binary=self.binary, session=1))
            worker.assert_not_called()
        with self.assertRaises(FileExistsError):
            study.validate(SimpleNamespace(study=self.output))

    def test_validation_refuses_missing_session_before_creating_evidence(self):
        with self.assertRaises(FileNotFoundError):
            study.validate(SimpleNamespace(study=self.output))
        self.assertFalse((self.output / "validation_summary.json").exists())
        self.assertFalse((self.output / "active.lock").exists())


class TimingBoundaryTests(unittest.TestCase):
    def test_residual_request_overhead_is_not_claimed_as_serialization(self):
        worker = SimpleNamespace(call=lambda request, timeout: {"status": "compiled", "parse_ns": 10, "compile_ns": 60})
        with patch.object(study.time, "perf_counter_ns", side_effect=[100, 200]):
            result = study.timed_request(worker, {}, 1)
        self.assertEqual(result["request_wall_ns"], 100)
        self.assertEqual(result["noncompile_nonparse_request_ns"], 30)
        self.assertNotIn("serialization_ns", result)

    def test_invalid_or_impossible_worker_times_fail_closed(self):
        for parse, compile_time in ((-1, 10), (0, 0), (0, 1.2), (60, 60)):
            worker = SimpleNamespace(call=lambda request, timeout: {"status": "compiled", "parse_ns": parse, "compile_ns": compile_time})
            with patch.object(study.time, "perf_counter_ns", side_effect=[100, 200]):
                result = study.timed_request(worker, {}, 1)
            self.assertEqual(result["status"], "error")


if __name__ == "__main__":
    unittest.main()
