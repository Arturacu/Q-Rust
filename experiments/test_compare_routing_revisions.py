"""Protocol regressions for immutable tasks and preserved revision failures."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import compare_routing_revisions as revisions
from test_comparison_study import case


class RevisionProtocolTests(unittest.TestCase):
    def setup_args(self, directory, request=None):
        directory = Path(directory)
        binary = directory / "binary"
        binary.write_bytes(b"fixed binary")
        variants = directory / "variants.json"
        variants.write_text(json.dumps({"before": {"binary": str(binary)},
                                       "after": {"binary": str(binary), "request": request or {}}}))
        return SimpleNamespace(variants=variants, output=directory / "results", case=None,
                               trials=2, timeout=1, machine_note="synthetic test")

    def test_variants_cannot_change_problem_or_paired_seed(self):
        for key in (*revisions.REQUEST_FIELDS, "seed"):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as temporary:
                args = self.setup_args(temporary, {key: "changed"})
                with self.assertRaisesRegex(ValueError, "must not alter"):
                    revisions.load_variants(args.variants)

    def test_binary_hash_and_request_controls_are_recorded(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.setup_args(temporary, {"lookahead_size": 20})
            variants = revisions.load_variants(args.variants)
        self.assertEqual(variants["before"]["binary_sha256"], revisions.study.digest(b"fixed binary"))
        self.assertEqual(variants["after"]["request"], {"lookahead_size": 20})

    def test_failed_trial_stays_in_results_and_replacement_worker_is_warmed(self):
        workers = []
        requests = []

        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                self.label = outdir.name
                self.calls = 0
                self.process = object()
                self.closed = False
                workers.append(self)

            def call(self, request, timeout):
                self.calls += 1
                requests.append((self.label, dict(request)))
                self.process = object()
                if self.label == "after" and self.calls == 2:
                    self.process = None
                    return {"status": "timeout"}
                return {"status": "compiled", "compile_ns": 10, "parse_ns": 1,
                        "out_qasm": request["qasm"], "initial_layout": [0], "final_layout": [0]}

            def close(self):
                self.closed = True

        with tempfile.TemporaryDirectory() as temporary, \
                patch.object(revisions.study, "Worker", FakeWorker), \
                patch.object(revisions.study, "build_cases", return_value=[case()]), \
                patch.object(revisions.platform, "platform", return_value="test"), \
                patch.object(revisions.subprocess, "check_output", return_value="commit"), \
                redirect_stdout(io.StringIO()):
            args = self.setup_args(temporary)
            revisions.timing(args)
            measured = [json.loads(l) for l in (args.output / "timed.jsonl").read_text().splitlines()]
            warmups = [json.loads(l) for l in (args.output / "warmup.jsonl").read_text().splitlines()]
        self.assertEqual(len(measured), 4)
        self.assertEqual(sum(r["status"] == "timeout" for r in measured), 1)
        self.assertEqual(len(warmups), 3)
        self.assertTrue(all(w.closed for w in workers))
        self.assertTrue(all(r["seed"] == revisions.SEED + r["trial"] for r in measured))
        self.assertTrue(all(req["qasm"] == case()["qasm"] for _, req in requests))
        self.assertEqual(next(w for w in workers if w.label == "after").calls, 4)

    def test_unknown_case_rejected_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temporary, patch.object(revisions.study, "build_cases", return_value=[case()]):
            args = self.setup_args(temporary)
            args.case = ["invented"]
            with self.assertRaisesRegex(ValueError, "unknown cases"):
                revisions.timing(args)
            self.assertFalse(args.output.exists())

    def test_already_started_workers_closed_on_later_start_failure(self):
        started = []

        class FakeWorker:
            def __init__(self, tool, binary, outdir):
                if outdir.name == "after":
                    raise RuntimeError("start failed")
                self.closed = False
                started.append(self)

            def close(self):
                self.closed = True

        with tempfile.TemporaryDirectory() as temporary, \
                patch.object(revisions.study, "Worker", FakeWorker), \
                patch.object(revisions.study, "build_cases", return_value=[case()]), \
                patch.object(revisions.platform, "platform", return_value="test"), \
                patch.object(revisions.subprocess, "check_output", return_value="commit"):
            with self.assertRaisesRegex(RuntimeError, "start failed"):
                revisions.timing(self.setup_args(temporary))
        self.assertEqual(len(started), 1)
        self.assertTrue(started[0].closed)


if __name__ == "__main__":
    unittest.main()
