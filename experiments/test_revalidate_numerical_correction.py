"""Guard evidence integrity and failure retention in correctness revalidation."""
import copy
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import revalidate_numerical_correction as study


class RevalidationEvidenceTests(unittest.TestCase):
    def setUp(self):
        qasm = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[1];'
        self.case = {"id": "one", "family": "test", "source_qasm": qasm, "qasm": qasm,
                     "source_sha256": study.digest(qasm), "normalized_sha256": study.digest(qasm),
                     "physical_width": 1, "n_logical": 1, "edges": [], "topology": "all_to_all"}
        self.row = {"id": "one__qrust__0", "case_id": "one", "tool": "qrust", "trial": 0,
                    "status": "compiled", "out_qasm": qasm, "output_sha256": study.digest(qasm),
                    "initial_layout": [0], "final_layout": [0]}
        self.protocol = {"cases": [self.case], "tools": ["qrust"], "trials": 1,
                         "source_sha256": {"experiments/comparison_validation.py": "frozen_oracle"}}

    def test_refuses_changed_oracle_and_changed_input_or_output(self):
        with self.assertRaisesRegex(ValueError, "oracle differs"):
            study.check_frozen_inputs(self.protocol, [self.row], "different_oracle")
        changed = copy.deepcopy(self.protocol)
        changed["cases"][0]["qasm"] += " x q[0];"
        with self.assertRaisesRegex(ValueError, "Fixture checksum"):
            study.check_frozen_inputs(changed, [self.row], "frozen_oracle")
        changed_row = {**self.row, "out_qasm": self.row["out_qasm"] + " x q[0];"}
        with self.assertRaisesRegex(ValueError, "Output checksum"):
            study.check_frozen_inputs(self.protocol, [changed_row], "frozen_oracle")

    def test_missing_duplicate_and_unexpected_trials_fail_closed(self):
        for rows in ([], [self.row, self.row], [{**self.row, "trial": 3}]):
            with self.assertRaises(ValueError):
                study.check_frozen_inputs(self.protocol, rows, "frozen_oracle")
        self.assertEqual(study.check_frozen_inputs(self.protocol, [self.row], "frozen_oracle"), {"one": self.case})

    def test_cache_does_not_conflate_different_layouts(self):
        case = {**self.case, "physical_width": 2, "edges": [[0, 1]]}
        left = {**self.row, "initial_layout": [0], "final_layout": [0]}
        right = {**self.row, "initial_layout": [1], "final_layout": [0]}
        self.assertNotEqual(study.artifact_key(case, left), study.artifact_key(case, right))
        right = {**left, "final_layout": [1]}
        self.assertNotEqual(study.artifact_key(case, left), study.artifact_key(case, right))

    def test_failed_and_inconclusive_results_are_retained(self):
        cache = {}
        failed = {"passed": False, "evidence_class": "full_isometry", "max_input_residual": 4.5e-7}
        with patch.object(study, "validate_output", return_value=failed):
            validation, key = study.validate_row(self.case, self.row, cache)
        rows = [{**self.row, "validation": validation, "validation_cache_key": key},
                {**self.row, "id": "unchecked", "validation": {"passed": None, "evidence_class": "unchecked_width"},
                 "validation_cache_key": "unchecked"}]
        summary = study.summarize(rows, cache)
        self.assertFalse(summary["all_passed"])
        self.assertEqual(summary["by_tool"]["qrust"]["failed"], 1)
        self.assertEqual(summary["by_tool"]["qrust"]["unchecked"], 1)
        self.assertEqual(summary["failed_or_unchecked_ids"], [self.row["id"], "unchecked"])

    def test_worker_error_does_not_get_a_successful_validation(self):
        with patch.object(study, "validate_output") as oracle:
            validation, key = study.validate_row(self.case, {"status": "timeout", "error": "expired"}, {})
        oracle.assert_not_called()
        self.assertFalse(validation["passed"])
        self.assertEqual(validation["error"], "expired")
        self.assertIsNone(key)

    def test_corrected_campaign_runs_from_source_export_without_git(self):
        """Exercise the full corrected path, not only a mocked provenance helper."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "source-export"
            adapter = root / "experiments/src/bin/comparison_runner.rs"
            adapter.parent.mkdir(parents=True)
            adapter.write_text("// exported worker adapter\n")
            untracked = root / "src/new_policy.rs"
            untracked.parent.mkdir()
            untracked.write_text("// source absent from any Git index\n")
            (root / "Cargo.toml").write_text('[package]\nname = "export"\n')
            (root / "Cargo.lock").write_text("version = 3\n")
            ignored = root / "target/generated.rs"
            ignored.parent.mkdir()
            ignored.write_text("// build artifact must not enter source manifest\n")
            binary = Path(directory) / "worker"
            binary.write_bytes(b"frozen-worker-test-double")
            args = SimpleNamespace(previous=root / "prior", output=Path(directory) / "results", binary=binary, timeout=1)
            prior = ({"master_seed": 123}, [self.row], {"one": self.case}, {}, {})
            worker_result = {k: self.row[k] for k in ("status", "out_qasm", "output_sha256", "initial_layout", "final_layout")}
            with patch.object(study, "HERE", root / "experiments"), \
                    patch.object(study, "load_previous", return_value=prior), \
                    patch.object(study, "Worker") as worker, \
                    patch.object(study.subprocess, "check_output", side_effect=FileNotFoundError("git unavailable")), \
                    redirect_stdout(io.StringIO()):
                worker.return_value.call.return_value = worker_result
                self.assertEqual(study.recompile_corrected(args), 0)
            report = json.loads((args.output / "summary.json").read_text())
            self.assertTrue(report["all_passed"])
            meta = report["provenance"]
            self.assertFalse(meta["git_available"])
            self.assertIsNone(meta["git_head"])
            self.assertIsNone(meta["git_diff_sha256"])
            self.assertEqual(meta["compiler_source_sha256"]["src/new_policy.rs"], study.digest(untracked.read_bytes()))
            self.assertEqual(set(meta["compiler_source_sha256"]), {
                "Cargo.toml", "Cargo.lock", "src/new_policy.rs", "experiments/src/bin/comparison_runner.rs"})
            self.assertEqual(meta["corrected_binary_sha256"], meta["corrected_binary_sha256_after"])
            worker.return_value.close.assert_called_once()

    def test_source_export_cannot_claim_enclosing_repository_commit(self):
        with tempfile.TemporaryDirectory() as directory:
            outer = Path(directory)
            exported = outer / "exported"
            exported.mkdir()
            (exported / "Cargo.toml").write_text("[workspace]\n")
            with patch.object(study.subprocess, "check_output", return_value=str(outer)) as git:
                meta = study.source_provenance(exported)
            self.assertFalse(meta["git_available"])
            self.assertIsNone(meta["git_head"])
            self.assertIn("enclosing", meta["git_unavailable_reason"])
            self.assertIn("Cargo.toml", meta["compiler_source_sha256"])
            git.assert_called_once()


if __name__ == "__main__":
    unittest.main()
