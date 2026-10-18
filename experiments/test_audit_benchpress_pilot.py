"""Small archive corruption controls for portable matched-pilot evidence replay.

No compiler or SDK runs. The fixture uses a two-qubit Bell preparation, three
timed outcomes and an independently computed Clifford-isometry check. Corruption
tests refresh outer receipts where appropriate so they exercise semantic binding,
not merely an incidental stale checksum.
"""
from collections import Counter
from contextlib import contextmanager
import copy
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import shutil
import statistics
import sys
import tempfile
import unittest

from audit_benchpress_pilot import audit


HAS_NUMPY = importlib.util.find_spec("numpy") is not None
HERE = Path(__file__).resolve().parent


def digest(data):
    return hashlib.sha256(data).hexdigest()


@contextmanager
def archived_modules(harness):
    """Use the archived modules, then restore any surrounding test imports."""
    names = ("comparison_validation", "scalable_clifford")
    previous = {name: sys.modules.pop(name, None) for name in names}
    original_path = list(sys.path)
    sys.path.insert(0, str(harness / "experiments"))
    importlib.invalidate_caches()
    try:
        yield
    finally:
        for name in names:
            sys.modules.pop(name, None)
            if previous[name] is not None:
                sys.modules[name] = previous[name]
        sys.path[:] = original_path


@unittest.skipUnless(HAS_NUMPY, "NumPy is required for independent evidence replay")
class BenchpressArchiveAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.evidence = root / "evidence"
        self.directory = self.evidence / "matched"
        self.directory.mkdir(parents=True)
        self.harness = root / "harness"
        (self.harness / "experiments").mkdir(parents=True)
        sources = {}
        for name in ("comparison_validation.py", "scalable_clifford.py"):
            destination = self.harness / "experiments" / name
            shutil.copyfile(HERE / name, destination)
            sources["/original/absolute/workspace/" + name] = digest(destination.read_bytes())
        source = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[2]; h q[0]; cx q[0],q[1];'
        output = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[2]; u3(pi/2,0,pi) q[0]; cx q[0],q[1];'
        case = {"id": "bell", "qasm": source, "qasm_sha256": digest(source.encode()),
                "n_logical": 2, "physical_width": 2, "edges": [[0, 1]], "topology": "linear"}
        schedule = [{"case_id": "bell", "tool": "qrust", "trial": trial,
                     "seed": 99 + max(0, trial)} for trial in (-1, 0, 1, 2)]
        self.protocol = {"source_sha256": sources, "cases": [case], "sessions": 1,
                         "tools": ["qrust"], "trials": 3, "warmups": 1,
                         "schedules": {"1": schedule}, "scope": "synthetic archive fixture",
                         "validation": "independent Clifford-isometry evidence"}
        self.raw = []
        for index, item in enumerate(schedule):
            scale = index or 1
            self.raw.append(dict(item, session=1, schedule_index=index,
                phase="warmup" if item["trial"] == -1 else "timed", status="compiled",
                out_qasm=output, output_sha256=digest(output.encode()),
                initial_layout=[0, 1], final_layout=[0, 1], output_qubits=2,
                compile_ns=100 * scale, parse_ns=10 * scale,
                request_wall_ns=200 * scale, sampled_worker_rss_bytes=1000 * scale))
        self.populate_validation()
        self.write_archive()

    def populate_validation(self):
        case = self.protocol["cases"][0]
        self.saved = copy.deepcopy(self.raw)
        with archived_modules(self.harness):
            from comparison_validation import common_metrics, parse_qasm
            from scalable_clifford import validate_clifford
            for row in self.saved:
                if row["status"] != "compiled":
                    continue
                check = validate_clifford(case["qasm"], row["out_qasm"], 2, 2,
                    row["initial_layout"], row["final_layout"], case["edges"])
                check["metrics"] = common_metrics(*parse_qasm(row["out_qasm"]))
                row["validation"] = check
                row["validation_cache_key"] = digest(json.dumps([
                    case["id"], row["output_sha256"], row["initial_layout"], row["final_layout"]]).encode())
        self.refresh_summary()

    def refresh_summary(self):
        counts = Counter()
        keys = set()
        for row in self.saved:
            status = ("passed" if row.get("validation", {}).get("passed") else
                      "unverified_output" if row["status"] == "compiled" else row["status"])
            counts[f'{row["phase"]}:{row["tool"]}:{status}'] += 1
            if row["status"] == "compiled":
                keys.add(row["validation_cache_key"])
        measured = [row for row in self.saved if row["phase"] == "timed"]
        warmups = [row for row in self.saved if row["phase"] != "timed"]
        eligible = len(measured) == 3 and bool(warmups) and all(
            row.get("validation", {}).get("passed") for row in [*measured, *warmups])
        cell = {"case_id": "bell", "tool": "qrust", "eligible": eligible,
                "status_counts": dict(Counter(row["status"] for row in measured)),
                "validation_errors": sorted({row["validation"].get("error", "") for row in measured
                    if "validation" in row and not row["validation"]["passed"]})}
        if eligible:
            values = {field: [row[field] for row in measured] for field in
                ("compile_ns", "parse_ns", "request_wall_ns", "sampled_worker_rss_bytes")}
            for field in ("cx_count", "gate_count", "depth", "two_q_depth"):
                values[field] = [row["validation"]["metrics"][field] for row in measured]
            cell["metrics"] = {field: {"median": statistics.median(items), "min": min(items),
                                      "max": max(items)} for field, items in values.items()}
        self.summary = {"counts": dict(counts), "unique_artifacts": len(keys),
                        "sessions": {"1": [cell]}, "scope": self.protocol["scope"],
                        "validation_limits": self.protocol["validation"]}

    def write_archive(self):
        protocol_path = self.directory / "protocol.json"
        protocol_path.write_text(json.dumps(self.protocol, sort_keys=True) + "\n")
        protocol_hash = digest(protocol_path.read_bytes())
        raw_path = self.directory / "session_1.jsonl"
        raw_path.write_text("".join(json.dumps(row) + "\n" for row in self.raw))
        (self.directory / "validated_session_1.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in self.saved))
        receipt = {"completed": True, "session": 1, "protocol_sha256": protocol_hash,
                   "rows_sha256": digest(raw_path.read_bytes()),
                   "scheduled_requests": len(self.protocol["schedules"]["1"])}
        (self.directory / "session_1_receipt.json").write_text(json.dumps(receipt))
        self.summary["protocol_sha256"] = protocol_hash
        (self.directory / "summary.json").write_text(json.dumps(self.summary))

    def run_audit(self):
        with archived_modules(self.harness):
            return audit(self.evidence, self.harness)

    def test_valid_archive_replays_complete_metrics(self):
        result = self.run_audit()
        self.assertTrue(result["passed"])
        self.assertEqual(result["rows_checked"], 4)
        self.assertEqual(result["unique_artifacts"], 1)
        cell = self.summary["sessions"]["1"][0]
        self.assertEqual(cell["metrics"]["compile_ns"], {"median": 200, "min": 100, "max": 300})
        self.assertEqual(cell["metrics"]["cx_count"], {"median": 1, "min": 1, "max": 1})

    def test_raw_file_change_breaks_original_receipt(self):
        with (self.directory / "session_1.jsonl").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "receipt"):
            self.run_audit()

    def test_saved_row_cannot_change_timing(self):
        self.saved[1]["compile_ns"] = 1
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "artifact"):
            self.run_audit()

    def test_summary_timing_median_is_recomputed(self):
        self.summary["sessions"]["1"][0]["metrics"]["compile_ns"]["median"] = 1
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "aggregates"):
            self.run_audit()

    def test_summary_quality_range_is_recomputed(self):
        self.summary["sessions"]["1"][0]["metrics"]["cx_count"]["max"] = 99
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "aggregates"):
            self.run_audit()

    def test_summary_eligibility_cannot_hide_a_cell(self):
        self.summary["sessions"]["1"][0]["eligible"] = False
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "aggregates"):
            self.run_audit()

    def test_scope_cannot_be_promoted_to_a_stronger_claim(self):
        self.summary["validation_limits"] = "formal proof on all possible programs"
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "scope|limits"):
            self.run_audit()

    def test_timed_row_cannot_be_reclassified_as_warmup(self):
        self.raw[1]["phase"] = self.saved[1]["phase"] = "warmup"
        self.refresh_summary()
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "Phase"):
            self.run_audit()

    def test_session_and_schedule_index_are_bound(self):
        for field, bad in (("session", 2), ("schedule_index", 2)):
            with self.subTest(field=field):
                before = self.raw[1][field]
                self.raw[1][field] = self.saved[1][field] = bad
                self.write_archive()
                with self.assertRaisesRegex(ValueError, "Schedule"):
                    self.run_audit()
                self.raw[1][field] = self.saved[1][field] = before

    def test_output_checksum_is_checked_after_receipt_binding(self):
        changed = self.raw[1]["out_qasm"].replace("cx q[0],q[1];", "")
        self.raw[1]["out_qasm"] = self.saved[1]["out_qasm"] = changed
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "Output checksum"):
            self.run_audit()

    def test_layout_is_part_of_the_validation_cache_binding(self):
        self.raw[1]["final_layout"] = self.saved[1]["final_layout"] = [1, 0]
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "validation replay"):
            self.run_audit()

    def test_forged_metric_in_saved_validation_is_rejected(self):
        self.saved[1]["validation"]["metrics"]["cx_count"] = 0
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "validation replay"):
            self.run_audit()

    def test_unverified_output_remains_visible_and_ineligible(self):
        row = self.raw[2]
        row["out_qasm"] = row["out_qasm"].replace("cx q[0],q[1];", "")
        row["output_sha256"] = digest(row["out_qasm"].encode())
        self.populate_validation()
        self.write_archive()
        result = self.run_audit()
        self.assertTrue(result["passed"], "Archive integrity can pass while an output is unverified")
        self.assertEqual(result["counts"]["timed:qrust:unverified_output"], 1)
        self.assertFalse(self.summary["sessions"]["1"][0]["eligible"])

    def test_failure_cannot_acquire_forged_passing_validation(self):
        keep = ("case_id", "tool", "trial", "seed", "session", "schedule_index", "phase")
        self.raw[2] = {key: self.raw[2][key] for key in keep}
        self.raw[2].update(status="timeout", error="synthetic resource limit")
        self.populate_validation()
        self.saved[2]["validation"] = {"passed": True}
        self.saved[2]["validation_cache_key"] = "fabricated"
        # Do not infer new aggregates from the forged validation: rejection
        # must happen before eligibility can dereference nonexistent timings.
        self.write_archive()
        with self.assertRaises(ValueError):
            self.run_audit()

    def test_source_hash_and_unique_filename_mapping_are_enforced(self):
        source = self.harness / "experiments" / "comparison_validation.py"
        with source.open("a") as stream:
            stream.write("\n# changed after protocol freeze\n")
        with self.assertRaisesRegex(ValueError, "source changed"):
            self.run_audit()

    def test_ambiguous_source_basenames_fail_closed(self):
        existing = next(value for key, value in self.protocol["source_sha256"].items()
                        if key.endswith("scalable_clifford.py"))
        self.protocol["source_sha256"]["/another/path/scalable_clifford.py"] = existing
        self.write_archive()
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            self.run_audit()

    def test_imported_oracle_must_come_from_archived_harness(self):
        with archived_modules(self.harness):
            oracle = importlib.import_module("comparison_validation")
            original = oracle.__file__
            try:
                oracle.__file__ = str(HERE / "comparison_validation.py")
                with self.assertRaisesRegex(ValueError, "oracle import"):
                    audit(self.evidence, self.harness)
            finally:
                oracle.__file__ = original


if __name__ == "__main__":
    unittest.main()
