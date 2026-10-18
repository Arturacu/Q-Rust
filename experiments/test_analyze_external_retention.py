"""Synthetic retention-study integrity, fixed denominators and target regressions."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HAS_NUMPY = importlib.util.find_spec("numpy") is not None
if HAS_NUMPY:
    from analyze_external_retention import (
        ALLOWED_COMPILER_CHANGES, FOCAL_WIDTHS, assert_comparable, check_preservation,
        compare, focal_case_ids, markdown,
    )
    from test_analyze_external_followup import synthetic_pair


def retention_pair():
    old, new, before, _ = synthetic_pair()
    old["source_sha256"] = {"experiments/comparison_cirq.py": "fixed", "oracle.py": "fixed"}
    old["compiler_provenance"] = {"compiler_source_sha256": {
        "src/transpiler/layout.rs": "old", "src/ir/gates.rs": "fixed"}}
    templates = copy.deepcopy(before)
    source_id = old["cases"][0]["id"]
    cases = []
    for fixture, width in FOCAL_WIDTHS.items():
        for topology in ("all_to_all", "linear"):
            case = copy.deepcopy(old["cases"][0])
            case.update(id=f"{fixture}__{topology}", fixture=fixture, n_logical=width,
                        physical_width=width, topology=topology,
                        family=fixture.split("__")[1].rsplit("_n", 1)[0],
                        suite="mqt" if fixture.startswith("mqt") else "qasmbench",
                        source_sha256=fixture + "source", normalized_sha256=fixture + "normalized")
            cases.append(case)
    old["cases"] = cases
    for session in (1, 2):
        before[session] = ([], [])
        for index, rows in enumerate(templates[session]):
            for case in cases:
                for source in rows:
                    if source["case_id"] != source_id:
                        continue
                    row = copy.deepcopy(source)
                    row.update(case_id=case["id"], id=f"{case['id']}__{source['tool']}__{source['trial']}",
                               artifact_sha256=case["id"] + source["tool"], output_sha256="same-hash",
                               initial_layout=list(range(case["n_logical"])),
                               final_layout=list(range(case["n_logical"])), output_qubits=case["n_logical"])
                    before[session][index].append(row)
    new = copy.deepcopy(old)
    new["compiler_provenance"]["compiler_source_sha256"]["src/transpiler/candidate_selection.rs"] = "new"
    new["qrust_binary_sha256"] = "new"
    after = copy.deepcopy(before)
    plan = {
        "baseline_protocol_sha256": "old-protocol", "baseline_binary_sha256": "old",
        "corpus_manifest_sha256": old["corpus_manifest_sha256"],
        "focal_cases": [{key: case[key] for key in ("id", "fixture", "n_logical", "topology",
                                                  "source_sha256", "normalized_sha256")} for case in cases],
        "predeclared_quality_ceilings": [],
    }
    for case in cases:
        recover = case["topology"] == "linear" and case["fixture"].startswith("qasmbench")
        metrics = {"cx_count": 9 if recover else 10, "depth": 9 if recover else 10,
                   "two_q_depth": 4 if recover else 5}
        plan["predeclared_quality_ceilings"].append({
            "case_id": case["id"], "purpose": "recover_previous_quality" if recover else "retain_current_quality",
            "reference_study": "compatibility" if recover else "routing", "reference_metrics": metrics,
        })
        if recover:
            for pair in after.values():
                for rows in pair:
                    for row in rows:
                        if row["case_id"] == case["id"] and row["tool"] == "qrust":
                            row["artifact_sha256"] += "recovered"
                            row["out_qasm"] = "recovered output"
                            if "validation" in row:
                                row["validation"]["metrics"].update(metrics)
    return old, new, before, after, plan


def analyze(pair):
    old, new, before, after, plan = pair
    return compare(old, new, before, after, plan=plan, baseline_protocol_sha256="old-protocol",
                   draws=2, expected_cases=12)


@unittest.skipUnless(HAS_NUMPY, "requires the pinned NumPy comparison environment")
class RetentionAnalysisTests(unittest.TestCase):
    def test_complete_grid_targets_identities_and_renderer_shape(self):
        pair = retention_pair()
        snapshot = copy.deepcopy(pair)
        result = analyze(pair)
        self.assertEqual(pair, snapshot)
        self.assertEqual(result["analysis"], "routing_retention_followup")
        self.assertEqual(result["non_qrust_identity"]["paired_trials"], 840)
        self.assertTrue(result["non_qrust_identity"]["all_identical"])
        for checks in result["same_seed_session_identity"].values():
            self.assertEqual(checks["identical_trials"], 504)
        acceptance = result["quality_acceptance"]
        self.assertTrue(acceptance["all_focal_targets_passed"])
        self.assertTrue(acceptance["full_cohort_nonworsening_passed"])
        self.assertEqual(acceptance["focal_session_cells"], 24)
        self.assertEqual(acceptance["focal_trials_checked"], 168)
        for info in result["sessions"].values():
            self.assertEqual(len(info["focal_cases"]), 12)
            self.assertEqual(len(info["case_changes"]["qrust"]["cells"]), 12)
            self.assertEqual(info["case_changes"]["qrust"]["outcomes"]["cx_count"],
                             {"improved": 3, "unchanged": 9, "regressed": 0})
            for profile in info["profiles"].values():
                self.assertEqual(profile["followup_on_original_common"]["pairwise"]["qiskit"]["paired_cells"], 12)
        self.assertIn("All twelve", markdown(result))
        json.dumps(result, allow_nan=False)

    def test_every_protocol_change_and_harness_edit_is_rejected(self):
        old, new, *_ = retention_pair()
        for key in ("corpus_manifest_sha256", "schedules", "recipes", "master_seed", "trials", "cases"):
            changed = copy.deepcopy(new)
            changed[key] = "changed"
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "settings differ"):
                assert_comparable(old, changed)
        new["source_sha256"]["oracle.py"] = "changed"
        with self.assertRaisesRegex(ValueError, "Unexpected source"):
            assert_comparable(old, new)

    def test_compiler_source_allowlist_is_explicit(self):
        old, new, *_ = retention_pair()
        for path in ALLOWED_COMPILER_CHANGES:
            new["compiler_provenance"]["compiler_source_sha256"][path] = "new"
        self.assertEqual(set(assert_comparable(old, new)["compiler_changed"]), ALLOWED_COMPILER_CHANGES)
        new["compiler_provenance"]["compiler_source_sha256"]["src/ir/gates.rs"] = "changed"
        with self.assertRaisesRegex(ValueError, "Unexpected source"):
            assert_comparable(old, new)

    def test_missing_duplicate_extra_or_failed_trial_does_not_shrink_cohort(self):
        for mode in ("missing", "duplicate", "extra", "failed"):
            pair = retention_pair()
            for rows in pair[3][1]:
                if mode == "missing":
                    rows.pop()
                elif mode in ("duplicate", "extra"):
                    row = copy.deepcopy(rows[-1])
                    if mode == "extra":
                        row["trial"] = 99
                        row["id"] += "extra"
                    rows.append(row)
            if mode == "failed":
                pair[3][1][1][0]["validation"]["passed"] = False
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                analyze(pair)

    def test_baseline_output_layout_hash_width_and_metrics_changes_are_rejected(self):
        mutations = {"out_qasm": "different", "initial_layout": [1, 0], "final_layout": [1, 0],
                     "artifact_sha256": "different", "output_sha256": "different", "output_qubits": 20}
        for field, value in mutations.items():
            pair = retention_pair()
            for rows in pair[3][1]:
                next(row for row in rows if row["tool"] == "qiskit")[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "compiler baseline changed"):
                analyze(pair)
        pair = retention_pair()
        next(row for row in pair[3][1][1] if row["tool"] == "cirq")["validation"]["metrics"]["depth"] += 1
        with self.assertRaisesRegex(ValueError, "compiler baseline changed"):
            analyze(pair)

    def test_plan_hash_metadata_focal_source_and_reference_mismatches_are_rejected(self):
        for field in ("baseline_protocol_sha256", "baseline_binary_sha256", "corpus_manifest_sha256"):
            pair = retention_pair()
            pair[4][field] = "changed"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "plan differs"):
                analyze(pair)
        pair = retention_pair()
        pair[4]["focal_cases"][0]["source_sha256"] = "changed"
        with self.assertRaisesRegex(ValueError, "focal source"):
            analyze(pair)
        pair = retention_pair()
        pair[4]["predeclared_quality_ceilings"][0]["reference_study"] = "compatibility"
        with self.assertRaisesRegex(ValueError, "ceiling reference"):
            analyze(pair)

    def test_missing_duplicate_or_altered_ceilings_are_rejected(self):
        for mode in ("missing", "duplicate", "metric", "retention", "recovery"):
            pair = retention_pair()
            caps = pair[4]["predeclared_quality_ceilings"]
            if mode == "missing":
                caps.pop()
            elif mode == "duplicate":
                caps[-1] = copy.deepcopy(caps[0])
            elif mode == "metric":
                caps[0]["reference_metrics"]["gate_count"] = 20
            elif mode == "retention":
                caps[0]["reference_metrics"]["depth"] = 100
            else:
                caps[-1]["reference_metrics"] = {"cx_count": 10, "depth": 10, "two_q_depth": 5}
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                analyze(pair)

    def test_one_bad_trial_cannot_hide_behind_successful_medians(self):
        pair = retention_pair()
        for rows in (pair[3][1][1], pair[3][2][1]):
            row = next(row for row in rows if row["tool"] == "qrust" and row["trial"] == 0)
            row["validation"]["metrics"]["two_q_depth"] += 1
        result = analyze(pair)
        acceptance = result["quality_acceptance"]
        self.assertFalse(acceptance["all_focal_targets_passed"])
        self.assertFalse(acceptance["full_cohort_nonworsening_passed"])
        self.assertEqual(len(acceptance["full_cohort_regressed_trials"]), 2)
        self.assertEqual(sum(not check["passed"] for check in acceptance["focal_targets"]), 2)
        self.assertEqual(result["sessions"]["1"]["case_changes"]["qrust"]["outcomes"]["two_q_depth"]["regressed"], 0)

    def test_regressions_and_zero_cost_transitions_are_retained(self):
        pair = retention_pair()
        case_id = pair[0]["cases"][0]["id"]
        pair[4]["predeclared_quality_ceilings"][0]["reference_metrics"]["cx_count"] = 0
        for sessions, value in ((pair[2], 0), (pair[3], 1)):
            for _, rows in sessions.values():
                for row in rows:
                    if row["case_id"] == case_id and row["tool"] == "qrust":
                        row["validation"]["metrics"]["cx_count"] = value
        result = analyze(pair)
        cell = result["sessions"]["1"]["focal_cases"][case_id]["qrust"]
        self.assertEqual(cell["metrics"]["cx_count"]["ratio_class"], "denominator_zero")
        self.assertIsNone(cell["metrics"]["cx_count"]["ratio"])
        self.assertEqual(result["sessions"]["1"]["case_changes"]["qrust"]["quality_adverse_cases"], [cell])
        json.dumps(result, allow_nan=False)

    def test_timings_drift_independently_of_identity_and_remain_separate(self):
        pair = retention_pair()
        for session, factor in ((1, 2), (2, 3)):
            for rows in pair[3][session]:
                for row in rows:
                    row["compile_ns"] *= factor
                    row["request_wall_ns"] *= factor
        result = analyze(pair)
        self.assertTrue(result["non_qrust_identity"]["all_identical"])
        for session, factor in (("1", 2), ("2", 3)):
            for info in result["sessions"][session]["recipes"].values():
                self.assertAlmostEqual(info["ratios"]["compile_ns"]["geomean"], factor)

    def test_cross_session_identity_includes_maps_and_retains_disagreements(self):
        pair = retention_pair()
        for rows in pair[3][2]:
            row = next(row for row in rows if row["tool"] == "qrust")
            row["final_layout"] = list(reversed(row["final_layout"]))
        result = analyze(pair)
        checks = result["same_seed_session_identity"]["followup"]
        self.assertFalse(checks["all_identical"])
        self.assertEqual(checks["identical_trials"], 503)
        self.assertEqual(checks["changed_trials"][0]["changed_fields"], ["final_layout"])

    def test_changed_seed_measured_field_session_or_focal_width_is_rejected(self):
        for field, value in (("seed", -1), ("session", 3)):
            pair = retention_pair()
            for rows in pair[3][1]:
                rows[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                analyze(pair)
        pair = retention_pair()
        pair[3][1][1][0]["request_wall_ns"] += 1
        with self.assertRaisesRegex(ValueError, "measured field"):
            analyze(pair)
        pair = retention_pair()
        pair[0]["cases"][0]["n_logical"] = 9
        with self.assertRaisesRegex(ValueError, "widths"):
            focal_case_ids(pair[0])

    def test_preservation_receipt_cannot_be_replaced(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "receipt.json"
            path.write_text('{"source_commit":"changed","sha256":{}}')
            with self.assertRaisesRegex(ValueError, "prespecified receipt"):
                check_preservation(Path(temporary), path)


if __name__ == "__main__":
    unittest.main()
