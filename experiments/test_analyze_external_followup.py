"""Synthetic before/after integrity and cohort-denominator regressions."""
import copy
import importlib.util
import unittest

HAS_NUMPY = importlib.util.find_spec("numpy") is not None
if HAS_NUMPY:
    from analyze_external_followup import assert_comparable, bind_validated_rows, compare, trial_changes
    from test_analyze_external_comparison import external_synthetic


def synthetic_pair():
    protocol, timed, checked = external_synthetic()
    protocol.update(sessions=2, master_seed=20, corpus_manifest_sha256="same-corpus", schedules={"1": "same", "2": "same"},
                    source_sha256={"experiments/comparison_cirq.py": "old", "oracle.py": "fixed"},
                    compiler_provenance={"compiler_source_sha256": {"src/ir/gates.rs": "old", "src/ir/gate_def.rs": "old"}},
                    qrust_binary_sha256="old", recipes={"qiskit": "fixed"})
    for rows in (timed, checked):
        for row in rows:
            row.update(seed=20 + row["trial"], session=1, request_wall_ns=100 + row["compile_ns"] + row["parse_ns"])
    sessions = {s: (copy.deepcopy(timed), copy.deepcopy(checked)) for s in (1, 2)}
    for s, pair in sessions.items():
        for rows in pair:
            for row in rows:
                row["session"] = s
    after_protocol = copy.deepcopy(protocol)
    after_protocol["source_sha256"]["experiments/comparison_cirq.py"] = "new"
    after_protocol["compiler_provenance"]["compiler_source_sha256"]["src/ir/gates.rs"] = "new"
    after_protocol["qrust_binary_sha256"] = "new"
    return protocol, after_protocol, sessions, copy.deepcopy(sessions)


@unittest.skipUnless(HAS_NUMPY, "requires pinned NumPy comparison environment")
class FollowupAnalysisTests(unittest.TestCase):
    def test_original_failure_excluded_without_shrinking_followup_cohort(self):
        old, new, before, after = synthetic_pair()
        bad = old["cases"][0]["id"]
        for pair in before.values():
            for row in pair[1]:
                if row["case_id"] == bad and row["tool"] == "qrust" and row["trial"] == 0:
                    row["validation"]["passed"] = False
        result = compare(old, new, before, after, draws=2)
        self.assertEqual(len(result["cohorts"]["original_common_case_ids"]), 2)
        self.assertEqual(len(result["cohorts"]["followup_common_case_ids"]), 3)
        self.assertEqual(result["cohorts"]["newly_common_case_ids"], [bad])
        for info in result["sessions"].values():
            self.assertEqual(info["recipes"]["qrust"]["before_eligible_cases"], 2)
            self.assertEqual(info["recipes"]["qrust"]["after_eligible_cases"], 3)
            self.assertEqual(info["recipes"]["qrust"]["paired_complete_cases"], 2)
            self.assertEqual(info["recipes"]["qiskit"]["paired_complete_cases"], 3)
            for stats in info["recipes_on_original_common"].values():
                self.assertEqual(stats["planned_cases"], 2)
                self.assertEqual(stats["paired_complete_cases"], 2)
            for profile in info["profiles"].values():
                self.assertEqual(profile["followup_on_original_common"]["pairwise"]["qiskit"]["paired_cells"], 2)
                self.assertEqual(profile["followup_on_all_common"]["pairwise"]["qiskit"]["paired_cells"], 3)
                self.assertEqual(profile["followup_on_newly_common"]["pairwise"]["qiskit"]["paired_cells"], 1)

    def test_same_seed_identity_pairing_is_required_even_when_rows_reordered(self):
        _, _, before, after = synthetic_pair()
        result = trial_changes(before[1][1], list(reversed(after[1][1])))
        self.assertEqual(result["by_recipe"]["qrust"]["paired_trials"], 21)
        after[1][1][0]["seed"] += 1
        with self.assertRaisesRegex(ValueError, "seed or session"):
            trial_changes(before[1][1], after[1][1])

    def test_corpus_schedule_recipe_and_unexpected_source_mismatches_rejected(self):
        old, new, _, _ = synthetic_pair()
        for key in ("corpus_manifest_sha256", "schedules", "master_seed", "recipes"):
            changed = copy.deepcopy(new)
            changed[key] = "changed"
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "settings differ"):
                assert_comparable(old, changed)
        new["source_sha256"]["oracle.py"] = "changed"
        with self.assertRaisesRegex(ValueError, "Unexpected source"):
            assert_comparable(old, new)

    def test_allowed_test_addition_and_cu1_changes_are_explicit(self):
        old, new, _, _ = synthetic_pair()
        new["compiler_provenance"]["compiler_source_sha256"]["tests/cu1_regressions.rs"] = "new"
        result = assert_comparable(old, new)
        self.assertEqual(result["compiler_changed"], ["src/ir/gates.rs", "tests/cu1_regressions.rs"])

    def test_timing_ratios_use_all_cell_trials_and_quality_changes_remain_visible(self):
        old, new, before, after = synthetic_pair()
        for pair in after.values():
            for rows in pair:
                for row in rows:
                    row["compile_ns"] *= 2
                    row["request_wall_ns"] *= 3
            row = pair[1][0]
            row["validation"]["metrics"]["cx_count"] += 1
            row["artifact_sha256"] += "changed"
        result = compare(old, new, before, after, draws=2)
        for info in result["sessions"].values():
            for stats in info["recipes"].values():
                self.assertAlmostEqual(stats["ratios"]["compile_ns"]["geomean"], 2)
                self.assertAlmostEqual(stats["ratios"]["request_wall_ns"]["geomean"], 3)
            self.assertEqual(info["trial_changes"]["by_recipe"]["qrust"]["quality_changed"], 1)
            self.assertEqual(info["trial_changes"]["by_recipe"]["qrust"]["artifact_changed"], 1)

    def test_validation_cannot_rewrite_error_status_request_time_or_seed(self):
        _, _, before, _ = synthetic_pair()
        timed, checked = before[1]
        bind_validated_rows(timed, checked)
        for key, value in (("status", "error"), ("request_wall_ns", 1), ("seed", -1)):
            changed = copy.deepcopy(checked)
            changed[0][key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "measured field"):
                bind_validated_rows(timed, changed)


if __name__ == "__main__":
    unittest.main()
