"""Fixed-cohort, integrity and adverse-case regressions for routing analysis."""
import copy
import importlib.util
import json
import unittest

HAS_NUMPY = importlib.util.find_spec("numpy") is not None
if HAS_NUMPY:
    from analyze_external_routing import (
        ALLOWED_COMPILER_CHANGES, FOCAL_FIXTURES, FOCAL_TOPOLOGIES,
        assert_comparable, compare, focal_case_ids, markdown, ratio_record,
    )
    from test_analyze_external_followup import synthetic_pair


def routing_pair():
    old, _, source_sessions, _ = synthetic_pair()
    source_cases = old["cases"]
    cases, case_sources = [], {}
    for index, fixture in enumerate(FOCAL_FIXTURES):
        for topology in FOCAL_TOPOLOGIES:
            case = copy.deepcopy(source_cases[index])
            case.update(id=fixture + "__" + topology, fixture=fixture,
                        topology=topology, n_logical=10, physical_width=10)
            cases.append(case)
            case_sources[case["id"]] = source_cases[index]["id"]
    old["cases"] = cases
    old["compiler_provenance"]["compiler_source_sha256"] = {
        "src/transpiler/layout.rs": "old", "src/ir/gates.rs": "unchanged"}
    new = copy.deepcopy(old)
    new["compiler_provenance"]["compiler_source_sha256"]["src/transpiler/layout.rs"] = "new"
    new["qrust_binary_sha256"] = "new"
    sessions = {}
    for session, pair in source_sessions.items():
        rebuilt = []
        for rows in pair:
            values = []
            for case_id, source_id in case_sources.items():
                for original in rows:
                    if original["case_id"] != source_id:
                        continue
                    row = copy.deepcopy(original)
                    row["case_id"] = case_id
                    row["id"] = f"{case_id}__{row['tool']}__{row['trial']}"
                    values.append(row)
            rebuilt.append(values)
        sessions[session] = tuple(rebuilt)
    return old, new, sessions, copy.deepcopy(sessions)


@unittest.skipUnless(HAS_NUMPY, "requires pinned NumPy comparison environment")
class RoutingAnalysisTests(unittest.TestCase):
    def test_fixed_cohort_profiles_and_focal_cases_are_complete_and_renderer_compatible(self):
        old, new, before, after = routing_pair()
        snapshot = copy.deepcopy((old, new, before, after))
        result = compare(old, new, before, after, draws=2, expected_cases=6)
        self.assertEqual((old, new, before, after), snapshot)
        self.assertEqual(result["cohorts"]["original_common_case_ids"], result["cohorts"]["followup_common_case_ids"])
        self.assertEqual(len(result["cohorts"]["focal_case_ids"]), 6)
        for session in result["sessions"].values():
            self.assertEqual(set(session["case_changes"]), set(old["tools"]))
            self.assertEqual(len(session["focal_cases"]), 6)
            for details in session["case_changes"].values():
                self.assertEqual(details["paired_cases"], 6)
                self.assertEqual(details["quality_adverse_cases"], [])
                self.assertEqual(details["outcomes"]["cx_count"], {"improved": 0, "unchanged": 6, "regressed": 0})
            for profile in session["profiles"].values():
                for key in ("original_on_original_common", "followup_on_original_common", "followup_on_all_common",
                            "original_on_focal", "followup_on_focal"):
                    self.assertEqual(profile[key]["complete_all_four_cases"], 6)
            for recipe in old["tools"]:
                self.assertEqual(session["recipes"][recipe]["paired_complete_cases"], 6)
                self.assertEqual(session["recipes_on_focal"][recipe]["paired_complete_cases"], 6)
        self.assertIn("CX depth", markdown(result))
        json.dumps(result, allow_nan=False)
        with self.assertRaisesRegex(ValueError, "Expected 100"):
            compare(old, new, before, after, draws=2)

    def test_frozen_harness_and_settings_cannot_change_even_for_previous_cirq_fix(self):
        old, new, _, _ = routing_pair()
        for name in ("experiments/comparison_cirq.py", "oracle.py"):
            changed = copy.deepcopy(new)
            changed["source_sha256"][name] = "changed"
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Unexpected source"):
                assert_comparable(old, changed)
        for field in ("cases", "schedules", "master_seed", "recipes", "corpus_manifest_sha256"):
            changed = copy.deepcopy(new)
            changed[field] = "changed"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "settings differ"):
                assert_comparable(old, changed)
        new["compiler_provenance"]["compiler_source_sha256"]["src/ir/gates.rs"] = "changed"
        with self.assertRaisesRegex(ValueError, "Unexpected source"):
            assert_comparable(old, new)

    def test_explicit_compiler_and_regression_paths_are_allowed(self):
        old, new, _, _ = routing_pair()
        for name in ALLOWED_COMPILER_CHANGES:
            new["compiler_provenance"]["compiler_source_sha256"][name] = "new"
        result = assert_comparable(old, new)
        self.assertEqual(result["compiler_changed"], sorted(ALLOWED_COMPILER_CHANGES))
        self.assertEqual(result["harness_changed"], [])

    def test_failed_light_recipe_or_missing_trial_cannot_shrink_the_cohort(self):
        old, new, before, after = routing_pair()
        for recipe in ("qrust", "cirq_one_pass", "tket_light"):
            changed = copy.deepcopy(after)
            row = next(r for r in changed[2][1] if r["tool"] == recipe)
            row["validation"]["passed"] = False
            with self.subTest(recipe=recipe), self.assertRaisesRegex(ValueError, "Full cohort lost"):
                compare(old, new, before, changed, draws=2, expected_cases=6)
        for rows in after[1]:
            rows.pop()
        with self.assertRaisesRegex(ValueError, "Full cohort lost"):
            compare(old, new, before, after, draws=2, expected_cases=6)

    def test_zero_cost_regression_and_different_depth_outcomes_remain_visible(self):
        old, new, before, after = routing_pair()
        case = old["cases"][0]["id"]
        for source, cx, depth, cx_depth in ((before, 0, 10, 4), (after, 2, 8, 5)):
            for pair in source.values():
                for row in pair[1]:
                    if row["case_id"] == case and row["tool"] == "qrust":
                        row["validation"]["metrics"].update(cx_count=cx, depth=depth, two_q_depth=cx_depth)
        result = compare(old, new, before, after, draws=2, expected_cases=6)
        for session in result["sessions"].values():
            metrics = session["focal_cases"][case]["qrust"]["metrics"]
            self.assertEqual(metrics["cx_count"]["ratio_class"], "denominator_zero")
            self.assertIsNone(metrics["cx_count"]["ratio"])
            self.assertEqual(metrics["depth"]["outcome"], "improved")
            self.assertEqual(metrics["two_q_depth"]["outcome"], "regressed")
            self.assertEqual([c["case_id"] for c in session["case_changes"]["qrust"]["quality_adverse_cases"]], [case])
            self.assertEqual(session["recipes"]["qrust"]["ratios"]["cx_count"]["zero_cases"]["denominator_zero"], 1)
        json.dumps(result, allow_nan=False)
        self.assertEqual(ratio_record(0, 0)["ratio_class"], "both_zero")
        self.assertEqual(ratio_record(2, 0)["ratio_class"], "numerator_zero")

    def test_sessions_seeds_validation_binding_and_focal_definitions_are_enforced(self):
        old, new, before, after = routing_pair()
        with self.assertRaisesRegex(ValueError, "two separate"):
            compare(old, new, {1: before[1]}, after, draws=2, expected_cases=6)
        altered = copy.deepcopy(after)
        for rows in altered[1]:
            rows[0]["seed"] += 1
        with self.assertRaisesRegex(ValueError, "seed or session"):
            compare(old, new, before, altered, draws=2, expected_cases=6)
        altered = copy.deepcopy(after)
        altered[1][1][0]["compile_ns"] += 1
        with self.assertRaisesRegex(ValueError, "measured field"):
            compare(old, new, before, altered, draws=2, expected_cases=6)
        old["cases"][0]["topology"] = "linear"
        with self.assertRaisesRegex(ValueError, "six predeclared"):
            focal_case_ids(old)

    def test_timing_changes_keep_sessions_separate_and_retain_unchanged_controls(self):
        old, new, before, after = routing_pair()
        for session, pair in after.items():
            for rows in pair:
                for row in rows:
                    if row["tool"] == "qrust":
                        row["compile_ns"] *= session + 1
                        row["request_wall_ns"] *= session + 2
        result = compare(old, new, before, after, draws=2, expected_cases=6)
        for session, details in result["sessions"].items():
            self.assertAlmostEqual(details["recipes"]["qrust"]["ratios"]["compile_ns"]["geomean"], int(session) + 1)
            self.assertAlmostEqual(details["recipes"]["qrust"]["ratios"]["request_wall_ns"]["geomean"], int(session) + 2)
            self.assertEqual(set(details["unchanged_recipe_drift_controls"]), set(old["tools"]) - {"qrust"})
            for control in details["unchanged_recipe_drift_controls"]:
                self.assertEqual(details["recipes"][control]["ratios"]["compile_ns"]["geomean"], 1)


if __name__ == "__main__":
    unittest.main()
