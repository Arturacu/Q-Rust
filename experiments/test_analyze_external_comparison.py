"""External-profile integrity tests with synthetic data, before measured results."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HAS_NUMPY = importlib.util.find_spec("numpy") is not None
if HAS_NUMPY:
    from analyze_external_comparison import analyze_profile, project_rows, recipe_sensitivity, session_agreement, suite_of, digest, verify_session_files
    from test_analyze_comparison_study import synthetic


def external_synthetic():
    protocol, timed, validated = synthetic()
    protocol["tools"] = ["qrust", "qiskit", "cirq", "tket", "cirq_one_pass", "tket_light"]
    for case in protocol["cases"]:
        case["suite"] = "qasmbench" if case["family"] == "a" else "mqt"
    for rows in (timed, validated):
        extra = []
        for row in rows:
            row["artifact_sha256"] = row["case_id"] + row["tool"]
            if row["tool"] in ("cirq", "tket"):
                other = copy.deepcopy(row)
                other["tool"] = "cirq_one_pass" if row["tool"] == "cirq" else "tket_light"
                other["id"] = f"{row['case_id']}__{other['tool']}__{row['trial']}"
                other["compile_ns"] //= 2
                extra.append(other)
        rows.extend(extra)
    return protocol, timed, validated


@unittest.skipUnless(HAS_NUMPY, "requires NumPy from the pinned comparison environment")
class ExternalAnalysisTests(unittest.TestCase):
    def test_changed_failure_status_or_warmup_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in ("timed", "warmup", "validated"):
                (directory / f"{name}.jsonl").write_text('{"status":"error"}\n')
            info = {"completed": True, "session": 1, "protocol_sha256": "frozen",
                    **{name + "_sha256": digest(directory / f"{name}.jsonl") for name in ("timed", "warmup")}}
            receipt = {"validated_sha256": digest(directory / "validated.jsonl")}
            (directory / "session.json").write_text(json.dumps(info))
            verify_session_files(directory, "frozen", 1, receipt)
            for name in ("timed", "warmup"):
                with self.subTest(name=name):
                    path = directory / f"{name}.jsonl"
                    original = path.read_text()
                    path.write_text('{"status":"compiled"}\n')
                    with self.assertRaisesRegex(ValueError, "changed after timing"):
                        verify_session_files(directory, "frozen", 1, receipt)
                    path.write_text(original)

    def test_real_community_metadata_and_fixture_names_are_recognized(self):
        for source, fixture, expected in (("MQT Bench", "mqtbench__ae_n5", "mqt"),
                                           ("QASMBench", "qasmbench__adder_n4", "qasmbench")):
            self.assertEqual(suite_of({"community_source": source, "fixture": fixture}), expected)
            self.assertEqual(suite_of({"fixture": fixture}), expected)

    def test_profiles_share_observations_without_mutating_or_duplicating_them(self):
        _, rows, _ = external_synthetic()
        snapshot = copy.deepcopy(rows)
        primary, light = project_rows(rows, "primary"), project_rows(rows, "light")
        self.assertEqual(len(primary), 84)
        self.assertEqual(len(light), 84)
        self.assertEqual([x for x in primary if x["tool"] == "qrust"],
                         [x for x in light if x["tool"] == "qrust"])
        self.assertEqual([x for x in primary if x["tool"] == "qiskit"],
                         [x for x in light if x["tool"] == "qiskit"])
        self.assertEqual(rows, snapshot)
        self.assertTrue(all(x["raw_tool"] == "cirq_one_pass" for x in light if x["tool"] == "cirq"))

    def test_failed_light_repeat_stays_in_denominator_and_out_of_paired_metrics(self):
        protocol, timed, validated = external_synthetic()
        row = next(x for x in validated if x["tool"] == "cirq_one_pass")
        row["validation"]["passed"] = False
        primary, _ = analyze_profile(protocol, timed, validated, "primary", 10)
        light, cells = analyze_profile(protocol, timed, validated, "light", 10)
        self.assertEqual(primary["coverage"]["excluded_cell_count"], 0)
        self.assertEqual(light["coverage"]["excluded_cell_count"], 1)
        self.assertEqual(light["coverage"]["planned_cases"], 3)
        self.assertEqual(light["pairwise"]["aggregate"]["cirq"]["paired_cells"], 2)
        bad = next(c for c in cells if not c["eligible"])
        self.assertEqual(bad["passed_trials"], 6)
        self.assertIsNone(bad["medians"])

    def test_sensitivity_uses_paired_within_session_values_and_retains_suite_strata(self):
        protocol, timed, validated = external_synthetic()
        primary, a = analyze_profile(protocol, timed, validated, "primary", 10)
        _, b = analyze_profile(protocol, timed, validated, "light", 10)
        sensitivity = recipe_sensitivity(a, b, 10)
        self.assertEqual(set(primary["pairwise"]["by_suite"]), {"qasmbench", "mqt"})
        for tool in ("cirq", "tket"):
            self.assertEqual(sensitivity[tool]["paired_cells"], 3)
            self.assertAlmostEqual(sensitivity[tool]["compile_ns"]["geomean"], 2)
            self.assertEqual(sensitivity[tool]["cx_count"]["geomean"], 1)

    def test_session_grid_mismatch_is_not_silently_pooled(self):
        protocol, _, validated = external_synthetic()
        with self.assertRaisesRegex(ValueError, "Session grids differ"):
            session_agreement(protocol, validated, validated[1:], {}, 10)

    def test_session_quality_changes_are_reported_separately(self):
        protocol, _, validated = external_synthetic()
        second = copy.deepcopy(validated)
        second[0]["validation"]["metrics"]["cx_count"] += 1
        result = session_agreement(protocol, validated, second, {}, 10)
        self.assertEqual(result["counts"]["quality_changed"], 1)
        self.assertEqual(len(result["changed_trials"]), 1)
        self.assertEqual(result["same_seed_comparisons"], 126)


if __name__ == "__main__":
    unittest.main()
