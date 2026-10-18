"""Synthetic-data tests; never use study outcomes to tune the analysis."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HAS_NUMPY = importlib.util.find_spec("numpy") is not None
if HAS_NUMPY:
    from analyze_comparison_study import (
        TOOLS, analyze, build_cells, family_balanced_ratio_summary,
        geometric_ratio_summary, narrative, write_csv,
    )


def synthetic():
    cases = [
        {"id": "fixture_a__linear", "fixture": "fixture_a", "family": "a", "topology": "linear", "n_logical": 3, "physical_width": 3},
        {"id": "fixture_a__all_to_all", "fixture": "fixture_a", "family": "a", "topology": "all_to_all", "n_logical": 3, "physical_width": 3},
        {"id": "fixture_b__linear", "fixture": "fixture_b", "family": "b", "topology": "linear", "n_logical": 3, "physical_width": 3},
    ]
    timed, validated = [], []
    for case in cases:
        for tool in TOOLS:
            scale = {"qrust": 1, "qiskit": 2, "cirq": 4, "tket": 8}[tool]
            for trial in range(7):
                row = {"id": f"{case['id']}__{tool}__{trial}", "case_id": case["id"], "tool": tool,
                       "trial": trial, "status": "compiled", "compile_ns": (100 + trial) * scale,
                       "parse_ns": 10 * scale, "out_qasm": "same output",
                       "initial_layout": [0, 1, 2], "final_layout": [0, 1, 2], "output_qubits": 3}
                timed.append(row)
                validated.append({**row, "validation": {"passed": True, "evidence_class": "synthetic_exact",
                                  "metrics": {"cx_count": 10, "gate_count": 20, "depth": 10, "two_q_depth": 5}}})
    return {"cases": cases, "trials": 7}, timed, validated


@unittest.skipUnless(HAS_NUMPY, "requires NumPy from the pinned comparison environment")
class StudyAnalysisTests(unittest.TestCase):
    def test_constant_ratios_and_complete_coverage(self):
        protocol, timed, validated = synthetic()
        summary, cells = analyze(protocol, timed, validated, draws=100)
        self.assertEqual(summary["coverage"]["excluded_cell_count"], 0)
        self.assertEqual(len(cells), 12)
        for tool, ratio in (("qiskit", 2), ("cirq", 4), ("tket", 8)):
            stat = summary["pairwise"]["aggregate"][tool]
            self.assertAlmostEqual(stat["compile_ns"]["geomean"], ratio)
            self.assertAlmostEqual(stat["parse_compile_ns"]["geomean"], ratio)
            self.assertEqual(stat["compile_ns"]["qrust_win_tie_loss"], {"win": 3, "tie": 0, "loss": 0})
            self.assertEqual(stat["cx_count"]["qrust_win_tie_loss"], {"win": 0, "tie": 3, "loss": 0})
        self.assertEqual(summary["complete_case_quality"]["pareto_counts"]["speed_cx"], {"qrust": 3, "qiskit": 0, "cirq": 0, "tket": 0})
        json.dumps(summary, allow_nan=False)

    def test_one_failed_repeat_excludes_whole_cell(self):
        protocol, timed, validated = synthetic()
        validated[0]["validation"]["passed"] = False
        summary, cells = analyze(protocol, timed, validated, draws=20)
        qrust = next(cell for cell in cells if cell["case_id"] == protocol["cases"][0]["id"] and cell["tool"] == "qrust")
        self.assertFalse(qrust["eligible"])
        self.assertIsNone(qrust["medians"])
        self.assertEqual(qrust["passed_trials"], 6)
        self.assertEqual(len(qrust["failed_trials"]), 1)
        self.assertEqual(summary["pairwise"]["aggregate"]["qiskit"]["paired_cells"], 2)
        self.assertEqual(summary["complete_case_quality"]["complete_all_four_cells"], 2)

    def test_missing_trial_and_missing_validation_remain_visible(self):
        protocol, timed, validated = synthetic()
        missing_id = timed[0]["id"]
        timed = timed[1:]
        validated = [row for row in validated if row["id"] not in {missing_id, timed[8]["id"]}]
        summary, _ = analyze(protocol, timed, validated, draws=20)
        self.assertEqual(summary["coverage"]["excluded_cell_count"], 2)
        self.assertTrue(any("received 6" in str(cell["exclusion_reasons"]) for cell in summary["excluded_cells"]))
        self.assertTrue(any("validation absent" in str(cell["failed_trials"]) for cell in summary["excluded_cells"]))

    def test_changed_measured_output_or_time_is_rejected(self):
        for key, value in (("out_qasm", "recompiled"), ("compile_ns", 99999),
                           ("initial_layout", [1, 0, 2]), ("final_layout", [0, 2, 1])):
            protocol, timed, validated = synthetic()
            validated[0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                build_cells(protocol, timed, validated)

    def test_duplicate_identifiers_or_trial_numbers_are_not_extra_repeats(self):
        protocol, timed, validated = synthetic()
        with self.assertRaises(ValueError):
            build_cells(protocol, timed + [copy.deepcopy(timed[0])], validated)
        timed[1]["trial"] = timed[0]["trial"]
        validated[1]["trial"] = timed[0]["trial"]
        cells = build_cells(protocol, timed, validated)
        self.assertFalse(cells[0]["eligible"])
        self.assertIn("duplicate trial number", cells[0]["exclusion_reasons"])

    def test_parse_compile_is_median_of_per_trial_sums(self):
        protocol, timed, validated = synthetic()
        parse = [90, 100, 10, 90, 100, 10, 100]
        compile_times = [100, 1, 101, 100, 1, 101, 1]
        for trial in range(7):
            for rows in (timed, validated):
                rows[trial]["parse_ns"] = parse[trial]
                rows[trial]["compile_ns"] = compile_times[trial]
        medians = build_cells(protocol, timed, validated)[0]["medians"]
        self.assertEqual(medians["parse_compile_ns"], 111)
        self.assertNotEqual(medians["parse_compile_ns"], medians["parse_ns"] + medians["compile_ns"])

    def test_zero_ratios_are_counted_without_fabricating_a_geomean(self):
        observations = [
            {"fixture": "a", "numerator": 0, "denominator": 0},
            {"fixture": "b", "numerator": 0, "denominator": 2},
            {"fixture": "c", "numerator": 2, "denominator": 0},
            {"fixture": "d", "numerator": 4, "denominator": 2},
        ]
        stat = geometric_ratio_summary(observations, label="zero", draws=20)
        self.assertEqual(stat["zero_cases"], {"both_zero": 1, "numerator_zero": 1, "denominator_zero": 1})
        self.assertEqual(stat["geomean"], 2)
        self.assertIsNone(stat["ci95"])

    def test_bootstrap_clusters_replicated_topologies_together(self):
        observations = [{"fixture": "a", "numerator": 4, "denominator": 1},
                        {"fixture": "b", "numerator": 16, "denominator": 1}]
        base = geometric_ratio_summary(observations, label="cluster", draws=100)
        replicated = geometric_ratio_summary(observations * 3, label="cluster", draws=100)
        self.assertEqual(base["bootstrap_fixture_clusters"], 2)
        for expected, actual in zip(base["ci95"], replicated["ci95"]):
            self.assertAlmostEqual(expected, actual)
        self.assertAlmostEqual(base["geomean"], 8)

    def test_zero_best_requires_zero_and_all_zero_is_explicit(self):
        protocol, timed, validated = synthetic()
        first = protocol["cases"][0]["id"]
        for row in validated:
            if row["case_id"] == first:
                row["validation"]["metrics"]["cx_count"] = 0 if row["tool"] != "qrust" else 1
                row["validation"]["metrics"]["two_q_depth"] = 0
        summary, _ = analyze(protocol, timed, validated, draws=20)
        cell = summary["complete_case_quality"]["cells"][0]
        self.assertFalse(cell["quality"]["cx_count"]["within_10_percent"])
        self.assertTrue(cell["quality"]["two_q_depth"]["all_zero"])
        self.assertTrue(cell["quality"]["two_q_depth"]["within_5_percent"])
        self.assertEqual(summary["worst_tails"]["qiskit"]["worst_cx_count"][0]["ratio_kind"], "positive numerator / zero baseline")
        json.dumps(summary, allow_nan=False)

    def test_family_balance_does_not_overweight_large_family(self):
        observations = [{"family": "large", "fixture": f"a{i}", "numerator": 4, "denominator": 1} for i in range(6)]
        observations += [{"family": "small", "fixture": "b", "numerator": 16, "denominator": 1},
                         {"family": "zero", "fixture": "c", "numerator": 0, "denominator": 1}]
        stat = family_balanced_ratio_summary(observations, label="families", draws=100)
        self.assertAlmostEqual(stat["geomean"], 8)
        self.assertEqual(stat["families"], 2)
        self.assertEqual(stat["families_without_positive_ratios"], ["zero"])

    def test_reusable_artifacts_include_all_cells_and_exclusions(self):
        protocol, timed, validated = synthetic()
        summary, cells = analyze(protocol, timed, validated, draws=20)
        report = narrative(summary)
        self.assertIn("not a proven optimum", report)
        self.assertIn("fixture-cluster", report)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cell_summary.csv"
            write_csv(path, cells)
            with path.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 12)
            self.assertEqual(rows[0]["median_compile_ns"], "103")


if __name__ == "__main__":
    unittest.main()
