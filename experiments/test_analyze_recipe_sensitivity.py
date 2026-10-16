"""Checks that sensitivity ratios preserve the primary convention and scope."""
import unittest

try:
    import analyze_recipe_sensitivity as analysis
except ModuleNotFoundError as error:
    if error.name != "numpy":
        raise
    analysis = None

from test_analyze_routing_revisions import records, protocol


def fixture():
    p = protocol()
    p["tools"] = ["qrust", "cirq_one_pass", "tket_two_qubit"]
    raw, validated = records(p)
    for r in raw + validated:
        r["recipe"] = r.pop("tool")
        r["kind"] = "measured"
        r.pop("id")
        if r["recipe"] != "qrust":
            r["compile_ns"] *= 4
    for r in validated:
        if r["recipe"] != "qrust":
            r["validation"]["metrics"]["cx_count"] = 2
        r["validation"]["evidence_class"] = "full_isometry"
        r["validation_cache_key"] = r["recipe"]
    return {"cases": p["cases"], "recipes": p["tools"], "measured_repeats": p["trials"]}, raw, validated


@unittest.skipIf(analysis is None, "optional NumPy statistics dependency not installed")
class SensitivityAnalysisTests(unittest.TestCase):
    def test_speed_and_quality_ratios_use_primary_convention(self):
        p, raw, checked = fixture()
        result = analysis.analyze(p, raw, checked, draws=0)
        for baseline in p["recipes"][1:]:
            aggregate = result["pairwise"]["aggregate"][baseline]
            self.assertEqual(aggregate["compile_ns"]["geomean"], 4)
            self.assertEqual(aggregate["cx_count"]["geomean"], 2)
        self.assertEqual(result["coverage"]["measured_rows"], 9)

    def test_near_best_explicitly_excludes_absent_primary_baseline(self):
        p, raw, checked = fixture()
        quality = analysis.analyze(p, raw, checked, draws=0)["selected_recipe_quality"]
        self.assertEqual(quality["baseline_recipes"], ["cirq_one_pass", "tket_two_qubit"])
        self.assertIn("excludes Qiskit", quality["definition"])
        self.assertEqual(quality["near_best_counts"]["cx_count"]["within_10_percent"], 0)

    def test_warmups_not_included_in_medians(self):
        p, raw, checked = fixture()
        warmup = dict(raw[0], kind="warmup", trial=-1, compile_ns=10**12)
        result = analysis.analyze(p, [*raw, warmup], checked, draws=0)
        self.assertEqual(result["coverage"]["measured_rows"], 9)
        self.assertEqual(result["coverage"]["all_rows_including_warmups_probes"], 10)
        self.assertEqual(result["pairwise"]["aggregate"]["cirq_one_pass"]["compile_ns"]["geomean"], 4)

    def test_failed_candidate_cell_does_not_contribute_partial_successes(self):
        p, raw, checked = fixture()
        checked[-1]["validation"]["passed"] = False
        result = analysis.analyze(p, raw, checked, draws=0)
        self.assertEqual(result["coverage"]["excluded_cells"], 1)
        self.assertEqual(result["pairwise"]["aggregate"]["tket_two_qubit"]["paired_cells"], 0)
        self.assertEqual(result["selected_recipe_quality"]["complete_cases"], 0)


if __name__ == "__main__":
    unittest.main()
