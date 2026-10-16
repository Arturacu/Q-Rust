"""Counterexamples for revision-study eligibility and ratio interpretation."""
import copy
import hashlib
import unittest

try:
    import analyze_routing_revisions as analysis
except ModuleNotFoundError as error:
    if error.name != "numpy":
        raise
    analysis = None


def protocol(cases=None, trials=3):
    return {"tools": ["before", "after"], "trials": trials, "cases": cases or [
        {"id": "a__linear", "fixture": "a", "family": "test", "topology": "linear", "n_logical": 2}]}


def records(spec, times=(12, 9, 15), counts=(4, 4, 4)):
    raw, validated = [], []
    output = 'OPENQASM 2.0; include "qelib1.inc"; qreg q[2];'
    for case in spec["cases"]:
        for tool in spec["tools"]:
            for i in range(spec["trials"]):
                row = {"id": f"{case['id']}__{tool}__{i}", "case_id": case["id"],
                       "tool": tool, "trial": i, "status": "compiled",
                       "compile_ns": times[i], "parse_ns": (1, 100, 1)[i],
                       "out_qasm": output, "output_sha256": hashlib.sha256(output.encode()).hexdigest()}
                raw.append(row)
                checked = dict(row)
                checked["validation"] = {"passed": True, "metrics": {
                    "cx_count": counts[i], "gate_count": counts[i]+2, "depth": counts[i]+1, "two_q_depth": counts[i]}}
                validated.append(checked)
    return raw, validated


@unittest.skipIf(analysis is None, "optional NumPy statistics dependency not installed")
class RevisionAnalysisTests(unittest.TestCase):
    def test_median_of_per_trial_sum_and_correct_ratio_orientation(self):
        p = protocol()
        raw, checked = records(p)
        for r in raw + checked:
            if r["tool"] == "after":
                r["compile_ns"] *= 2
        for r in checked:
            if r["tool"] == "after":
                r["validation"]["metrics"]["cx_count"] = 2
        result = analysis.analyze(p, raw, checked, "before", draws=0)
        aggregate = result["comparisons"]["after"]["aggregate"]
        self.assertAlmostEqual(aggregate["compile_ns"]["geomean"], 0.5)
        self.assertAlmostEqual(aggregate["cx_count"]["geomean"], 0.5)
        self.assertEqual(aggregate["compile_ns"]["candidate_win_tie_loss"], {"win": 0, "tie": 0, "loss": 1})
        self.assertEqual(aggregate["cx_count"]["candidate_win_tie_loss"], {"win": 1, "tie": 0, "loss": 0})
        before = next(c for c in result["cells"] if c["tool"] == "before")
        self.assertEqual(before["medians"]["parse_compile_ns"], 16)
        self.assertNotEqual(before["medians"]["parse_compile_ns"], before["medians"]["compile_ns"]+1)

    def test_failed_trial_excludes_entire_pair_without_cherry_picking(self):
        p = protocol()
        raw, checked = records(p)
        checked[-1]["validation"]["passed"] = False
        result = analysis.analyze(p, raw, checked, "before", draws=0)
        self.assertEqual(result["coverage"]["excluded_cell_count"], 1)
        self.assertEqual(result["comparisons"]["after"]["aggregate"]["paired_cases"], 0)
        self.assertEqual(result["cells"][-1]["medians"], None)

    def test_missing_and_unexpected_trial_ids_are_not_complete(self):
        for replacement in (None, 99):
            with self.subTest(replacement=replacement):
                p = protocol()
                raw, checked = records(p)
                if replacement is None:
                    raw.pop(); checked.pop()
                else:
                    raw[-1]["trial"] = checked[-1]["trial"] = replacement
                cells = analysis.build_cells(p, raw, checked)
                self.assertFalse(cells[-1]["eligible"])

    def test_changed_validated_timing_or_output_is_rejected(self):
        for field, value in (("compile_ns", 999), ("out_qasm", "changed")):
            p = protocol(); raw, checked = records(p)
            checked[0][field] = value
            with self.assertRaisesRegex(ValueError, "changed measured record"):
                analysis.build_cells(p, raw, checked)

    def test_stale_output_digest_rejected_even_if_validator_says_passed(self):
        p = protocol(); raw, checked = records(p)
        raw[0]["out_qasm"] = checked[0]["out_qasm"] = "changed"
        with self.assertRaisesRegex(ValueError, "output hash mismatch"):
            analysis.build_cells(p, raw, checked)

    def test_duplicate_record_and_unmeasured_validation_rejected(self):
        p = protocol(); raw, checked = records(p)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            analysis.build_cells(p, raw + [copy.deepcopy(raw[0])], checked)
        checked[-1]["id"] = "not-measured"
        with self.assertRaisesRegex(ValueError, "unmeasured"):
            analysis.build_cells(p, raw, checked)

    def test_zero_quality_not_logged_or_hidden(self):
        p = protocol(); raw, checked = records(p, counts=(0, 0, 0))
        for r in checked:
            if r["tool"] == "after":
                r["validation"]["metrics"]["cx_count"] = 2
        result = analysis.analyze(p, raw, checked, "before", draws=0)
        quality = result["comparisons"]["after"]["aggregate"]["cx_count"]
        self.assertIsNone(quality["geomean"])
        self.assertEqual(quality["zero_cases"]["denominator_zero"], 1)
        self.assertEqual(quality["candidate_win_tie_loss"]["loss"], 1)

    def test_one_fixture_many_topologies_has_no_spurious_interval(self):
        cases = [{"id": f"a__{t}", "fixture": "a", "family": "test", "topology": t, "n_logical": 2}
                 for t in ("linear", "all_to_all")]
        p = protocol(cases); raw, checked = records(p)
        result = analysis.analyze(p, raw, checked, "before", draws=50)
        stats = result["comparisons"]["after"]["aggregate"]["compile_ns"]
        self.assertEqual(stats["paired_cells"], 2)
        self.assertEqual(stats["bootstrap_fixture_clusters"], 1)
        self.assertIsNone(stats["ci95"])

    def test_two_fixture_bootstrap_is_reproducible(self):
        cases = [{"id": f"{f}__linear", "fixture": f, "family": "test", "topology": "linear", "n_logical": 2}
                 for f in ("a", "b")]
        p = protocol(cases); raw, checked = records(p)
        for r in raw + checked:
            if r["case_id"].startswith("b") and r["tool"] == "after":
                r["compile_ns"] *= 4
        a = analysis.analyze(p, raw, checked, "before", draws=100)
        b = analysis.analyze(p, raw, checked, "before", draws=100)
        stats = a["comparisons"]["after"]["aggregate"]["compile_ns"]
        self.assertEqual(a, b)
        self.assertAlmostEqual(stats["geomean"], 0.5)
        self.assertEqual(stats["ci95"], [0.25, 1.0])


if __name__ == "__main__":
    unittest.main()
