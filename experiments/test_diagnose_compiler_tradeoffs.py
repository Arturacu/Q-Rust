"""Semantic and baseline-parity checks for separately labeled recipe ablations."""
from __future__ import annotations

import importlib.util
import unittest

from diagnose_compiler_tradeoffs import RECIPES, cirq_trial, tket_trial, compile_recipe, validate_record_grid

HAS_DEPS = all(importlib.util.find_spec(name) is not None
               for name in ("cirq", "pytket", "numpy", "qiskit"))


def example(topology="linear", spare=False):
    source = ('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[3];\n'
              'h q[0];\nry(0.37) q[1];\ncx q[0],q[2];\n'
              'rz(0.21) q[2];\ncx q[1],q[0];\nswap q[0],q[2];\n')
    width = 4 if spare else 3
    return {"qasm": source, "source_qasm": source, "n_logical": 3,
            "physical_width": width, "edges": [(i, i + 1) for i in range(width - 1)],
            "topology": topology}


class NativeWorkerDispatchTests(unittest.TestCase):
    def test_native_timer_and_pipeline_are_preserved(self):
        from unittest.mock import Mock
        worker = Mock()
        worker.call.return_value = {"status": "compiled", "compile_ns": 1234,
                                    "pipeline": "canonical", "seed_used": False}
        result = compile_recipe(example(), "qrust", worker, seed=99)
        self.assertEqual(result["compile_ns"], 1234)
        self.assertEqual(result["pipeline"], "canonical")
        self.assertFalse(result["seed_used"])
        request, timeout = worker.call.call_args.args
        self.assertEqual(request["seed"], 99)
        self.assertEqual(timeout, 60)
        self.assertNotIn("source_qasm", request)

    def test_dispatch_rejects_missing_worker_and_unknown_recipe(self):
        with self.assertRaisesRegex(ValueError, "persistent native worker"):
            compile_recipe(example(), "qrust")
        with self.assertRaisesRegex(ValueError, "Unknown recipe"):
            compile_recipe(example(), "other")


class GridIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.protocol = {"case_ids": ["a", "b"], "recipes": ["qrust", "tket_two_qubit"],
                         "measured_repeats": 2, "stage_probe_repeats": 1, "warmups": 1}
        self.rows = []
        for case in self.protocol["case_ids"]:
            for recipe in self.protocol["recipes"]:
                schedule = [(-1, "warmup"), (0, "measured"), (1, "measured")]
                if recipe.startswith("tket_"):
                    schedule += [(-1, "stage_warmup"), (0, "stage_probe")]
                self.rows += [{"case_id": case, "recipe": recipe, "trial": trial, "kind": kind}
                              for trial, kind in schedule]

    def test_complete_grid_accepts_any_execution_order(self):
        result = validate_record_grid(self.protocol, list(reversed(self.rows)))
        self.assertEqual(result["expected_records"], 16)
        self.assertEqual(result["actual_records"], 16)

    def test_missing_warmup_and_missing_probe_are_rejected(self):
        for index in (0, 7):
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, "1 missing"):
                validate_record_grid(self.protocol, self.rows[:index] + self.rows[index + 1:])

    def test_duplicate_and_extra_trial_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "1 duplicate"):
            validate_record_grid(self.protocol, self.rows + [self.rows[0]])
        extra = dict(self.rows[0], trial=99, kind="measured")
        with self.assertRaisesRegex(ValueError, "1 unexpected"):
            validate_record_grid(self.protocol, self.rows + [extra])

    def test_boolean_trial_index_is_not_an_integer_trial(self):
        rows = [dict(r) for r in self.rows]
        rows[1]["trial"] = False
        with self.assertRaisesRegex(ValueError, "integer"):
            validate_record_grid(self.protocol, rows)


@unittest.skipUnless(HAS_DEPS, "Optional comparison dependencies not installed")
class TradeoffDiagnosticTests(unittest.TestCase):
    def test_existing_campaign_directory_is_never_overwritten(self):
        from argparse import Namespace
        from pathlib import Path
        from tempfile import TemporaryDirectory
        from unittest.mock import patch
        from diagnose_compiler_tradeoffs import timing
        with TemporaryDirectory() as directory:
            output = Path(directory)
            saved = output / "tradeoffs_protocol.json"
            saved.write_text("preserve prior artifact")
            args = Namespace(protocol=Path("unused"), all_cases=False, recipes=None, output=output)
            with patch("diagnose_compiler_tradeoffs.load_cases", return_value=[]):
                with self.assertRaises(FileExistsError):
                    timing(args)
            self.assertEqual(saved.read_text(), "preserve prior artifact")

    def validate(self, case, row):
        from comparison_validation import validate_output
        edges = None if case["topology"] == "all_to_all" else [e for a, b in case["edges"] for e in [(a, b), (b, a)]]
        result = validate_output(case["source_qasm"], row["out_qasm"], case["n_logical"],
                                 case["physical_width"], edges,
                                 row["initial_layout"], row["final_layout"])
        self.assertTrue(result["passed"], result)
        self.assertEqual(result["evidence_class"], "full_isometry")
        self.assertGreater(row["compile_ns"], 0)
        return result

    def test_all_six_recipes_preserve_routed_arbitrary_input_with_spare_wire(self):
        case = example(spare=True)
        for recipe in RECIPES:
            with self.subTest(recipe=recipe):
                row = cirq_trial(case, recipe) if recipe.startswith("cirq_") else tket_trial(case, recipe)
                self.validate(case, row)

    def test_cirq_original_matches_frozen_adapter(self):
        from comparison_cirq import compile_one
        case = example()
        reference = compile_one(case["qasm"], case["n_logical"], case["physical_width"],
                                case["edges"], case["topology"], seed=7)
        actual = cirq_trial(case, "cirq_original")
        for field in ("out_qasm", "initial_layout", "final_layout"):
            self.assertEqual(actual[field], reference[field])

    def test_both_cirq_synthesis_recipes_close_imported_full_cz_over_cx(self):
        for recipe in ("cirq_original", "cirq_one_pass"):
            for topology in ("all_to_all", "linear"):
                with self.subTest(recipe=recipe, topology=topology):
                    case = example(topology)
                    case["qasm"] = ('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[3];\n'
                                    'cu1(pi) q[2],q[0];\nh q[1];\ncx q[0],q[1];\n')
                    case["source_qasm"] = case["qasm"]
                    row = cirq_trial(case, recipe)
                    self.validate(case, row)
                    self.assertNotIn("cu1", row["out_qasm"])
                    self.assertIn("cx", row["out_qasm"])

    def test_tket_original_matches_frozen_adapter(self):
        from comparison_qiskit_tket import compile_one
        case = example()
        reference = compile_one("tket", case["qasm"], case["n_logical"], case["physical_width"],
                                case["edges"], case["topology"], seed=7)
        actual = tket_trial(case, "tket_original")
        for field in ("out_qasm", "initial_layout", "final_layout"):
            self.assertEqual(actual[field], reference[field])

    def test_tket_stage_probes_match_sequencepass_outputs(self):
        case = example()
        for recipe in (r for r in RECIPES if r.startswith("tket_")):
            with self.subTest(recipe=recipe):
                sequence = tket_trial(case, recipe)
                stages = tket_trial(case, recipe, stage_probe=True)
                self.validate(case, stages)
                for field in ("out_qasm", "initial_layout", "final_layout"):
                    self.assertEqual(sequence[field], stages[field])

    def test_native_cx_rejects_unsupported_high_arity_input(self):
        case = example("all_to_all")
        case["qasm"] += 'ccx q[0],q[1],q[2];\n'
        with self.assertRaisesRegex(ValueError, "outside 1q/CX/SWAP"):
            cirq_trial(case, "cirq_native_cx")

    def test_full_corpus_sensitivity_recipes_lower_ccx_inside_compilation(self):
        case = example()
        case["qasm"] += 'ccx q[0],q[1],q[2];\n'
        case["source_qasm"] = case["qasm"]
        self.validate(case, cirq_trial(case, "cirq_one_pass"))
        self.validate(case, tket_trial(case, "tket_two_qubit"))


if __name__ == "__main__":
    unittest.main()
