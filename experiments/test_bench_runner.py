"""Prevent benchmark aggregation from losing failed or unvalidated outputs."""
import copy
import unittest
from unittest.mock import patch
import bench_runner


def record(cx, total, parse, status='ok'):
    return {'status': status, 'post': {'cx': cx, 'gates': cx, 'depth': cx,
            'two_q': cx, 'swaps': 0}, 'timing_us': {'parse': parse,
            'transpile': total-parse, 'total': total, 'logic_opt': 0,
            'routing': 0, 'basis_decomp': 0, 'cleanup': 0},
            'verification': {'verified': status == 'ok', 'verdict_class': 'Numerical'},
            'initial_layout': [0], 'final_layout': [0], 'out_qasm': f'output {cx}'}


class BenchmarkAggregationTests(unittest.TestCase):
    def test_all_trials_are_preserved_and_not_mutated(self):
        runs = [record(3, 100, 90), record(4, 150, 10), record(8, 200, 180)]
        original = copy.deepcopy(runs)
        with patch.object(bench_runner, '_run_once', side_effect=runs):
            result = bench_runner.run_one('runner', {}, 'linear', 3, True, False, reps=3)
        self.assertEqual(runs, original)
        self.assertEqual(result['trials'], original)
        self.assertEqual(result['timing_us']['transpile'], 20)
        self.assertNotIn('out_qasm', result)
        self.assertNotIn('final_layout', result)

    def test_failed_repeat_is_visible_and_prevents_success_status(self):
        runs = [record(3, 100, 10), {'status': 'error', 'error': 'routing failed'}]
        with patch.object(bench_runner, '_run_once', side_effect=runs):
            result = bench_runner.run_one('runner', {}, 'linear', 3, True, False, reps=2)
        self.assertNotEqual(result['status'], 'ok')
        self.assertEqual(result['trial_status_counts'], {'ok': 1, 'error': 1})
        self.assertEqual(len(result['trials']), 2)
        self.assertFalse(result['verification']['verified'])

    def test_no_success_does_not_discard_later_errors(self):
        runs = [{'status': 'error', 'error': 'a'}, {'status': 'error', 'error': 'b'}]
        with patch.object(bench_runner, '_run_once', side_effect=runs):
            result = bench_runner.run_one('runner', {}, 'linear', 3, True, False, reps=2)
        self.assertEqual(result['trials'], runs)

    def test_zero_repetitions_rejected(self):
        with self.assertRaises(ValueError):
            bench_runner.run_one('runner', {}, 'linear', 3, True, False, reps=0)


if __name__ == '__main__':
    unittest.main()
