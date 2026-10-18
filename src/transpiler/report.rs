//! Per-stage transpilation report.
//!
//! Tracks gate count, two-qubit gate count, depth, and (when routing runs)
//! the initial/final logical-to-physical layouts.
//!
//! Stage names are stable strings (`"1. parsed"`, `"2. optimized"`,
//! `"3. routed+decomposed"` or `"3. decomposed"`) so the CLI smoke
//! tests and downstream consumers can pin against them.

use crate::ir::{Circuit, GateType, Operation};

/// A snapshot of circuit metrics at a single pipeline stage.
#[derive(Debug, Clone, PartialEq)]
pub struct StageSnapshot {
    pub stage: String,
    pub num_ops: usize,
    pub num_2q_gates: usize,
    pub num_swaps: usize,
    pub depth: usize,
}

impl StageSnapshot {
    pub fn capture(stage: impl Into<String>, circuit: &Circuit) -> Self {
        let mut num_2q_gates = 0usize;
        let mut num_swaps = 0usize;
        for op in &circuit.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                if qubits.len() >= 2 {
                    num_2q_gates += 1;
                    if matches!(name, GateType::SWAP) {
                        num_swaps += 1;
                    }
                }
            }
        }
        Self {
            stage: stage.into(),
            num_ops: circuit.operations.len(),
            num_2q_gates,
            num_swaps,
            depth: circuit.depth(),
        }
    }
}

/// Optional per-pass evidence. Timing excludes the following metric capture.
#[derive(Debug, Clone, PartialEq)]
pub struct PassDiagnostic {
    pub elapsed_ns: u64,
    pub layout_route_evaluations: Option<usize>,
    pub output: StageSnapshot,
}

/// Final-output protection for Clifford rewrites before or after routing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CliffordRetentionStats {
    /// Costs are [two-qubit gates, total depth, two-qubit depth] after all
    /// lowering, direction repair and terminal readout restoration.
    /// This is the winner of the original routing portfolio, which can itself
    /// retain an unrewritten local routing candidate.
    pub rewritten_output_cost: [usize; 3],
    pub retained_output_cost: Option<[usize; 3]>,
    pub selected_output_cost: [usize; 3],
    pub selected_retained_input: bool,
    /// An unavailable optional continuation cannot certify nonregression.
    pub retained_compilation_failed: bool,
    /// Exact rewrites from the original portfolio winner, or the first changed
    /// candidate if the winner had none. Includes work later discarded.
    pub attempted_simplification: super::clifford_simplification::CliffordSimplificationStats,
}

/// Full transpilation report.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TranspilationReport {
    /// Exact encoded-SWAP recognition counters, when enabled.
    pub encoded_swap_stats: Option<super::encoded_swaps::EncodedSwapStats>,
    /// Exact cross-block Clifford simplification counters, when enabled.
    pub clifford_simplification_stats:
        Option<super::clifford_simplification::CliffordSimplificationStats>,
    /// Final-output comparison with an independently lowered pre-Clifford input.
    pub clifford_retention_stats: Option<CliffordRetentionStats>,
    /// Shared numerical replacement ledger; not a proof of equivalence.
    pub numerical_budget: super::numerical_policy::NumericalBudget,
    pub stages: Vec<StageSnapshot>,
    pub pass_diagnostics: Vec<PassDiagnostic>,
    /// Physical location of each logical input qubit, when routing ran.
    pub initial_layout: Option<Vec<usize>>,
    /// Physical location of each logical output qubit, when routing ran.
    pub final_layout: Option<Vec<usize>>,
    /// Cheap block-synthesis counters, when that cleanup pass ran. These do
    /// not enable per-pass clocks or diagnostic circuit captures.
    pub two_qubit_block_stats: Option<super::block_synthesis::TwoQubitBlockStats>,
    /// Exact CX-RZ-CX patterns recovered before routing, when enabled.
    pub phase_gadgets_recovered: Option<usize>,
    /// Numerical local-fusion counters, without detailed profiling clocks.
    pub post_block_fusion_stats: Option<super::local_cleanup::PostBlockFusionStats>,
}

impl TranspilationReport {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, snapshot: StageSnapshot) {
        self.stages.push(snapshot);
    }

    /// Renders the report as a list of human-readable lines (one per stage,
    /// plus a header row).
    pub fn format_lines(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(self.stages.len() + 1);
        out.push(format!(
            "{:<28} {:>8} {:>8} {:>8} {:>8}",
            "stage", "ops", "2q", "swaps", "depth"
        ));
        for s in &self.stages {
            out.push(format!(
                "{:<28} {:>8} {:>8} {:>8} {:>8}",
                s.stage, s.num_ops, s.num_2q_gates, s.num_swaps, s.depth
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snapshot_captures_2q_gates_and_swaps() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        let snap = StageSnapshot::capture("test", &c);
        assert_eq!(snap.num_ops, 3);
        assert_eq!(snap.num_2q_gates, 2);
        assert_eq!(snap.num_swaps, 1);
    }

    #[test]
    fn test_format_lines_has_header_plus_one_per_stage() {
        let mut r = TranspilationReport::new();
        let c = Circuit::new(1, 0);
        r.push(StageSnapshot::capture("parsed", &c));
        r.push(StageSnapshot::capture("optimized", &c));
        let lines = r.format_lines();
        assert_eq!(lines.len(), 3);
    }
}
