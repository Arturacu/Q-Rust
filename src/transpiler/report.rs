//! [E2E-NEW-FEATURE] Per-stage transpilation report.
//!
//! Tracks gate count, two-qubit gate count, and depth as the circuit
//! flows through the pipeline, mirroring thesis Figure 4.1's dashed
//! per-stage annotations ("Fewer gates", "SWAPs inserted", "{U,CX} only").
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

/// Full transpilation report.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TranspilationReport {
    pub stages: Vec<StageSnapshot>,
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