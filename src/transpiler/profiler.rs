use crate::ir::{Circuit, GateType};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashMap;

/// A heuristic report exposing algebraic and structural depth properties ahead of time.
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileReport {
    /// Histogram of standard physical gate operations natively present.
    pub gate_counts: HashMap<GateType, usize>,
    /// Exact count of purely self-inverse operators present (e.g. H, CX).
    pub self_inverse_count: usize,
    /// Exact count of continuous parametric generators (e.g. RZ, RX).
    pub continuous_rotation_count: usize,
    /// Absolute total boundary size.
    pub total_gates: usize,
}

impl ProfileReport {
    pub fn can_bypass_inverse_cancellation(&self) -> bool {
        if self.self_inverse_count == 0 && self.continuous_rotation_count < 2 {
            return true;
        }
        false
    }
}

/// A lightweight $O(N)$ pre-flight structural scanner acting as an Analysis Pass.
pub struct CircuitProfilerPass;

impl Pass for CircuitProfilerPass {
    fn name(&self) -> &str {
        "CircuitProfilerPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        let mut counts = HashMap::new();
        let mut self_inverse = 0;
        let mut continuous_rotation = 0;

        for op in &circuit.operations {
            if let crate::ir::Operation::Gate { name, .. } = op {
                *counts.entry(name.clone()).or_insert(0) += 1;

                match name {
                    GateType::H
                    | GateType::X
                    | GateType::Y
                    | GateType::Z
                    | GateType::CX
                    | GateType::CY
                    | GateType::CZ
                    | GateType::SWAP
                    | GateType::CCX => {
                        self_inverse += 1;
                    }
                    GateType::RX
                    | GateType::RY
                    | GateType::RZ
                    | GateType::U
                    | GateType::CRX
                    | GateType::CRY
                    | GateType::CRZ => {
                        continuous_rotation += 1;
                    }
                    _ => {}
                }
            }
        }

        let report = ProfileReport {
            gate_counts: counts,
            self_inverse_count: self_inverse,
            continuous_rotation_count: continuous_rotation,
            total_gates: circuit.operations.len(),
        };

        property_set.insert("profile_report", report);
        circuit.clone()
    }
}
