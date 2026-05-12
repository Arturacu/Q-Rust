//! Analysis-only profiler pass.
//!
//! Emits a [`ProfileReport`] into the [`PropertySet`] under the key
//! `"profile_report"`. Downstream passes can use this to decide whether
//! expensive analyses are worth running.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashMap;

/// Summary of a circuit's algebraic and structural properties.
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileReport {
    /// Histogram of gate-type occurrences.
    pub gate_counts: HashMap<GateType, usize>,
    /// Count of self-inverse gates present (e.g. H, CX, CCX).
    pub self_inverse_count: usize,
    /// Count of parametric rotation gates present.
    pub continuous_rotation_count: usize,
    /// Total number of operations in the circuit.
    pub total_gates: usize,
}

impl ProfileReport {
    /// Returns `true` if the circuit has no algebraic opportunity for
    /// inverse cancellation.
    pub fn can_bypass_inverse_cancellation(&self) -> bool {
        self.self_inverse_count == 0 && self.continuous_rotation_count < 2
    }
}

/// Cheap O(N) structural scan that populates [`ProfileReport`].
#[derive(Debug, Clone, Copy)]
pub struct CircuitProfilerPass;

impl Pass for CircuitProfilerPass {
    fn name(&self) -> &str {
        "CircuitProfilerPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        let mut counts: HashMap<GateType, usize> = HashMap::new();
        let mut self_inverse = 0usize;
        let mut continuous_rotation = 0usize;

        for op in &circuit.operations {
            if let Operation::Gate { name, .. } = op {
                *counts.entry(name.clone()).or_default() += 1;
                match name {
                    GateType::H
                    | GateType::X
                    | GateType::Y
                    | GateType::Z
                    | GateType::CX
                    | GateType::CY
                    | GateType::CZ
                    | GateType::SWAP
                    | GateType::CCX => self_inverse += 1,
                    GateType::RX
                    | GateType::RY
                    | GateType::RZ
                    | GateType::U
                    | GateType::CRX
                    | GateType::CRY
                    | GateType::CRZ => continuous_rotation += 1,
                    _ => {}
                }
            }
        }

        property_set.insert(
            "profile_report",
            ProfileReport {
                gate_counts: counts,
                self_inverse_count: self_inverse,
                continuous_rotation_count: continuous_rotation,
                total_gates: circuit.operations.len(),
            },
        );

        circuit.clone()
    }
}
