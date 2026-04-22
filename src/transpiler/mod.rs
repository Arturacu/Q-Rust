pub mod dag;
pub mod decomposition;
pub mod layout;
pub mod optimization;
pub mod pass;
pub mod pauli_tracker;
pub mod profiler;
pub mod property_set;
pub mod routing;
pub mod synthesis;

use crate::ir::Circuit;
use pass::PassManager;

/// Configuration for the transpiler.
#[derive(Clone, Debug)]
pub struct TranspilerConfig {
    /// If true, decomposes all gates to the basis set (U, CX).
    pub decompose_basis: bool,
    /// Optimization level (0-3). Currently unused but reserved.
    pub optimization_level: u8,
    /// Explicit target hardware configuration topology.
    /// If None, safely defaults strictly to an All-to-All implicit software simulation bound.
    pub backend: Option<crate::backend::Backend>,
}

impl Default for TranspilerConfig {
    fn default() -> Self {
        Self {
            decompose_basis: true,
            optimization_level: 1,
            backend: None,
        }
    }
}

/// Transpiles a circuit according to the given configuration.
///
/// This is the main entry point for the transpiler. It constructs a
/// PassManager mapping exactly to standard execution stages:
/// 1. Type Promotion (Crystallization)
/// 2. Algebraic Optimization (Pauli Tracking, Commutation, Rotations)
/// 3. Basis Unrolling (Decomposition)
/// 4. Hardware routing (Pending Phase 5 features)
pub fn transpile(circuit: &Circuit, config: Option<TranspilerConfig>) -> Circuit {
    let config = config.unwrap_or_default();
    let mut pm = PassManager::new();

    // Stage 1: Mathematical Type Promotion
    if config.optimization_level >= 1 {
        pm.add_pass(Box::new(optimization::GateCrystallizationPass {
            epsilon: 1e-9,
        }));
    }

    // Stage 2: Pre-Routing Algebraic Simplification
    if config.optimization_level >= 2 {
        pm.add_pass(Box::new(optimization::RotationMergePass {}));
        pm.add_pass(Box::new(optimization::CrossConjugationPass {}));
        pm.add_pass(Box::new(optimization::InverseCancellationPass {}));
        pm.add_pass(Box::new(optimization::CommutationCancellationPass {}));
        // Note: PauliTrackerPass requires explicit invocation / integration mechanics for bounding contexts
        // pm.add_pass(Box::new(pauli_tracker::PauliTrackerPass {}));
    }

    // Stage 2.5: Hardware Routing (only when targeting a specific backend)
    if let Some(ref backend) = config.backend {
        // Stage 2.5a: Layout discovery via iterative SABRE (opt level ≥ 2)
        if config.optimization_level >= 2 {
            let (num_trials, num_iterations) = match config.optimization_level {
                2 => (10, 3), // Moderate search
                _ => (50, 5), // Full search with many stochastic restarts
            };
            pm.add_pass(Box::new(layout::SabreLayoutPass {
                backend: backend.clone(),
                num_trials,
                num_iterations,
            }));
        }

        // Stage 2.5b: BeamSABRE routing (uses layout from above if available)
        let (beam_width, branch_factor, bidir_iters) = match config.optimization_level {
            0 | 1 => (1, 1, 1), // Greedy (= LightSABRE)
            2 => (4, 3, 2),     // Moderate beam search
            _ => (8, 5, 4),     // Full BeamSABRE
        };
        pm.add_pass(Box::new(routing::BeamSabrePass {
            backend: backend.clone(),
            beam_width,
            branch_factor,
            bidir_iterations: bidir_iters,
        }));
    }

    // Stage 3: Hardware Unrolling (Basis Decomposition)
    if config.decompose_basis {
        pm.add_pass(Box::new(decomposition::BasisDecompositionPass {}));
    }

    // Stage 4: Post-Routing Graph Fusions
    if config.optimization_level >= 3 {
        pm.add_pass(Box::new(optimization::GateFusionPass {}));
        pm.add_pass(Box::new(optimization::SwapSimplificationPass {}));
        // Inverse again just in case fusion generated new identities natively
        pm.add_pass(Box::new(optimization::InverseCancellationPass {}));
    }

    pm.run(circuit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::Backend;
    use crate::ir::{GateType, Operation};
    use crate::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};

    #[test]
    fn test_transpile_default() {
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });

        // Default config should decompose H -> U
        let new_circuit = transpile(&circuit, None);

        assert_eq!(new_circuit.operations.len(), 1);
        match &new_circuit.operations[0] {
            Operation::Gate { name, .. } => {
                if let GateType::U = name {
                    // OK
                } else {
                    panic!("Expected U gate, got {:?}", name);
                }
            }
            _ => panic!("Expected Gate operation"),
        }
    }

    #[test]
    fn test_transpile_no_decompose() {
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });

        let config = TranspilerConfig {
            decompose_basis: false,
            optimization_level: 0,
            backend: None,
        };

        // Should keep H gate
        let new_circuit = transpile(&circuit, Some(config));

        assert_eq!(new_circuit.operations.len(), 1);
        match &new_circuit.operations[0] {
            Operation::Gate { name, .. } => assert_eq!(*name, GateType::H),
            _ => panic!("Expected Gate operation"),
        }
    }

    #[test]
    fn test_pipeline_e2e_fidelity() {
        // GHZ-3 circuit on linear(3) backend
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 2],
            params: vec![],
        });

        let backend = Backend::linear(3);
        let mut pm = PassManager::new();

        // Add layout and routing
        pm.add_pass(Box::new(layout::SabreLayoutPass {
            backend: backend.clone(),
            num_trials: 10,
            num_iterations: 2,
        }));
        pm.add_pass(Box::new(routing::BeamSabrePass {
            backend: backend.clone(),
            beam_width: 2,
            branch_factor: 2,
            bidir_iterations: 1,
        }));

        let routed = pm.run(&circuit);

        let u_orig = circuit_to_unitary(&circuit);
        let u_phys = circuit_to_unitary(&routed);

        let initial = pm.property_set.get::<Vec<usize>>("initial_layout").unwrap();
        let final_l = pm.property_set.get::<Vec<usize>>("final_layout").unwrap();

        let u_log = extract_logical_unitary(&u_phys, 3, initial, final_l);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);

        // Final connectivity check
        for op in &routed.operations {
            if let Operation::Gate { qubits, .. } = op {
                if qubits.len() == 2 {
                    assert!(backend.is_adjacent(qubits[0], qubits[1]));
                }
            }
        }
    }
}
