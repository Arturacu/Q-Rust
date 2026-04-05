pub mod dag;
pub mod decomposition;
pub mod optimization;
pub mod pass;
pub mod synthesis;

use crate::ir::Circuit;
use pass::PassManager;

/// Configuration for the transpiler.
pub struct TranspilerConfig {
    /// If true, decomposes all gates to the basis set (U, CX).
    pub decompose_basis: bool,
    /// Optimization level (0-3). Currently unused but reserved.
    pub optimization_level: u8,
}

impl Default for TranspilerConfig {
    fn default() -> Self {
        Self {
            decompose_basis: true,
            optimization_level: 1,
        }
    }
}

/// Transpiles a circuit according to the given configuration.
///
/// This is the main entry point for the transpiler. It constructs a
/// PassManager based on the config and runs it on the circuit.
pub fn transpile(circuit: &Circuit, config: Option<TranspilerConfig>) -> Circuit {
    let config = config.unwrap_or_default();
    let mut pm = PassManager::new();

    if config.decompose_basis {
        pm.add_pass(Box::new(decomposition::BasisDecompositionPass));
    }

    pm.run(circuit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{GateType, Operation};

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
        };

        // Should keep H gate
        let new_circuit = transpile(&circuit, Some(config));

        assert_eq!(new_circuit.operations.len(), 1);
        match &new_circuit.operations[0] {
            Operation::Gate { name, .. } => assert_eq!(*name, GateType::H),
            _ => panic!("Expected Gate operation"),
        }
    }
}
