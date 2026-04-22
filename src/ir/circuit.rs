use super::operations::Operation;
use super::registry::GateRegistry;

/// Intermediate Representation of a Quantum Circuit.
///
/// A `Circuit` consists of a sequence of operations and metadata about the
/// number of qubits and classical bits required.
#[derive(Debug, Clone, Default)]
pub struct Circuit {
    /// Total number of qubits in the circuit.
    pub num_qubits: usize,
    /// Total number of classical bits in the circuit.
    pub num_cbits: usize,
    /// Sequence of operations (gates, measurements, etc.).
    pub operations: Vec<Operation>,
    /// Registry of user-defined custom gates.
    pub custom_gates: GateRegistry,
}

impl Circuit {
    /// Creates a new empty circuit.
    ///
    /// # Arguments
    ///
    /// * `num_qubits` - The number of qubits to allocate.
    /// * `num_cbits` - The number of classical bits to allocate.
    pub fn new(num_qubits: usize, num_cbits: usize) -> Self {
        Self {
            num_qubits,
            num_cbits,
            operations: Vec::new(),
            custom_gates: GateRegistry::new(),
        }
    }

    /// Adds a gate definition to the circuit's registry.
    pub fn register_custom_gate(
        &mut self,
        name: String,
        params: Vec<String>,
        qubits: Vec<String>,
        body: Vec<crate::ir::ast::ParsedStatement>,
    ) {
        self.custom_gates.register(name, params, qubits, body);
    }

    /// Adds an operation to the circuit.
    pub fn add_op(&mut self, op: Operation) {
        self.operations.push(op);
    }

    /// Validates the circuit and returns a list of warnings.
    ///
    /// Checks:
    /// - Presence of at least one measurement.
    pub fn validate(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        let has_measurement = self
            .operations
            .iter()
            .any(|op| matches!(op, Operation::Measure { .. }));

        if !has_measurement {
            warnings.push("Warning: No measurements found. The circuit will not produce classical output on hardware.".to_string());
        }

        warnings
    }

    /// Exports the Circuit as an OpenQASM 2.0 string.
    ///
    /// If `property_set` is provided, includes `initial_layout` and `final_layout` mapping
    /// arrays as metadata comments so downstream tools can track hardware routing results.
    pub fn to_qasm(
        &self,
        property_set: Option<&crate::transpiler::property_set::PropertySet>,
    ) -> String {
        let mut qasm = String::new();
        qasm.push_str("OPENQASM 2.0;\n");
        qasm.push_str("include \"qelib1.inc\";\n\n");

        // Write custom metadata
        if let Some(props) = property_set {
            if let Some(initial) = props.get::<Vec<usize>>("initial_layout") {
                qasm.push_str(&format!("// qrust_initial_layout: {:?}\n", initial));
            }
            if let Some(final_l) = props.get::<Vec<usize>>("final_layout") {
                qasm.push_str(&format!("// qrust_final_layout: {:?}\n", final_l));
            }
        }

        qasm.push_str(&format!("qreg q[{}];\n", self.num_qubits));
        if self.num_cbits > 0 {
            qasm.push_str(&format!("creg c[{}];\n", self.num_cbits));
        }

        qasm.push('\n');

        for op in &self.operations {
            qasm.push_str(&op.to_qasm());
            qasm.push('\n');
        }

        qasm
    }

    /// Calculates the quantum depth of the circuit.
    ///
    /// Depth is defined as the length of the longest path in the circuit's
    /// dependency graph (critical path).
    pub fn depth(&self) -> usize {
        let mut qubit_depths = vec![0usize; self.num_qubits];

        for op in &self.operations {
            match op {
                Operation::Gate { qubits, .. } => {
                    let mut max_d = 0;
                    for &q in qubits {
                        if q < qubit_depths.len() {
                            max_d = max_d.max(qubit_depths[q]);
                        }
                    }
                    let new_d = max_d + 1;
                    for &q in qubits {
                        if q < qubit_depths.len() {
                            qubit_depths[q] = new_d;
                        }
                    }
                }
                Operation::Measure { qubit, .. } => {
                    if *qubit < qubit_depths.len() {
                        qubit_depths[*qubit] += 1;
                    }
                }
                Operation::Reset { qubit } => {
                    if *qubit < qubit_depths.len() {
                        qubit_depths[*qubit] += 1;
                    }
                }
                Operation::Barrier { qubits } => {
                    let mut max_d = 0;
                    for &q in qubits {
                        if q < qubit_depths.len() {
                            max_d = max_d.max(qubit_depths[q]);
                        }
                    }
                    for &q in qubits {
                        if q < qubit_depths.len() {
                            qubit_depths[q] = max_d;
                        }
                    }
                }
            }
        }

        *qubit_depths.iter().max().unwrap_or(&0)
    }

    /// Returns the total number of gate operations in the circuit.
    pub fn gate_count(&self) -> usize {
        self.operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { .. }))
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::gates::GateType;

    #[test]
    fn test_circuit_creation() {
        let circuit = Circuit::new(2, 2);
        assert_eq!(circuit.num_qubits, 2);
        assert_eq!(circuit.num_cbits, 2);
        assert!(circuit.operations.is_empty());
    }

    #[test]
    fn test_add_op() {
        let mut circuit = Circuit::new(1, 0);
        let op = Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        };
        circuit.add_op(op.clone());
        assert_eq!(circuit.operations.len(), 1);
        assert_eq!(circuit.operations[0], op);
    }

    #[test]
    fn test_validation_no_measurements() {
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let warnings = circuit.validate();
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("No measurements found"));
    }

    #[test]
    fn test_validation_with_measurements() {
        let mut circuit = Circuit::new(1, 1);
        circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        let warnings = circuit.validate();
        assert!(warnings.is_empty());
    }

    #[test]
    fn test_to_qasm() {
        let mut circuit = Circuit::new(2, 2);
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
        circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });

        let qasm = circuit.to_qasm(None);
        assert!(qasm.contains("OPENQASM 2.0;"));
        assert!(qasm.contains("qreg q[2];"));
        assert!(qasm.contains("creg c[2];"));
        assert!(qasm.contains("h q[0];"));
        assert!(qasm.contains("cx q[0], q[1];"));
        assert!(qasm.contains("measure q[0] -> c[0];"));
    }

    #[test]
    fn test_to_qasm_with_layout() {
        use crate::transpiler::property_set::PropertySet;
        let circuit = Circuit::new(2, 0);
        let mut props = PropertySet::new();
        props.insert("initial_layout", vec![1usize, 0usize]);
        props.insert("final_layout", vec![0usize, 1usize]);

        let qasm = circuit.to_qasm(Some(&props));
        assert!(qasm.contains("// qrust_initial_layout: [1, 0]"));
        assert!(qasm.contains("// qrust_final_layout: [0, 1]"));
    }

    #[test]
    fn test_gate_count_mixed_ops() {
        let mut circuit = Circuit::new(2, 2);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Barrier { qubits: vec![0, 1] });
        circuit.add_op(Operation::Reset { qubit: 0 });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });

        // Gates: H, CX. Measurements, Barriers, and Resets are excluded.
        assert_eq!(circuit.gate_count(), 2);
    }

    #[test]
    fn test_circuit_depth_complex_dependencies() {
        let mut circuit = Circuit::new(4, 0);

        // Parallel Layer 1
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![2],
            params: vec![],
        });

        // Dependent Layer 2
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![2, 3],
            params: vec![],
        });

        // Cross-Dependency Layer 3
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 2],
            params: vec![],
        });

        // Final Depth Layer 4
        circuit.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![2],
            params: vec![],
        });

        assert_eq!(circuit.depth(), 4);
    }

    #[test]
    fn test_depth_synchronization_at_barrier() {
        let mut circuit = Circuit::new(3, 0);

        // Q0 reaches Depth 10
        for _ in 0..10 {
            circuit.add_op(Operation::Gate {
                name: GateType::X,
                qubits: vec![0],
                params: vec![],
            });
        }

        // Barrier synchronizes Q1 and Q2 (at Depth 0) up to Depth 10
        circuit.add_op(Operation::Barrier {
            qubits: vec![0, 1, 2],
        });

        // Next operation on Q1 starts from Depth 10 + 1
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });

        assert_eq!(circuit.depth(), 11);
    }
}
