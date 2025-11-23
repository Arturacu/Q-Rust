use super::operations::Operation;

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
        }
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
}
