use super::gates::GateType;

/// Represents a single operation in the quantum circuit.
///
/// Operations can be quantum gates, measurements, resets, or barriers.
#[derive(Debug, Clone, PartialEq)]
pub enum Operation {
    /// A quantum gate application.
    Gate {
        /// Type of the gate (e.g., H, CX).
        name: GateType,
        /// Indices of the qubits involved.
        qubits: Vec<usize>,
        /// Parameters for the gate (if any).
        params: Vec<f64>,
    },
    /// A measurement operation.
    Measure {
        /// Index of the qubit to measure.
        qubit: usize,
        /// Index of the classical bit to store the result.
        cbit: usize,
    },
    /// Reset a qubit to the |0> state.
    Reset {
        /// Index of the qubit to reset.
        qubit: usize,
    },
    /// A barrier to prevent optimizations across a boundary.
    Barrier {
        /// Indices of the qubits involved in the barrier.
        qubits: Vec<usize>,
    },
}
