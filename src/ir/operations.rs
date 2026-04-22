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

impl Operation {
    /// Returns the OpenQASM 2.0 string representation of the operation.
    pub fn to_qasm(&self) -> String {
        match self {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                let mut s = name.to_qasm_name();
                if !params.is_empty() {
                    let p_strs: Vec<String> = params.iter().map(|p| format!("{:.6}", p)).collect();
                    s.push_str(&format!("({})", p_strs.join(", ")));
                }
                let q_strs: Vec<String> = qubits.iter().map(|q| format!("q[{}]", q)).collect();
                s.push_str(&format!(" {};", q_strs.join(", ")));
                s
            }
            Operation::Measure { qubit, cbit } => {
                format!("measure q[{}] -> c[{}];", qubit, cbit)
            }
            Operation::Reset { qubit } => {
                format!("reset q[{}];", qubit)
            }
            Operation::Barrier { qubits } => {
                let q_strs: Vec<String> = qubits.iter().map(|q| format!("q[{}]", q)).collect();
                format!("barrier {};", q_strs.join(", "))
            }
        }
    }
}
