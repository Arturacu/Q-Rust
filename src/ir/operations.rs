//! Circuit operations (gates, measurements, resets, barriers, conditionals).

use super::gates::GateType;
use std::fmt;

/// A classical condition used to guard a conditional operation.
///
/// Represents `if (creg == value) op` from QASM 2.0.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde-ir", derive(serde::Serialize, serde::Deserialize))]
pub struct ClassicalCondition {
    /// Classical register name.
    pub creg: String,
    /// Required decimal value for the register contents.
    pub value: u64,
}

/// A single operation in a [`crate::ir::Circuit`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde-ir", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Operation {
    Gate {
        name: GateType,
        qubits: Vec<usize>,
        params: Vec<f64>,
    },
    Measure {
        qubit: usize,
        cbit: usize,
    },
    Reset {
        qubit: usize,
    },
    Barrier {
        qubits: Vec<usize>,
    },
    /// Classically-conditioned inner operation (QASM 2.0 `if (c==v) op`).
    Conditional {
        condition: ClassicalCondition,
        op: Box<Operation>,
    },
}

impl Operation {
    fn write_qasm<W: fmt::Write>(&self, w: &mut W) -> fmt::Result {
        match self {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                w.write_str(name.to_qasm_name())?;
                if !params.is_empty() {
                    w.write_char('(')?;
                    for (i, p) in params.iter().enumerate() {
                        if i > 0 {
                            w.write_str(", ")?;
                        }
                        write!(w, "{:.10}", p)?;
                    }
                    w.write_char(')')?;
                }
                w.write_char(' ')?;
                for (i, q) in qubits.iter().enumerate() {
                    if i > 0 {
                        w.write_str(", ")?;
                    }
                    write!(w, "q[{}]", q)?;
                }
                w.write_char(';')
            }
            Operation::Measure { qubit, cbit } => {
                write!(w, "measure q[{}] -> c[{}];", qubit, cbit)
            }
            Operation::Reset { qubit } => write!(w, "reset q[{}];", qubit),
            Operation::Barrier { qubits } => {
                w.write_str("barrier ")?;
                for (i, q) in qubits.iter().enumerate() {
                    if i > 0 {
                        w.write_str(", ")?;
                    }
                    write!(w, "q[{}]", q)?;
                }
                w.write_char(';')
            }
            Operation::Conditional { condition, op } => {
                write!(w, "if({}=={}) ", condition.creg, condition.value)?;
                op.write_qasm(w)
            }
        }
    }

    pub fn to_qasm(&self) -> String {
        let mut s = String::with_capacity(32);
        // Writing into a String is infallible.
        let _ = self.write_qasm(&mut s);
        s
    }

    pub fn qubits(&self) -> &[usize] {
        match self {
            Operation::Gate { qubits, .. } | Operation::Barrier { qubits } => qubits,
            Operation::Measure { qubit, .. } | Operation::Reset { qubit } => {
                std::slice::from_ref(qubit)
            }
            Operation::Conditional { op, .. } => op.qubits(),
        }
    }

    /// Returns true if this operation is a barrier. Optimization passes
    /// must not reorder/combine operations across barriers.
    #[inline]
    pub fn is_barrier(&self) -> bool {
        matches!(self, Operation::Barrier { .. })
    }

    /// Returns true if this operation is classically conditioned.
    #[inline]
    pub fn is_conditional(&self) -> bool {
        matches!(self, Operation::Conditional { .. })
    }
}

impl fmt::Display for Operation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_qasm(f)
    }
}
