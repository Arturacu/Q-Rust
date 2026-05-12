//! Intermediate representation of a quantum circuit.

use super::ast::ParsedStatement;
use super::gates::GateType;
use super::operations::Operation;
use super::registry::GateRegistry;
use std::collections::HashMap;
use std::fmt;

#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde-ir", derive(serde::Serialize, serde::Deserialize))]
pub struct Circuit {
    pub num_qubits: usize,
    pub num_cbits: usize,
    pub operations: Vec<Operation>,
    pub custom_gates: GateRegistry,
}

impl Circuit {
    pub fn new(num_qubits: usize, num_cbits: usize) -> Self {
        Self {
            num_qubits,
            num_cbits,
            operations: Vec::new(),
            custom_gates: GateRegistry::new(),
        }
    }

    pub fn register_custom_gate(
        &mut self,
        name: String,
        params: Vec<String>,
        qubits: Vec<String>,
        body: Vec<ParsedStatement>,
    ) {
        self.custom_gates.register(name, params, qubits, body);
    }

    #[inline]
    pub fn add_op(&mut self, op: Operation) {
        self.operations.push(op);
    }

    pub fn validate(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        let has_measurement = self
            .operations
            .iter()
            .any(|op| matches!(op, Operation::Measure { .. }));
        if !has_measurement {
            warnings.push(
                "Warning: No measurements found. The circuit will not produce \
                 classical output on hardware."
                    .to_string(),
            );
        }
        warnings
    }

    pub fn to_qasm(
        &self,
        property_set: Option<&crate::transpiler::property_set::PropertySet>,
    ) -> String {
        let mut qasm = String::with_capacity(256 + self.operations.len() * 24);
        qasm.push_str("OPENQASM 2.0;\n");
        qasm.push_str("include \"qelib1.inc\";\n\n");

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

    /// Returns the depth (critical-path length). Barriers bump each affected
    /// wire's depth tracker to the max so subsequent ops line up after the
    /// barrier, but barriers themselves don't add depth.
    pub fn depth(&self) -> usize {
        let mut qd = vec![0usize; self.num_qubits];
        for op in &self.operations {
            match op {
                Operation::Gate { qubits, .. } => {
                    let m = qubits
                        .iter()
                        .filter_map(|&q| qd.get(q).copied())
                        .max()
                        .unwrap_or(0);
                    let new_d = m + 1;
                    for &q in qubits {
                        if let Some(slot) = qd.get_mut(q) {
                            *slot = new_d;
                        }
                    }
                }
                Operation::Measure { qubit, .. } | Operation::Reset { qubit } => {
                    if let Some(slot) = qd.get_mut(*qubit) {
                        *slot += 1;
                    }
                }
                Operation::Barrier { qubits } => {
                    let m = qubits
                        .iter()
                        .filter_map(|&q| qd.get(q).copied())
                        .max()
                        .unwrap_or(0);
                    for &q in qubits {
                        if let Some(slot) = qd.get_mut(q) {
                            *slot = m;
                        }
                    }
                }
                Operation::Conditional { op, .. } => {
                    let qubits = op.qubits();
                    let m = qubits
                        .iter()
                        .filter_map(|&q| qd.get(q).copied())
                        .max()
                        .unwrap_or(0);
                    let new_d = m + 1;
                    for &q in qubits {
                        if let Some(slot) = qd.get_mut(q) {
                            *slot = new_d;
                        }
                    }
                }
            }
        }
        qd.into_iter().max().unwrap_or(0)
    }

    /// Returns the number of [`Operation::Gate`] operations.
    pub fn gate_count(&self) -> usize {
        self.operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { .. }))
            .count()
    }

    /// Counts occurrences of each gate kind in the circuit.
    ///
    /// - [`Operation::Gate`] is tallied under its `GateType`.
    /// - [`Operation::Barrier`] is tallied under [`GateType::Barrier`].
    /// - [`Operation::Conditional`] is recursed into.
    /// - Measurement/Reset are not counted.
    pub fn count_ops(&self) -> HashMap<GateType, usize> {
        let mut counts: HashMap<GateType, usize> = HashMap::new();
        fn walk(op: &Operation, counts: &mut HashMap<GateType, usize>) {
            match op {
                Operation::Gate { name, .. } => {
                    *counts.entry(name.clone()).or_insert(0) += 1;
                }
                Operation::Barrier { .. } => {
                    *counts.entry(GateType::Barrier).or_insert(0) += 1;
                }
                Operation::Conditional { op, .. } => walk(op, counts),
                _ => {}
            }
        }
        for op in &self.operations {
            walk(op, &mut counts);
        }
        counts
    }
}

impl fmt::Display for Circuit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Circuit({} qubits, {} cbits, {} ops, depth={})",
            self.num_qubits,
            self.num_cbits,
            self.operations.len(),
            self.depth()
        )?;
        for (i, op) in self.operations.iter().enumerate() {
            writeln!(f, "  {:>3}: {}", i, op)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::gates::GateType;

    #[test]
    fn test_circuit_creation() {
        let c = Circuit::new(2, 2);
        assert_eq!(c.num_qubits, 2);
        assert_eq!(c.num_cbits, 2);
        assert!(c.operations.is_empty());
    }

    #[test]
    fn test_add_op() {
        let mut c = Circuit::new(1, 0);
        let op = Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        };
        c.add_op(op.clone());
        assert_eq!(c.operations, vec![op]);
    }

    #[test]
    fn test_depth_simple() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        assert_eq!(c.depth(), 2);
    }

    #[test]
    fn test_count_ops() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Barrier { qubits: vec![0, 1] });
        let counts = c.count_ops();
        assert_eq!(counts.get(&GateType::H), Some(&2));
        assert_eq!(counts.get(&GateType::CX), Some(&1));
        assert_eq!(counts.get(&GateType::Barrier), Some(&1));
    }

    #[test]
    fn test_validation_with_measurements() {
        let mut c = Circuit::new(1, 1);
        c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        assert!(c.validate().is_empty());
    }

    #[test]
    fn test_to_qasm_with_barrier() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Barrier { qubits: vec![0, 1] });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let qasm = c.to_qasm(None);
        assert!(qasm.contains("barrier q[0], q[1];"));
    }
}
