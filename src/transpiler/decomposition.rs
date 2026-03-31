use crate::ir::{Circuit, GateDefinition, GateType, Operation};

/// Decomposes a circuit into the basis gate set: U, CX.
///
/// All other gates are expanded into sequences of these basis gates
/// using the decomposition rules defined in `GateDefinition::decompose()`.
pub fn decompose_basis(circuit: &Circuit) -> Circuit {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);

    for op in &circuit.operations {
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                expand_gate(&mut result, name, qubits, params);
            }
            Operation::Measure { qubit, cbit } => {
                result.add_op(Operation::Measure {
                    qubit: *qubit,
                    cbit: *cbit,
                });
            }
            _ => result.add_op(op.clone()),
        }
    }
    result
}

/// Expands a single gate into a sequence of basis gates (U, CX) and adds them to the circuit.
///
/// Uses `GateDefinition::decompose()` as the single source of truth for
/// decomposition rules. Results are recursively expanded until only
/// basis gates remain.
fn expand_gate(circuit: &mut Circuit, name: &GateType, qubits: &[usize], params: &[f64]) {
    match name.decompose(qubits, params) {
        None => {
            // Basis gate or non-decomposable: keep as-is
            circuit.add_op(Operation::Gate {
                name: name.clone(),
                qubits: qubits.to_vec(),
                params: params.to_vec(),
            });
        }
        Some(ops) => {
            // Recursively expand each resulting operation
            for op in ops {
                if let Operation::Gate {
                    name: sub_name,
                    qubits: sub_qubits,
                    params: sub_params,
                } = op
                {
                    if sub_name.is_basis() {
                        // Already a basis gate, add directly
                        circuit.add_op(Operation::Gate {
                            name: sub_name,
                            qubits: sub_qubits,
                            params: sub_params,
                        });
                    } else {
                        // Needs further decomposition
                        expand_gate(circuit, &sub_name, &sub_qubits, &sub_params);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::gates::GateType;
    use std::f64::consts::PI;

    #[test]
    fn test_decompose_h() {
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });

        let decomposed = decompose_basis(&circuit);
        assert_eq!(decomposed.operations.len(), 1);

        match &decomposed.operations[0] {
            Operation::Gate { name, params, .. } => {
                // H -> U(pi/2, 0, pi)
                if let GateType::U = name {
                    assert!((params[0] - PI / 2.0).abs() < 1e-6);
                    assert!(params[1].abs() < 1e-6);
                    assert!((params[2] - PI).abs() < 1e-6);
                } else {
                    panic!("Expected U gate");
                }
            }
            _ => panic!("Expected Gate operation"),
        }
    }

    #[test]
    fn test_decompose_swap() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let decomposed = decompose_basis(&circuit);
        assert_eq!(decomposed.operations.len(), 3);

        // Should be 3 CX gates
        for op in &decomposed.operations {
            match op {
                Operation::Gate { name, .. } => {
                    assert_eq!(*name, GateType::CX);
                }
                _ => panic!("Expected Gate operation"),
            }
        }
    }

    #[test]
    fn test_decompose_ccx() {
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CCX,
            qubits: vec![0, 1, 2],
            params: vec![],
        });

        let decomposed = decompose_basis(&circuit);

        // CCX decomposes into a sequence of H, CX, T, Tdg.
        // These are then further decomposed into U and CX.
        // The exact count depends on the decomposition, but it should be > 1
        // and contain only U and CX gates.
        assert!(decomposed.operations.len() > 1);

        for op in &decomposed.operations {
            match op {
                Operation::Gate { name, .. } => match name {
                    GateType::U | GateType::CX => {} // OK
                    _ => panic!("Found non-basis gate in CCX decomposition: {:?}", name),
                },
                _ => panic!("Expected Gate operation"),
            }
        }

        // Standard decomposition has 15 gates (6 CX + single qubit gates)
        // 6 CX + 9 single qubit gates (H, T, Tdg) -> 15 ops
        assert_eq!(decomposed.operations.len(), 15);
    }
}

/// A pass that decomposes all gates into the basis set (U, CX).
pub struct BasisDecompositionPass;

impl crate::transpiler::pass::Pass for BasisDecompositionPass {
    fn name(&self) -> &str {
        "BasisDecompositionPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        decompose_basis(circuit)
    }
}
