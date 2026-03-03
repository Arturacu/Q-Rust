use crate::ir::{Circuit, GateType, Operation};
use std::f64::consts::PI;

/// Decomposes a circuit into the basis gate set: u1, u2, u3, CX.
///
/// All other gates are expanded into sequences of these basis gates.
pub fn decompose_basis(circuit: &Circuit) -> Circuit {
    let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);

    for op in &circuit.operations {
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                expand_gate(&mut new_circuit, name, qubits, params);
            }
            // Preserve other operations (Measure, Reset, Barrier)
            _ => new_circuit.add_op(op.clone()),
        }
    }

    new_circuit
}

/// Expands a single gate into a sequence of basis gates (U, CX) and adds them to the circuit.
///
/// # Arguments
///
/// * `circuit` - The circuit to add the expanded gates to.
/// * `name` - The type of the gate to expand.
/// * `qubits` - The indices of the qubits the gate acts on.
/// * `params` - The parameters of the gate (if any).
fn expand_gate(circuit: &mut Circuit, name: &GateType, qubits: &[usize], params: &[f64]) {
    match name {
        // Basis gates - keep as is
        GateType::U | GateType::CX => {
            circuit.add_op(Operation::Gate {
                name: name.clone(),
                qubits: qubits.to_vec(),
                params: params.to_vec(),
            });
        }

        // Identity
        GateType::ID => {
            // id q -> u1(0) q
            circuit.add_op(Operation::Gate {
                name: GateType::U, // u1(0) is U(0, 0, 0) effectively or U(0, 0, lambda) with lambda=0
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, 0.0],
            });
        }

        // Pauli Gates
        GateType::X => {
            // x q -> u3(pi, 0, pi) q
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![PI, 0.0, PI],
            });
        }
        GateType::Y => {
            // y q -> u3(pi, pi/2, pi/2) q
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![PI, PI / 2.0, PI / 2.0],
            });
        }
        GateType::Z => {
            // z q -> u1(pi) q -> U(0, 0, pi)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, PI],
            });
        }

        // Hadamard
        GateType::H => {
            // h q -> u2(0, pi) q -> U(pi/2, 0, pi)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![PI / 2.0, 0.0, PI],
            });
        }

        // Phase Gates
        GateType::S => {
            // s q -> u1(pi/2) q -> U(0, 0, pi/2)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, PI / 2.0],
            });
        }
        GateType::Sdg => {
            // sdg q -> u1(-pi/2) q -> U(0, 0, -pi/2)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, -PI / 2.0],
            });
        }
        GateType::T => {
            // t q -> u1(pi/4) q -> U(0, 0, pi/4)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, PI / 4.0],
            });
        }
        GateType::Tdg => {
            // tdg q -> u1(-pi/4) q -> U(0, 0, -pi/4)
            circuit.add_op(Operation::Gate {
                name: GateType::U,
                qubits: qubits.to_vec(),
                params: vec![0.0, 0.0, -PI / 4.0],
            });
        }

        // Rotation Gates
        GateType::RX => {
            // rx(theta) q -> u3(theta, -pi/2, pi/2) q
            if let Some(&theta) = params.first() {
                circuit.add_op(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![theta, -PI / 2.0, PI / 2.0],
                });
            }
        }
        GateType::RY => {
            // ry(theta) q -> u3(theta, 0, 0) q
            if let Some(&theta) = params.first() {
                circuit.add_op(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![theta, 0.0, 0.0],
                });
            }
        }
        GateType::RZ => {
            // rz(phi) q -> u1(phi) q -> U(0, 0, phi)
            if let Some(&phi) = params.first() {
                circuit.add_op(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, phi],
                });
            }
        }

        // Multi-Qubit Gates
        GateType::SWAP => {
            // swap a, b -> cx a,b; cx b,a; cx a,b;
            let a = qubits[0];
            let b = qubits[1];

            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, b],
                params: vec![],
            });
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![b, a],
                params: vec![],
            });
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, b],
                params: vec![],
            });
        }

        // Toffoli (CCX) decomposition
        GateType::CCX => {
            // Standard decomposition into H, CX, T, Tdg:
            // h c; cx b,c; tdg c; cx a,c; t c; cx b,c; tdg c; cx a,c; t b; t c; h c; cx a,b; t a; tdg b; cx a,b;
            let a = qubits[0];
            let b = qubits[1];
            let c = qubits[2];

            // h c
            expand_gate(circuit, &GateType::H, &[c], &[]);
            // cx b,c
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![b, c],
                params: vec![],
            });
            // tdg c
            expand_gate(circuit, &GateType::Tdg, &[c], &[]);
            // cx a,c
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, c],
                params: vec![],
            });
            // t c
            expand_gate(circuit, &GateType::T, &[c], &[]);
            // cx b,c
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![b, c],
                params: vec![],
            });
            // tdg c
            expand_gate(circuit, &GateType::Tdg, &[c], &[]);
            // cx a,c
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, c],
                params: vec![],
            });
            // t b
            expand_gate(circuit, &GateType::T, &[b], &[]);
            // t c
            expand_gate(circuit, &GateType::T, &[c], &[]);
            // h c
            expand_gate(circuit, &GateType::H, &[c], &[]);
            // cx a,b
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, b],
                params: vec![],
            });
            // t a
            expand_gate(circuit, &GateType::T, &[a], &[]);
            // tdg b
            expand_gate(circuit, &GateType::Tdg, &[b], &[]);
            // cx a,b
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, b],
                params: vec![],
            });
        }

        GateType::Custom(_) => {
            // Custom gates should have been expanded during parsing/IR generation ideally.
            // If they reach here, we can't decompose them without their definition.
            // For now, we keep them.
            circuit.add_op(Operation::Gate {
                name: name.clone(),
                qubits: qubits.to_vec(),
                params: params.to_vec(),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::gates::GateType;

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
