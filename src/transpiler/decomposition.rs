use crate::ir::{Circuit, GateDefinition, GateType, Operation};

/// Decomposes a circuit into the basis gate set: U, CX.
///
/// All other gates are expanded into sequences of these basis gates
/// using the decomposition rules defined in `GateDefinition::decompose()`.
pub fn decompose_basis(circuit: &Circuit) -> Circuit {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    result.custom_gates = circuit.custom_gates.clone();

    for op in &circuit.operations {
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                expand_gate(&mut result, &circuit.custom_gates, name, qubits, params);
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

/// Unrolls ONLY custom gates, leaving standard gates untouched.
pub fn unroll_custom_gates(circuit: &Circuit) -> Circuit {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    result.custom_gates = circuit.custom_gates.clone();

    for op in &circuit.operations {
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                if matches!(name, GateType::Custom(_)) {
                    expand_gate_custom_only(
                        &mut result,
                        &circuit.custom_gates,
                        name,
                        qubits,
                        params,
                    );
                } else {
                    result.add_op(op.clone());
                }
            }
            _ => result.add_op(op.clone()),
        }
    }
    result
}

fn expand_gate_custom_only(
    circuit: &mut Circuit,
    registry: &crate::ir::registry::GateRegistry,
    name: &GateType,
    qubits: &[usize],
    params: &[f64],
) {
    if let GateType::Custom(ref custom_name) = name {
        if let Some(def) = registry.get(custom_name) {
            let mut scope_params = std::collections::HashMap::new();
            for (p_name, &val) in def.params.iter().zip(params.iter()) {
                scope_params.insert(p_name.clone(), val);
            }

            let mut scope_qubits = std::collections::HashMap::new();
            for (q_name, &idx) in def.qubits.iter().zip(qubits.iter()) {
                scope_qubits.insert(q_name.clone(), idx);
            }

            for stmt in &def.body {
                if let crate::ir::ast::ParsedStatement::Gate(
                    inner_name_str,
                    inner_qubits,
                    inner_params,
                ) = stmt
                {
                    let mut resolved_qubits = Vec::new();
                    for (q_reg, _) in inner_qubits {
                        if let Some(&resolved) = scope_qubits.get(q_reg) {
                            resolved_qubits.push(resolved);
                        } else {
                            panic!("Unknown qubit {} in custom gate {}", q_reg, custom_name);
                        }
                    }

                    let mut resolved_params = Vec::new();
                    for p in inner_params {
                        resolved_params.push(
                            p.evaluate_with_scope(&scope_params)
                                .expect("Failed to evaluate parameter"),
                        );
                    }

                    let inner_gatetype = inner_name_str
                        .parse::<GateType>()
                        .unwrap_or_else(|_| GateType::Custom(inner_name_str.clone()));

                    if matches!(inner_gatetype, GateType::Custom(_)) {
                        expand_gate_custom_only(
                            circuit,
                            registry,
                            &inner_gatetype,
                            &resolved_qubits,
                            &resolved_params,
                        );
                    } else {
                        circuit.add_op(Operation::Gate {
                            name: inner_gatetype,
                            qubits: resolved_qubits,
                            params: resolved_params,
                        });
                    }
                }
            }
        } else {
            panic!("Custom gate {} not found in registry", custom_name);
        }
    }
}

/// Expands a single gate into a sequence of basis gates (U, CX) and adds them to the circuit.
///
/// Uses `GateDefinition::decompose()` as the single source of truth for
/// decomposition rules. Results are recursively expanded until only
/// basis gates remain.
fn expand_gate(
    circuit: &mut Circuit,
    registry: &crate::ir::registry::GateRegistry,
    name: &GateType,
    qubits: &[usize],
    params: &[f64],
) {
    if let GateType::Custom(ref custom_name) = name {
        if let Some(def) = registry.get(custom_name) {
            let mut scope_params = std::collections::HashMap::new();
            for (p_name, &val) in def.params.iter().zip(params.iter()) {
                scope_params.insert(p_name.clone(), val);
            }

            let mut scope_qubits = std::collections::HashMap::new();
            for (q_name, &idx) in def.qubits.iter().zip(qubits.iter()) {
                scope_qubits.insert(q_name.clone(), idx);
            }

            for stmt in &def.body {
                if let crate::ir::ast::ParsedStatement::Gate(
                    inner_name_str,
                    inner_qubits,
                    inner_params,
                ) = stmt
                {
                    let mut resolved_qubits = Vec::new();
                    for (q_reg, _) in inner_qubits {
                        if let Some(&resolved) = scope_qubits.get(q_reg) {
                            resolved_qubits.push(resolved);
                        } else {
                            panic!("Unknown qubit {} in custom gate {}", q_reg, custom_name);
                        }
                    }

                    let mut resolved_params = Vec::new();
                    for p in inner_params {
                        resolved_params.push(
                            p.evaluate_with_scope(&scope_params)
                                .expect("Failed to evaluate parameter"),
                        );
                    }

                    let inner_gatetype = inner_name_str
                        .parse::<GateType>()
                        .unwrap_or_else(|_| GateType::Custom(inner_name_str.clone()));

                    expand_gate(
                        circuit,
                        registry,
                        &inner_gatetype,
                        &resolved_qubits,
                        &resolved_params,
                    );
                }
            }
        } else {
            panic!("Custom gate {} not found in registry", custom_name);
        }
    } else {
        match name.decompose(qubits, params) {
            None => {
                circuit.add_op(Operation::Gate {
                    name: name.clone(),
                    qubits: qubits.to_vec(),
                    params: params.to_vec(),
                });
            }
            Some(ops) => {
                for op in ops {
                    if let Operation::Gate {
                        name: sub_name,
                        qubits: sub_qubits,
                        params: sub_params,
                    } = op
                    {
                        if sub_name.is_basis() {
                            circuit.add_op(Operation::Gate {
                                name: sub_name,
                                qubits: sub_qubits,
                                params: sub_params,
                            });
                        } else {
                            expand_gate(circuit, registry, &sub_name, &sub_qubits, &sub_params);
                        }
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

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        decompose_basis(circuit)
    }
}
