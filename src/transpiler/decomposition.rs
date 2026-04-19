//! Basis decomposition pass with caching of custom gate expansions.

use crate::error::{QRustError, Result};
use crate::ir::ast::ParsedStatement;
use crate::ir::registry::GateRegistry;
use crate::ir::{Circuit, GateDefinition, GateType, Operation};
use std::collections::HashMap;

/// Cached expansion of a custom gate — stored as a symbolic template in
/// terms of formal parameter names (indices). To avoid retaining string
/// keys for every expansion, we currently cache per *(name, arity_pattern)*.
///
/// For the common case where a custom gate is invoked many times with
/// different concrete parameter values, we still need to re-evaluate the
/// expressions per invocation. The cache therefore stores the *pre-parsed,
/// symbolically resolved* sub-gate templates so we skip the recursive
/// walk through gate-def bodies on repeated invocations.
#[derive(Clone)]
struct CachedTemplate {
    /// Each entry: (GateType, qubit-name indices into def's formal qubit
    /// list, param AST references cloned from the body).
    steps: Vec<(GateType, Vec<usize>, Vec<crate::ir::ast::Expr>)>,
}

pub fn try_decompose_basis(circuit: &Circuit) -> Result<Circuit> {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    result.custom_gates = circuit.custom_gates.clone();
    let mut cache: HashMap<String, CachedTemplate> = HashMap::new();

    for op in &circuit.operations {
        expand_op(&mut result, &circuit.custom_gates, op, &mut cache)?;
    }
    Ok(result)
}

pub fn decompose_basis(circuit: &Circuit) -> Circuit {
    try_decompose_basis(circuit).expect("decompose_basis failed")
}

pub fn try_unroll_custom_gates(circuit: &Circuit) -> Result<Circuit> {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    result.custom_gates = circuit.custom_gates.clone();
    let mut cache: HashMap<String, CachedTemplate> = HashMap::new();

    for op in &circuit.operations {
        match op {
            Operation::Gate { name, qubits, params } if matches!(name, GateType::Custom(_)) => {
                expand_gate_custom_only(
                    &mut result,
                    &circuit.custom_gates,
                    name,
                    qubits,
                    params,
                    &mut cache,
                )?;
            }
            other => result.add_op(other.clone()),
        }
    }
    Ok(result)
}

pub fn unroll_custom_gates(circuit: &Circuit) -> Circuit {
    try_unroll_custom_gates(circuit).unwrap_or_else(|_| circuit.clone())
}

fn expand_op(
    circuit: &mut Circuit,
    registry: &GateRegistry,
    op: &Operation,
    cache: &mut HashMap<String, CachedTemplate>,
) -> Result<()> {
    match op {
        Operation::Gate { name, qubits, params } => {
            expand_gate(circuit, registry, name, qubits, params, cache)
        }
        Operation::Conditional { condition, op: inner } => {
            // Expand the inner op into a fresh sub-circuit and re-wrap each
            // resulting basis gate under the same condition.
            let mut tmp = Circuit::new(circuit.num_qubits, circuit.num_cbits);
            expand_op(&mut tmp, registry, inner, cache)?;
            for sub in tmp.operations {
                circuit.add_op(Operation::Conditional {
                    condition: condition.clone(),
                    op: Box::new(sub),
                });
            }
            Ok(())
        }
        other => {
            circuit.add_op(other.clone());
            Ok(())
        }
    }
}

fn get_or_build_template(
    registry: &GateRegistry,
    name: &str,
    cache: &mut HashMap<String, CachedTemplate>,
) -> Result<CachedTemplate> {
    if let Some(t) = cache.get(name) {
        return Ok(t.clone());
    }
    let def = registry.get(name).ok_or_else(|| {
        QRustError::Decomposition(format!("custom gate '{}' not in registry", name))
    })?;
    let qubit_index: HashMap<&str, usize> = def
        .qubits
        .iter()
        .enumerate()
        .map(|(i, q)| (q.as_str(), i))
        .collect();

    let mut steps = Vec::new();
    for stmt in &def.body {
        if let ParsedStatement::Gate(inner_name_str, inner_qubits, inner_params) = stmt {
            let mut qidx = Vec::with_capacity(inner_qubits.len());
            for (q_reg, _) in inner_qubits {
                match qubit_index.get(q_reg.as_str()) {
                    Some(&i) => qidx.push(i),
                    None => {
                        return Err(QRustError::Decomposition(format!(
                            "unknown qubit '{}' in custom gate '{}'",
                            q_reg, name
                        )))
                    }
                }
            }
            let gt = inner_name_str
                .parse::<GateType>()
                .unwrap_or_else(|_| GateType::Custom(inner_name_str.clone()));
            steps.push((gt, qidx, inner_params.clone()));
        }
    }
    let t = CachedTemplate { steps };
    cache.insert(name.to_string(), t.clone());
    Ok(t)
}

fn expand_gate_custom_only(
    circuit: &mut Circuit,
    registry: &GateRegistry,
    name: &GateType,
    qubits: &[usize],
    params: &[f64],
    cache: &mut HashMap<String, CachedTemplate>,
) -> Result<()> {
    let GateType::Custom(ref custom_name) = name else {
        return Ok(());
    };
    let def = registry.get(custom_name).ok_or_else(|| {
        QRustError::Decomposition(format!("custom gate '{}' not in registry", custom_name))
    })?;
    let param_names: Vec<String> = def.params.clone();

    let template = get_or_build_template(registry, custom_name, cache)?;

    let mut scope_params = HashMap::with_capacity(param_names.len());
    for (p_name, &val) in param_names.iter().zip(params.iter()) {
        scope_params.insert(p_name.clone(), val);
    }

    for (gate_type, qidx, params_expr) in &template.steps {
        let resolved_qubits: Vec<usize> = qidx.iter().map(|&i| qubits[i]).collect();
        let mut resolved_params = Vec::with_capacity(params_expr.len());
        for p in params_expr {
            resolved_params.push(p.evaluate_with_scope(&scope_params)?);
        }
        if matches!(gate_type, GateType::Custom(_)) {
            expand_gate_custom_only(circuit, registry, gate_type, &resolved_qubits, &resolved_params, cache)?;
        } else {
            circuit.add_op(Operation::Gate {
                name: gate_type.clone(),
                qubits: resolved_qubits,
                params: resolved_params,
            });
        }
    }
    Ok(())
}

fn expand_gate(
    circuit: &mut Circuit,
    registry: &GateRegistry,
    name: &GateType,
    qubits: &[usize],
    params: &[f64],
    cache: &mut HashMap<String, CachedTemplate>,
) -> Result<()> {
    if let GateType::Custom(ref custom_name) = name {
        let def = registry.get(custom_name).ok_or_else(|| {
            QRustError::Decomposition(format!("custom gate '{}' not in registry", custom_name))
        })?;
        let param_names: Vec<String> = def.params.clone();
        let template = get_or_build_template(registry, custom_name, cache)?;

        let mut scope_params = HashMap::with_capacity(param_names.len());
        for (p_name, &val) in param_names.iter().zip(params.iter()) {
            scope_params.insert(p_name.clone(), val);
        }

        for (gate_type, qidx, params_expr) in &template.steps {
            let resolved_qubits: Vec<usize> = qidx.iter().map(|&i| qubits[i]).collect();
            let mut resolved_params = Vec::with_capacity(params_expr.len());
            for p in params_expr {
                resolved_params.push(p.evaluate_with_scope(&scope_params)?);
            }
            expand_gate(
                circuit,
                registry,
                gate_type,
                &resolved_qubits,
                &resolved_params,
                cache,
            )?;
        }
        return Ok(());
    }

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
                        expand_gate(circuit, registry, &sub_name, &sub_qubits, &sub_params, cache)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// Rewrites a CX whose (control, target) direction violates the backend
/// coupling map — but the reverse edge exists — by sandwiching with H on
/// both wires: `CX(a,b) = (H⊗H) · CX(b,a) · (H⊗H)`.
pub fn reorient_cx_for_coupling(
    circuit: &Circuit,
    backend: &crate::backend::Backend,
) -> Circuit {
    let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    out.custom_gates = circuit.custom_gates.clone();

    for op in &circuit.operations {
        match op {
            Operation::Gate { name: GateType::CX, qubits, params } if qubits.len() == 2 => {
                let (c, t) = (qubits[0], qubits[1]);
                if backend.has_directed_edge(c, t) {
                    out.add_op(op.clone());
                } else if backend.has_directed_edge(t, c) {
                    // Insert H·H CX(t,c) H·H
                    out.add_op(Operation::Gate { name: GateType::H, qubits: vec![c], params: vec![] });
                    out.add_op(Operation::Gate { name: GateType::H, qubits: vec![t], params: vec![] });
                    out.add_op(Operation::Gate { name: GateType::CX, qubits: vec![t, c], params: params.clone() });
                    out.add_op(Operation::Gate { name: GateType::H, qubits: vec![c], params: vec![] });
                    out.add_op(Operation::Gate { name: GateType::H, qubits: vec![t], params: vec![] });
                } else {
                    // Neither direction available; leave as-is (routing will handle).
                    out.add_op(op.clone());
                }
            }
            other => out.add_op(other.clone()),
        }
    }
    out
}

#[derive(Debug, Clone, Copy)]
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
        try_decompose_basis(circuit).unwrap_or_else(|_| circuit.clone())
    }
}

/// Pass that rewrites wrong-direction CX gates relative to a coupling map.
#[derive(Debug, Clone)]
pub struct CxDirectionPass {
    pub backend: crate::backend::Backend,
}

impl crate::transpiler::pass::Pass for CxDirectionPass {
    fn name(&self) -> &str {
        "CxDirectionPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        reorient_cx_for_coupling(circuit, &self.backend)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::gates::GateType;
    use std::f64::consts::PI;

    #[test]
    fn test_decompose_h() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let d = decompose_basis(&c);
        assert_eq!(d.operations.len(), 1);
        match &d.operations[0] {
            Operation::Gate { name: GateType::U, params, .. } => {
                assert!((params[0] - PI / 2.0).abs() < 1e-6);
                assert!(params[1].abs() < 1e-6);
                assert!((params[2] - PI).abs() < 1e-6);
            }
            _ => panic!("expected U gate"),
        }
    }

    #[test]
    fn test_decompose_swap() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        let d = decompose_basis(&c);
        assert_eq!(d.operations.len(), 3);
    }

    #[test]
    fn test_missing_custom_gate_returns_err() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::Custom("nope".into()),
            qubits: vec![0],
            params: vec![],
        });
        assert!(try_decompose_basis(&c).is_err());
    }

    #[test]
    fn test_preserves_barrier() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
        c.add_op(Operation::Barrier { qubits: vec![0, 1] });
        c.add_op(Operation::Gate { name: GateType::X, qubits: vec![1], params: vec![] });
        let d = decompose_basis(&c);
        assert!(d.operations.iter().any(|op| matches!(op, Operation::Barrier { .. })));
    }
}