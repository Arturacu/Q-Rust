//! OpenQASM 2.0 parser.

pub mod rules;

use self::rules::{comment, creg, gate_call, include, measure, openqasm_version, qreg};
use crate::error::{QRustError, Result};
use crate::ir::ast::{Expr, ParsedStatement};
use crate::ir::{ClassicalCondition, Circuit, GateType, Operation};
use nom::{branch::alt, character::complete::multispace0};
use std::collections::HashMap;
use std::f64::consts::PI;

fn resolve_argument(
    arg: &(String, Option<usize>),
    qregs: &HashMap<String, (usize, usize)>,
    scope_qubits: &HashMap<String, usize>,
) -> Result<Vec<usize>> {
    let (name, idx) = arg;
    if let Some(mapped_idx) = scope_qubits.get(name) {
        if idx.is_some() {
            return Err(QRustError::ParseError(format!(
                "cannot index qubit argument '{}' inside gate definition",
                name
            )));
        }
        Ok(vec![*mapped_idx])
    } else if let Some(&(start, size)) = qregs.get(name) {
        if let Some(i) = idx {
            if *i < size {
                Ok(vec![start + i])
            } else {
                Err(QRustError::IndexOutOfBounds {
                    name: name.clone(),
                    index: *i,
                    size,
                })
            }
        } else {
            Ok((0..size).map(|i| start + i).collect())
        }
    } else {
        Err(QRustError::Undefined(format!(
            "Undefined quantum register or argument: {}",
            name
        )))
    }
}

#[derive(Default)]
struct ParseContext {
    qregs: HashMap<String, (usize, usize)>,
    cregs: HashMap<String, (usize, usize)>,
    gate_defs: HashMap<String, (Vec<String>, Vec<String>, Vec<ParsedStatement>)>,
}

pub fn parse_qasm(input: &str) -> Result<Circuit> {
    let mut circuit = Circuit::new(0, 0);
    let mut ctx = ParseContext::default();
    let mut total_qubits = 0;
    let mut total_cbits = 0;

    let mut current = input;

    loop {
        let (rem, _) = multispace0::<&str, nom::error::Error<&str>>(current)
            .map_err(|e| QRustError::ParseError(e.to_string()))?;
        current = rem;
        if current.is_empty() {
            return Err(QRustError::ParseError(
                "Empty source or missing OPENQASM header".into(),
            ));
        }
        if let Ok((rem, _)) = comment(current) {
            current = rem;
            continue;
        }
        break;
    }

    let (rem, version) = openqasm_version(current).map_err(|_| {
        QRustError::ParseError(
            "Missing or invalid OPENQASM header. File must start with 'OPENQASM 2.0;'".into(),
        )
    })?;
    if version != "2.0" {
        return Err(QRustError::Unsupported(format!(
            "OpenQASM version '{}' (only '2.0' is supported)",
            version
        )));
    }
    current = rem;

    loop {
        let (rem, _) = multispace0::<&str, nom::error::Error<&str>>(current)
            .map_err(|e| QRustError::ParseError(e.to_string()))?;
        current = rem;
        if current.is_empty() {
            break;
        }
        if let Ok((rem, _)) = comment(current) {
            current = rem;
            continue;
        }

        let (rem, stmt) = alt((
            include,
            qreg,
            creg,
            measure,
            rules::barrier,
            rules::reset,
            rules::gate_def,
            rules::if_stmt,
            gate_call,
        ))(current)
        .map_err(|_| {
            let snippet: String = current.chars().take(60).collect();
            QRustError::ParseError(format!("parse error near: {}", snippet))
        })?;
        current = rem;

        handle_statement(&mut circuit, &mut ctx, &mut total_qubits, &mut total_cbits, stmt)?;
    }

    circuit.num_qubits = total_qubits;
    circuit.num_cbits = total_cbits;
    Ok(circuit)
}

fn handle_statement(
    circuit: &mut Circuit,
    ctx: &mut ParseContext,
    total_qubits: &mut usize,
    total_cbits: &mut usize,
    stmt: ParsedStatement,
) -> Result<()> {
    match stmt {
        ParsedStatement::Ignore => {}
        ParsedStatement::Include(filename) => {
            if filename != "qelib1.inc" {
                return Err(QRustError::Unsupported(format!(
                    "Includes are not supported. Please resolve all imports before parsing. \
                     Found: 'include \"{}\"'",
                    filename
                )));
            }
        }
        ParsedStatement::QReg(name, size) => {
            ctx.qregs.insert(name, (*total_qubits, size));
            *total_qubits += size;
        }
        ParsedStatement::CReg(name, size) => {
            ctx.cregs.insert(name, (*total_cbits, size));
            *total_cbits += size;
        }
        ParsedStatement::GateDef(name, params, qubits, body) => {
            circuit.register_custom_gate(
                name.clone(),
                params.clone(),
                qubits.clone(),
                body.clone(),
            );
            ctx.gate_defs.insert(name, (params, qubits, body));
        }
        ParsedStatement::Gate(name, qubits, params) => {
            // Intercept the reset sentinel from rules::reset.
            if name == "__reset__" {
                emit_reset(circuit, ctx, &qubits, None)?;
            } else {
                emit_resolved_gate_call(circuit, ctx, &name, &qubits, &params, None)?;
            }
        }
        ParsedStatement::Measure((q_name, q_idx), (c_name, c_idx)) => {
            emit_measure(circuit, ctx, &q_name, q_idx, &c_name, c_idx, None)?;
        }
        ParsedStatement::Barrier(args) => {
            let mut qubits: Vec<usize> = Vec::new();
            if args.is_empty() {
                qubits = (0..*total_qubits).collect();
            } else {
                for arg in &args {
                    let indices = resolve_argument(
                        &(arg.0.clone(), arg.1),
                        &ctx.qregs,
                        &HashMap::new(),
                    )?;
                    qubits.extend(indices);
                }
            }
            circuit.add_op(Operation::Barrier { qubits });
        }
        ParsedStatement::If(creg, value, inner) => {
            if !ctx.cregs.contains_key(&creg) {
                return Err(QRustError::Undefined(format!(
                    "Undefined classical register in `if`: {}",
                    creg
                )));
            }
            let condition = ClassicalCondition {
                creg,
                value: value as u64,
            };
            match *inner {
                ParsedStatement::Gate(name, qubits, params) => {
                    if name == "__reset__" {
                        emit_reset(circuit, ctx, &qubits, Some(condition))?;
                    } else {
                        emit_resolved_gate_call(
                            circuit,
                            ctx,
                            &name,
                            &qubits,
                            &params,
                            Some(condition),
                        )?;
                    }
                }
                ParsedStatement::Measure((q_name, q_idx), (c_name, c_idx)) => {
                    emit_measure(
                        circuit, ctx, &q_name, q_idx, &c_name, c_idx,
                        Some(condition),
                    )?;
                }
                ParsedStatement::Barrier(_) => {
                    return Err(QRustError::Unsupported(
                        "`if` over barrier is not supported".into(),
                    ));
                }
                other => {
                    return Err(QRustError::Unsupported(format!(
                        "`if` over {:?} is not supported",
                        other
                    )));
                }
            }
        }
    }
    Ok(())
}

fn emit_reset(
    circuit: &mut Circuit,
    ctx: &ParseContext,
    qubits: &[(String, Option<usize>)],
    condition: Option<ClassicalCondition>,
) -> Result<()> {
    let mut all: Vec<usize> = Vec::new();
    for arg in qubits {
        let indices =
            resolve_argument(&(arg.0.clone(), arg.1), &ctx.qregs, &HashMap::new())?;
        all.extend(indices);
    }
    for q in all {
        let op = Operation::Reset { qubit: q };
        match &condition {
            Some(cond) => circuit.add_op(Operation::Conditional {
                condition: cond.clone(),
                op: Box::new(op),
            }),
            None => circuit.add_op(op),
        }
    }
    Ok(())
}

fn emit_measure(
    circuit: &mut Circuit,
    ctx: &ParseContext,
    q_name: &str,
    q_idx: Option<usize>,
    c_name: &str,
    c_idx: Option<usize>,
    condition: Option<ClassicalCondition>,
) -> Result<()> {
    let q_indices =
        resolve_argument(&(q_name.to_string(), q_idx), &ctx.qregs, &HashMap::new())?;
    let c_indices = if let Some(&(start, size)) = ctx.cregs.get(c_name) {
        if let Some(i) = c_idx {
            if i < size {
                vec![start + i]
            } else {
                return Err(QRustError::IndexOutOfBounds {
                    name: c_name.to_string(),
                    index: i,
                    size,
                });
            }
        } else {
            (0..size).map(|i| start + i).collect()
        }
    } else {
        return Err(QRustError::Undefined(format!(
            "Undefined classical register: {}",
            c_name
        )));
    };
    if q_indices.len() != c_indices.len() {
        return Err(QRustError::SizeMismatch(
            "measure register size mismatch".into(),
        ));
    }
    for (q, c) in q_indices.into_iter().zip(c_indices) {
        let op = Operation::Measure { qubit: q, cbit: c };
        match &condition {
            Some(cond) => circuit.add_op(Operation::Conditional {
                condition: cond.clone(),
                op: Box::new(op),
            }),
            None => circuit.add_op(op),
        }
    }
    Ok(())
}

fn emit_resolved_gate_call(
    circuit: &mut Circuit,
    ctx: &ParseContext,
    name: &str,
    qubits: &[(String, Option<usize>)],
    params: &[Expr],
    condition: Option<ClassicalCondition>,
) -> Result<()> {
    let mut args_indices = Vec::with_capacity(qubits.len());
    let mut max_len = 1;
    for q_arg in qubits {
        let indices = resolve_argument(q_arg, &ctx.qregs, &HashMap::new())?;
        max_len = max_len.max(indices.len());
        args_indices.push(indices);
    }
    for indices in &args_indices {
        if indices.len() != 1 && indices.len() != max_len {
            return Err(QRustError::SizeMismatch(
                "register size mismatch in gate call".into(),
            ));
        }
    }
    for i in 0..max_len {
        let current_qubits: Vec<usize> = args_indices
            .iter()
            .map(|indices| if indices.len() == 1 { indices[0] } else { indices[i] })
            .collect();
        emit_gate(circuit, ctx, name, params, &current_qubits, condition.clone())?;
    }
    Ok(())
}

fn emit_gate(
    circuit: &mut Circuit,
    ctx: &ParseContext,
    name: &str,
    params: &[Expr],
    qubits: &[usize],
    condition: Option<ClassicalCondition>,
) -> Result<()> {
    let empty_scope = HashMap::new();
    let mut eval_params = Vec::with_capacity(params.len());
    for p in params {
        eval_params.push(p.evaluate_with_scope(&empty_scope)?);
    }

    let gate_type = name
        .parse::<GateType>()
        .unwrap_or_else(|_| GateType::Custom(name.to_string()));

    let emit = |circuit: &mut Circuit, op: Operation| match &condition {
        Some(cond) => circuit.add_op(Operation::Conditional {
            condition: cond.clone(),
            op: Box::new(op),
        }),
        None => circuit.add_op(op),
    };

    if !matches!(gate_type, GateType::Custom(_)) {
        let final_params = match name {
            "u2" if eval_params.len() == 2 => vec![PI / 2.0, eval_params[0], eval_params[1]],
            _ => eval_params,
        };
        emit(
            circuit,
            Operation::Gate {
                name: gate_type,
                qubits: qubits.to_vec(),
                params: final_params,
            },
        );
        return Ok(());
    }

    if let GateType::Custom(ref n) = gate_type {
        if let Some((def_params, def_qubits, _)) = ctx.gate_defs.get(n) {
            if eval_params.len() != def_params.len() {
                return Err(QRustError::ParseError(format!(
                    "gate '{}' expects {} params, got {}",
                    n,
                    def_params.len(),
                    eval_params.len()
                )));
            }
            if qubits.len() != def_qubits.len() {
                return Err(QRustError::ParseError(format!(
                    "gate '{}' expects {} qubits, got {}",
                    n,
                    def_qubits.len(),
                    qubits.len()
                )));
            }
            emit(
                circuit,
                Operation::Gate {
                    name: gate_type,
                    qubits: qubits.to_vec(),
                    params: eval_params,
                },
            );
            return Ok(());
        }
    }
    Err(QRustError::UnknownGate(name.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_gates() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[3];
            u3(0.1, 0.2, 0.3) q[0];
            id q[0];
            s q[1];
            sdg q[1];
            t q[2];
            tdg q[2];
            ccx q[0], q[1], q[2];
        "#;
        let circuit = parse_qasm(qasm).expect("parse");
        assert_eq!(circuit.operations.len(), 7);
    }

    #[test]
    fn test_header() {
        assert_eq!(openqasm_version("OPENQASM 2.0;"), Ok(("", "2.0".to_string())));
    }

    #[test]
    fn test_valid_file() {
        assert!(parse_qasm("OPENQASM 2.0; qreg q[1];").is_ok());
    }

    #[test]
    fn test_invalid_version() {
        let err = parse_qasm("OPENQASM 3.0; qreg q[1];").unwrap_err();
        assert!(matches!(err, QRustError::Unsupported(_)));
    }

    #[test]
    fn test_missing_header() {
        let err = parse_qasm("qreg q[1];").unwrap_err();
        assert!(matches!(err, QRustError::ParseError(_)));
    }

    #[test]
    fn test_garbage() {
        let err = parse_qasm("NOT A QASM FILE").unwrap_err();
        assert!(matches!(err, QRustError::ParseError(_)));
    }

    #[test]
    fn test_barrier_parses() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[2];
            h q[0];
            barrier q[0], q[1];
            cx q[0], q[1];
        "#;
        let c = parse_qasm(qasm).unwrap();
        assert!(c.operations.iter().any(|op| matches!(op, Operation::Barrier { .. })));
    }

    #[test]
    fn test_reset_parses() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[2];
            reset q[0];
            reset q;
        "#;
        let c = parse_qasm(qasm).unwrap();
        let resets = c
            .operations
            .iter()
            .filter(|op| matches!(op, Operation::Reset { .. }))
            .count();
        // 1 indexed reset + 2 from register-wide.
        assert_eq!(resets, 3);
    }

    #[test]
    fn test_conditional_gate() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[1];
            creg c[1];
            if(c==1) x q[0];
        "#;
        let circ = parse_qasm(qasm).unwrap();
        assert!(circ.operations.iter().any(|op| matches!(op, Operation::Conditional { .. })));
    }

    #[test]
    fn test_conditional_measure() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[1];
            creg c[1];
            if(c==0) measure q[0] -> c[0];
        "#;
        let circ = parse_qasm(qasm).unwrap();
        assert_eq!(circ.operations.len(), 1);
        assert!(matches!(circ.operations[0], Operation::Conditional { .. }));
    }
}