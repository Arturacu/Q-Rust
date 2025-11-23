pub mod ast;
pub mod rules;

use self::ast::ParsedStatement;
use self::rules::{comment, creg, gate_call, include, measure, openqasm_version, qreg};
use crate::ir::{Circuit, GateType, Operation};
use nom::{branch::alt, character::complete::multispace0};
use std::collections::HashMap;

// --- Resolution & Mapping ---

fn map_gate_type(name: &str, params: &[f64]) -> GateType {
    match name {
        "h" => GateType::H,
        "x" => GateType::X,
        "y" => GateType::Y,
        "z" => GateType::Z,
        "cx" => GateType::CX,
        "rx" => GateType::RX(params.first().cloned().unwrap_or(0.0)),
        "ry" => GateType::RY(params.first().cloned().unwrap_or(0.0)),
        "rz" => GateType::RZ(params.first().cloned().unwrap_or(0.0)),
        "u1" => GateType::RZ(params.first().cloned().unwrap_or(0.0)), // u1(lambda) = RZ(lambda)
        "u2" => GateType::U(
            std::f64::consts::FRAC_PI_2,
            params.first().cloned().unwrap_or(0.0),
            params.get(1).cloned().unwrap_or(0.0),
        ), // u2(phi, lambda) = U(pi/2, phi, lambda)
        "u3" | "U" => GateType::U(
            params.first().cloned().unwrap_or(0.0),
            params.get(1).cloned().unwrap_or(0.0),
            params.get(2).cloned().unwrap_or(0.0),
        ),
        "id" => GateType::ID,
        "s" => GateType::S,
        "sdg" => GateType::Sdg,
        "t" => GateType::T,
        "tdg" => GateType::Tdg,
        "swap" => GateType::SWAP,
        "ccx" => GateType::CCX,
        _ => GateType::Custom(name.to_string()),
    }
}

use self::ast::Expr;

fn evaluate_expr(expr: &Expr, params: &HashMap<String, f64>) -> Result<f64, String> {
    match expr {
        Expr::Float(val) => Ok(*val),
        Expr::Var(name) => {
            if name == "pi" {
                Ok(std::f64::consts::PI)
            } else {
                params
                    .get(name)
                    .cloned()
                    .ok_or_else(|| format!("Undefined parameter: {}", name))
            }
        }
        Expr::Add(lhs, rhs) => Ok(evaluate_expr(lhs, params)? + evaluate_expr(rhs, params)?),
        Expr::Sub(lhs, rhs) => Ok(evaluate_expr(lhs, params)? - evaluate_expr(rhs, params)?),
        Expr::Mul(lhs, rhs) => Ok(evaluate_expr(lhs, params)? * evaluate_expr(rhs, params)?),
        Expr::Div(lhs, rhs) => Ok(evaluate_expr(lhs, params)? / evaluate_expr(rhs, params)?),
    }
}

fn resolve_argument(
    arg: &(String, Option<usize>),
    qregs: &HashMap<String, (usize, usize)>,
    scope_qubits: &HashMap<String, usize>, // For inside gate defs (arg name -> resolved index)
) -> Result<Vec<usize>, String> {
    let (name, idx) = arg;
    if let Some(mapped_idx) = scope_qubits.get(name) {
        // It's a qubit argument in a gate definition
        if idx.is_some() {
            return Err(format!("Cannot index a qubit argument: {}", name));
        }
        Ok(vec![*mapped_idx])
    } else if let Some(&(start, size)) = qregs.get(name) {
        // It's a global register
        if let Some(i) = idx {
            if *i < size {
                Ok(vec![start + i])
            } else {
                Err(format!("Qubit index out of bounds: {}[{}]", name, i))
            }
        } else {
            // Broadcasting: return all qubits in register
            Ok((0..size).map(|i| start + i).collect())
        }
    } else {
        Err(format!("Undefined quantum register or argument: {}", name))
    }
}

struct ParseContext {
    qregs: HashMap<String, (usize, usize)>,
    cregs: HashMap<String, (usize, usize)>,
    gate_defs: HashMap<String, (Vec<String>, Vec<String>, Vec<ParsedStatement>)>,
}

pub fn parse_qasm(input: &str) -> Result<Circuit, String> {
    let mut circuit = Circuit::new(0, 0);
    let mut ctx = ParseContext {
        qregs: HashMap::new(),
        cregs: HashMap::new(),
        gate_defs: HashMap::new(),
    };
    let mut total_qubits = 0;
    let mut total_cbits = 0;

    let mut current_input = input;

    // 1. Skip initial comments/whitespace and parse Header
    loop {
        // Consume whitespace
        let (rem, _) = multispace0::<&str, nom::error::Error<&str>>(current_input)
            .map_err(|e| e.to_string())?;
        current_input = rem;

        if current_input.is_empty() {
            return Err("Empty file or missing OPENQASM header".to_string());
        }

        // Check for comment
        if let Ok((rem, _)) = comment(current_input) {
            current_input = rem;
            continue;
        }
        break;
    }

    let (rem, version) = openqasm_version(current_input).map_err(|_| {
        "Missing or invalid OPENQASM header. File must start with 'OPENQASM 2.0;'".to_string()
    })?;

    if version != "2.0" {
        return Err(format!(
            "Unsupported OpenQASM version: '{}'. Only '2.0' is supported.",
            version
        ));
    }
    current_input = rem;

    // 2. Parse remaining statements
    loop {
        // Consume whitespace
        let (rem, _) = multispace0::<&str, nom::error::Error<&str>>(current_input)
            .map_err(|e| e.to_string())?;
        current_input = rem;

        if current_input.is_empty() {
            break;
        }

        // Check for comment
        if let Ok((rem, _)) = comment(current_input) {
            current_input = rem;
            continue;
        }

        let (rem, stmt) = alt((
            include,
            qreg,
            creg,
            measure,
            rules::barrier,
            rules::gate_def,
            rules::if_stmt,
            gate_call,
        ))(current_input)
        .map_err(|_e| format!("Parse error at: {}", current_input))?;

        current_input = rem;

        match stmt {
            ParsedStatement::Ignore => {}
            ParsedStatement::Include(filename) => {
                return Err(format!(
                    "Includes are not supported. Please resolve all imports before parsing. Found: 'include \"{}\"'",
                    filename
                ));
            }
            ParsedStatement::QReg(name, size) => {
                ctx.qregs.insert(name, (total_qubits, size));
                total_qubits += size;
            }
            ParsedStatement::CReg(name, size) => {
                ctx.cregs.insert(name, (total_cbits, size));
                total_cbits += size;
            }
            ParsedStatement::GateDef(name, params, qubits, body) => {
                ctx.gate_defs.insert(name, (params, qubits, body));
            }
            ParsedStatement::Gate(name, qubits, params) => {
                // Top-level gate call
                let mut args_indices = Vec::new();
                let mut max_len = 1;

                for q_arg in &qubits {
                    let indices = resolve_argument(q_arg, &ctx.qregs, &HashMap::new())?;
                    if indices.len() > max_len {
                        max_len = indices.len();
                    }
                    args_indices.push(indices);
                }

                // Validate broadcasting
                for indices in &args_indices {
                    if indices.len() != 1 && indices.len() != max_len {
                        return Err("Register size mismatch in gate call".to_string());
                    }
                }

                // Iterate max_len times
                for i in 0..max_len {
                    let mut current_qubits = Vec::new();
                    for indices in &args_indices {
                        if indices.len() == 1 {
                            current_qubits.push(indices[0]);
                        } else {
                            current_qubits.push(indices[i]);
                        }
                    }

                    expand_gate(
                        &mut circuit,
                        &ctx,
                        &name,
                        &params,
                        &current_qubits,
                        &HashMap::new(),
                        &HashMap::new(),
                    )?;
                }
            }
            ParsedStatement::Measure((q_name, q_idx), (c_name, c_idx)) => {
                let q_indices = resolve_argument(&(q_name, q_idx), &ctx.qregs, &HashMap::new())?;
                let c_indices = if let Some(&(start, size)) = ctx.cregs.get(&c_name) {
                    if let Some(i) = c_idx {
                        if i < size {
                            vec![start + i]
                        } else {
                            return Err(format!("Index out of bounds"));
                        }
                    } else {
                        (0..size).map(|i| start + i).collect()
                    }
                } else {
                    return Err(format!("Undefined creg"));
                };

                if q_indices.len() != c_indices.len() {
                    return Err("Measure register size mismatch".to_string());
                }

                for (q, c) in q_indices.iter().zip(c_indices.iter()) {
                    circuit.add_op(Operation::Measure {
                        qubit: *q,
                        cbit: *c,
                    });
                }
            }
            ParsedStatement::Barrier(_) => {} // Ignore top level barrier
            ParsedStatement::If(_, _, _) => {
                return Err("Conditional operations not yet supported".to_string());
            }
        }
    }

    circuit.num_qubits = total_qubits;
    circuit.num_cbits = total_cbits;
    Ok(circuit)
}

fn expand_gate(
    circuit: &mut Circuit,
    ctx: &ParseContext,
    name: &str,
    params: &[Expr],
    qubits: &[usize],
    scope_params: &HashMap<String, f64>,
    scope_qubits: &HashMap<String, usize>,
) -> Result<(), String> {
    // 1. Evaluate parameters
    let mut eval_params = Vec::new();
    for p in params {
        eval_params.push(evaluate_expr(p, scope_params)?);
    }

    // 2. Check standard gate
    let gate_type = map_gate_type(name, &eval_params);
    if !matches!(gate_type, GateType::Custom(_)) {
        circuit.add_op(Operation::Gate {
            name: gate_type,
            qubits: qubits.to_vec(),
            params: eval_params,
        });
        return Ok(());
    }

    // 3. Custom gate expansion
    if let GateType::Custom(ref n) = gate_type {
        if let Some((def_params, def_qubits, body)) = ctx.gate_defs.get(n) {
            if eval_params.len() != def_params.len() {
                return Err(format!(
                    "Gate {} expects {} params, got {}",
                    n,
                    def_params.len(),
                    eval_params.len()
                ));
            }
            if qubits.len() != def_qubits.len() {
                return Err(format!(
                    "Gate {} expects {} qubits, got {}",
                    n,
                    def_qubits.len(),
                    qubits.len()
                ));
            }

            let mut new_scope_params = scope_params.clone();
            for (j, param_name) in def_params.iter().enumerate() {
                new_scope_params.insert(param_name.clone(), eval_params[j]);
            }
            let mut new_scope_qubits = scope_qubits.clone();
            for (j, qubit_name) in def_qubits.iter().enumerate() {
                new_scope_qubits.insert(qubit_name.clone(), qubits[j]);
            }

            for stmt in body {
                match stmt {
                    ParsedStatement::Gate(g_name, g_qubits, g_params) => {
                        // Resolve qubits in body
                        let mut g_resolved_qubits = Vec::new();
                        for q_arg in g_qubits {
                            let indices = resolve_argument(q_arg, &ctx.qregs, &new_scope_qubits)?;
                            if indices.len() != 1 {
                                return Err(format!("Broadcasting not allowed in gate definition body for argument {:?}", q_arg));
                            }
                            g_resolved_qubits.push(indices[0]);
                        }
                        expand_gate(
                            circuit,
                            ctx,
                            g_name,
                            g_params,
                            &g_resolved_qubits,
                            &new_scope_params,
                            &new_scope_qubits,
                        )?;
                    }
                    ParsedStatement::Barrier(_) => {}
                    _ => {}
                }
            }
            return Ok(());
        }
    }

    Err(format!("Unknown gate: {}", name))
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
        let circuit = parse_qasm(qasm).expect("Failed to parse standard gates");
        assert_eq!(circuit.operations.len(), 7);
        // Check first gate is U
        if let Operation::Gate { name, .. } = &circuit.operations[0] {
            assert!(matches!(name, GateType::U(..)));
        } else {
            panic!("Expected U gate");
        }
    }

    #[test]
    fn test_header() {
        assert_eq!(
            openqasm_version("OPENQASM 2.0;"),
            Ok(("", "2.0".to_string()))
        );
    }

    #[test]
    fn test_qreg() {
        assert_eq!(
            qreg("qreg q[2];"),
            Ok(("", ParsedStatement::QReg("q".to_string(), 2)))
        );
    }

    #[test]
    fn test_creg() {
        assert_eq!(
            creg("creg c[3];"),
            Ok(("", ParsedStatement::CReg("c".to_string(), 3)))
        );
    }

    #[test]
    fn test_gate_call_no_params() {
        assert_eq!(
            gate_call("h q[0];"),
            Ok((
                "",
                ParsedStatement::Gate("h".to_string(), vec![("q".to_string(), Some(0))], vec![])
            ))
        );
    }

    #[test]
    fn test_gate_call_with_params() {
        assert_eq!(
            gate_call("rx(1.57) q[0];"),
            Ok((
                "",
                ParsedStatement::Gate(
                    "rx".to_string(),
                    vec![("q".to_string(), Some(0))],
                    vec![ast::Expr::Float(1.57)]
                )
            ))
        );
    }

    #[test]
    fn test_measure() {
        assert_eq!(
            measure("measure q[0] -> c[0];"),
            Ok((
                "",
                ParsedStatement::Measure(("q".to_string(), Some(0)), ("c".to_string(), Some(0)))
            ))
        );
    }

    #[test]
    fn test_valid_file() {
        let qasm = "OPENQASM 2.0; qreg q[1];";
        assert!(parse_qasm(qasm).is_ok());
    }

    #[test]
    fn test_invalid_version() {
        let qasm = "OPENQASM 3.0; qreg q[1];";
        let err = parse_qasm(qasm).unwrap_err();
        assert!(err.contains("Unsupported OpenQASM version"));
    }

    #[test]
    fn test_missing_header() {
        let qasm = "qreg q[1];";
        let err = parse_qasm(qasm).unwrap_err();
        assert!(err.contains("Missing or invalid OPENQASM header"));
    }

    #[test]
    fn test_gate_def_simple() {
        let input = "gate mygate q { h q; }";
        let result = rules::gate_def(input);
        assert!(
            result.is_ok(),
            "Failed to parse simple gate def: {:?}",
            result
        );
    }

    #[test]
    fn test_gate_def_with_params() {
        let input = "gate u2(phi,lambda) q { h q; }";
        let result = rules::gate_def(input);
        assert!(
            result.is_ok(),
            "Failed to parse gate def with params: {:?}",
            result
        );
    }

    #[test]
    fn test_gate_def_with_uppercase_call() {
        let input = "gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }";
        let result = rules::gate_def(input);
        assert!(
            result.is_ok(),
            "Failed to parse gate def with uppercase call: {:?}",
            result
        );
    }

    #[test]
    fn test_garbage() {
        let qasm = "NOT A QASM FILE";
        let err = parse_qasm(qasm).unwrap_err();
        assert!(err.contains("Missing or invalid OPENQASM header"));
    }
}
