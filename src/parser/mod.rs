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

pub fn parse_qasm(input: &str) -> Result<Circuit, String> {
    let mut circuit = Circuit::new(0, 0);
    let mut qregs: HashMap<String, (usize, usize)> = HashMap::new(); // name -> (start_index, size)
    let mut cregs: HashMap<String, (usize, usize)> = HashMap::new(); // name -> (start_index, size)
    let mut total_qubits = 0;
    let mut total_cbits = 0;

    let mut current_input = input;

    // 1. Skip initial comments/whitespace and parse Header
    loop {
        current_input = multispace0::<&str, nom::error::Error<&str>>(current_input)
            .map_err(|e| e.to_string())?
            .0;
        if current_input.is_empty() {
            return Err("Empty file or missing OPENQASM header".to_string());
        }

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
        current_input = multispace0::<&str, nom::error::Error<&str>>(current_input)
            .map_err(|e| e.to_string())?
            .0;
        if current_input.is_empty() {
            break;
        }

        if let Ok((rem, _)) = comment(current_input) {
            current_input = rem;
            continue;
        }

        let (rem, stmt) = alt((include, qreg, creg, measure, gate_call))(current_input)
            .map_err(|_e| format!("Parse error at: {}", current_input))?;

        current_input = rem;

        match stmt {
            ParsedStatement::Ignore => {}
            ParsedStatement::QReg(name, size) => {
                qregs.insert(name, (total_qubits, size));
                total_qubits += size;
            }
            ParsedStatement::CReg(name, size) => {
                cregs.insert(name, (total_cbits, size));
                total_cbits += size;
            }
            ParsedStatement::Gate(name, qubits, params) => {
                let mut resolved_qubits = Vec::new();
                for (q_name, q_idx) in qubits {
                    if let Some(&(start, size)) = qregs.get(&q_name) {
                        if q_idx < size {
                            resolved_qubits.push(start + q_idx);
                        } else {
                            return Err(format!(
                                "Qubit index out of bounds: {}[{}]",
                                q_name, q_idx
                            ));
                        }
                    } else {
                        return Err(format!("Undefined quantum register: {}", q_name));
                    }
                }

                let gate_type = map_gate_type(&name, &params);
                circuit.add_op(Operation::Gate {
                    name: gate_type,
                    qubits: resolved_qubits,
                    params,
                });
            }
            ParsedStatement::Measure((q_name, q_idx), (c_name, c_idx)) => {
                let resolved_q = if let Some(&(start, size)) = qregs.get(&q_name) {
                    if q_idx < size {
                        start + q_idx
                    } else {
                        return Err(format!("Qubit index out of bounds: {}[{}]", q_name, q_idx));
                    }
                } else {
                    return Err(format!("Undefined quantum register: {}", q_name));
                };

                let resolved_c = if let Some(&(start, size)) = cregs.get(&c_name) {
                    if c_idx < size {
                        start + c_idx
                    } else {
                        return Err(format!(
                            "Classical bit index out of bounds: {}[{}]",
                            c_name, c_idx
                        ));
                    }
                } else {
                    return Err(format!("Undefined classical register: {}", c_name));
                };

                circuit.add_op(Operation::Measure {
                    qubit: resolved_q,
                    cbit: resolved_c,
                });
            }
            _ => {}
        }
    }

    circuit.num_qubits = total_qubits;
    circuit.num_cbits = total_cbits;
    Ok(circuit)
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
                ParsedStatement::Gate("h".to_string(), vec![("q".to_string(), 0)], vec![])
            ))
        );
    }

    #[test]
    fn test_gate_call_with_params() {
        assert_eq!(
            gate_call("rx(1.57) q[0];"),
            Ok((
                "",
                ParsedStatement::Gate("rx".to_string(), vec![("q".to_string(), 0)], vec![1.57])
            ))
        );
    }

    #[test]
    fn test_measure() {
        assert_eq!(
            measure("measure q[0] -> c[0];"),
            Ok((
                "",
                ParsedStatement::Measure(("q".to_string(), 0), ("c".to_string(), 0))
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
    fn test_garbage() {
        let qasm = "NOT A QASM FILE";
        let err = parse_qasm(qasm).unwrap_err();
        assert!(err.contains("Missing or invalid OPENQASM header"));
    }
}
