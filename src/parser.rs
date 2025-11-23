use crate::ir::{Circuit, GateType, Operation};
use nom::{
    IResult,
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{alpha1, alphanumeric1, char, digit1, multispace0, space0, space1},
    combinator::{map, map_res, opt, recognize, value},
    multi::{many0, separated_list0},
    sequence::{delimited, pair, tuple},
};
use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone)]
#[allow(dead_code)]
enum ParsedStatement {
    QReg(String, usize),
    CReg(String, usize),
    Gate(String, Vec<(String, usize)>, Vec<f64>),
    Measure((String, usize), (String, usize)),
    Barrier(Vec<(String, usize)>),
    Reset((String, usize)),
    Ignore,
}

// --- Helpers ---

fn identifier(input: &str) -> IResult<&str, String> {
    map(
        recognize(pair(
            alt((alpha1, tag("_"))),
            many0(alt((alphanumeric1, tag("_")))),
        )),
        |s: &str| s.to_string(),
    )(input)
}

fn usize_parser(input: &str) -> IResult<&str, usize> {
    map_res(digit1, |s: &str| s.parse::<usize>())(input)
}

fn float_parser(input: &str) -> IResult<&str, f64> {
    map_res(
        recognize(tuple((
            opt(tag("-")),
            digit1,
            opt(tuple((tag("."), digit1))),
        ))),
        |s: &str| s.parse::<f64>(),
    )(input)
}

fn comment(input: &str) -> IResult<&str, ()> {
    value((), pair(tag("//"), take_while(|c| c != '\n')))(input)
}

// --- QASM Parsers ---

fn header(input: &str) -> IResult<&str, ParsedStatement> {
    value(
        ParsedStatement::Ignore,
        tuple((tag("OPENQASM"), space1, tag("2.0"), tag(";"))),
    )(input)
}

fn include(input: &str) -> IResult<&str, ParsedStatement> {
    value(
        ParsedStatement::Ignore,
        tuple((
            tag("include"),
            space1,
            delimited(char('"'), take_while1(|c| c != '"'), char('"')),
            tag(";"),
        )),
    )(input)
}

fn qreg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("qreg"),
            space1,
            identifier,
            delimited(char('['), usize_parser, char(']')),
            tag(";"),
        )),
        |(_, _, name, size, _)| ParsedStatement::QReg(name, size),
    )(input)
}

fn creg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("creg"),
            space1,
            identifier,
            delimited(char('['), usize_parser, char(']')),
            tag(";"),
        )),
        |(_, _, name, size, _)| ParsedStatement::CReg(name, size),
    )(input)
}

fn qubit_ref(input: &str) -> IResult<&str, (String, usize)> {
    pair(identifier, delimited(char('['), usize_parser, char(']')))(input)
}

fn gate_call(input: &str) -> IResult<&str, ParsedStatement> {
    let (input, name) = identifier(input)?;
    let (input, params) = opt(delimited(
        char('('),
        separated_list0(tuple((space0, char(','), space0)), float_parser),
        char(')'),
    ))(input)?;
    let (input, _) = space1(input)?;
    let (input, qubits) = separated_list0(tuple((space0, char(','), space0)), qubit_ref)(input)?;
    let (input, _) = tag(";")(input)?;

    Ok((
        input,
        ParsedStatement::Gate(name, qubits, params.unwrap_or_default()),
    ))
}

fn measure(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("measure"),
            space1,
            qubit_ref,
            space0,
            tag("->"),
            space0,
            qubit_ref,
            tag(";"),
        )),
        |(_, _, q, _, _, _, c, _)| ParsedStatement::Measure(q, c),
    )(input)
}

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

        let (rem, stmt) = alt((header, include, qreg, creg, measure, gate_call))(current_input)
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
    fn test_header() {
        assert_eq!(header("OPENQASM 2.0;"), Ok(("", ParsedStatement::Ignore)));
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
}
