#[cfg(test)]
mod tests {
    use q_rust::ir::{GateType, Operation};
    use q_rust::parser::parse_qasm;

    #[test]
    fn test_simple_circuit() {
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[2];
            creg c[2];
            h q[0];
            cx q[0], q[1];
            measure q[0] -> c[0];
            measure q[1] -> c[1];
        "#;

        let circuit = parse_qasm(qasm).expect("Failed to parse QASM");

        assert_eq!(circuit.num_qubits, 2);
        assert_eq!(circuit.num_cbits, 2);
        assert_eq!(circuit.operations.len(), 4);

        match &circuit.operations[0] {
            Operation::Gate { name, qubits, .. } => {
                assert_eq!(*name, GateType::H);
                assert_eq!(*qubits, vec![0]);
            }
            _ => panic!("Expected H gate"),
        }

        match &circuit.operations[1] {
            Operation::Gate { name, qubits, .. } => {
                assert_eq!(*name, GateType::CX);
                assert_eq!(*qubits, vec![0, 1]);
            }
            _ => panic!("Expected CX gate"),
        }
    }

    #[test]
    fn test_custom_gate_expansion() {
        let qasm = r#"
OPENQASM 2.0;

gate my_rotation(theta) q { U(theta, 0, pi/2) q; }

qreg q[1];
my_rotation(1.57) q[0];
"#;

        let circuit = parse_qasm(qasm).expect("Failed to parse");
        assert_eq!(circuit.num_qubits, 1);
        assert_eq!(circuit.operations.len(), 1);

        match &circuit.operations[0] {
            Operation::Gate { name, qubits, params } => {
                assert!(matches!(name, GateType::U(_, _, _)));
                assert_eq!(*qubits, vec![0]);
                assert_eq!(params.len(), 3);
                assert!((params[0] - 1.57).abs() < 1e-10);
                assert!((params[1] - 0.0).abs() < 1e-10);
                assert!((params[2] - std::f64::consts::PI / 2.0).abs() < 1e-10);
            }
            _ => panic!("Expected U gate"),
        }
    }
}
