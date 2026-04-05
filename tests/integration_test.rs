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
        use q_rust::transpiler::decomposition::decompose_basis;
        let qasm = r#"
OPENQASM 2.0;

gate my_rotation(theta) q { U(theta, 0, pi/2) q; }

qreg q[1];
my_rotation(1.57) q[0];
"#;

        let parsed_circuit = parse_qasm(qasm).expect("Failed to parse");
        assert_eq!(parsed_circuit.num_qubits, 1);
        assert_eq!(parsed_circuit.operations.len(), 1);

        // The parser should defer expansion, emitting a Custom gate
        if let Operation::Gate { name, .. } = &parsed_circuit.operations[0] {
            assert_eq!(*name, GateType::Custom("my_rotation".to_string()));
        } else {
            panic!("Expected Custom gate from parser");
        }

        // Transpiler should unroll it
        let circuit = decompose_basis(&parsed_circuit);
        assert_eq!(circuit.operations.len(), 1);

        match &circuit.operations[0] {
            Operation::Gate {
                name,
                qubits,
                params,
            } => {
                assert_eq!(*name, GateType::U);
                assert_eq!(*qubits, vec![0]);
                assert_eq!(params.len(), 3);
                assert!((params[0] - 1.57).abs() < 1e-10);
                assert!((params[1] - 0.0).abs() < 1e-10);
                assert!((params[2] - std::f64::consts::PI / 2.0).abs() < 1e-10);
            }
            _ => panic!("Expected U gate"),
        }
    }

    #[test]
    fn test_decomposition_integration() {
        use q_rust::transpiler::decomposition::decompose_basis;

        // Circuit with various gates: H, CX, X, Y, Z, RX
        let qasm = r#"
            OPENQASM 2.0;
            qreg q[2];
            h q[0];
            x q[1];
            cx q[0], q[1];
            y q[0];
            z q[1];
            rx(1.57) q[0];
        "#;

        let circuit = parse_qasm(qasm).expect("Failed to parse");
        let decomposed = decompose_basis(&circuit);

        // Verify all gates are in basis set (U or CX)
        for op in &decomposed.operations {
            match op {
                Operation::Gate { name, .. } => match name {
                    GateType::U | GateType::CX => {} // OK
                    _ => panic!("Found non-basis gate in CCX decomposition: {:?}", name),
                },
                _ => {} // Ignore non-gate ops
            }
        }

        // Verify count:
        // H -> 1 U
        // X -> 1 U
        // CX -> 1 CX
        // Y -> 1 U
        // Z -> 1 U
        // RX -> 1 U
        // Total: 6 gates
        assert_eq!(decomposed.operations.len(), 6);
    }
}
