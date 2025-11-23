#[cfg(test)]
mod tests {
    use q_rust::ir::{GateType, Operation};
    use q_rust::parser::parse_qasm;

    #[test]
    fn test_simple_circuit() {
        let qasm = r#"
            OPENQASM 2.0;
            include "qelib1.inc";
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
}
