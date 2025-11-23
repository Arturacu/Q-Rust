use q_rust::parser::parse_qasm;

#[test]
fn test_teleportation_circuit() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[3];
        creg c0[1];
        creg c1[1];
        creg c2[1];

        // Bell pair
        h q[1];
        cx q[1], q[2];

        // Prepare payload
        rx(0.5) q[0];

        // Teleportation
        cx q[0], q[1];
        h q[0];
        measure q[0] -> c0[0];
        measure q[1] -> c1[0];

        // Correction (conditional logic not yet supported, but gates are)
        z q[2]; 
        x q[2];
    "#;

    let circuit = parse_qasm(qasm).expect("Failed to parse teleportation circuit");
    assert_eq!(circuit.num_qubits, 3);
    assert_eq!(circuit.num_cbits, 3);

    // Verify operation count (H, CX, RX, CX, H, M, M, Z, X) = 9 operations
    assert_eq!(circuit.operations.len(), 9);
}

#[test]
fn test_whitespace_tolerance() {
    let qasm = r#"
        OPENQASM 2.0;

        qreg    q[2]   ;
          creg  c[2];

        h   q[0]  ; // Comment
        
        cx q[0] ,  q[1];
    "#;
    let circuit = parse_qasm(qasm).expect("Failed to parse whitespace");
    assert_eq!(circuit.num_qubits, 2);
    assert_eq!(circuit.operations.len(), 2);
}

#[test]
fn test_undeclared_qubit() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[1];
        x q[10]; // Index out of bounds
    "#;
    let err = parse_qasm(qasm).unwrap_err();
    assert!(err.contains("Qubit index out of bounds"));
}

#[test]
fn test_unknown_register() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[1];
        x r[0]; // Unknown register
    "#;
    let err = parse_qasm(qasm).unwrap_err();
    assert!(err.contains("Undefined quantum register"));
}

#[test]
fn test_missing_semicolon() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[1]
        x q[0];
    "#;
    let err = parse_qasm(qasm).unwrap_err();
    assert!(err.contains("Parse error"));
}

#[test]
fn test_parameterized_gate() {
    let qasm = "OPENQASM 2.0; qreg q[1]; rx(0.5) q[0];";
    let circuit = parse_qasm(qasm).expect("Failed to parse parameterized gate");
    assert_eq!(circuit.operations.len(), 1);
}

#[test]
fn test_custom_include_error() {
    let qasm = r#"
        OPENQASM 2.0;
        include "custom.inc";
    "#;
    let err = parse_qasm(qasm).unwrap_err();
    assert!(err.contains("Includes are not supported"));
}

#[test]
fn test_qelib1_include_rejected() {
    let qasm = r#"
        OPENQASM 2.0;
        include "qelib1.inc";
    "#;
    let err = parse_qasm(qasm).unwrap_err();
    assert!(err.contains("Includes are not supported"));
}

#[test]
fn test_inverse_qft() {
    let qasm = r#"
        OPENQASM 2.0;

        // 2-parameter single-qubit gate
        gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }

        // 1-parameter single-qubit gate
        gate u1(lambda) q { U(0,0,lambda) q; }

        // Hadamard
        gate h a { u2(0,pi) a; }

        // --- End minimal qelib1.inc subset ---

        // QFT and measure
        qreg q[4];
        creg c[4];

        h q;
        barrier q;

        h q[0];
        measure q[0] -> c[0];

        h q[1];
        measure q[1] -> c[1];

        h q[2];
        measure q[2] -> c[2];

        h q[3];
        measure q[3] -> c[3];
    "#;

    let circuit = parse_qasm(qasm).expect("Failed to parse inverse QFT");
    assert_eq!(circuit.num_qubits, 4);
    assert_eq!(circuit.num_cbits, 4);
    assert_eq!(circuit.operations.len(), 12);
}
