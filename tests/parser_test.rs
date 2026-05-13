//! Parser integration tests.

use q_rust::parser::parse_qasm;

#[test]
fn test_teleportation_circuit() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[3];
        creg c0[1];
        creg c1[1];
        creg c2[1];

        h q[1];
        cx q[1], q[2];

        rx(0.5) q[0];

        cx q[0], q[1];
        h q[0];
        measure q[0] -> c0[0];
        measure q[1] -> c1[0];

        z q[2];
        x q[2];
    "#;
    let circuit = parse_qasm(qasm).expect("parse");
    assert_eq!(circuit.num_qubits, 3);
    assert_eq!(circuit.num_cbits, 3);
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
    let circuit = parse_qasm(qasm).expect("parse");
    assert_eq!(circuit.num_qubits, 2);
    assert_eq!(circuit.operations.len(), 2);
}

#[test]
fn test_undeclared_qubit() {
    let err = parse_qasm("OPENQASM 2.0; qreg q[1]; x q[10];").unwrap_err();
    assert!(format!("{err}").contains("Index out of bounds"));
}

#[test]
fn test_unknown_register() {
    let err = parse_qasm("OPENQASM 2.0; qreg q[1]; x r[0];").unwrap_err();
    assert!(format!("{err}").contains("Undefined"));
}

#[test]
fn test_missing_semicolon() {
    let err = parse_qasm(
        r#"OPENQASM 2.0;
qreg q[1]
x q[0];"#,
    )
    .unwrap_err();
    assert!(format!("{err}").to_lowercase().contains("parse error"));
}

#[test]
fn test_parameterized_gate() {
    let c = parse_qasm("OPENQASM 2.0; qreg q[1]; rx(0.5) q[0];").expect("parse");
    assert_eq!(c.operations.len(), 1);
}

#[test]
fn test_custom_include_error() {
    let err = parse_qasm("OPENQASM 2.0; include \"custom.inc\";").unwrap_err();
    assert!(format!("{err}").contains("Includes are not supported"));
}

#[test]
fn test_qelib1_include_accepted() {
    let qasm = r#"
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[1];
        h q[0];
    "#;
    let c = parse_qasm(qasm).expect("parse");
    assert_eq!(c.num_qubits, 1);
    assert_eq!(c.operations.len(), 1);
}

#[test]
fn test_qasm_round_trip() {
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
    let p1 = parse_qasm(qasm).expect("parse 1");
    let gen = p1.to_qasm(None);
    let p2 = parse_qasm(&gen).expect("parse 2");
    assert_eq!(p1.num_qubits, p2.num_qubits);
    assert_eq!(p1.num_cbits, p2.num_cbits);
    assert_eq!(p1.operations, p2.operations);
}

#[test]
fn test_inverse_qft() {
    let qasm = r#"
        OPENQASM 2.0;

        gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }
        gate u1(lambda) q { U(0,0,lambda) q; }
        gate h a { u2(0,pi) a; }

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
    let c = parse_qasm(qasm).expect("parse");
    assert_eq!(c.num_qubits, 4);
    assert_eq!(c.num_cbits, 4);
    assert_eq!(c.operations.len(), 13);
}
