//! End-to-end transpilation demo.

use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[3];

        // Redundant H pair
        h q[0];
        h q[0];

        // GHZ preparation
        h q[0];
        cx q[0], q[1];
        cx q[1], q[2];

        // Inverse-pair controlled gates
        cz q[0], q[1];
        cz q[0], q[1];

        // Rotation mergeable pair
        crz(0.1) q[1], q[2];
        crz(0.2) q[1], q[2];

        // Commuting RZs
        rz(0.1) q[0];
        rz(0.2) q[0];

        // Final measurement
        creg c[3];
        measure q -> c;
    "#;

    println!("--- Q-Rust E2E Transpiler Demo ---");
    let circuit = parse_qasm(qasm).expect("parse failure");
    let original = circuit.operations.len();
    println!("Original gate count: {original}");

    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .build();

    let transpiled = transpile(&circuit, Some(config)).expect("transpile failed");
    let final_count = transpiled.operations.len();
    println!("Transpiled gate count: {final_count}");
    println!(
        "Reduction: {:.1}%",
        (1.0 - (final_count as f64 / original as f64)) * 100.0
    );

    let u_orig = circuit_to_unitary(&circuit);
    let u_trans = circuit_to_unitary(&transpiled);
    let fid = unitary_fidelity(&u_orig, &u_trans);
    println!("Unitary fidelity: {fid:.10}");
    if fid > 0.999_999 {
        println!("✅ Transpilation is mathematically equivalent.");
    } else {
        println!("❌ Fidelity error!");
    }

    println!("\nFinal gate sequence:");
    for (i, op) in transpiled.operations.iter().enumerate() {
        println!("  {i:02}: {op:?}");
    }
}
