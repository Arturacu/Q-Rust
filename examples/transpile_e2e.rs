use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() {
    // 1. Define a source QASM circuit (GHZ-3 State with some redundant gates for optimization)
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[3];
        
        // Redundant H gates that should cancel
        h q[0];
        h q[0];
        
        // GHZ State preparation
        h q[0];
        cx q[0], q[1];
        cx q[1], q[2];
        
        // New: Controlled gate optimizations
        // Two CZ gates should cancel
        cz q[0], q[1];
        cz q[0], q[1];
        
        // Two CRZ gates should merge
        crz(0.1) q[1], q[2];
        crz(0.2) q[1], q[2];
        
        // Commuting RZ gates that should merge
        rz(0.1) q[0];
        rz(0.2) q[0];
        
        // Final measurement
        creg c[3];
        measure q -> c;
    "#;

    println!("--- Q-Rust E2E Transpiler Demo ---");

    // 2. Parse the QASM into an IR Circuit
    let circuit = parse_qasm(qasm).expect("Failed to parse QASM");
    let original_count = circuit.operations.len();
    println!("Original gate count: {}", original_count);

    // 3. Configure the transpiler for Level 3 Optimization (Max)
    // - optimization_level: 3 (Enables Fusion, Commutation, Cancellation)
    // - decompose_basis: true (Unrolls to U, CX)
    // - backend: None (Defaults to All-to-All coupling)
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: true,
        backend: None,
    };

    // 4. Run the Transpilation Pipeline
    println!("Running transpiler with Level 3 Optimization...");
    let transpiled_circuit = transpile(&circuit, Some(config));
    let final_count = transpiled_circuit.operations.len();

    println!("Transpiled gate count: {}", final_count);
    println!(
        "Optimization reduction: {:.1}%",
        (1.0 - (final_count as f64 / original_count as f64)) * 100.0
    );

    // 5. Verify physical correctness via Unitary Simulation
    let u_orig = circuit_to_unitary(&circuit);
    let u_trans = circuit_to_unitary(&transpiled_circuit);
    let fidelity = unitary_fidelity(&u_orig, &u_trans);

    println!("Unitary Fidelity: {:.10}", fidelity);
    if fidelity > 0.999999 {
        println!("✅ Transpilation is mathematically equivalent (Equivalence up to Global Phase).");
    } else {
        println!("❌ Transpilation fidelity error!");
    }

    // 6. Inspect the resulting gates
    println!("\nFinal Optimized Gate Sequence:");
    for (i, op) in transpiled_circuit.operations.iter().enumerate() {
        println!("  {:02}: {:?}", i, op);
    }
}
