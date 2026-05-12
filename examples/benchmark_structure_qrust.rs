use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{transpile, TranspilerConfig};
use std::time::Instant;

fn main() {
    println!("--- Structured Algorithm Benchmark ---");
    let b = Backend::linear(10);
    let mut c = Circuit::new(10, 0);

    // Build an inverse QFT-like structure
    for i in 0..10 {
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![i],
            params: vec![],
        });
        for j in i + 1..10 {
            c.add_op(Operation::Gate {
                name: GateType::CRZ,
                qubits: vec![j, i],
                params: vec![std::f64::consts::PI / (2_i32.pow((j - i) as u32) as f64)],
            });
        }
    }

    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .backend(b)
        .build();

    let start = Instant::now();
    let r = transpile(&c, Some(cfg));
    let duration = start.elapsed();

    println!("Original QFT-like ops: {}", c.operations.len());
    println!("Transpiled into standard basis: {}", r.operations.len());
    println!("Time taken: {:?}", duration);
}
