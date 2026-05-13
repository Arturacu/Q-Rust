use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{transpile, TranspilerConfig};
use std::time::Instant;

fn build_massive(num_qubits: usize, depth: usize) -> Circuit {
    let mut c = Circuit::new(num_qubits, 0);
    for _ in 0..depth {
        for i in 0..num_qubits {
            c.add_op(Operation::Gate {
                name: GateType::H,
                qubits: vec![i],
                params: vec![],
            });
        }
        for i in 0..num_qubits - 1 {
            if i % 2 == 0 {
                c.add_op(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![i, i + 1],
                    params: vec![],
                });
            }
        }
    }
    c
}

fn main() {
    println!("--- Massive 400-Qubit Stress Test ---");
    let num_qubits = 400;
    let b = Backend::grid(20, 20);

    let c = build_massive(num_qubits, 5); // depth 5

    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .backend(b.clone())
        .build();

    let start = Instant::now();
    let r = transpile(&c, Some(cfg)).expect("transpile failed");
    let duration = start.elapsed();

    println!("Original ops: {}", c.operations.len());
    println!("Routed ops: {}", r.operations.len());
    println!("Time taken: {:?}", duration);
}
