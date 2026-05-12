use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{transpile, TranspilerConfig};
use std::time::Instant;

fn build_chain(num_qubits: usize) -> Circuit {
    let mut c = Circuit::new(num_qubits, 0);
    for i in 0..num_qubits - 1 {
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![i, i + 1],
            params: vec![],
        });
    }
    for i in (0..num_qubits - 1).rev() {
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![i, i + 1],
            params: vec![],
        });
    }
    c
}

fn main() {
    println!("--- Linear Chain Routing Benchmark ---");
    let num_qubits = 50;
    let _b = Backend::linear(num_qubits); // Actually we want a bad topology
    let bad_b = Backend::grid(10, 5); // 50 qubits

    let c = build_chain(num_qubits);

    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .backend(bad_b.clone())
        .build();

    let start = Instant::now();
    let r = transpile(&c, Some(cfg));
    let duration = start.elapsed();

    println!("Original CX count: {}", c.operations.len());
    let new_ops = r.operations.len();
    println!("Routed Ops (including SWAPs): {}", new_ops);
    println!(
        "Overhead ratio: {:.2}",
        new_ops as f64 / c.operations.len() as f64
    );
    println!("Time taken: {:?}", duration);
}
