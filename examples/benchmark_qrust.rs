use q_rust::backend::{Backend, BackendConfig};
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{transpile, TranspilerConfig};
use std::time::Instant;

fn load_backend(name: &str) -> Backend {
    let path = format!("tests/fixtures/{name}");
    let json = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {path}"));
    let cfg: BackendConfig = serde_json::from_str(&json).unwrap();
    Backend::from_config(cfg)
}

fn build_circuit(num_qubits: usize, gates: &[(GateType, Vec<usize>, Vec<f64>)]) -> Circuit {
    let mut c = Circuit::new(num_qubits, 0);
    for (name, qubits, params) in gates {
        c.add_op(Operation::Gate {
            name: name.clone(),
            qubits: qubits.clone(),
            params: params.clone(),
        });
    }
    c
}

fn main() {
    println!("--- Basic Execution Benchmark ---");
    let b = load_backend("ibm_heavy_hex_16q.json");

    // Generate 10-qubit chain
    let mut gates = vec![];
    for i in 0..10 {
        gates.push((GateType::H, vec![i], vec![]));
    }
    for i in 0..9 {
        gates.push((GateType::CX, vec![i, i + 1], vec![]));
    }
    for i in 0..10 {
        gates.push((GateType::RZ, vec![i], vec![0.5]));
    }

    let c = build_circuit(10, &gates);

    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .backend(b.clone())
        .build();

    let start = Instant::now();
    let r = transpile(&c, Some(cfg)).expect("transpile failed");
    let duration = start.elapsed();

    println!("Original ops: {}", c.operations.len());
    println!("Routed ops: {}", r.operations.len());
    println!("Time taken: {:?}", duration);
}
