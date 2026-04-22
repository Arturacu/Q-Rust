use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn build_ghz(n: usize) -> Circuit {
    let mut c = Circuit::new(n, 0);
    c.add_op(Operation::Gate {
        name: GateType::H,
        qubits: vec![0],
        params: vec![],
    });
    for i in 0..n - 1 {
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![i, i + 1],
            params: vec![],
        });
    }
    c
}

fn build_reverse_chain(n: usize) -> Circuit {
    let mut c = Circuit::new(n, 0);
    for i in (0..n - 1).rev() {
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![i + 1, i],
            params: vec![],
        });
    }
    c
}

fn build_qft(n: usize) -> Circuit {
    let mut c = Circuit::new(n, 0);
    let mut angle = std::f64::consts::PI;
    for target in 0..n {
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![target],
            params: vec![],
        });
        for control in target + 1..n {
            angle /= 2.0;
            // Standard QFT uses controlled-phase, but we'll approximate with generic 2Q gates for routing stress test
            c.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![control, target],
                params: vec![],
            });
        }
    }
    c
}

fn count_cx(c: &Circuit) -> usize {
    c.operations
        .iter()
        .filter(|op| {
            if let Operation::Gate { name, .. } = op {
                *name == GateType::CX
            } else {
                false
            }
        })
        .count()
}

fn main() {
    // Topologies
    let quito_json = std::fs::read_to_string("tests/fixtures/ibm_quito_5q.json").unwrap();
    let backend = Backend::from_config(serde_json::from_str(&quito_json).unwrap());

    let workloads = vec![
        ("QFT (5Q)", build_qft(5)),
        ("GHZ (5Q)", build_ghz(5)),
        ("Reverse Chain (5Q)", build_reverse_chain(5)),
    ];

    println!("=== Q-RUST TRANSPILER (Optimization Level 3) ===");

    for (name, circ) in workloads {
        let config = TranspilerConfig {
            optimization_level: 3,
            decompose_basis: true,
            backend: Some(backend.clone()),
        };
        let transpiled = transpile(&circ, Some(config));
        let cx_count = count_cx(&transpiled);

        println!("\n--- {} ---", name);
        println!("CX count:    {}", cx_count);
        // Note: Q-Rust depth is currently evaluated loosely, so we track raw gate metrics
        println!("Total gates: {}", transpiled.operations.len());
    }
}
