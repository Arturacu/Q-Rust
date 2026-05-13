//! Comparative transpilation across workloads.

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
    for target in 0..n {
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![target],
            params: vec![],
        });
        for control in target + 1..n {
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
            matches!(
                op,
                Operation::Gate {
                    name: GateType::CX,
                    ..
                }
            )
        })
        .count()
}

fn main() {
    let quito = std::fs::read_to_string("tests/fixtures/ibm_quito_5q.json").unwrap();
    let backend = Backend::from_config(serde_json::from_str(&quito).unwrap());

    let workloads = vec![
        ("QFT (5Q)", build_qft(5)),
        ("GHZ (5Q)", build_ghz(5)),
        ("Reverse Chain (5Q)", build_reverse_chain(5)),
    ];

    println!("=== Q-Rust transpiler (opt 3) ===");
    for (name, circ) in workloads {
        let cfg = TranspilerConfig::builder()
            .optimization_level(3)
            .decompose_basis(true)
            .backend(backend.clone())
            .build();
        let out = transpile(&circ, Some(cfg)).expect("transpile failed");
        println!("\n--- {name} ---");
        println!(
            "CX: {}   total gates: {}   depth: {}",
            count_cx(&out),
            out.operations.len(),
            out.depth()
        );
    }
}
