//! BeamSABRE routing demo.

use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn count_swaps(circuit: &Circuit) -> usize {
    circuit
        .operations
        .iter()
        .filter(|op| {
            matches!(
                op,
                Operation::Gate {
                    name: GateType::SWAP,
                    ..
                }
            )
        })
        .count()
}

fn check_adjacency(circuit: &Circuit, backend: &Backend) {
    for op in &circuit.operations {
        if let Operation::Gate { name, qubits, .. } = op {
            if qubits.len() == 2 {
                assert!(
                    backend.is_adjacent(qubits[0], qubits[1]),
                    "❌ Non-adjacent: {:?} on {:?}",
                    name,
                    qubits
                );
            }
        }
    }
    println!("  ✅ all 2-qubit gates on adjacent physical qubits");
}

fn main() {
    println!("=== BeamSABRE Routing Demo ===\n");

    println!("--- Test 1: GHZ-3 on Linear(3) ---");
    let ghz_src = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[1];
        cx q[1], q[2];
    "#;
    let ghz = parse_qasm(ghz_src).expect("parse");
    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .backend(Backend::linear(3))
        .build();
    let out = transpile(&ghz, Some(cfg)).expect("transpile failed");
    println!(
        "  gates: {}  swaps: {}",
        out.operations.len(),
        count_swaps(&out)
    );
    check_adjacency(&out, &Backend::linear(3));
    let fid = unitary_fidelity(&circuit_to_unitary(&ghz), &circuit_to_unitary(&out));
    println!("  fidelity (no permutation assumed): {fid:.10}");

    println!("\n--- Test 2: non-adjacent CX on Linear(3) ---");
    let nonadj_src = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[2];
    "#;
    let nonadj = parse_qasm(nonadj_src).expect("parse");
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .backend(Backend::linear(3))
        .build();
    let out = transpile(&nonadj, Some(cfg)).expect("transpile failed");
    println!(
        "  gates: {}  swaps: {}",
        out.operations.len(),
        count_swaps(&out)
    );
    check_adjacency(&out, &Backend::linear(3));

    println!("\n--- Test 3: star pattern on Grid(2,2) ---");
    let grid_src = r#"
        OPENQASM 2.0;
        qreg q[4];
        h q[0];
        cx q[0], q[1];
        cx q[0], q[2];
        cx q[0], q[3];
    "#;
    let grid = parse_qasm(grid_src).expect("parse");
    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(false)
        .backend(Backend::grid(2, 2))
        .build();
    let out = transpile(&grid, Some(cfg)).expect("transpile failed");
    println!(
        "  gates: {}  swaps: {}",
        out.operations.len(),
        count_swaps(&out)
    );
    check_adjacency(&out, &Backend::grid(2, 2));
}
