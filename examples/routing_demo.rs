use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() {
    println!("=== BeamSABRE Routing Demo ===\n");

    // ─── Test 1: Linear(3) — all CX already adjacent ─────────────────
    println!("--- Test 1: GHZ-3 on Linear(3) (already routable) ---");
    let qasm_ghz = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[1];
        cx q[1], q[2];
    "#;
    let circuit_ghz = parse_qasm(qasm_ghz).expect("Failed to parse");
    let config_lin = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: true,
        backend: Some(Backend::linear(3)),
    };
    let result_lin = transpile(&circuit_ghz, Some(config_lin));
    let swap_count_1 = count_swaps(&result_lin);
    println!(
        "  Gate count: {}  |  SWAPs: {}",
        result_lin.operations.len(),
        swap_count_1
    );
    check_adjacency(&result_lin, &Backend::linear(3));

    // Fidelity works here because layout is trivial (no permutation)
    let fid = unitary_fidelity(
        &circuit_to_unitary(&circuit_ghz),
        &circuit_to_unitary(&result_lin),
    );
    println!("  Unitary Fidelity: {:.10}", fid);

    // ─── Test 2: Non-adjacent CX on Linear(3) ─────────────────────────
    println!("\n--- Test 2: CX(0,2) on Linear(3) (needs routing) ---");
    let qasm_nonadj = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[2];
    "#;
    let circuit_nonadj = parse_qasm(qasm_nonadj).expect("Failed to parse");

    // With bidir_iterations > 1, the router may find a layout that avoids SWAPs
    let config_nonadj = TranspilerConfig {
        optimization_level: 1, // beam_width=1 (greedy), bidir=1 → forced SWAP
        decompose_basis: true,
        backend: Some(Backend::linear(3)),
    };
    let result_nonadj = transpile(&circuit_nonadj, Some(config_nonadj));
    let swap_count_2 = count_swaps(&result_nonadj);
    println!(
        "  Gate count: {}  |  SWAPs: {}",
        result_nonadj.operations.len(),
        swap_count_2
    );
    check_adjacency(&result_nonadj, &Backend::linear(3));

    // ─── Test 3: 4-qubit circuit on Grid(2,2) ─────────────────────────
    println!("\n--- Test 3: Star pattern on Grid(2,2) ---");
    // On grid(2,2):  0-1
    //                | |
    //                2-3
    // CX(0,3) uses the diagonal — NOT adjacent. Routing needed.
    let qasm_grid = r#"
        OPENQASM 2.0;
        qreg q[4];
        h q[0];
        cx q[0], q[1];
        cx q[0], q[2];
        cx q[0], q[3];
    "#;
    let circuit_grid = parse_qasm(qasm_grid).expect("Failed to parse");
    let config_grid = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: false, // Keep gates readable
        backend: Some(Backend::grid(2, 2)),
    };
    let result_grid = transpile(&circuit_grid, Some(config_grid));
    let swap_count_3 = count_swaps(&result_grid);
    println!(
        "  Gate count: {}  |  SWAPs: {}",
        result_grid.operations.len(),
        swap_count_3
    );
    check_adjacency(&result_grid, &Backend::grid(2, 2));

    println!("\n  Routed circuit:");
    for (i, op) in result_grid.operations.iter().enumerate() {
        println!("    {:02}: {:?}", i, op);
    }

    // ─── Summary ──────────────────────────────────────────────────────
    println!("\n=== Summary ===");
    println!(
        "  Test 1 (already routable): {} gates, {} SWAPs",
        result_lin.operations.len(),
        swap_count_1
    );
    println!(
        "  Test 2 (forced routing):   {} gates, {} SWAPs",
        result_nonadj.operations.len(),
        swap_count_2
    );
    println!(
        "  Test 3 (grid routing):     {} gates, {} SWAPs",
        result_grid.operations.len(),
        swap_count_3
    );
    println!("  All 2-qubit gates respect hardware adjacency ✅");
}

fn count_swaps(circuit: &q_rust::ir::Circuit) -> usize {
    circuit
        .operations
        .iter()
        .filter(|op| {
            matches!(
                op,
                q_rust::ir::Operation::Gate {
                    name: q_rust::ir::GateType::SWAP,
                    ..
                }
            )
        })
        .count()
}

fn check_adjacency(circuit: &q_rust::ir::Circuit, backend: &Backend) {
    for op in &circuit.operations {
        if let q_rust::ir::Operation::Gate { name, qubits, .. } = op {
            if qubits.len() == 2 {
                assert!(
                    backend.is_adjacent(qubits[0], qubits[1]),
                    "  ❌ Non-adjacent: {:?} on {:?}",
                    name,
                    qubits
                );
            }
        }
    }
    println!("  ✅ All 2-qubit gates on adjacent physical qubits");
}
