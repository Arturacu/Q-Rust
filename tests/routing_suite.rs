//! Comprehensive routing test suite.
//!
//! Tests BeamSABRE routing across:
//! - Real IBM hardware topologies (Quito 5q, Nairobi 7q, Heavy-Hex 16q)
//! - Synthetic topologies (linear, grid, all-to-all)
//! - Diverse circuit patterns (Bell, GHZ, QFT, random CX, controlled rotations)
//! - Edge cases (single gate, already routable, more qubits than needed)
//! - Different optimization levels / beam widths

use q_rust::backend::{Backend, BackendConfig};
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::{BeamSabrePass, Layout};
use q_rust::transpiler::{transpile, TranspilerConfig};

// ─── Helpers ───────────────────────────────────────────────────────────────

/// Load a backend from a JSON fixture file.
fn load_backend(filename: &str) -> Backend {
    let path = format!("tests/fixtures/{}", filename);
    let json = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Failed to read fixture: {}", path));
    let config: BackendConfig = serde_json::from_str(&json).unwrap();
    Backend::from_config(config)
}

/// Assert that every 2-qubit gate in the circuit acts on adjacent physical qubits.
fn assert_all_adjacent(circuit: &Circuit, backend: &Backend) {
    for (i, op) in circuit.operations.iter().enumerate() {
        if let Operation::Gate { name, qubits, .. } = op {
            if qubits.len() >= 2 {
                assert!(
                    backend.is_adjacent(qubits[0], qubits[1]),
                    "Gate #{} {:?} on qubits {:?} violates adjacency on backend '{}'",
                    i,
                    name,
                    qubits,
                    backend.name
                );
            }
        }
    }
}

/// Count the number of SWAP gates in a circuit.
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

/// Count the number of 2-qubit gates (excluding SWAP) in a circuit.
fn count_2q_gates(circuit: &Circuit) -> usize {
    circuit
        .operations
        .iter()
        .filter(|op| {
            if let Operation::Gate { name, qubits, .. } = op {
                qubits.len() >= 2 && *name != GateType::SWAP
            } else {
                false
            }
        })
        .count()
}

/// Route a circuit through BeamSABRE with the given backend and parameters.
fn route(circuit: &Circuit, backend: &Backend, beam_width: usize, bidir: usize) -> Circuit {
    let (routed, _) = route_with_props(circuit, backend, beam_width, bidir);
    routed
}

/// Route a circuit and also return the resulting property set (which contains the final layout).
fn route_with_props(
    circuit: &Circuit,
    backend: &Backend,
    beam_width: usize,
    bidir: usize,
) -> (Circuit, PropertySet) {
    let pass = BeamSabrePass {
        backend: backend.clone(),
        beam_width,
        branch_factor: beam_width.max(1),
        bidir_iterations: bidir,
    };
    let mut props = PropertySet::new();
    let routed = pass.run(circuit, &mut props);
    (routed, props)
}

/// Build a circuit programmatically from a list of (gate, qubits, params).
fn build_circuit(num_qubits: usize, gates: &[(GateType, Vec<usize>, Vec<f64>)]) -> Circuit {
    let mut circuit = Circuit::new(num_qubits, 0);
    for (name, qubits, params) in gates {
        circuit.add_op(Operation::Gate {
            name: name.clone(),
            qubits: qubits.clone(),
            params: params.clone(),
        });
    }
    circuit
}

// ═══════════════════════════════════════════════════════════════════════════
// 1. LAYOUT TESTS
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_layout_roundtrip_10_swaps() {
    let mut layout = Layout::trivial(8, 8);
    let swaps = vec![
        (0, 1),
        (3, 7),
        (2, 5),
        (6, 4),
        (1, 3),
        (0, 7),
        (5, 6),
        (4, 2),
        (7, 0),
        (3, 1),
    ];
    for (a, b) in &swaps {
        layout.swap_physical(*a, *b);
        // invariant: l2p and p2l are consistent
        for l in 0..8 {
            let p = layout.l2p(l);
            assert_eq!(
                layout.physical_to_logical[p], l,
                "Inconsistent after swap({},{}): l2p[{}]={} but p2l[{}]={}",
                a, b, l, p, p, layout.physical_to_logical[p]
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// 2. ADJACENCY-ONLY ROUTING TESTS (no fidelity — just topology correctness)
// ═══════════════════════════════════════════════════════════════════════════

// ─── Bell State ────────────────────────────────────────────────────────────

#[test]
fn test_bell_state_linear_3() {
    let c = build_circuit(
        2,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_swaps(&routed), 0); // Already adjacent
}

#[test]
fn test_bell_state_ibm_quito() {
    let c = build_circuit(
        2,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ─── GHZ States ────────────────────────────────────────────────────────────

#[test]
fn test_ghz3_linear_3() {
    let c = build_circuit(
        3,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_swaps(&routed), 0); // Chain pattern matches linear
}

#[test]
fn test_ghz4_linear_4() {
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_swaps(&routed), 0);
}

#[test]
fn test_ghz5_ibm_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
            (GateType::CX, vec![3, 4], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    // Quito is T-shaped, not linear, so CX(2,3) may not be adjacent — expect SWAPs
}

#[test]
fn test_ghz5_ibm_nairobi() {
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
            (GateType::CX, vec![3, 4], vec![]),
        ],
    );
    let b = load_backend("ibm_nairobi_7q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_ghz3_grid_2x2() {
    let c = build_circuit(
        3,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
        ],
    );
    let b = Backend::grid(2, 2);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ─── Star Pattern (hub qubit talks to all others) ──────────────────────────

#[test]
fn test_star_4q_grid_2x2() {
    // CX(0,1), CX(0,2), CX(0,3) — qubit 0 connects to everything
    // On grid(2,2), 0-3 diagonal is NOT adjacent
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
        ],
    );
    let b = Backend::grid(2, 2);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_star_5q_ibm_quito() {
    // Hub qubit 0 connects to 1,2,3,4
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::CX, vec![0, 4], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
}

// ─── Reverse-order CX chain (worst case for trivial layout) ────────────────

#[test]
fn test_reverse_chain_linear_5() {
    // CX(4,3), CX(3,2), CX(2,1), CX(1,0) — reversed direction
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![4, 3], vec![]),
            (GateType::CX, vec![3, 2], vec![]),
            (GateType::CX, vec![2, 1], vec![]),
            (GateType::CX, vec![1, 0], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_swaps(&routed), 0); // Reversed but still all adjacent
}

// ─── Non-contiguous qubit pairs (skip connections) ─────────────────────────

#[test]
fn test_skip_connections_linear_5() {
    // CX(0,2), CX(1,3), CX(2,4) — all skip one qubit
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::CX, vec![2, 4], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_long_range_cx_linear_5() {
    // CX(0, 4) — maximum distance on linear(5)
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 4], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ─── Controlled Rotation Gates ─────────────────────────────────────────────

#[test]
fn test_crz_routing_linear_3() {
    let c = build_circuit(3, &[(GateType::CRZ, vec![0, 2], vec![0.5])]);
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_crx_cry_crz_ibm_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::CRX, vec![0, 2], vec![0.3]),
            (GateType::CRY, vec![1, 4], vec![0.7]),
            (GateType::CRZ, vec![3, 0], vec![1.1]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ─── CZ and CY gates ──────────────────────────────────────────────────────

#[test]
fn test_cz_cy_routing_linear_4() {
    let c = build_circuit(
        4,
        &[
            (GateType::CZ, vec![0, 2], vec![]),
            (GateType::CY, vec![1, 3], vec![]),
            (GateType::CZ, vec![0, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ─── SWAP gate in input (should be treated as a 2q gate) ───────────────────

#[test]
fn test_swap_gate_routing() {
    let c = build_circuit(3, &[(GateType::SWAP, vec![0, 2], vec![])]);
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 3. HEAVY-HEX TOPOLOGY TESTS (16 qubits — realistic IBM topology)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_ghz8_heavy_hex_16q() {
    let mut gates = vec![(GateType::H, vec![0], vec![])];
    for i in 0..7 {
        gates.push((GateType::CX, vec![i, i + 1], vec![]));
    }
    let c = build_circuit(8, &gates);
    let b = load_backend("ibm_heavy_hex_16q.json");
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_all_pairs_4q_heavy_hex() {
    // Every pair of 4 qubits needs a CX — maximally entangling
    let c = build_circuit(
        4,
        &[
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
        ],
    );
    let b = load_backend("ibm_heavy_hex_16q.json");
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_2q_gates(&routed), 6); // All 6 CX preserved
}

#[test]
fn test_long_range_heavy_hex() {
    // CX(0, 14) — opposite corners of the heavy-hex
    let c = build_circuit(16, &[(GateType::CX, vec![0, 14], vec![])]);
    let b = load_backend("ibm_heavy_hex_16q.json");
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 4. EDGE CASES
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_single_qubit_only_circuit() {
    let c = build_circuit(
        3,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::RZ, vec![1], vec![0.5]),
            (GateType::X, vec![2], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    // No 2q gates → pass-through
    assert_eq!(routed.operations.len(), 3);
    assert_eq!(count_swaps(&routed), 0);
}

#[test]
fn test_empty_circuit() {
    let c = Circuit::new(3, 0);
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_eq!(routed.operations.len(), 0);
}

#[test]
fn test_fewer_logical_than_physical() {
    // 2 logical qubits on a 7-qubit backend
    let c = build_circuit(2, &[(GateType::CX, vec![0, 1], vec![])]);
    let b = load_backend("ibm_nairobi_7q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

#[test]
fn test_circuit_with_measurements() {
    let mut c = Circuit::new(3, 3);
    c.add_op(Operation::Gate {
        name: GateType::H,
        qubits: vec![0],
        params: vec![],
    });
    c.add_op(Operation::Gate {
        name: GateType::CX,
        qubits: vec![0, 2],
        params: vec![],
    });
    c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    c.add_op(Operation::Measure { qubit: 2, cbit: 2 });

    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);

    // Measurements should still be present
    let meas_count = routed
        .operations
        .iter()
        .filter(|op| matches!(op, Operation::Measure { .. }))
        .count();
    assert_eq!(meas_count, 2);
}

#[test]
fn test_circuit_with_barriers() {
    let mut c = Circuit::new(3, 0);
    c.add_op(Operation::Gate {
        name: GateType::CX,
        qubits: vec![0, 1],
        params: vec![],
    });
    c.add_op(Operation::Barrier {
        qubits: vec![0, 1, 2],
    });
    c.add_op(Operation::Gate {
        name: GateType::CX,
        qubits: vec![1, 2],
        params: vec![],
    });

    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 5. BEAM WIDTH COMPARISON TESTS
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_beam_width_1_vs_8_quality() {
    // A harder circuit: all-pairs on 5 qubits on linear(5)
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![0, 4], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![2, 4], vec![]),
            (GateType::CX, vec![1, 4], vec![]),
        ],
    );
    let b = Backend::linear(5);

    let routed_greedy = route(&c, &b, 1, 1);
    let routed_beam = route(&c, &b, 8, 3);

    assert_all_adjacent(&routed_greedy, &b);
    assert_all_adjacent(&routed_beam, &b);

    // Beam search should use <= the number of SWAPs as greedy
    let swaps_greedy = count_swaps(&routed_greedy);
    let swaps_beam = count_swaps(&routed_beam);
    assert!(
        swaps_beam <= swaps_greedy,
        "Beam search ({} SWAPs) should not be worse than greedy ({} SWAPs)",
        swaps_beam,
        swaps_greedy
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// 6. END-TO-END PIPELINE TESTS (transpile() with backend)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_e2e_pipeline_linear_3() {
    let c = build_circuit(
        3,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: true,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

#[test]
fn test_e2e_pipeline_ibm_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::RZ, vec![2], vec![0.5]),
            (GateType::CX, vec![2, 4], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let config = TranspilerConfig {
        optimization_level: 2,
        decompose_basis: true,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

#[test]
fn test_e2e_pipeline_ibm_nairobi() {
    let c = build_circuit(
        7,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![3, 5], vec![]),
            (GateType::CX, vec![4, 6], vec![]),
            (GateType::CX, vec![0, 6], vec![]), // Long-range!
        ],
    );
    let b = load_backend("ibm_nairobi_7q.json");
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: true,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

#[test]
fn test_e2e_pipeline_heavy_hex() {
    let c = build_circuit(
        8,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![3, 7], vec![]),
            (GateType::CX, vec![4, 5], vec![]),
            (GateType::CX, vec![0, 7], vec![]), // Long-range
        ],
    );
    let b = load_backend("ibm_heavy_hex_16q.json");
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: false,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

#[test]
fn test_e2e_pipeline_no_backend_unchanged() {
    // Without a backend, routing should not be invoked
    let c = build_circuit(3, &[(GateType::CX, vec![0, 2], vec![])]);
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: false,
        backend: None,
    };
    let result = transpile(&c, Some(config));
    // CX(0,2) should remain — no routing, no adjacency enforcement
    assert_eq!(count_swaps(&result), 0);
}

// ═══════════════════════════════════════════════════════════════════════════
// 7. QFT-LIKE PATTERNS (all-pairs within each layer)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_qft_pattern_4q_linear() {
    // QFT-like: each qubit CX-connects to all qubits after it
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::H, vec![1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::H, vec![2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_2q_gates(&routed), 6); // All 6 CX preserved
}

#[test]
fn test_qft_pattern_5q_ibm_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::CX, vec![0, 4], vec![]),
            (GateType::H, vec![1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::CX, vec![1, 4], vec![]),
            (GateType::H, vec![2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
            (GateType::CX, vec![2, 4], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 8, 3);
    assert_all_adjacent(&routed, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 8. MIXED GATE SETS (single-qubit + multi-qubit interleaved)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_mixed_gates_ibm_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::RZ, vec![0], vec![0.3]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::RX, vec![1], vec![0.7]),
            (GateType::CZ, vec![1, 3], vec![]),
            (GateType::H, vec![2], vec![]),
            (GateType::CRZ, vec![2, 4], vec![1.2]),
            (GateType::Y, vec![3], vec![]),
            (GateType::CX, vec![3, 0], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);

    // All single-qubit gates should be preserved
    let sq_count = routed
        .operations
        .iter()
        .filter(|op| {
            if let Operation::Gate { qubits, .. } = op {
                qubits.len() == 1
            } else {
                false
            }
        })
        .count();
    assert!(
        sq_count >= 5,
        "Expected at least 5 single-qubit gates, got {}",
        sq_count
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// 9. REPEATED GATE PATTERNS
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_repeated_cx_same_pair() {
    let c = build_circuit(
        3,
        &[
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
    assert_eq!(count_2q_gates(&routed), 3);
}

#[test]
fn test_alternating_cx_pairs() {
    // CX(0,1) CX(2,3) CX(0,1) CX(2,3) — two independent pairs
    let c = build_circuit(
        4,
        &[
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let routed = route(&c, &b, 4, 2);
    assert_all_adjacent(&routed, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 10. OPTIMIZATION LEVEL TESTS (full pipeline)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_opt_level_0_skips_routing() {
    let c = build_circuit(3, &[(GateType::CX, vec![0, 2], vec![])]);
    let b = Backend::linear(3);
    let config = TranspilerConfig {
        optimization_level: 0,
        decompose_basis: false,
        backend: Some(b),
    };
    let result = transpile(&c, Some(config));
    // Even at level 0, routing should still work (it's hardware-correctness, not optimization)
    // The CX(0,2) should be routed to adjacent qubits
    assert_all_adjacent(&result, &Backend::linear(3));
}

#[test]
fn test_opt_level_1_greedy_routing() {
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![0, 4], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let config = TranspilerConfig {
        optimization_level: 1,
        decompose_basis: false,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

#[test]
fn test_opt_level_3_full_beam() {
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![0, 4], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
            (GateType::CX, vec![2, 0], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let config = TranspilerConfig {
        optimization_level: 3,
        decompose_basis: false,
        backend: Some(b.clone()),
    };
    let result = transpile(&c, Some(config));
    assert_all_adjacent(&result, &b);
}

// ═══════════════════════════════════════════════════════════════════════════
// 11. CORRECTNESS / FIDELITY TESTS
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_routing_fidelity_bell_state_quito() {
    let c = build_circuit(
        2,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");

    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);

    let final_layout: Vec<usize> = props.get::<Vec<usize>>("final_layout").unwrap().clone();
    let initial_layout: Vec<usize> = props.get::<Vec<usize>>("initial_layout").unwrap().clone();

    let mut padded_routed = routed.clone();
    padded_routed.num_qubits = b.num_qubits;
    let u_routed_physical = circuit_to_unitary(&padded_routed);

    let u_routed_logical = extract_logical_unitary(
        &u_routed_physical,
        c.num_qubits,
        &initial_layout,
        &final_layout,
    );

    let fid = unitary_fidelity(&u_orig, &u_routed_logical);
    assert!(
        (fid - 1.0).abs() < 1e-9,
        "Fidelity is {}, expected 1.0",
        fid
    );
}

#[test]
fn test_routing_fidelity_ghz4_linear() {
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![2, 3], vec![]),
        ],
    );
    // A linear 4 backend, but let's reverse the circuit so SWAPs are needed
    let c_reversed = build_circuit(
        4,
        &[
            (GateType::H, vec![3], vec![]),
            (GateType::CX, vec![3, 0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
        ],
    );
    let b = Backend::linear(4);

    let u_orig = circuit_to_unitary(&c_reversed);
    let (routed, props) = route_with_props(&c_reversed, &b, 4, 3);

    let final_layout: Vec<usize> = props.get::<Vec<usize>>("final_layout").unwrap().clone();
    let initial_layout: Vec<usize> = props.get::<Vec<usize>>("initial_layout").unwrap().clone();

    let mut padded_routed = routed.clone();
    padded_routed.num_qubits = b.num_qubits;
    let u_routed_physical = circuit_to_unitary(&padded_routed);

    let u_routed_logical = extract_logical_unitary(
        &u_routed_physical,
        c_reversed.num_qubits,
        &initial_layout,
        &final_layout,
    );

    let fid = unitary_fidelity(&u_orig, &u_routed_logical);
    assert!(
        (fid - 1.0).abs() < 1e-9,
        "Fidelity is {}, expected 1.0",
        fid
    );
}

#[test]
fn test_routing_fidelity_qft_pattern_nairobi() {
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
            (GateType::H, vec![1], vec![]),
            (GateType::CX, vec![1, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
        ],
    );
    let b = load_backend("ibm_nairobi_7q.json"); // 7q backend
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 7;
            r
        }),
        4,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_star_grid2x2() {
    let c = build_circuit(
        4,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![0, 3], vec![]),
        ],
    );
    let b = Backend::grid(2, 2); // 4q
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 4;
            r
        }),
        4,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_reverse_chain_linear5() {
    let c = build_circuit(
        5,
        &[
            (GateType::CX, vec![4, 3], vec![]),
            (GateType::CX, vec![3, 2], vec![]),
            (GateType::CX, vec![2, 1], vec![]),
            (GateType::CX, vec![1, 0], vec![]),
        ],
    );
    let b = Backend::linear(5);
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 5;
            r
        }),
        5,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_skip_connection_linear4() {
    let c = build_circuit(
        4,
        &[
            (GateType::CX, vec![0, 2], vec![]),
            (GateType::CX, vec![1, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 4;
            r
        }),
        4,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_mixed_gates_quito() {
    let c = build_circuit(
        5,
        &[
            (GateType::CRX, vec![0, 2], vec![0.3]),
            (GateType::CRY, vec![1, 4], vec![0.7]),
            (GateType::CRZ, vec![3, 0], vec![1.1]),
        ],
    );
    let b = load_backend("ibm_quito_5q.json");
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 5;
            r
        }),
        5,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_cz_cy_linear4() {
    let c = build_circuit(
        4,
        &[
            (GateType::CZ, vec![0, 2], vec![]),
            (GateType::CY, vec![1, 3], vec![]),
        ],
    );
    let b = Backend::linear(4);
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 8, 3);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 4;
            r
        }),
        4,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_single_qubits_linear3() {
    let c = build_circuit(
        3,
        &[
            (GateType::H, vec![0], vec![]),
            (GateType::RZ, vec![1], vec![0.5]),
            (GateType::X, vec![2], vec![]),
        ],
    );
    let b = Backend::linear(3);
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);

    // If there are no 2q gates, BeamSabre skips and leaves layouts empty. Default to trivial.
    let initial_layout = props
        .get::<Vec<usize>>("initial_layout")
        .cloned()
        .unwrap_or_else(|| vec![0usize, 1, 2]);
    let final_layout = props
        .get::<Vec<usize>>("final_layout")
        .cloned()
        .unwrap_or_else(|| vec![0usize, 1, 2]);

    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 3;
            r
        }),
        3,
        &initial_layout,
        &final_layout,
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}

#[test]
fn test_routing_fidelity_swap_in_input_grid2x2() {
    let c = build_circuit(
        4,
        &[
            (GateType::CX, vec![0, 1], vec![]),
            (GateType::SWAP, vec![0, 3], vec![]), // explicit swap in original
            (GateType::CX, vec![1, 2], vec![]),
        ],
    );
    let b = Backend::grid(2, 2);
    let u_orig = circuit_to_unitary(&c);
    let (routed, props) = route_with_props(&c, &b, 4, 2);
    let u_routed_logical = extract_logical_unitary(
        &circuit_to_unitary(&{
            let mut r = routed.clone();
            r.num_qubits = 4;
            r
        }),
        4,
        &props.get::<Vec<usize>>("initial_layout").unwrap(),
        &props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_routed_logical) - 1.0).abs() < 1e-9);
}
