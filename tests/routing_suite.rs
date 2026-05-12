//! Routing correctness test suite.

use q_rust::backend::{Backend, BackendConfig};
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::{BeamSabrePass, Layout};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn load_backend(name: &str) -> Backend {
    let path = format!("tests/fixtures/{name}");
    let json = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {path}"));
    let cfg: BackendConfig = serde_json::from_str(&json).unwrap();
    Backend::from_config(cfg)
}

fn assert_all_adjacent(circuit: &Circuit, backend: &Backend) {
    for (i, op) in circuit.operations.iter().enumerate() {
        if let Operation::Gate { name, qubits, .. } = op {
            if qubits.len() >= 2 {
                assert!(
                    backend.is_adjacent(qubits[0], qubits[1]),
                    "gate #{i} {name:?} on {qubits:?} violates adjacency on {}",
                    backend.name
                );
            }
        }
    }
}

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

fn count_2q(circuit: &Circuit) -> usize {
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

fn route(circuit: &Circuit, backend: &Backend, beam: usize, bidir: usize) -> Circuit {
    route_with_props(circuit, backend, beam, bidir).0
}

fn route_with_props(
    circuit: &Circuit,
    backend: &Backend,
    beam: usize,
    bidir: usize,
) -> (Circuit, PropertySet) {
    let pass = BeamSabrePass {
        backend: backend.clone(),
        beam_width: beam,
        branch_factor: beam.max(1),
        bidir_iterations: bidir,
    };
    let mut props = PropertySet::new();
    let routed = pass.run(circuit, &mut props);
    (routed, props)
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
        for l in 0..8 {
            let p = layout.l2p(l);
            assert_eq!(
                layout.physical_to_logical[p], l,
                "inconsistent after swap({a},{b})"
            );
        }
    }
}

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
    let r = route(&c, &b, 4, 2);
    assert_all_adjacent(&r, &b);
    assert_eq!(count_swaps(&r), 0);
}

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
    let r = route(&c, &b, 4, 2);
    assert_all_adjacent(&r, &b);
    assert_eq!(count_swaps(&r), 0);
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
    let r = route(&c, &b, 4, 2);
    assert_all_adjacent(&r, &b);
}

#[test]
fn test_star_4q_grid_2x2() {
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
    let r = route(&c, &b, 4, 2);
    assert_all_adjacent(&r, &b);
}

#[test]
fn test_all_pairs_4q_heavy_hex() {
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
    let r = route(&c, &b, 8, 3);
    assert_all_adjacent(&r, &b);
    assert_eq!(count_2q(&r), 6);
}

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
    let r = route(&c, &b, 4, 2);
    assert_eq!(r.operations.len(), 3);
    assert_eq!(count_swaps(&r), 0);
}

#[test]
fn test_empty_circuit() {
    let c = Circuit::new(3, 0);
    let b = Backend::linear(3);
    let r = route(&c, &b, 4, 2);
    assert_eq!(r.operations.len(), 0);
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
    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .decompose_basis(true)
        .backend(b.clone())
        .build();
    let r = transpile(&c, Some(cfg)).expect("transpile failed");
    assert_all_adjacent(&r, &b);
}

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
    let mut padded = routed.clone();
    padded.num_qubits = b.num_qubits;
    let u_phys = circuit_to_unitary(&padded);
    let u_log = extract_logical_unitary(&u_phys, c.num_qubits, &initial_layout, &final_layout);
    let fid = unitary_fidelity(&u_orig, &u_log);
    assert!((fid - 1.0).abs() < 1e-9, "fidelity = {fid}");
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
    let mut padded = routed.clone();
    padded.num_qubits = 5;
    let u_log = extract_logical_unitary(
        &circuit_to_unitary(&padded),
        5,
        props.get::<Vec<usize>>("initial_layout").unwrap(),
        props.get::<Vec<usize>>("final_layout").unwrap(),
    );
    assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
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
    let initial = props
        .get::<Vec<usize>>("initial_layout")
        .cloned()
        .unwrap_or_else(|| vec![0usize, 1, 2]);
    let final_l = props
        .get::<Vec<usize>>("final_layout")
        .cloned()
        .unwrap_or_else(|| vec![0usize, 1, 2]);
    let mut padded = routed.clone();
    padded.num_qubits = 3;
    let u_log = extract_logical_unitary(&circuit_to_unitary(&padded), 3, &initial, &final_l);
    assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
}
