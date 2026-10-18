//! Candidate retention checked with independent analytic one-qubit, CX, CU1
//! and SWAP evolution.
//! No Q-Rust gate matrices, decomposer, simulator or verifier are used below.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, transpile_with_report, RoutingSearchConfig, TranspilerConfig};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

fn local_blocks() -> Circuit {
    let mut qasm = "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[3];".to_string();
    for (a, b) in [(0, 1), (1, 2), (0, 1)] {
        qasm += &format!(
            "u(0.21,0,0) q[{a}]; u(0,0,0.17) q[{b}]; cx q[{a}],q[{b}];
             u(0.31,0,0) q[{a}]; cx q[{b}],q[{a}]; u(0,0,0.47) q[{b}];
             cx q[{a}],q[{b}]; swap q[{a}],q[{b}];"
        );
    }
    parse_qasm(&qasm).unwrap()
}

fn evolve(circuit: &Circuit, mut state: Vec<C>) -> Vec<C> {
    for op in &circuit.operations {
        let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        else {
            assert!(matches!(op, Operation::Barrier { .. }));
            continue;
        };
        match name {
            GateType::U | GateType::H | GateType::X | GateType::RZ => {
                let phase = |theta| C::from_polar(1.0, theta);
                let zero = C::new(0.0, 0.0);
                let one = C::new(1.0, 0.0);
                let m = match name {
                    GateType::U => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [C::new(c, 0.0), -phase(params[2]) * s],
                            [phase(params[1]) * s, phase(params[1] + params[2]) * c],
                        ]
                    }
                    GateType::H => [
                        [one * FRAC_1_SQRT_2, one * FRAC_1_SQRT_2],
                        [one * FRAC_1_SQRT_2, -one * FRAC_1_SQRT_2],
                    ],
                    GateType::X => [[zero, one], [one, zero]],
                    GateType::RZ => [
                        [phase(-params[0] / 2.0), zero],
                        [zero, phase(params[0] / 2.0)],
                    ],
                    _ => unreachable!(),
                };
                for low in 0..state.len() {
                    if low & (1 << qubits[0]) != 0 {
                        continue;
                    }
                    let high = low | (1 << qubits[0]);
                    let (a, b) = (state[low], state[high]);
                    state[low] = m[0][0] * a + m[0][1] * b;
                    state[high] = m[1][0] * a + m[1][1] * b;
                }
            }
            GateType::CU1 => {
                for (basis, amplitude) in state.iter_mut().enumerate() {
                    if qubits.iter().all(|q| basis & (1 << q) != 0) {
                        *amplitude *= C::from_polar(1.0, params[0]);
                    }
                }
            }
            GateType::CX => {
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) != 0 && basis & (1 << qubits[1]) == 0 {
                        state.swap(basis, basis | (1 << qubits[1]));
                    }
                }
            }
            GateType::SWAP => {
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) == 0 && basis & (1 << qubits[1]) != 0 {
                        state.swap(basis, basis ^ (1 << qubits[0]) ^ (1 << qubits[1]));
                    }
                }
            }
            other => panic!("unsupported analytic oracle gate {other:?}"),
        }
    }
    state
}

fn verify_all_inputs(source: &Circuit, output: &Circuit, initial: &[usize], final_map: &[usize]) {
    let embed = |basis: usize, layout: &[usize]| -> usize {
        layout
            .iter()
            .enumerate()
            .map(|(q, p)| ((basis >> q) & 1) << p)
            .sum()
    };
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    for basis in 0..1 << source.num_qubits {
        let mut input = vec![C::new(0.0, 0.0); 1 << source.num_qubits];
        input[basis] = C::new(1.0, 0.0);
        let mut column = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        for (row, value) in evolve(source, input).into_iter().enumerate() {
            column[embed(row, final_map)] = value;
        }
        expected.extend(column);
        let mut input = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        input[embed(basis, initial)] = C::new(1.0, 0.0);
        actual.extend(evolve(output, input));
    }
    let overlap: C = expected
        .iter()
        .zip(&actual)
        .map(|(a, b)| a.conj() * b)
        .sum();
    assert!(overlap.norm() > 0.0);
    let phase = overlap / overlap.norm();
    let residual = expected
        .iter()
        .zip(actual)
        .map(|(a, b)| (b - phase * a).norm_sqr())
        .sum::<f64>()
        .sqrt();
    assert!(
        residual < 1e-9,
        "shared-phase full-isometry residual {residual:e}"
    );
}

fn costs(circuit: &Circuit) -> (usize, usize, usize) {
    let mut count = 0;
    let mut depth = vec![0usize; circuit.num_qubits];
    for op in &circuit.operations {
        if let Operation::Gate { qubits, .. } = op {
            if qubits.len() == 2 {
                count += 1;
                let level = 1 + depth[qubits[0]].max(depth[qubits[1]]);
                depth[qubits[0]] = level;
                depth[qubits[1]] = level;
            }
        }
    }
    (count, circuit.depth(), depth.into_iter().max().unwrap_or(0))
}

#[test]
fn retaining_local_swap_blocks_recovers_lowered_quality_without_changing_semantics() {
    let source = local_blocks();
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(3))
        .target_basis(["u", "cx"])
        .build();
    let mut no_retention = config.clone();
    no_retention.retain_routing_candidates = false;
    let mut explicit_swaps = config.clone();
    explicit_swaps.routing_search = RoutingSearchConfig {
        elide_output_permutations: false,
        ..explicit_swaps.routing_search
    };
    let (current, current_report) = transpile_with_report(&source, Some(no_retention)).unwrap();
    let (retained, retained_report) = transpile_with_report(&source, Some(explicit_swaps)).unwrap();
    let (selected, selected_report) = transpile_with_report(&source, Some(config.clone())).unwrap();
    assert_eq!(costs(&current), (12, 21, 12));
    assert_eq!(costs(&retained), (7, 14, 7));
    assert_eq!(costs(&selected), costs(&retained));
    assert_eq!(selected, transpile(&source, Some(config)).unwrap());
    for (output, report) in [
        (current, current_report),
        (retained, retained_report),
        (selected, selected_report),
    ] {
        verify_all_inputs(
            &source,
            &output,
            &report.initial_layout.unwrap(),
            &report.final_layout.unwrap(),
        );
        assert!(report.numerical_budget.charged_error <= 1e-9);
    }
}

#[test]
fn candidate_maps_preserve_idle_physical_wires_and_terminal_readout() {
    let source = local_blocks();
    let mut measured = source.clone();
    measured.num_cbits = 2;
    measured.add_op(Operation::Measure { qubit: 2, cbit: 0 });
    measured.add_op(Operation::Measure { qubit: 0, cbit: 1 });
    let mut directed = Backend::new("one-way line", 4);
    directed.set_coupling_map([(0, 1), (1, 2), (2, 3)]);
    for backend in [Backend::linear(4), directed] {
        let config = TranspilerConfig::builder()
            .optimization_level(3)
            .backend(backend.clone())
            .target_basis(["u", "cx"])
            .build();
        let (output, report) = transpile_with_report(&measured, Some(config.clone())).unwrap();
        assert_eq!(output, transpile(&measured, Some(config)).unwrap());
        let initial = report.initial_layout.unwrap();
        let final_map = report.final_layout.unwrap();
        let readout: Vec<_> = output
            .operations
            .iter()
            .filter_map(|op| match op {
                Operation::Measure { qubit, cbit } => Some((*qubit, *cbit)),
                _ => None,
            })
            .collect();
        assert_eq!(readout, vec![(final_map[2], 0), (final_map[0], 1)]);
        let mut unitary = output.clone();
        unitary
            .operations
            .retain(|op| !matches!(op, Operation::Measure { .. }));
        verify_all_inputs(&source, &unitary, &initial, &final_map);
        for op in &unitary.operations {
            if let Operation::Gate {
                name: GateType::CX,
                qubits,
                ..
            } = op
            {
                assert!(backend.has_directed_edge(qubits[0], qubits[1]));
            }
        }
    }
}

fn assert_cost_bounds(actual: (usize, usize, usize), bound: (usize, usize, usize)) {
    assert!(
        actual.0 <= bound.0 && actual.1 <= bound.1 && actual.2 <= bound.2,
        "cost {actual:?} exceeds componentwise bound {bound:?}"
    );
}

fn check_phase_regression(source: Circuit, historical_bound: (usize, usize, usize)) {
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(source.num_qubits))
        .target_basis(["u", "cx"])
        .build();
    let mut disabled = config.clone();
    disabled.retain_routing_candidates = false;
    let (current, current_report) = transpile_with_report(&source, Some(disabled)).unwrap();
    let (selected, report) = transpile_with_report(&source, Some(config.clone())).unwrap();
    assert_cost_bounds(costs(&selected), historical_bound);
    assert_cost_bounds(costs(&selected), costs(&current));
    assert_eq!(selected, transpile(&source, Some(config)).unwrap());
    verify_all_inputs(
        &source,
        &current,
        &current_report.initial_layout.unwrap(),
        &current_report.final_layout.unwrap(),
    );
    verify_all_inputs(
        &source,
        &selected,
        &report.initial_layout.unwrap(),
        &report.final_layout.unwrap(),
    );
}

#[test]
fn retained_qft_phase_schedule_recovers_pre_regression_costs() {
    // Procedural four-wire QFT with the preparation used by QASMBench qft_n4.
    let mut qasm = "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[4]; x q[0]; x q[2];".to_string();
    for target in 0..4 {
        for prior in 0..target {
            let theta = PI / (1usize << (target - prior)) as f64;
            qasm += &format!("cu1({theta}) q[{target}],q[{prior}];");
        }
        qasm += &format!("h q[{target}];");
    }
    check_phase_regression(parse_qasm(&qasm).unwrap(), (17, 29, 15));
}

#[test]
fn retained_phase_estimation_schedule_preserves_parity_rotation_opportunities() {
    // QASMBench pea_n5: powers of a phase unitary followed by inverse QFT.
    // Keeping the RZ/CX pairs matters: recovery produces commuting RZZ gates,
    // so disabling all diagonal commutation is not an equivalent candidate.
    let mut qasm = "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[5];".to_string();
    for q in 0..4 {
        qasm += &format!("h q[{q}];");
    }
    let theta = 3.0 * PI / 8.0;
    for control in (0..4).rev() {
        for _ in 0..1 << (3 - control) {
            qasm += &format!(
                "u1({}) q[4]; cx q[{control}],q[4]; u1({theta}) q[4]; cx q[{control}],q[4];",
                -theta
            );
        }
    }
    for target in 0..4 {
        for prior in 0..target {
            let theta = -PI / (1usize << (target - prior)) as f64;
            qasm += &format!("cu1({theta}) q[{prior}],q[{target}];");
        }
        qasm += &format!("h q[{target}];");
    }
    check_phase_regression(parse_qasm(&qasm).unwrap(), (27, 37, 20));
}
