//! Independent all-input checks for encoded permutations and Clifford rewrites.
//!
//! This oracle uses literal matrices and bit-indexed amplitudes, not Q-Rust's
//! gate definitions, simulator, synthesizer, or verifier. Every comparison uses
//! one shared phase across the entire clean-ancilla isometry. Consequently it
//! checks relative phases and arbitrary superpositions, including entanglement
//! with an external reference, rather than just computational-basis outcomes.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, ClassicalCondition, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::clifford_simplification::CliffordSimplificationPass;
use q_rust::transpiler::encoded_swaps::EncodedSwapRecognitionPass;
use q_rust::transpiler::numerical_policy::CIRCUIT_ERROR_ALLOWANCE;
use q_rust::transpiler::optimization::{
    ApplyOutputPermutationPass, ComposeOutputPermutationPass, ExtractOutputPermutationPass,
};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::{transpile, transpile_with_report, RoutingSearchConfig, TranspilerConfig};
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};

fn gate(name: GateType, wires: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: wires.to_vec(),
        params: params.to_vec(),
    }
}

fn source(width: usize, body: &str) -> Circuit {
    parse_qasm(&format!(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[{width}]; {body}"
    ))
    .unwrap()
}

fn phase(theta: f64) -> C {
    C::from_polar(1.0, theta)
}

fn one_qubit(state: &mut [C], q: usize, control: Option<usize>, u: [[C; 2]; 2]) {
    for low in 0..state.len() {
        if low & (1 << q) != 0 || control.is_some_and(|c| low & (1 << c) == 0) {
            continue;
        }
        let high = low | (1 << q);
        let (a, b) = (state[low], state[high]);
        state[low] = u[0][0] * a + u[0][1] * b;
        state[high] = u[1][0] * a + u[1][1] * b;
    }
}

fn evolve(circuit: &Circuit, mut state: Vec<C>) -> Vec<C> {
    let zero = C::new(0.0, 0.0);
    let one = C::new(1.0, 0.0);
    let i = C::new(0.0, 1.0);
    assert_eq!(state.len(), 1 << circuit.num_qubits);
    for op in &circuit.operations {
        let Operation::Gate {
            name,
            params,
            qubits,
        } = op
        else {
            assert!(matches!(op, Operation::Barrier { .. }));
            continue;
        };
        assert!(params.iter().all(|x| x.is_finite()));
        match name {
            GateType::CX => {
                for low in 0..state.len() {
                    if low & (1 << qubits[0]) != 0 && low & (1 << qubits[1]) == 0 {
                        state.swap(low, low | (1 << qubits[1]));
                    }
                }
            }
            GateType::CZ => {
                for (basis, amplitude) in state.iter_mut().enumerate() {
                    if qubits.iter().all(|q| basis & (1 << q) != 0) {
                        *amplitude = -*amplitude;
                    }
                }
            }
            GateType::SWAP => {
                for low in 0..state.len() {
                    if low & (1 << qubits[0]) == 0 && low & (1 << qubits[1]) != 0 {
                        state.swap(low, low ^ (1 << qubits[0]) ^ (1 << qubits[1]));
                    }
                }
            }
            _ => {
                let u = match name {
                    GateType::ID => [[one, zero], [zero, one]],
                    GateType::H => [
                        [one * FRAC_1_SQRT_2, one * FRAC_1_SQRT_2],
                        [one * FRAC_1_SQRT_2, -one * FRAC_1_SQRT_2],
                    ],
                    GateType::X => [[zero, one], [one, zero]],
                    GateType::Y | GateType::CY => [[zero, -i], [i, zero]],
                    GateType::Z => [[one, zero], [zero, -one]],
                    GateType::S => [[one, zero], [zero, i]],
                    GateType::Sdg => [[one, zero], [zero, -i]],
                    GateType::T => [[one, zero], [zero, phase(PI / 4.0)]],
                    GateType::Tdg => [[one, zero], [zero, phase(-PI / 4.0)]],
                    GateType::SX => [
                        [(one + i) / 2.0, (one - i) / 2.0],
                        [(one - i) / 2.0, (one + i) / 2.0],
                    ],
                    GateType::RX => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [[one * c, -i * s], [-i * s, one * c]]
                    }
                    GateType::RY => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [[one * c, -one * s], [one * s, one * c]]
                    }
                    GateType::RZ => [
                        [phase(-params[0] / 2.0), zero],
                        [zero, phase(params[0] / 2.0)],
                    ],
                    GateType::U => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [one * c, -phase(params[2]) * s],
                            [phase(params[1]) * s, phase(params[1] + params[2]) * c],
                        ]
                    }
                    _ => panic!("unsupported independent oracle gate: {name:?}"),
                };
                let (q, control) = if matches!(name, GateType::CY) {
                    (qubits[1], Some(qubits[0]))
                } else {
                    (qubits[0], None)
                };
                one_qubit(&mut state, q, control, u);
            }
        }
    }
    state
}

fn embed(value: usize, layout: &[usize]) -> usize {
    layout
        .iter()
        .enumerate()
        .map(|(l, p)| ((value >> l) & 1) << p)
        .sum()
}

fn isometry_error(
    source: &Circuit,
    output: &Circuit,
    initial: &[usize],
    final_map: &[usize],
) -> f64 {
    assert_eq!(initial.len(), source.num_qubits);
    assert_eq!(final_map.len(), source.num_qubits);
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    for column in 0..1 << source.num_qubits {
        let mut logical_input = vec![C::new(0.0, 0.0); 1 << source.num_qubits];
        logical_input[column] = C::new(1.0, 0.0);
        let logical = evolve(source, logical_input);
        let mut physical = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        physical[embed(column, initial)] = C::new(1.0, 0.0);
        let observed = evolve(output, physical);
        let mut wanted = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        for (row, amplitude) in logical.into_iter().enumerate() {
            wanted[embed(row, final_map)] = amplitude;
        }
        expected.extend(wanted);
        actual.extend(observed);
    }
    let overlap: C = expected
        .iter()
        .zip(&actual)
        .map(|(e, a)| e.conj() * a)
        .sum();
    let global = if overlap.norm() > 1e-15 {
        overlap / overlap.norm()
    } else {
        C::new(1.0, 0.0)
    };
    // Frobenius rather than separately normalized columns: relative phases and
    // all leakage onto supposedly clean physical spare wires count as errors.
    expected
        .iter()
        .zip(&actual)
        .map(|(e, a)| (*a - global * e).norm_sqr())
        .sum::<f64>()
        .sqrt()
}

fn assert_same(source: &Circuit, actual: &Circuit) {
    let identity = (0..source.num_qubits).collect::<Vec<_>>();
    let error = isometry_error(source, actual, &identity, &identity);
    assert!(
        error < 2e-11,
        "full-operator shared-phase residual {error:e}"
    );
}

fn clifford(source: &Circuit, allow_swaps: bool) -> Circuit {
    CliffordSimplificationPass {
        allow_swap_replacements: allow_swaps,
    }
    .run(source, &mut PropertySet::new())
}

fn cx_count(circuit: &Circuit) -> usize {
    circuit
        .operations
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

#[test]
fn independent_oracle_rejects_relative_phase_and_dirty_spares() {
    let id = source(2, "");
    assert!(isometry_error(&id, &source(2, "z q[0];"), &[0, 1], &[0, 1]) > 1.0);
    assert!(isometry_error(&id, &source(3, "x q[2];"), &[0, 1], &[0, 1]) > 1.0);
    // X Y Z differs from identity only by a common phase.
    assert_same(&id, &source(2, "x q[0]; y q[0]; z q[0];"));
    // A wrong three-cycle cannot hide behind an inverse-equals-forward swap.
    let permutations = source(3, "swap q[0],q[1]; swap q[1],q[2];");
    assert!(isometry_error(&permutations, &source(3, ""), &[0, 1, 2], &[2, 0, 1]) > 1.0);
}

#[test]
fn cross_block_pair_reduction_preserves_phases_and_both_wire_orders() {
    for (a, b, other) in [(0, 2, 1), (2, 0, 1), (1, 0, 2)] {
        for local in ["s", "sdg"] {
            let circuit = source(
                3,
                &format!(
                    "cx q[{a}],q[{b}]; {local} q[{a}]; cx q[{a}],q[{other}]; cx q[{a}],q[{b}];"
                ),
            );
            for allow_swaps in [false, true] {
                let output = clifford(&circuit, allow_swaps);
                assert_same(&circuit, &output);
                assert!(
                    cx_count(&output) < cx_count(&circuit),
                    "cross-block pair was not reduced"
                );
            }
            // Transporting a target phase creates Y axes and their signs.
            // Its best CX count may trade off against local-gate depth, so
            // acceptance is deliberately not required for this second case.
            let y_axis = source(
                3,
                &format!(
                    "cx q[{a}],q[{b}]; {local} q[{b}]; cx q[{a}],q[{other}]; cx q[{a}],q[{b}];"
                ),
            );
            assert_same(&y_axis, &clifford(&y_axis, true));
        }
    }
}

#[test]
fn noncommuting_cross_wire_interactions_do_not_authorize_false_cancellation() {
    for middle in ["cx q[1],q[2];", "cy q[1],q[2];", "cz q[0],q[2];"] {
        for local in ["h q[0];", "s q[1];", "sdg q[0]; h q[1];"] {
            let circuit = source(3, &format!("cx q[0],q[1]; {local} {middle} cx q[0],q[1];"));
            assert_same(&circuit, &clifford(&circuit, true));
        }
    }
}

#[test]
fn every_signed_controlled_axis_pair_matches_an_independent_dense_operator() {
    // Each word conjugates the first endpoint's Z_control or X_target onto
    // one of ±X, ±Y, ±Z. This tests all 72 signed-axis/orientation products
    // independently of the production Clifford table's own consistency tests.
    let from_z = ["", "x", "h", "h z", "h s", "h sdg"];
    let from_x = ["", "z", "h", "h x", "s", "sdg"];
    let mut changed = 0;
    for control_frame in from_z {
        for target_frame in from_x {
            for reverse in [false, true] {
                let mut circuit = source(3, "cx q[2],q[0];");
                for (word, wire) in [(control_frame, 2), (target_frame, 0)] {
                    for name in word.split_whitespace() {
                        circuit
                            .operations
                            .extend(source(3, &format!("{name} q[{wire}];")).operations);
                    }
                }
                circuit.add_op(gate(
                    GateType::CX,
                    if reverse { &[0, 2] } else { &[2, 0] },
                    &[],
                ));
                for allow_swaps in [false, true] {
                    let output = clifford(&circuit, allow_swaps);
                    assert_same(&circuit, &output);
                    changed += usize::from(output.operations != circuit.operations);
                }
            }
        }
    }
    assert!(
        changed > 0,
        "axis enumeration must exercise accepted rewrites"
    );
}

#[test]
fn arbitrary_clifford_networks_preserve_all_columns_with_a_shared_phase() {
    // A small deterministic generator avoids any dependency on a testing RNG.
    // Its cases include reversed/nonadjacent wire orders and disjoint gates.
    let mut random = 0x4b8a_61cf_2e7d_0935u64;
    let mut next = || {
        random = random
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (random >> 32) as usize
    };
    for width in 2..=5 {
        for _case in 0..12 {
            let mut circuit = Circuit::new(width, 0);
            for _ in 0..45 {
                let a = next() % width;
                let choice = next() % 13;
                let b = (a + 1 + next() % (width - 1)) % width;
                let (name, wires) = match choice {
                    0 => (GateType::H, vec![a]),
                    1 => (GateType::S, vec![a]),
                    2 => (GateType::Sdg, vec![a]),
                    3 => (GateType::X, vec![a]),
                    4 => (GateType::Y, vec![a]),
                    5 => (GateType::Z, vec![a]),
                    6 => (GateType::SX, vec![a]),
                    7 => (GateType::CY, vec![a, b]),
                    8 => (GateType::CZ, vec![a, b]),
                    _ => (GateType::CX, vec![a, b]),
                };
                circuit.add_op(gate(name, &wires, &[]));
            }
            for allow_swaps in [false, true] {
                assert_same(&circuit, &clifford(&circuit, allow_swaps));
            }
        }
    }
}

#[test]
fn nonclifford_and_near_clifford_angles_are_not_silently_rounded() {
    for angle in [
        FRAC_PI_2 + 9e-10,
        FRAC_PI_2 - 9e-10,
        FRAC_PI_2 + 1e-13,
        9e-10,
        0.37,
    ] {
        let circuit = source(
            3,
            &format!("cx q[0],q[2]; rz({angle}) q[2]; cx q[0],q[1]; cx q[0],q[2];"),
        );
        let output = clifford(&circuit, true);
        assert_same(&circuit, &output);
        assert!(output.operations.iter().any(|op| matches!(op, Operation::Gate { name: GateType::RZ, params, .. } if params[0].to_bits() == angle.to_bits())));
    }
}

#[test]
fn barriers_and_dynamic_boundaries_are_never_crossed() {
    let endpoint = gate(GateType::CX, &[0, 1], &[]);
    for boundary in [
        Operation::Barrier { qubits: vec![0, 1] },
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Reset { qubit: 1 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 1,
            },
            op: Box::new(gate(GateType::X, &[1], &[])),
        },
    ] {
        let mut circuit = Circuit::new(2, 1);
        circuit.operations = vec![endpoint.clone(), boundary, endpoint.clone()];
        assert_eq!(clifford(&circuit, true).operations, circuit.operations);
    }
}

fn encoded_swap(a: usize, b: usize) -> String {
    format!("cx q[{a}],q[{b}]; h q[{a}]; h q[{b}]; cx q[{a}],q[{b}]; h q[{a}]; h q[{b}]; cx q[{a}],q[{b}];")
}

#[test]
fn encoded_swap_recognition_preserves_global_phase_and_disjoint_operations() {
    for (a, b) in [(0, 2), (2, 0)] {
        let circuit = source(
            4,
            &format!(
                "x q[{a}]; y q[{a}]; z q[{a}];
             cx q[{a}],q[{b}]; h q[{b}]; u(0.4,-0.3,0.7) q[1]; h q[{a}];
             cx q[{a}],q[{b}]; h q[{a}]; h q[{b}]; s q[3]; cx q[{a}],q[{b}];
             u(-0.5,0.6,0.2) q[{a}]; cx q[{b}],q[1];"
            ),
        );
        let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
        assert_same(&circuit, &output);
        assert!(output.operations.iter().any(|op| matches!(
            op,
            Operation::Gate {
                name: GateType::SWAP,
                ..
            }
        )));
        assert!(cx_count(&output) <= cx_count(&circuit) - 3);
    }
    for (a, b) in [(0, 1), (1, 0)] {
        let circuit = source(
            2,
            &format!("cx q[{a}],q[{b}]; cx q[{b}],q[{a}]; cx q[{a}],q[{b}];"),
        );
        let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
        assert_same(&circuit, &output);
        assert_eq!(output.operations.len(), 1);
        assert!(matches!(
            &output.operations[0],
            Operation::Gate {
                name: GateType::SWAP,
                ..
            }
        ));
    }
}

#[test]
fn encoded_swap_recognition_does_not_cross_fences_or_approximate_rotations() {
    let circuit = source(4, &encoded_swap(0, 2));
    for boundary in [
        Operation::Barrier { qubits: vec![3] },
        Operation::Measure { qubit: 3, cbit: 0 },
        Operation::Reset { qubit: 3 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 1,
            },
            op: Box::new(gate(GateType::X, &[3], &[])),
        },
    ] {
        let mut fenced = circuit.clone();
        fenced.num_cbits = 1;
        fenced.operations.insert(2, boundary);
        let output = EncodedSwapRecognitionPass.run(&fenced, &mut PropertySet::new());
        assert_eq!(output.operations, fenced.operations);
    }
    for epsilon in [9e-10, -9e-10, 1e-13] {
        let mut almost = circuit.clone();
        // U(pi/2, 0, pi) is H, but this has a real relative-phase offset.
        almost.operations[1] = gate(GateType::U, &[0], &[FRAC_PI_2, epsilon, PI]);
        let output = EncodedSwapRecognitionPass.run(&almost, &mut PropertySet::new());
        assert_eq!(output.operations, almost.operations);
        assert_same(&almost, &output);
    }
}

#[test]
fn encoded_swap_recognition_keeps_relative_pauli_signs_and_third_wire_dependencies() {
    for prefix in ["z q[0];", "x q[2];", "y q[0];", "cz q[0],q[2];"] {
        let circuit = source(3, &format!("{prefix} {}", encoded_swap(0, 2)));
        let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
        assert_same(&circuit, &output);
    }
    for bridge in [
        gate(GateType::CX, &[0, 1], &[]),
        gate(GateType::CY, &[1, 2], &[]),
    ] {
        let mut circuit = source(3, &encoded_swap(0, 2));
        circuit.operations.insert(3, bridge);
        let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
        assert_eq!(output.operations, circuit.operations);
        assert_same(&circuit, &output);
    }
}

#[test]
fn independent_interleaved_swap_windows_keep_their_original_wire_dependencies() {
    let a = source(4, &encoded_swap(0, 2));
    let b = source(4, &encoded_swap(3, 1));
    let mut circuit = Circuit::new(4, 0);
    for (left, right) in a.operations.into_iter().zip(b.operations) {
        circuit.operations.extend([left, right]);
    }
    let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
    assert_same(&circuit, &output);
    assert_eq!(output.operations.len(), 2);
    assert!(output.operations.iter().all(|op| matches!(
        op,
        Operation::Gate {
            name: GateType::SWAP,
            ..
        }
    )));
}

#[test]
fn generated_clifford_prefixes_and_inverse_words_do_not_corrupt_swap_recognition() {
    let mut random = 0x17cb_043a_728e_d659u64;
    let mut next = || {
        random = random
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (random >> 32) as usize
    };
    let mut recognized = 0;
    for case in 0..36 {
        let pair = if case % 2 == 0 { [0, 3] } else { [3, 0] };
        let mut forward = Vec::new();
        let mut inverse_words = Vec::new();
        for _ in 0..9 {
            let wire = pair[next() % 2];
            let (name, wires) = match next() % 9 {
                0 => (GateType::H, vec![wire]),
                1 => (GateType::S, vec![wire]),
                2 => (GateType::Sdg, vec![wire]),
                3 => (GateType::SX, vec![wire]),
                4 => (GateType::Y, vec![wire]),
                5 => (GateType::X, vec![wire]),
                6 => (GateType::CY, pair.to_vec()),
                7 => (GateType::CZ, pair.to_vec()),
                _ => (GateType::CX, pair.to_vec()),
            };
            let inverse = match name {
                GateType::S => vec![gate(GateType::Sdg, &wires, &[])],
                GateType::Sdg => vec![gate(GateType::S, &wires, &[])],
                GateType::SX => vec![gate(GateType::SX, &wires, &[]); 3],
                _ => vec![gate(name.clone(), &wires, &[])],
            };
            forward.push(gate(name, &wires, &[]));
            inverse_words.push(inverse);
        }
        let mut circuit = Circuit::new(4, 0);
        circuit.operations.extend(forward);
        for word in inverse_words.into_iter().rev() {
            circuit.operations.extend(word);
        }
        circuit
            .operations
            .extend(source(4, &encoded_swap(pair[0], pair[1])).operations);
        // An unsupported local rotation and an entangler on wholly disjoint
        // wires must retain their order while the selected window contracts.
        circuit
            .operations
            .insert(5, gate(GateType::U, &[1], &[0.3, -0.7, 0.2]));
        circuit
            .operations
            .insert(13, gate(GateType::CX, &[1, 2], &[]));
        let output = EncodedSwapRecognitionPass.run(&circuit, &mut PropertySet::new());
        assert_same(&circuit, &output);
        recognized += usize::from(output.operations.iter().any(|op| {
            matches!(
                op,
                Operation::Gate {
                    name: GateType::SWAP,
                    ..
                }
            )
        }));
    }
    assert_eq!(recognized, 36);
}

#[test]
fn repeated_near_hadamards_do_not_accumulate_unreported_rounding_error() {
    let mut circuit = Circuit::new(2, 0);
    for _ in 0..96 {
        let mut almost = source(2, &encoded_swap(0, 1));
        almost.operations[1] = gate(GateType::U, &[0], &[FRAC_PI_2, 9e-10, PI]);
        circuit.operations.extend(almost.operations);
    }
    let config = TranspilerConfig::builder()
        .backend(Backend::all_to_all(2))
        .target_basis(["u", "cx"])
        .optimization_level(3)
        .retain_routing_candidates(false)
        .post_routing_cleanup(false)
        .build();
    let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
    let error = isometry_error(
        &circuit,
        &output,
        report.initial_layout.as_ref().unwrap(),
        report.final_layout.as_ref().unwrap(),
    );
    assert!(error < 1e-10, "accumulated near-H residual {error:e}");
    assert!(report.numerical_budget.charged_error <= CIRCUIT_ERROR_ALLOWANCE);
}

#[test]
fn independent_optimization_switches_preserve_the_same_routed_isometry() {
    let circuit = source(
        3,
        &format!(
            "h q[1]; cx q[1],q[2]; {} cx q[0],q[2]; s q[2]; cx q[0],q[1]; cx q[0],q[2];",
            encoded_swap(2, 0)
        ),
    );
    for recognize in [false, true] {
        for simplify in [false, true] {
            for elide in [false, true] {
                let config = TranspilerConfig::builder()
                    .backend(Backend::linear(4))
                    .target_basis(["u", "cx"])
                    .optimization_level(3)
                    .encoded_swap_recognition(recognize)
                    .clifford_simplification(simplify)
                    .routing_search(RoutingSearchConfig {
                        elide_output_permutations: elide,
                        ..Default::default()
                    })
                    .build();
                let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
                let error = isometry_error(
                    &circuit,
                    &output,
                    &report.initial_layout.as_ref().unwrap()[..3],
                    &report.final_layout.as_ref().unwrap()[..3],
                );
                assert!(
                    error < 2e-8,
                    "recognize {recognize}, simplify {simplify}, elide {elide}: {error:e}"
                );
            }
        }
    }
}

#[test]
fn disabling_interior_permutation_elision_keeps_the_encoded_swap_physical() {
    // Gates on both output wires ensure this is not a trailing-SWAP case;
    // the pre-existing trailing-elision option has a separate contract.
    let circuit = source(
        2,
        &format!(
            "{} u(0.37,-0.2,0.6) q[0]; u(-0.41,0.3,0.8) q[1];",
            encoded_swap(0, 1)
        ),
    );
    for elide in [false, true] {
        let config = TranspilerConfig::builder()
            .backend(Backend::all_to_all(2))
            .target_basis(["u", "cx"])
            .optimization_level(3)
            .encoded_swap_recognition(true)
            .clifford_simplification(false)
            .post_routing_cleanup(false)
            .retain_routing_candidates(false)
            .routing_search(RoutingSearchConfig {
                elide_output_permutations: elide,
                ..Default::default()
            })
            .build();
        let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
        assert_eq!(cx_count(&output), if elide { 0 } else { 3 });
        if !elide {
            assert_eq!(report.initial_layout, report.final_layout);
        }
        let error = isometry_error(
            &circuit,
            &output,
            report.initial_layout.as_ref().unwrap(),
            report.final_layout.as_ref().unwrap(),
        );
        assert!(error < 2e-11);
    }
}

#[test]
fn two_noninvolutive_extractions_compose_in_logical_to_physical_order() {
    let circuit = source(
        4,
        &format!(
            "swap q[0],q[1]; swap q[1],q[2]; u(0.3,-0.2,0.7) q[0];
         {} {} u(-0.4,0.6,0.2) q[3]; cx q[0],q[2];",
            encoded_swap(1, 2),
            encoded_swap(2, 3)
        ),
    );
    let mut properties = PropertySet::new();
    let first = ExtractOutputPermutationPass.run(&circuit, &mut properties);
    let p1 = properties
        .get::<Vec<usize>>("logical_output_permutation")
        .unwrap()
        .clone();
    assert_eq!(p1, vec![1, 2, 0, 3]);
    let recognized = EncodedSwapRecognitionPass.run(&first, &mut properties);
    let mut isolated = PropertySet::new();
    ExtractOutputPermutationPass.run(&recognized, &mut isolated);
    let p2 = isolated
        .get::<Vec<usize>>("logical_output_permutation")
        .unwrap();
    assert_eq!(p2, &vec![3, 1, 0, 2]);
    // Both maps are genuine three-cycles and they do not commute.
    assert_ne!(
        p1.iter().map(|&q| p1[q]).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    assert_ne!(
        p2.iter().map(|&q| p2[q]).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    let composed = p1.iter().map(|&q| p2[q]).collect::<Vec<_>>();
    assert_ne!(composed, p2.iter().map(|&q| p1[q]).collect::<Vec<_>>());
    let output = ComposeOutputPermutationPass.run(&recognized, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("logical_output_permutation"),
        Some(&composed)
    );
    assert!(isometry_error(&circuit, &output, &[0, 1, 2, 3], &composed) < 2e-11);
    let physical = vec![4usize, 0, 3, 1];
    properties.insert("final_layout", physical.clone());
    ApplyOutputPermutationPass.run(&output, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&composed.iter().map(|&q| physical[q]).collect::<Vec<_>>())
    );
}

#[test]
fn composition_preserves_prior_maps_without_new_swaps_and_rejects_invalid_maps() {
    let unchanged = source(3, "h q[0]; cx q[0],q[2]; u(0.2,0.3,-0.4) q[1];");
    let previous = vec![1usize, 2, 0];
    let mut properties = PropertySet::new();
    properties.insert("logical_output_permutation", previous.clone());
    assert_eq!(
        ComposeOutputPermutationPass
            .run(&unchanged, &mut properties)
            .operations,
        unchanged.operations
    );
    assert_eq!(
        properties.get::<Vec<usize>>("logical_output_permutation"),
        Some(&previous)
    );
    for invalid in [
        vec![0usize, 0, 2],
        vec![1, 2],
        vec![0, 1, 3],
        vec![0, 1, 2, 3],
    ] {
        let circuit = source(3, "swap q[0],q[2]; h q[0];");
        properties.insert("logical_output_permutation", invalid.clone());
        let output = ComposeOutputPermutationPass.run(&circuit, &mut properties);
        assert_eq!(output.operations, circuit.operations);
        assert_eq!(
            properties.get::<Vec<usize>>("logical_output_permutation"),
            Some(&invalid)
        );
    }
    for suffix in [
        Operation::Measure { qubit: 2, cbit: 0 },
        Operation::Reset { qubit: 0 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 1,
            },
            op: Box::new(gate(GateType::X, &[1], &[])),
        },
    ] {
        let mut circuit = source(3, "swap q[0],q[2]; h q[0];");
        circuit.num_cbits = 1;
        circuit.add_op(suffix);
        properties.insert("logical_output_permutation", previous.clone());
        assert_eq!(
            ComposeOutputPermutationPass
                .run(&circuit, &mut properties)
                .operations,
            circuit.operations
        );
        assert_eq!(
            properties.get::<Vec<usize>>("logical_output_permutation"),
            Some(&previous)
        );
    }
}

#[test]
fn combined_w_transfer_pipeline_preserves_arbitrary_inputs_and_clean_spares() {
    for width in [3, 4] {
        let mut circuit = Circuit::new(width, 0);
        circuit.add_op(gate(GateType::X, &[0], &[]));
        for target in 1..width {
            let control = target - 1;
            let angle = 2.0 * (1.0 / ((width - target + 1) as f64).sqrt()).acos();
            circuit.operations.extend([
                gate(GateType::RY, &[target], &[angle / 2.0]),
                gate(GateType::CX, &[control, target], &[]),
                gate(GateType::RY, &[target], &[-angle / 2.0]),
                gate(GateType::CX, &[control, target], &[]),
                gate(GateType::CX, &[target, control], &[]),
            ]);
        }
        let mut directed = Backend::new("backwards_line", width + 1);
        directed.set_coupling_map((0..width).map(|q| (q + 1, q)));
        for backend in [
            Backend::all_to_all(width),
            Backend::linear(width + 1),
            directed,
        ] {
            let config = TranspilerConfig::builder()
                .optimization_level(3)
                .backend(backend.clone())
                .target_basis(["u", "cx"])
                .build();
            let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
            let error = isometry_error(
                &circuit,
                &output,
                &report.initial_layout.as_ref().unwrap()[..width],
                &report.final_layout.as_ref().unwrap()[..width],
            );
            assert!(
                error < 2e-8,
                "W-{width} combined pipeline isometry residual {error:e}"
            );
            for op in &output.operations {
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
}

#[test]
fn new_switches_have_no_effect_below_optimization_level_three() {
    let circuit = source(
        3,
        &format!(
            "{} cx q[0],q[1]; s q[0]; cx q[0],q[2]; cx q[0],q[1];",
            encoded_swap(0, 2)
        ),
    );
    for level in 0..=2 {
        let config = |enabled| {
            TranspilerConfig::builder()
                .optimization_level(level)
                .backend(Backend::linear(3))
                .encoded_swap_recognition(enabled)
                .clifford_simplification(enabled)
                .target_basis(["u", "cx"])
                .build()
        };
        let (off, off_report) = transpile_with_report(&circuit, Some(config(false))).unwrap();
        let (on, on_report) = transpile_with_report(&circuit, Some(config(true))).unwrap();
        assert_eq!(off.operations, on.operations);
        assert_eq!(off_report, on_report);
        assert!(on_report.encoded_swap_stats.is_none());
        assert!(on_report.clifford_simplification_stats.is_none());
    }
}

#[test]
fn encoded_three_cycle_preserves_routed_isometry_and_plain_report_parity() {
    let circuit = source(
        3,
        &format!(
            "u(0.3,-0.2,0.7) q[0]; {} u(-0.7,0.9,0.2) q[0]; {} cx q[2],q[0]; u(0.2,0.6,-0.8) q[1];",
            encoded_swap(0, 2),
            encoded_swap(2, 1)
        ),
    );
    let mut directed = Backend::new("forward_line_with_spare", 4);
    directed.set_coupling_map([(0, 1), (1, 2), (2, 3)]);
    for backend in [Backend::all_to_all(3), Backend::linear(4), directed] {
        for level in 0..=3 {
            for elide in [false, true] {
                let config = TranspilerConfig::builder()
                    .backend(backend.clone())
                    .target_basis(["u", "cx"])
                    .optimization_level(level)
                    .routing_search(RoutingSearchConfig {
                        elide_output_permutations: elide,
                        ..Default::default()
                    })
                    .build();
                let plain = transpile(&circuit, Some(config.clone())).unwrap();
                let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
                assert_eq!(plain.operations, output.operations);
                assert!(report.numerical_budget.charged_error <= CIRCUIT_ERROR_ALLOWANCE);
                for op in &output.operations {
                    if let Operation::Gate {
                        name: GateType::CX,
                        qubits,
                        ..
                    } = op
                    {
                        assert!(backend.has_directed_edge(qubits[0], qubits[1]));
                    }
                }
                let error = isometry_error(
                    &circuit,
                    &output,
                    &report.initial_layout.as_ref().unwrap()[..3],
                    &report.final_layout.as_ref().unwrap()[..3],
                );
                assert!(
                    error < 2e-8,
                    "level {level}, elide {elide}, routed isometry residual {error:e}"
                );
            }
        }
    }
}

#[test]
fn encoded_swaps_keep_terminal_readout_labels_and_repeated_write_order() {
    let mut circuit = source(
        3,
        &format!(
            "h q[0]; cx q[0],q[1]; {} {}",
            encoded_swap(0, 2),
            encoded_swap(2, 1)
        ),
    );
    circuit.num_cbits = 2;
    circuit.operations.extend([
        Operation::Measure { qubit: 0, cbit: 1 },
        Operation::Measure { qubit: 2, cbit: 1 },
        Operation::Measure { qubit: 1, cbit: 0 },
    ]);
    for elide in [false, true] {
        let config = TranspilerConfig::builder()
            .backend(Backend::linear(4))
            .target_basis(["u", "cx"])
            .optimization_level(3)
            .routing_search(RoutingSearchConfig {
                elide_output_permutations: elide,
                ..Default::default()
            })
            .build();
        let (mut output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
        let map = report.final_layout.as_ref().unwrap();
        let measurements = output
            .operations
            .iter()
            .filter_map(|op| match op {
                Operation::Measure { qubit, cbit } => Some((*qubit, *cbit)),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(measurements, vec![(map[0], 1), (map[2], 1), (map[1], 0)]);
        let mut unitary = circuit.clone();
        unitary
            .operations
            .retain(|op| !matches!(op, Operation::Measure { .. }));
        output
            .operations
            .retain(|op| !matches!(op, Operation::Measure { .. }));
        let error = isometry_error(
            &unitary,
            &output,
            &report.initial_layout.as_ref().unwrap()[..3],
            &map[..3],
        );
        assert!(error < 2e-8);
    }
}

#[test]
fn optimization_cannot_erase_unsupported_dynamics_before_validation() {
    for body in [
        "cx q[0],q[1]; measure q[0] -> c[0]; cx q[0],q[1];",
        "cx q[0],q[1]; reset q[0]; cx q[0],q[1];",
        "cx q[0],q[1]; if(c==1) x q[1]; cx q[0],q[1];",
    ] {
        let circuit = parse_qasm(&format!(
            "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[2]; creg c[1]; {body}"
        ))
        .unwrap();
        let config = TranspilerConfig::builder()
            .optimization_level(3)
            .backend(Backend::all_to_all(2))
            .build();
        assert!(transpile(&circuit, Some(config.clone())).is_err());
        assert!(transpile_with_report(&circuit, Some(config)).is_err());
    }
}
