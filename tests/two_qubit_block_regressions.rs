//! Independent analytic-oracle regressions for two-qubit block synthesis.
//!
//! The oracle deliberately does not call GateDefinition, decompose(), the
//! production simulator, or either synthesizer. Bit-indexed amplitudes and
//! textbook gate matrices give an independent operand/phase convention.
use nalgebra::DMatrix;
use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, ClassicalCondition, GateType, Operation};
use q_rust::transpiler::block_synthesis::{TwoQubitBlockStats, TwoQubitBlockSynthesisPass};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::synthesis::{kak::KakSynthesizer, Synthesizer};
use q_rust::transpiler::{transpile_with_report, TranspilerConfig};
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI};

fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    }
}

fn real(x: f64) -> C {
    C::new(x, 0.0)
}

fn phase(x: f64) -> C {
    C::from_polar(1.0, x)
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
    assert_eq!(state.len(), 1 << circuit.num_qubits);
    let z = real(0.0);
    let o = real(1.0);
    let i = C::new(0.0, 1.0);
    for operation in &circuit.operations {
        let Operation::Gate {
            name,
            qubits,
            params,
        } = operation
        else {
            assert!(
                matches!(operation, Operation::Barrier { .. }),
                "nonunitary oracle input"
            );
            continue;
        };
        assert!(params.iter().all(|p| p.is_finite()));
        let q = qubits[0];
        match name {
            GateType::CX | GateType::CCX => {
                let target = *qubits.last().unwrap();
                for low in 0..state.len() {
                    if low & (1 << target) == 0
                        && qubits[..qubits.len() - 1]
                            .iter()
                            .all(|c| low & (1 << c) != 0)
                    {
                        state.swap(low, low | (1 << target));
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
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) == 0 && basis & (1 << qubits[1]) != 0 {
                        state.swap(basis, basis ^ (1 << qubits[0]) ^ (1 << qubits[1]));
                    }
                }
            }
            GateType::RXX | GateType::RYY | GateType::RZZ => {
                let before = state.clone();
                let (s, c) = (params[0] / 2.0).sin_cos();
                state.fill(z);
                for (basis, amplitude) in before.into_iter().enumerate() {
                    state[basis] += c * amplitude;
                    let mut mapped = basis;
                    let mut pauli_phase = o;
                    for &wire in qubits {
                        match name {
                            GateType::RXX => mapped ^= 1 << wire,
                            GateType::RYY => {
                                pauli_phase *= if basis & (1 << wire) == 0 { i } else { -i };
                                mapped ^= 1 << wire;
                            }
                            GateType::RZZ => {
                                if basis & (1 << wire) != 0 {
                                    pauli_phase = -pauli_phase;
                                }
                            }
                            _ => unreachable!(),
                        }
                    }
                    state[mapped] += -i * s * pauli_phase * amplitude;
                }
            }
            _ => {
                let (target, control) = if matches!(
                    name,
                    GateType::CRX | GateType::CRY | GateType::CRZ | GateType::CH | GateType::CY
                ) {
                    (qubits[1], Some(q))
                } else {
                    (q, None)
                };
                let matrix = match name {
                    GateType::ID => [[o, z], [z, o]],
                    GateType::H | GateType::CH => [
                        [real(FRAC_1_SQRT_2), real(FRAC_1_SQRT_2)],
                        [real(FRAC_1_SQRT_2), real(-FRAC_1_SQRT_2)],
                    ],
                    GateType::X => [[z, o], [o, z]],
                    GateType::Y | GateType::CY => [[z, -i], [i, z]],
                    GateType::Z => [[o, z], [z, -o]],
                    GateType::SX => [
                        [(o + i) / 2.0, (o - i) / 2.0],
                        [(o - i) / 2.0, (o + i) / 2.0],
                    ],
                    GateType::S => [[o, z], [z, i]],
                    GateType::Sdg => [[o, z], [z, -i]],
                    GateType::T => [[o, z], [z, phase(PI / 4.0)]],
                    GateType::Tdg => [[o, z], [z, phase(-PI / 4.0)]],
                    GateType::RX | GateType::CRX => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [[real(c), -i * s], [-i * s, real(c)]]
                    }
                    GateType::RY | GateType::CRY => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [[real(c), real(-s)], [real(s), real(c)]]
                    }
                    GateType::RZ | GateType::CRZ => {
                        [[phase(-params[0] / 2.0), z], [z, phase(params[0] / 2.0)]]
                    }
                    GateType::U => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [real(c), -phase(params[2]) * s],
                            [phase(params[1]) * s, phase(params[1] + params[2]) * c],
                        ]
                    }
                    _ => panic!("unsupported independent oracle gate {name:?}"),
                };
                one_qubit(&mut state, target, control, matrix);
            }
        }
    }
    state
}

fn unitary(circuit: &Circuit) -> DMatrix<C> {
    let dimension = 1 << circuit.num_qubits;
    let columns: Vec<_> = (0..dimension)
        .map(|column| {
            let mut input = vec![real(0.0); dimension];
            input[column] = real(1.0);
            evolve(circuit, input)
        })
        .collect();
    DMatrix::from_fn(dimension, dimension, |r, c| columns[c][r])
}

/// One phase for the entire operator, never a separate phase per input.
fn residual(expected: &DMatrix<C>, actual: &DMatrix<C>) -> f64 {
    assert_eq!(expected.shape(), actual.shape());
    assert!(actual.iter().all(|x| x.re.is_finite() && x.im.is_finite()));
    let overlap: C = expected
        .iter()
        .zip(actual.iter())
        .map(|(e, a)| e.conj() * a)
        .sum();
    let factor = if overlap.norm() > 1e-15 {
        overlap / overlap.norm()
    } else {
        real(1.0)
    };
    (0..expected.ncols())
        .map(|column| {
            (0..expected.nrows())
                .map(|row| (actual[(row, column)] - factor * expected[(row, column)]).norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max)
}

fn assert_equivalent(source: &Circuit, actual: &Circuit, tolerance: f64) {
    let error = residual(&unitary(source), &unitary(actual));
    assert!(
        error <= tolerance,
        "shared-phase operator residual {error:e} > {tolerance:e}"
    );
}

fn assert_full_isometry(
    source: &Circuit,
    actual: &Circuit,
    initial: &[usize],
    final_layout: &[usize],
) {
    let logical = unitary(source);
    let physical = unitary(actual);
    let embed = |basis: usize, layout: &[usize]| -> usize {
        layout
            .iter()
            .enumerate()
            .map(|(logical, physical)| ((basis >> logical) & 1) << physical)
            .sum()
    };
    let mut expected = DMatrix::zeros(1 << actual.num_qubits, 1 << source.num_qubits);
    let mut observed = expected.clone();
    for column in 0..logical.ncols() {
        for row in 0..logical.nrows() {
            expected[(embed(row, final_layout), column)] = logical[(row, column)];
        }
        for row in 0..physical.nrows() {
            observed[(row, column)] = physical[(row, embed(column, initial))];
        }
    }
    let error = residual(&expected, &observed);
    assert!(error <= 1e-8, "full physical isometry residual {error:e}");
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

fn transfer(circuit: &mut Circuit, control: usize, target: usize, angle: f64) {
    circuit.add_op(gate(GateType::RY, &[target], &[angle / 2.0]));
    circuit.add_op(gate(GateType::CX, &[control, target], &[]));
    circuit.add_op(gate(GateType::RY, &[target], &[-angle / 2.0]));
    circuit.add_op(gate(GateType::CX, &[control, target], &[]));
    circuit.add_op(gate(GateType::CX, &[target, control], &[]));
}

fn local_wrappers(circuit: &mut Circuit, seed: u64) {
    let f = seed as f64;
    circuit.add_op(gate(
        GateType::U,
        &[0],
        &[0.17 + f * 0.031, -0.23 + f * 0.017, 0.61 - f * 0.011],
    ));
    circuit.add_op(gate(
        GateType::U,
        &[1],
        &[-0.39 + f * 0.019, 0.31 - f * 0.013, -0.79 + f * 0.023],
    ));
}

#[test]
fn oracle_detects_relative_phase_even_when_basis_probabilities_match() {
    let source = Circuit::new(2, 0);
    let mut phase_error = source.clone();
    phase_error.add_op(gate(GateType::Z, &[0], &[]));
    assert!(residual(&unitary(&source), &unitary(&phase_error)) > 1.0);
    let superposition = vec![real(0.5); 4];
    assert_ne!(
        evolve(&source, superposition.clone()),
        evolve(&phase_error, superposition)
    );
}

#[test]
#[should_panic(expected = "full physical isometry residual")]
fn oracle_detects_dirty_ancillas_without_projecting_or_renormalizing() {
    let source = Circuit::new(1, 0);
    let mut dirty = Circuit::new(2, 0);
    dirty.add_op(gate(GateType::H, &[1], &[]));
    assert_full_isometry(&source, &dirty, &[0], &[0]);
}

#[test]
fn kak_emits_zero_one_two_or_three_cx_for_known_local_equivalence_classes() {
    for (class, bound) in [(0, 0), (1, 1), (2, 2), (3, 3)] {
        for seed in 0..8 {
            let mut source = Circuit::new(2, 0);
            local_wrappers(&mut source, seed);
            match class {
                0 => {}
                1 => source.add_op(gate(GateType::CX, &[1, 0], &[])),
                2 => transfer(&mut source, 1, 0, 0.31 + seed as f64 * 0.19),
                3 => {
                    source.add_op(gate(GateType::RXX, &[0, 1], &[0.71]));
                    source.add_op(gate(GateType::RYY, &[1, 0], &[0.43]));
                    source.add_op(gate(GateType::RZZ, &[0, 1], &[0.29]));
                }
                _ => unreachable!(),
            }
            local_wrappers(&mut source, 13 + seed);
            let input = unitary(&source) * phase(0.13 * seed as f64);
            let output = KakSynthesizer
                .synthesize(&input, &[])
                .expect("well-conditioned 2q unitary");
            assert!(cx_count(&output) <= bound, "class {class}, seed {seed}");
            assert!(
                residual(&input, &unitary(&output)) <= 1e-9,
                "class {class}, seed {seed}"
            );
        }
    }
}

#[test]
fn kak_preserves_seeded_generic_operators_and_operand_directions() {
    let mut rng = 0xc0ffee1234567890u64;
    for sample in 0..24 {
        let mut source = Circuit::new(2, 0);
        for layer in 0..6 {
            for q in 0..2 {
                let mut angles = [0.0; 3];
                for angle in &mut angles {
                    rng = rng
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    *angle = ((rng >> 11) as f64 / ((1u64 << 53) as f64) - 0.5) * 2.0 * PI;
                }
                source.add_op(gate(GateType::U, &[q], &angles));
            }
            source.add_op(gate(GateType::CX, &[layer % 2, 1 - layer % 2], &[]));
        }
        let output = KakSynthesizer
            .synthesize(&unitary(&source), &[])
            .expect("seeded generic unitary");
        assert!(cx_count(&output) <= 3, "sample {sample}");
        assert_equivalent(&source, &output, 1e-9);
    }
}

#[test]
fn kak_rejects_nonfinite_nonunitary_and_wrong_shape_inputs() {
    for size in [0, 2, 3, 8] {
        assert!(KakSynthesizer
            .synthesize(&DMatrix::identity(size, size), &[])
            .is_none());
    }
    assert!(KakSynthesizer
        .synthesize(&DMatrix::zeros(4, 3), &[])
        .is_none());
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.9, 1.0 + 5e-9] {
        let mut input = DMatrix::identity(4, 4);
        input[(0, 0)] = C::new(invalid, 0.0);
        assert!(KakSynthesizer.synthesize(&input, &[]).is_none());
    }
    let mut imaginary_nan = DMatrix::identity(4, 4);
    imaginary_nan[(1, 2)] = C::new(0.0, f64::NAN);
    assert!(KakSynthesizer.synthesize(&imaginary_nan, &[]).is_none());
}

#[test]
fn kak_degenerate_and_tiny_interactions_are_retained_or_safely_rejected() {
    for angle in [0.0, 1e-12, 2e-8, 9e-8, 1e-6, FRAC_PI_2, PI, PI + 1e-9] {
        let mut source = Circuit::new(2, 0);
        source.add_op(gate(GateType::RZZ, &[1, 0], &[angle]));
        local_wrappers(&mut source, 7);
        if let Some(output) = KakSynthesizer.synthesize(&unitary(&source), &[]) {
            assert_equivalent(&source, &output, 1e-10);
        }
    }
}

fn run_blocks(source: &Circuit) -> Circuit {
    TwoQubitBlockSynthesisPass.run(source, &mut PropertySet::new())
}

/// Exact RY=U(theta,0,0) input preparation, independent of production lowering.
fn transfer_ucx(circuit: &mut Circuit, control: usize, target: usize, angle: f64) {
    circuit.add_op(gate(GateType::U, &[target], &[angle / 2.0, 0.0, 0.0]));
    circuit.add_op(gate(GateType::CX, &[control, target], &[]));
    circuit.add_op(gate(GateType::U, &[target], &[-angle / 2.0, 0.0, 0.0]));
    circuit.add_op(gate(GateType::CX, &[control, target], &[]));
    circuit.add_op(gate(GateType::CX, &[target, control], &[]));
}

fn w_circuit(width: usize) -> Circuit {
    let mut source = Circuit::new(width, 0);
    source.add_op(gate(GateType::X, &[0], &[]));
    for target in 1..width {
        let angle = 2.0 * (1.0 / ((width - target + 1) as f64).sqrt()).acos();
        transfer(&mut source, target - 1, target, angle);
    }
    source
}

fn qft(width: usize) -> Circuit {
    let mut source = Circuit::new(width, 0);
    for target in (0..width).rev() {
        source.add_op(gate(GateType::H, &[target], &[]));
        for control in (0..target).rev() {
            let angle = PI / (1 << (target - control)) as f64;
            source.add_op(gate(GateType::U, &[control], &[0.0, 0.0, angle / 2.0]));
            source.add_op(gate(GateType::CRZ, &[control, target], &[angle]));
        }
    }
    for q in 0..width / 2 {
        source.add_op(gate(GateType::SWAP, &[q, width - 1 - q], &[]));
    }
    source
}

#[test]
fn w_transfer_reduces_three_cx_to_two_for_every_input_and_reversed_nonadjacent_wires() {
    for pair in [[0, 2], [2, 0]] {
        for angle in [0.37, 1.234, 2.0 * (1.0 / 3.0_f64.sqrt()).acos()] {
            let mut source = Circuit::new(4, 0);
            transfer_ucx(&mut source, pair[0], pair[1], angle);
            let mut properties = PropertySet::new();
            let output = TwoQubitBlockSynthesisPass.run(&source, &mut properties);
            assert_eq!(cx_count(&source), 3);
            assert_eq!(
                cx_count(&output),
                2,
                "angle={angle}, pair={pair:?}, stats={:?}",
                properties.get::<TwoQubitBlockStats>("two_qubit_block_stats")
            );
            assert_equivalent(&source, &output, 1e-9);
            // The untouched qubits are arbitrary inputs, not assumed clean.
            assert!(output
                .operations
                .iter()
                .all(|op| op.qubits().iter().all(|q| pair.contains(q))));
        }
    }
}

#[test]
fn collection_may_cross_disjoint_operations_but_preserves_their_order_and_unitary() {
    let mut source = Circuit::new(4, 0);
    transfer_ucx(&mut source, 2, 0, 0.73);
    source
        .operations
        .insert(2, gate(GateType::U, &[1], &[0.51, -0.42, 0.17]));
    source
        .operations
        .insert(4, gate(GateType::CX, &[1, 3], &[]));
    source.add_op(gate(GateType::U, &[3], &[-0.38, 0.29, 0.61]));
    transfer_ucx(&mut source, 0, 2, 1.07);
    let output = run_blocks(&source);
    assert!(cx_count(&output) < cx_count(&source));
    assert_equivalent(&source, &output, 1e-9);
    let disjoint = |c: &Circuit| {
        c.operations
            .iter()
            .filter(|op| op.qubits().iter().all(|q| [1, 3].contains(q)))
            .cloned()
            .collect::<Vec<_>>()
    };
    assert_eq!(disjoint(&source), disjoint(&output));
}

#[test]
fn a_third_wire_interaction_terminates_both_dependency_chains() {
    for seed in 0..12 {
        let mut source = Circuit::new(3, 0);
        transfer_ucx(&mut source, 0, 2, 0.4 + seed as f64 * 0.13);
        source
            .operations
            .insert(3, gate(GateType::CX, &[2, 1], &[]));
        source
            .operations
            .insert(4, gate(GateType::U, &[1], &[0.41, 0.57, -0.61]));
        source.add_op(gate(GateType::CX, &[1, 0], &[]));
        transfer_ucx(&mut source, 2, 0, 1.2 - seed as f64 * 0.07);
        let output = run_blocks(&source);
        assert!(cx_count(&output) <= cx_count(&source));
        assert_equivalent(&source, &output, 1e-9);
    }
}

#[test]
fn barriers_readout_reset_conditionals_and_unknown_gates_are_global_collection_boundaries() {
    // Without the middle boundary these gates cancel completely. Crossing any
    // boundary below would change the execution contract even if the unitary
    // prefix/suffix happened to simplify together.
    let boundaries = [
        Operation::Barrier { qubits: vec![2] },
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Reset { qubit: 0 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".to_owned(),
                value: 1,
            },
            op: Box::new(gate(GateType::X, &[2], &[])),
        },
        gate(GateType::Custom("opaque".to_owned()), &[2], &[]),
        gate(GateType::CCX, &[0, 1, 2], &[]),
    ];
    for boundary in boundaries {
        let mut source = Circuit::new(3, 1);
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(gate(GateType::U, &[1], &[0.29, 0.0, 0.0]));
        source.add_op(boundary);
        source.add_op(gate(GateType::U, &[1], &[-0.29, 0.0, 0.0]));
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        assert_eq!(run_blocks(&source).operations, source.operations);
    }
}

#[test]
fn pipeline_w_state_improvement_preserves_arbitrary_inputs_not_only_preparation() {
    for width in [3, 4] {
        let source = w_circuit(width);
        let config = |enabled| {
            TranspilerConfig::builder()
                .optimization_level(3)
                .backend(Backend::all_to_all(width))
                .target_basis(["u", "cx"])
                .two_qubit_block_synthesis(enabled)
                .build()
        };
        let (before, before_report) = transpile_with_report(&source, Some(config(false))).unwrap();
        let (after, after_report) = transpile_with_report(&source, Some(config(true))).unwrap();
        assert_eq!(cx_count(&before), 3 * (width - 1));
        assert_eq!(cx_count(&after), 2 * (width - 1));
        for (output, report) in [(&before, &before_report), (&after, &after_report)] {
            assert_full_isometry(
                &source,
                output,
                report.initial_layout.as_ref().unwrap(),
                report.final_layout.as_ref().unwrap(),
            );
        }
    }
}

#[test]
fn qft_pipeline_preserves_fourier_columns_with_clean_ancillas_on_a_directed_line() {
    for width in [2, 3, 4] {
        let source = qft(width);
        let exact = DMatrix::from_fn(1 << width, 1 << width, |r, c| {
            phase(2.0 * PI * (r * c) as f64 / (1 << width) as f64) / ((1 << width) as f64).sqrt()
        });
        assert!(
            residual(&exact, &unitary(&source)) < 1e-12,
            "fixture itself must equal Fourier matrix"
        );
        let mut backend = Backend::new("one_way_line", width + 1);
        backend.set_coupling_map((0..width).map(|q| (q + 1, q)));
        let config = TranspilerConfig::builder()
            .optimization_level(3)
            .backend(backend.clone())
            .target_basis(["u", "cx"])
            .two_qubit_block_synthesis(true)
            .build();
        let (output, report) = transpile_with_report(&source, Some(config)).unwrap();
        assert_full_isometry(
            &source,
            &output,
            report.initial_layout.as_ref().unwrap(),
            report.final_layout.as_ref().unwrap(),
        );
        for op in &output.operations {
            match op {
                Operation::Gate {
                    name: GateType::CX,
                    qubits,
                    ..
                } => assert!(backend.has_directed_edge(qubits[0], qubits[1])),
                Operation::Gate {
                    name: GateType::U, ..
                } => {}
                _ => panic!("output is not in the requested U/CX basis"),
            }
        }
    }
}

#[test]
fn terminal_readout_keeps_repeated_classical_writes_in_order_after_block_synthesis() {
    let mut source = w_circuit(3);
    source.num_cbits = 2;
    source.add_op(gate(GateType::SWAP, &[0, 2], &[]));
    source.add_op(gate(GateType::SWAP, &[1, 2], &[]));
    let prefix = source.clone();
    let readout = [(2, 0), (0, 0), (1, 1)];
    for (qubit, cbit) in readout {
        source.add_op(Operation::Measure { qubit, cbit });
    }
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(4))
        .target_basis(["u", "cx"])
        .two_qubit_block_synthesis(true)
        .build();
    let (mut output, report) = transpile_with_report(&source, Some(config)).unwrap();
    let final_layout = report.final_layout.as_ref().unwrap();
    let expected: Vec<_> = readout
        .iter()
        .map(|&(q, cbit)| Operation::Measure {
            qubit: final_layout[q],
            cbit,
        })
        .collect();
    assert_eq!(
        &output.operations[output.operations.len() - 3..],
        expected.as_slice()
    );
    output.operations.truncate(output.operations.len() - 3);
    assert_full_isometry(
        &prefix,
        &output,
        report.initial_layout.as_ref().unwrap(),
        final_layout,
    );
}

#[test]
fn many_small_coherent_interactions_cannot_be_erased_by_a_local_fidelity_threshold() {
    let mut source = Circuit::new(2, 0);
    for _ in 0..96 {
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(gate(GateType::U, &[1], &[0.0, 0.0, 8e-8]));
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(Operation::Barrier { qubits: vec![0, 1] });
    }
    let output = run_blocks(&source);
    assert_equivalent(&source, &output, 1e-9);
    assert_eq!(
        cx_count(&output),
        cx_count(&source),
        "small but nonzero interactions are not local unitaries"
    );
    assert_eq!(
        output
            .operations
            .iter()
            .filter(|op| op.is_barrier())
            .count(),
        96
    );
}

#[test]
fn repeated_accepted_blocks_obey_a_whole_circuit_coherent_error_budget() {
    let mut source = Circuit::new(2, 0);
    for index in 0..64 {
        transfer_ucx(
            &mut source,
            index % 2,
            1 - index % 2,
            0.31 + index as f64 * 0.017,
        );
        source.add_op(Operation::Barrier { qubits: vec![0, 1] });
    }
    let output = run_blocks(&source);
    assert!(cx_count(&output) < cx_count(&source));
    assert_equivalent(&source, &output, 1e-9);
    assert_eq!(
        output
            .operations
            .iter()
            .filter(|op| op.is_barrier())
            .count(),
        64
    );
}

#[test]
fn invalid_block_parameters_are_retained_or_rejected_without_panicking() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut source = Circuit::new(2, 0);
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(gate(GateType::U, &[1], &[value, 0.1, 0.2]));
        source.add_op(gate(GateType::CX, &[1, 0], &[]));
        let pass_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            TwoQubitBlockSynthesisPass.try_run(&source, &mut PropertySet::new())
        }))
        .expect("invalid numerical input must not panic");
        if let Ok(retained) = pass_result {
            // NaN is not equal to itself, hence compare its preserved encoding.
            assert_eq!(source.to_qasm(None), retained.to_qasm(None));
        }
        let config = TranspilerConfig::builder()
            .optimization_level(3)
            .two_qubit_block_synthesis(true)
            .build();
        assert!(transpile_with_report(&source, Some(config)).is_err());
    }
}

#[test]
fn finite_parameters_with_unstable_matrix_arithmetic_retain_the_original_block() {
    let mut source = Circuit::new(2, 0);
    source.add_op(gate(GateType::CX, &[0, 1], &[]));
    // Each parameter is finite, but adding the two phases overflows. Numerical
    // failure must not silently remove the block or emit nonfinite replacements.
    source.add_op(gate(GateType::U, &[1], &[0.33, f64::MAX, f64::MAX]));
    source.add_op(gate(GateType::CX, &[1, 0], &[]));
    let output = run_blocks(&source);
    assert_eq!(source.operations, output.operations);
}

#[test]
fn malformed_lowered_operations_are_safe_global_boundaries() {
    for malformed in [
        gate(GateType::U, &[], &[0.1, 0.2, 0.3]),
        gate(GateType::U, &[0], &[0.1, 0.2]),
        gate(GateType::U, &[0, 1], &[0.1, 0.2, 0.3]),
        gate(GateType::CX, &[0, 0], &[]),
        gate(GateType::CX, &[0, 2], &[]),
        gate(GateType::CX, &[0, 1], &[0.5]),
    ] {
        let mut source = Circuit::new(2, 0);
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(malformed);
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        assert_eq!(run_blocks(&source).operations, source.operations);
    }
}

#[test]
fn pass_metadata_is_per_run_and_cannot_leak_an_old_residual_budget() {
    let mut source = Circuit::new(2, 0);
    transfer_ucx(&mut source, 0, 1, 1.23);
    let mut properties = PropertySet::new();
    properties.insert("initial_layout", vec![1usize, 0]);
    properties.insert("final_layout", vec![0usize, 1]);
    let output = TwoQubitBlockSynthesisPass.run(&source, &mut properties);
    assert_equivalent(&source, &output, 1e-9);
    let stats = properties
        .get::<TwoQubitBlockStats>("two_qubit_block_stats")
        .unwrap();
    assert_eq!(stats.accepted, 1);
    assert_eq!(stats.cx_removed, 1);
    assert!(stats.max_residual <= 1e-10 && stats.residual_sum <= 1e-9);
    let empty = Circuit::new(2, 0);
    let rerun = TwoQubitBlockSynthesisPass.run(&empty, &mut properties);
    assert!(rerun.operations.is_empty());
    assert_eq!(
        properties.get::<TwoQubitBlockStats>("two_qubit_block_stats"),
        Some(&TwoQubitBlockStats::default())
    );
    assert_eq!(
        properties.get::<Vec<usize>>("initial_layout"),
        Some(&vec![1, 0])
    );
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&vec![0, 1])
    );
}

#[test]
fn absorbed_leading_locals_must_not_move_replacement_before_an_intervening_interaction() {
    let mut source = Circuit::new(3, 0);
    source.add_op(gate(GateType::U, &[0], &[0.43, -0.29, 0.61]));
    let intervening = gate(GateType::CX, &[1, 2], &[]);
    source.add_op(intervening.clone());
    transfer_ucx(&mut source, 0, 1, 1.19);
    let output = run_blocks(&source);
    assert_eq!(cx_count(&output), 3);
    assert_equivalent(&source, &output, 1e-9);
    // The first leading U acts only on q0 and may move right across CX(1,2).
    // The replacement acts on q1 too and must remain after CX(1,2).
    assert_eq!(output.operations.first(), Some(&intervening));
}

#[test]
fn predecessor_absorption_cannot_reuse_gates_owned_by_an_earlier_block() {
    for first_block_is_reducible in [false, true] {
        let mut source = Circuit::new(3, 0);
        if first_block_is_reducible {
            transfer_ucx(&mut source, 0, 1, 0.73);
        } else {
            // A generic two-CX class has no strict CX improvement. The trailing
            // U is still owned by this processed, retained block.
            source.add_op(gate(GateType::CX, &[0, 1], &[]));
            source.add_op(gate(GateType::U, &[0], &[0.37, -0.29, 0.17]));
            source.add_op(gate(GateType::U, &[1], &[0.91, 0.23, -0.51]));
            source.add_op(gate(GateType::CX, &[1, 0], &[]));
        }
        source.add_op(gate(GateType::U, &[1], &[0.57, 0.12, -0.31]));
        let prefix_length = source.operations.len();
        transfer_ucx(&mut source, 1, 2, 1.07);
        let output = run_blocks(&source);
        assert_eq!(cx_count(&output), 4);
        assert_equivalent(&source, &output, 1e-9);
        if !first_block_is_reducible {
            assert_eq!(
                output.operations[..prefix_length],
                source.operations[..prefix_length]
            );
        }
    }
}

#[test]
fn a_global_boundary_also_stops_absorbing_predecessor_locals() {
    let mut source = Circuit::new(3, 0);
    source.add_op(gate(GateType::U, &[0], &[0.27, -0.39, 0.11]));
    source.add_op(gate(GateType::U, &[1], &[-0.37, 0.19, 0.71]));
    source.add_op(Operation::Barrier { qubits: vec![2] });
    transfer_ucx(&mut source, 0, 1, 0.93);
    let output = run_blocks(&source);
    assert_eq!(cx_count(&output), 2);
    assert_eq!(output.operations[..3], source.operations[..3]);
    assert_equivalent(&source, &output, 1e-9);
}
