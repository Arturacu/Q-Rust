//! Full-operator checks with independent analytic matrices. No production
//! gate definition, decomposition or simulator is used as the expected oracle.
use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, ClassicalCondition, GateType, Operation};
use q_rust::transpiler::local_cleanup::{PostBlockFusionPass, PostBlockFusionStats};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::phase_gadgets::RecoverPhaseGadgetsPass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::BeamSabrePass;
use q_rust::transpiler::{transpile_with_report, RoutingSearchConfig, TranspilerConfig};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

fn gate(name: GateType, wires: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: wires.to_vec(),
        params: params.to_vec(),
    }
}

fn phase(angle: f64) -> C {
    C::from_polar(1.0, angle)
}

fn evolve(circuit: &Circuit, mut state: Vec<C>) -> Vec<C> {
    let zero = C::new(0.0, 0.0);
    let one = C::new(1.0, 0.0);
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
        assert!(params.iter().all(|p| p.is_finite()));
        match name {
            GateType::CX => {
                for index in 0..state.len() {
                    if index & (1 << qubits[0]) != 0 && index & (1 << qubits[1]) == 0 {
                        state.swap(index, index | (1 << qubits[1]));
                    }
                }
            }
            GateType::SWAP => {
                for index in 0..state.len() {
                    if index & (1 << qubits[0]) == 0 && index & (1 << qubits[1]) != 0 {
                        state.swap(index, index ^ (1 << qubits[0]) ^ (1 << qubits[1]));
                    }
                }
            }
            GateType::RZZ => {
                for (index, value) in state.iter_mut().enumerate() {
                    let different = ((index >> qubits[0]) ^ (index >> qubits[1])) & 1;
                    *value *= phase(if different == 0 {
                        -params[0] / 2.0
                    } else {
                        params[0] / 2.0
                    });
                }
            }
            GateType::CZ => {
                for (index, value) in state.iter_mut().enumerate() {
                    if qubits.iter().all(|q| index & (1 << q) != 0) {
                        *value = -*value;
                    }
                }
            }
            _ => {
                let (wire, control) = if matches!(name, GateType::CRZ | GateType::CRY) {
                    (qubits[1], Some(qubits[0]))
                } else {
                    (qubits[0], None)
                };
                let matrix = match name {
                    GateType::U => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [C::new(c, 0.0), -phase(params[2]) * s],
                            [phase(params[1]) * s, phase(params[1] + params[2]) * c],
                        ]
                    }
                    GateType::H => [
                        [C::new(FRAC_1_SQRT_2, 0.0), C::new(FRAC_1_SQRT_2, 0.0)],
                        [C::new(FRAC_1_SQRT_2, 0.0), C::new(-FRAC_1_SQRT_2, 0.0)],
                    ],
                    GateType::X => [[zero, one], [one, zero]],
                    GateType::Z => [[one, zero], [zero, -one]],
                    GateType::RX => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [C::new(c, 0.0), C::new(0.0, -s)],
                            [C::new(0.0, -s), C::new(c, 0.0)],
                        ]
                    }
                    GateType::RY | GateType::CRY => {
                        let (s, c) = (params[0] / 2.0).sin_cos();
                        [
                            [C::new(c, 0.0), C::new(-s, 0.0)],
                            [C::new(s, 0.0), C::new(c, 0.0)],
                        ]
                    }
                    GateType::RZ | GateType::CRZ => [
                        [phase(-params[0] / 2.0), zero],
                        [zero, phase(params[0] / 2.0)],
                    ],
                    _ => panic!("unsupported independent oracle gate {name:?}"),
                };
                for low in 0..state.len() {
                    if low & (1 << wire) != 0 || control.is_some_and(|q| low & (1 << q) == 0) {
                        continue;
                    }
                    let high = low | (1 << wire);
                    let (a, b) = (state[low], state[high]);
                    state[low] = matrix[0][0] * a + matrix[0][1] * b;
                    state[high] = matrix[1][0] * a + matrix[1][1] * b;
                }
            }
        }
    }
    state
}

fn columns(circuit: &Circuit) -> Vec<Vec<C>> {
    let size = 1 << circuit.num_qubits;
    (0..size)
        .map(|column| {
            let mut input = vec![C::new(0.0, 0.0); size];
            input[column] = C::new(1.0, 0.0);
            evolve(circuit, input)
        })
        .collect()
}

fn error(expected: &[Vec<C>], actual: &[Vec<C>]) -> f64 {
    assert_eq!(expected.len(), actual.len());
    assert!(actual
        .iter()
        .flatten()
        .all(|z| z.re.is_finite() && z.im.is_finite()));
    let overlap: C = expected
        .iter()
        .flatten()
        .zip(actual.iter().flatten())
        .map(|(a, b)| a.conj() * b)
        .sum();
    let shared_phase = if overlap.norm() > 1e-15 {
        overlap / overlap.norm()
    } else {
        C::new(1.0, 0.0)
    };
    expected
        .iter()
        .zip(actual)
        .map(|(a, b)| {
            assert_eq!(a.len(), b.len());
            a.iter()
                .zip(b)
                .map(|(a, b)| (*b - shared_phase * a).norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max)
}

fn equivalent(source: &Circuit, output: &Circuit, tolerance: f64) {
    let residual = error(&columns(source), &columns(output));
    assert!(
        residual <= tolerance,
        "shared-phase full-operator residual {residual:e}"
    );
}

fn full_isometry(source: &Circuit, output: &Circuit, initial: &[usize], final_map: &[usize]) {
    let embed = |index: usize, layout: &[usize]| -> usize {
        layout
            .iter()
            .enumerate()
            .map(|(q, p)| ((index >> q) & 1) << p)
            .sum()
    };
    let source_columns = columns(source);
    let physical_size = 1 << output.num_qubits;
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    for (column, source_column) in source_columns.iter().enumerate() {
        let mut wanted = vec![C::new(0.0, 0.0); physical_size];
        for (row, value) in source_column.iter().enumerate() {
            wanted[embed(row, final_map)] = *value;
        }
        expected.push(wanted);
        let mut input = vec![C::new(0.0, 0.0); physical_size];
        input[embed(column, initial)] = C::new(1.0, 0.0);
        actual.push(evolve(output, input));
    }
    let residual = error(&expected, &actual);
    assert!(
        residual < 1e-8,
        "full physical isometry residual {residual:e}"
    );
}

/// Independent unit-cost ASAP depths: the CX-only graph omits all U vertices.
fn metrics(circuit: &Circuit) -> (usize, usize, usize, usize) {
    let (mut gates, mut cx) = (0, 0);
    let mut depth = vec![0; circuit.num_qubits];
    let mut entangling = depth.clone();
    for op in &circuit.operations {
        if let Operation::Gate { name, qubits, .. } = op {
            gates += 1;
            let layer = qubits.iter().map(|q| depth[*q]).max().unwrap() + 1;
            for q in qubits {
                depth[*q] = layer;
            }
            if *name == GateType::CX {
                cx += 1;
                let layer = qubits.iter().map(|q| entangling[*q]).max().unwrap() + 1;
                for q in qubits {
                    entangling[*q] = layer;
                }
            }
        }
    }
    (
        gates,
        cx,
        *depth.iter().max().unwrap(),
        *entangling.iter().max().unwrap(),
    )
}

fn w_circuit(width: usize) -> Circuit {
    let mut source = Circuit::new(width, 0);
    source.add_op(gate(GateType::X, &[0], &[]));
    for target in 1..width {
        let angle = 2.0 * (1.0 / ((width - target + 1) as f64).sqrt()).acos();
        source.add_op(gate(GateType::CRY, &[target - 1, target], &[angle]));
        source.add_op(gate(GateType::CX, &[target, target - 1], &[]));
    }
    source
}

fn qft(width: usize) -> Circuit {
    let mut source = Circuit::new(width, 0);
    for target in (0..width).rev() {
        source.add_op(gate(GateType::H, &[target], &[]));
        for control in (0..target).rev() {
            let angle = PI / (1usize << (target - control)) as f64;
            source.add_op(gate(GateType::U, &[control], &[0.0, 0.0, angle / 2.0]));
            // Analytic CRZ expansion followed by the control phase above is
            // CP, not CRZ alone. This matches the common normalized input.
            source.add_op(gate(GateType::RZ, &[target], &[angle / 2.0]));
            source.add_op(gate(GateType::CX, &[control, target], &[]));
            source.add_op(gate(GateType::RZ, &[target], &[-angle / 2.0]));
            source.add_op(gate(GateType::CX, &[control, target], &[]));
        }
    }
    for q in 0..width / 2 {
        source.add_op(gate(GateType::SWAP, &[q, width - 1 - q], &[]));
    }
    source
}

fn fusion(source: &Circuit) -> Circuit {
    PostBlockFusionPass.run(source, &mut PropertySet::new())
}

#[test]
fn local_fusion_preserves_multiplication_order_with_entangled_inputs_and_disjoint_gates() {
    let mut source = Circuit::new(3, 0);
    source.operations = vec![
        gate(GateType::CX, &[0, 2], &[]),
        gate(GateType::U, &[2], &[0.31, -0.47, 0.91]),
        gate(GateType::U, &[1], &[-0.21, 0.37, -0.11]),
        gate(GateType::CX, &[0, 1], &[]),
        gate(GateType::U, &[2], &[-0.73, 0.19, -0.52]),
        gate(GateType::CX, &[2, 0], &[]),
    ];
    let output = fusion(&source);
    equivalent(&source, &output, 1e-11);
    assert_eq!(metrics(&output).0, metrics(&source).0 - 1);
    assert_eq!(metrics(&output).1, metrics(&source).1);
    let mut wrong = source.clone();
    wrong.operations.swap(1, 4);
    assert!(
        error(&columns(&source), &columns(&wrong)) > 0.1,
        "noncommuting local order must be observable"
    );
}

#[test]
fn local_fusion_preserves_tiny_coherent_rotations_and_both_euler_endpoints() {
    for theta in [0.0, 1e-13, -1e-11, 1e-8, PI - 1e-13, PI, PI + 1e-11] {
        let mut source = Circuit::new(2, 0);
        source.add_op(gate(GateType::CX, &[1, 0], &[]));
        source.add_op(gate(GateType::U, &[0], &[theta, -0.81, 1.23]));
        source.add_op(gate(GateType::U, &[0], &[-theta, 0.39, -0.47]));
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        equivalent(&source, &fusion(&source), 1e-11);
    }
    let mut source = Circuit::new(1, 0);
    for _ in 0..512 {
        source.add_op(gate(GateType::U, &[0], &[8e-9, 0.0, 0.0]));
    }
    let output = fusion(&source);
    assert_eq!(metrics(&output).0, 1);
    equivalent(&source, &output, 1e-11);
    assert!(error(&columns(&source), &columns(&Circuit::new(1, 0))) > 1e-6);
}

#[test]
fn local_fusion_handles_global_scalar_phases_but_cannot_erase_relative_phase() {
    let mut scalar = Circuit::new(2, 0);
    scalar.add_op(gate(GateType::U, &[1], &[2.0 * PI, 0.0, 0.0]));
    scalar.add_op(gate(GateType::U, &[1], &[0.0, 0.0, 0.0]));
    equivalent(&scalar, &fusion(&scalar), 1e-12);
    let mut relative = Circuit::new(2, 0);
    relative.add_op(gate(GateType::U, &[1], &[0.0, 0.0, 0.47]));
    relative.add_op(gate(GateType::U, &[1], &[0.0, 0.0, 0.23]));
    let output = fusion(&relative);
    assert_eq!(metrics(&output).0, 1);
    equivalent(&relative, &output, 1e-12);
    assert!(error(&columns(&relative), &columns(&Circuit::new(2, 0))) > 0.3);
}

#[test]
fn local_fusion_global_boundaries_retain_operation_order_and_classical_writes() {
    for boundary in [
        Operation::Barrier { qubits: vec![2] },
        Operation::Measure { qubit: 2, cbit: 0 },
        Operation::Reset { qubit: 2 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 1,
            },
            op: Box::new(gate(GateType::X, &[2], &[])),
        },
        gate(GateType::Custom("opaque".into()), &[2], &[]),
    ] {
        let mut source = Circuit::new(3, 1);
        source.operations = vec![
            gate(GateType::U, &[0], &[0.4, 0.1, 0.2]),
            boundary,
            gate(GateType::U, &[0], &[-0.4, -0.2, -0.1]),
        ];
        assert_eq!(fusion(&source).operations, source.operations);
    }
    let mut source = Circuit::new(2, 1);
    source.operations = vec![
        gate(GateType::U, &[0], &[0.3, 0.0, 0.0]),
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Measure { qubit: 1, cbit: 0 },
        gate(GateType::U, &[0], &[-0.3, 0.0, 0.0]),
    ];
    assert_eq!(fusion(&source).operations, source.operations);
}

#[test]
fn local_fusion_retains_malformed_nonfinite_and_overflowing_runs_without_panicking() {
    for invalid in [
        gate(GateType::U, &[0], &[f64::NAN, 0.0, 0.0]),
        gate(GateType::U, &[0], &[f64::INFINITY, 0.0, 0.0]),
        gate(GateType::U, &[0], &[0.3, f64::MAX, f64::MAX]),
        gate(GateType::U, &[0], &[0.1, 0.2]),
        gate(GateType::U, &[0, 1], &[0.1, 0.2, 0.3]),
        gate(GateType::U, &[2], &[0.1, 0.2, 0.3]),
        gate(GateType::CX, &[0, 0], &[]),
    ] {
        let mut source = Circuit::new(2, 0);
        source.operations = vec![
            gate(GateType::U, &[0], &[0.3, 0.1, 0.7]),
            invalid,
            gate(GateType::U, &[0], &[-0.2, 0.4, -0.5]),
        ];
        let output = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| fusion(&source)))
            .expect("bad runs must remain safe");
        assert_eq!(output.to_qasm(None), source.to_qasm(None));
    }
}

#[test]
fn local_fusion_resets_its_budget_and_preserves_layout_metadata() {
    let mut source = Circuit::new(2, 0);
    for index in 0..64 {
        source.add_op(gate(GateType::U, &[index % 2], &[1e-8, 0.13, -0.29]));
        source.add_op(gate(GateType::U, &[index % 2], &[-2e-8, -0.47, 0.61]));
        source.add_op(Operation::Barrier { qubits: vec![0, 1] });
    }
    let mut properties = PropertySet::new();
    properties.insert("initial_layout", vec![2usize, 0]);
    properties.insert("final_layout", vec![1usize, 2]);
    let output = PostBlockFusionPass.run(&source, &mut properties);
    equivalent(&source, &output, 1e-10);
    let stats = properties
        .get::<PostBlockFusionStats>("post_block_fusion_stats")
        .unwrap();
    assert_eq!(stats.accepted, 64);
    assert_eq!(stats.u_gates_removed, 64);
    assert!(stats.max_residual <= 1e-12 && stats.residual_sum <= 1e-10);
    PostBlockFusionPass.run(&Circuit::new(2, 0), &mut properties);
    assert_eq!(
        properties.get::<PostBlockFusionStats>("post_block_fusion_stats"),
        Some(&PostBlockFusionStats::default())
    );
    assert_eq!(
        properties.get::<Vec<usize>>("initial_layout"),
        Some(&vec![2, 0])
    );
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&vec![1, 2])
    );
}

#[test]
fn parity_recovery_is_exact_for_reversed_nonadjacent_operands_and_noncanonical_angles() {
    for pair in [[0, 2], [2, 0]] {
        for angle in [-9.1, -PI, 0.0, 1e-12, 0.43, 4.0 * PI + 0.2] {
            let mut source = Circuit::new(4, 0);
            source.operations = vec![
                gate(GateType::CX, &pair, &[]),
                gate(GateType::U, &[1], &[0.39, -0.7, 0.81]),
                gate(GateType::RZ, &[pair[1]], &[angle]),
                gate(GateType::CX, &[1, 3], &[]),
                gate(GateType::CX, &pair, &[]),
            ];
            let mut properties = PropertySet::new();
            let output = RecoverPhaseGadgetsPass
                .try_run(&source, &mut properties)
                .unwrap();
            assert_eq!(properties.get::<usize>("phase_gadgets_recovered"), Some(&1));
            assert_eq!(output.operations[0], gate(GateType::RZZ, &pair, &[angle]));
            assert_eq!(
                output.operations[1..],
                [source.operations[1].clone(), source.operations[3].clone()]
            );
            equivalent(&source, &output, 1e-12);
            RecoverPhaseGadgetsPass
                .try_run(&Circuit::new(4, 0), &mut properties)
                .unwrap();
            assert_eq!(properties.get::<usize>("phase_gadgets_recovered"), Some(&0));
        }
    }
}

#[test]
fn parity_recovery_rejects_wrong_control_reverse_cx_and_shared_wire_fences() {
    for middle in [
        gate(GateType::RZ, &[0], &[0.37]),
        gate(GateType::RX, &[1], &[0.37]),
        gate(GateType::H, &[0], &[]),
        gate(GateType::CX, &[1, 2], &[]),
        gate(GateType::SWAP, &[0, 2], &[]),
    ] {
        let mut source = Circuit::new(3, 0);
        source.operations = vec![
            gate(GateType::CX, &[0, 1], &[]),
            middle,
            gate(GateType::RZ, &[1], &[0.51]),
            gate(GateType::CX, &[0, 1], &[]),
        ];
        let output = RecoverPhaseGadgetsPass
            .try_run(&source, &mut PropertySet::new())
            .unwrap();
        assert_eq!(output.operations, source.operations);
        equivalent(&source, &output, 1e-12);
    }
    let mut source = Circuit::new(2, 0);
    source.operations = vec![
        gate(GateType::CX, &[0, 1], &[]),
        gate(GateType::RZ, &[1], &[0.51]),
        gate(GateType::CX, &[1, 0], &[]),
    ];
    let output = RecoverPhaseGadgetsPass
        .try_run(&source, &mut PropertySet::new())
        .unwrap();
    assert_eq!(output.operations, source.operations);
    let mut wrong = Circuit::new(2, 0);
    wrong.add_op(gate(GateType::RZZ, &[0, 1], &[0.51]));
    assert!(error(&columns(&source), &columns(&wrong)) > 1.0);
}

#[test]
fn recovery_cannot_cross_a_global_barrier_or_readout_and_rejects_malformed_inputs() {
    for boundary in [
        Operation::Barrier { qubits: vec![2] },
        Operation::Measure { qubit: 2, cbit: 0 },
    ] {
        let mut source = Circuit::new(3, 1);
        source.operations = vec![
            gate(GateType::CX, &[0, 1], &[]),
            gate(GateType::RZ, &[1], &[0.51]),
            boundary,
            gate(GateType::CX, &[0, 1], &[]),
        ];
        assert_eq!(
            RecoverPhaseGadgetsPass
                .try_run(&source, &mut PropertySet::new())
                .unwrap()
                .operations,
            source.operations
        );
    }
    for malformed in [
        gate(GateType::CX, &[0, 0], &[]),
        gate(GateType::CX, &[0, 3], &[]),
        gate(GateType::CX, &[0, 1], &[0.1]),
        gate(GateType::RZ, &[1], &[f64::NAN]),
        gate(GateType::RZ, &[1], &[]),
    ] {
        let mut source = Circuit::new(3, 0);
        source.operations = vec![
            gate(GateType::CX, &[0, 1], &[]),
            malformed,
            gate(GateType::CX, &[0, 1], &[]),
        ];
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            RecoverPhaseGadgetsPass.try_run(&source, &mut PropertySet::new())
        }))
        .expect("invalid input must fail without panicking");
        assert!(result.is_err());
    }
}

fn route(
    source: &Circuit,
    layout: Vec<usize>,
    iterations: usize,
    commute: bool,
) -> (Circuit, PropertySet) {
    let mut properties = PropertySet::new();
    properties.insert("sabre_initial_layout", layout);
    properties.insert("routing_commute_diagonal", commute);
    let router = BeamSabrePass {
        backend: Backend::linear(source.num_qubits + 1),
        bidir_iterations: iterations,
        ..Default::default()
    };
    let output = router.try_run(source, &mut properties).unwrap();
    full_isometry(
        source,
        &output,
        properties.get::<Vec<usize>>("initial_layout").unwrap(),
        properties.get::<Vec<usize>>("final_layout").unwrap(),
    );
    (output, properties)
}

#[test]
fn diagonal_routing_preserves_all_predecessors_before_local_and_shared_third_wire_fences() {
    for sample in 0..8 {
        let mut source = Circuit::new(4, 0);
        let a = 0.31 + sample as f64 * 0.17;
        source.operations = vec![
            gate(GateType::RZZ, &[0, 3], &[a]),
            gate(GateType::RZZ, &[0, 1], &[-2.1 * a]),
            gate(GateType::RZ, &[0], &[0.23]),
            gate(GateType::RZZ, &[2, 0], &[3.7 * a]),
            gate(GateType::H, &[0], &[]),
            gate(GateType::RZZ, &[0, 2], &[-a]),
            gate(GateType::CRZ, &[1, 3], &[0.61]),
            gate(GateType::RX, &[3], &[0.47]),
            gate(GateType::RZZ, &[0, 3], &[-0.83]),
            gate(GateType::CX, &[0, 1], &[]),
            gate(GateType::RZZ, &[1, 2], &[a]),
            gate(GateType::H, &[1], &[]),
        ];
        for iterations in [1, 3] {
            route(&source, vec![0, 2, 4, 1], iterations, true);
        }
    }
}

#[test]
fn infinitesimal_nondiagonal_local_gate_is_a_structural_fence_not_a_tolerance_decision() {
    let mut source = Circuit::new(4, 0);
    source.operations = vec![
        gate(GateType::RZZ, &[0, 3], &[0.37]),
        gate(GateType::RZZ, &[0, 1], &[0.61]),
        gate(GateType::U, &[0], &[1e-15, 0.0, 0.0]),
        gate(GateType::RZZ, &[0, 2], &[0.83]),
    ];
    let (output, properties) = route(&source, vec![0, 1, 2, 3], 3, true);
    let mut current = properties
        .get::<Vec<usize>>("initial_layout")
        .unwrap()
        .clone();
    let mut preceding = Vec::new();
    let mut saw_fence = false;
    for op in &output.operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            if *name == GateType::SWAP {
                for physical in &mut current {
                    if *physical == qubits[0] {
                        *physical = qubits[1];
                    } else if *physical == qubits[1] {
                        *physical = qubits[0];
                    }
                }
            } else if *name == GateType::RZZ {
                if params[0] == 0.83 {
                    assert!(saw_fence, "post-fence parity phase executed early");
                } else {
                    assert!(!saw_fence, "pre-fence parity phase executed late");
                    preceding.push(params[0]);
                }
            } else if *name == GateType::U {
                assert_eq!(qubits, &vec![current[0]]);
                assert_eq!(params, &vec![1e-15, 0.0, 0.0]);
                assert_eq!(
                    preceding.len(),
                    2,
                    "the gate must wait for every diagonal predecessor"
                );
                saw_fence = true;
            }
        }
    }
    assert!(saw_fence);
}

#[test]
fn long_diagonal_epochs_preserve_full_isometry_across_every_noncommuting_fence() {
    let mut source = Circuit::new(4, 0);
    for layer in 0..3 {
        // Cross both the 64-interaction cap and its next rollover. Distinct
        // phases make incorrect duplication or dropping observable.
        for index in 0..131 {
            let target = 1 + index % 3;
            let angle = (index as f64 + 1.0) * (0.0013 + layer as f64 * 0.0007);
            source.add_op(gate(GateType::RZZ, &[0, target], &[angle]));
            if index % 17 == 0 {
                source.add_op(gate(GateType::RZ, &[target], &[-0.091]));
            }
        }
        source.add_op(gate(GateType::H, &[0], &[]));
        source.add_op(gate(GateType::CX, &[2, 1], &[]));
    }
    for iterations in [1, 3] {
        route(&source, vec![0, 2, 4, 1], iterations, true);
    }
}

#[test]
fn phase_transport_through_swap_requires_relabeling_and_noncommuting_reorderings_are_detected() {
    let mut source = Circuit::new(3, 0);
    source.operations = vec![
        gate(GateType::RZZ, &[0, 2], &[0.73]),
        gate(GateType::SWAP, &[0, 1], &[]),
    ];
    let mut correct = Circuit::new(3, 0);
    correct.operations = vec![
        gate(GateType::SWAP, &[0, 1], &[]),
        gate(GateType::RZZ, &[1, 2], &[0.73]),
    ];
    equivalent(&source, &correct, 1e-12);
    let mut wrong = correct.clone();
    wrong.operations[1] = gate(GateType::RZZ, &[0, 2], &[0.73]);
    assert!(error(&columns(&source), &columns(&wrong)) > 0.5);
    for fence in [
        gate(GateType::H, &[0], &[]),
        gate(GateType::CX, &[1, 0], &[]),
    ] {
        source.operations = vec![gate(GateType::RZZ, &[0, 2], &[0.73]), fence];
        wrong = source.clone();
        wrong.operations.swap(0, 1);
        assert!(error(&columns(&source), &columns(&wrong)) > 0.4);
    }
}

fn config(backend: Backend, recover: bool, commute: bool, fuse: bool) -> TranspilerConfig {
    TranspilerConfig::builder()
        .optimization_level(3)
        .backend(backend)
        .target_basis(["u", "cx"])
        .recover_phase_gadgets(recover)
        .post_block_fusion(fuse)
        .routing_search(RoutingSearchConfig {
            layout_trials: Some(2),
            layout_iterations: Some(2),
            commute_diagonal: commute,
            ..Default::default()
        })
        .build()
}

#[test]
fn all_refinement_switches_preserve_fourier_columns_and_directed_clean_ancilla_isometries() {
    for width in [3, 4] {
        let source = qft(width);
        let dimension = 1 << width;
        let exact: Vec<Vec<C>> = (0..dimension)
            .map(|column| {
                (0..dimension)
                    .map(|row| {
                        phase(2.0 * PI * (row * column) as f64 / dimension as f64)
                            / (dimension as f64).sqrt()
                    })
                    .collect()
            })
            .collect();
        assert!(
            error(&exact, &columns(&source)) < 1e-12,
            "fixture must itself be Fourier"
        );
        for recover in [false, true] {
            for commute in [false, true] {
                for fuse in [false, true] {
                    let mut backend = Backend::new("directed_line", width + 1);
                    backend.set_coupling_map((0..width).map(|q| (q + 1, q)));
                    let (output, report) = transpile_with_report(
                        &source,
                        Some(config(backend.clone(), recover, commute, fuse)),
                    )
                    .unwrap();
                    full_isometry(
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
                            _ => panic!("compiled output escaped requested U/CX basis"),
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn normalized_qft14_routing_improves_resources_and_preserves_fourier_samples() {
    let width = 14;
    let source = qft(width);
    let build_config = |enabled| {
        TranspilerConfig::builder()
            .optimization_level(3)
            .backend(Backend::linear(width))
            .target_basis(["u", "cx"])
            .recover_phase_gadgets(enabled)
            .post_block_fusion(false)
            .routing_search(RoutingSearchConfig {
                commute_diagonal: enabled,
                ..Default::default()
            })
            .build()
    };
    // Use the common harness's normalized RZ/CX input and unchanged default
    // search budget, keeping post-block fusion off to isolate QFT routing.
    let (before, _) = transpile_with_report(&source, Some(build_config(false))).unwrap();
    let (after, report) = transpile_with_report(&source, Some(build_config(true))).unwrap();
    let before_metrics = metrics(&before);
    let after_metrics = metrics(&after);
    assert!(
        after_metrics.1 <= 280 && after_metrics.1 * 100 < before_metrics.1 * 85,
        "CX benefit regressed: {before_metrics:?} -> {after_metrics:?}"
    );
    assert!(
        after_metrics.2 <= 165 && after_metrics.2 < before_metrics.2,
        "depth benefit regressed: {before_metrics:?} -> {after_metrics:?}"
    );

    // Complete Fourier checks above cover small circuits. At width14, check
    // five widely separated Fourier columns and one coherent superposition
    // directly against the closed form, using one shared phase for all six.
    let dimension = 1usize << width;
    let initial = report.initial_layout.as_ref().unwrap();
    let final_map = report.final_layout.as_ref().unwrap();
    let embed = |index: usize, layout: &[usize]| -> usize {
        layout
            .iter()
            .enumerate()
            .map(|(q, p)| ((index >> q) & 1) << p)
            .sum()
    };
    let mut inputs: Vec<Vec<(usize, C)>> = [0, 1, 37, 8193, 16383]
        .into_iter()
        .map(|index| vec![(index, C::new(1.0, 0.0))])
        .collect();
    inputs.push(vec![
        (19, C::new(0.5, 0.0)),
        (4097, C::new(0.0, 0.5)),
        (12288, C::new(FRAC_1_SQRT_2, 0.0)),
    ]);
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    for terms in inputs {
        let mut input = vec![C::new(0.0, 0.0); dimension];
        for (index, amplitude) in &terms {
            input[embed(*index, initial)] = *amplitude;
        }
        let mut wanted = vec![C::new(0.0, 0.0); dimension];
        for row in 0..dimension {
            wanted[embed(row, final_map)] = terms
                .iter()
                .map(|(column, amplitude)| {
                    // Reduce the integer phase exactly before evaluating sin/cos.
                    amplitude
                        * phase(2.0 * PI * ((row * column) % dimension) as f64 / dimension as f64)
                        / (dimension as f64).sqrt()
                })
                .sum();
        }
        expected.push(wanted);
        actual.push(evolve(&after, input));
    }
    assert!(error(&expected, &actual) < 1e-8);
}

#[test]
fn w_cleanup_reduces_depth_without_changing_cx_layout_or_arbitrary_input_behavior() {
    for width in [3, 4, 14] {
        let source = w_circuit(width);
        for backend in [Backend::all_to_all(width), Backend::linear(width)] {
            let (before, before_report) =
                transpile_with_report(&source, Some(config(backend.clone(), true, true, false)))
                    .unwrap();
            let (after, after_report) =
                transpile_with_report(&source, Some(config(backend, true, true, true))).unwrap();
            let before_metrics = metrics(&before);
            let after_metrics = metrics(&after);
            assert_eq!(before_report.initial_layout, after_report.initial_layout);
            assert_eq!(before_report.final_layout, after_report.final_layout);
            assert_eq!(after_metrics.1, 2 * (width - 1));
            assert_eq!(after_metrics.1, before_metrics.1);
            assert_eq!(after_metrics.3, before_metrics.3);
            assert!(after_metrics.0 < before_metrics.0 && after_metrics.2 < before_metrics.2);
            if width == 14 {
                // This large case asserts the resource claim only; the small
                // widths above receive complete arbitrary-input validation.
                assert_eq!(before_metrics, (104, 26, 65, 26));
                assert_eq!(after_metrics, (92, 26, 53, 26));
            } else {
                for (output, report) in [(&before, &before_report), (&after, &after_report)] {
                    full_isometry(
                        &source,
                        output,
                        report.initial_layout.as_ref().unwrap(),
                        report.final_layout.as_ref().unwrap(),
                    );
                }
            }
        }
    }
}

#[test]
fn refinement_pipeline_preserves_global_barriers_and_noninvolutory_output_readout() {
    let mut source = qft(3);
    source.num_cbits = 2;
    source
        .operations
        .insert(5, Operation::Barrier { qubits: vec![2] });
    source.add_op(gate(GateType::SWAP, &[0, 1], &[]));
    source.add_op(gate(GateType::SWAP, &[1, 2], &[]));
    let prefix = source.clone();
    let readout = [
        Operation::Measure { qubit: 2, cbit: 0 },
        Operation::Barrier { qubits: vec![0, 2] },
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Measure { qubit: 1, cbit: 1 },
    ];
    source.operations.extend(readout.iter().cloned());
    for recover in [false, true] {
        for commute in [false, true] {
            let (mut output, report) = transpile_with_report(
                &source,
                Some(config(Backend::linear(4), recover, commute, true)),
            )
            .unwrap();
            let final_map = report.final_layout.as_ref().unwrap();
            let expected: Vec<_> = readout
                .iter()
                .map(|op| match op {
                    Operation::Measure { qubit, cbit } => Operation::Measure {
                        qubit: final_map[*qubit],
                        cbit: *cbit,
                    },
                    Operation::Barrier { qubits } => Operation::Barrier {
                        qubits: qubits.iter().map(|q| final_map[*q]).collect(),
                    },
                    _ => unreachable!(),
                })
                .collect();
            assert_eq!(
                &output.operations[output.operations.len() - 4..],
                expected.as_slice()
            );
            output.operations.truncate(output.operations.len() - 4);
            assert_eq!(
                output
                    .operations
                    .iter()
                    .filter(|op| op.is_barrier())
                    .count(),
                1
            );
            full_isometry(
                &prefix,
                &output,
                report.initial_layout.as_ref().unwrap(),
                final_map,
            );
        }
    }
}
