//! Controlled-phase routing checks against an independent bitwise oracle.
//! Every input column is checked with one shared phase, including idle wires.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::BeamSabrePass;
use q_rust::transpiler::{transpile_with_report, TranspilerConfig};

fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    }
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
            GateType::U => {
                let (s, c) = (params[0] / 2.0).sin_cos();
                let matrix = [
                    [C::new(c, 0.0), -C::from_polar(s, params[2])],
                    [
                        C::from_polar(s, params[1]),
                        C::from_polar(c, params[1] + params[2]),
                    ],
                ];
                for low in 0..state.len() {
                    if low & (1 << qubits[0]) != 0 {
                        continue;
                    }
                    let high = low | (1 << qubits[0]);
                    let (a, b) = (state[low], state[high]);
                    state[low] = matrix[0][0] * a + matrix[0][1] * b;
                    state[high] = matrix[1][0] * a + matrix[1][1] * b;
                }
            }
            other => panic!("unsupported oracle gate {other:?}"),
        }
    }
    state
}

fn verify(source: &Circuit, output: &Circuit, initial: &[usize], final_map: &[usize]) {
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
        expected.push(column);
        let mut input = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        input[embed(basis, initial)] = C::new(1.0, 0.0);
        actual.push(evolve(output, input));
    }
    let overlap: C = expected
        .iter()
        .flatten()
        .zip(actual.iter().flatten())
        .map(|(a, b)| a.conj() * b)
        .sum();
    assert!(overlap.norm() > 0.0);
    let phase = overlap / overlap.norm();
    let residual = expected
        .iter()
        .zip(&actual)
        .map(|(a, b)| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (*b - phase * a).norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max);
    assert!(
        residual < 1e-9,
        "shared-phase full-isometry residual {residual:e}"
    );
    for op in &output.operations {
        if let Operation::Gate { qubits, .. } = op {
            if qubits.len() == 2 {
                assert_eq!(qubits[0].abs_diff(qubits[1]), 1);
            }
        }
    }
}

fn route(source: &Circuit, physical: usize, commute: bool) -> (Circuit, PropertySet) {
    let pass = BeamSabrePass {
        backend: Backend::linear(physical),
        beam_width: 8,
        branch_factor: 5,
        bidir_iterations: 1,
        ..BeamSabrePass::default()
    };
    let mut properties = PropertySet::new();
    properties.insert("routing_commute_diagonal", commute);
    let output = pass.try_run(source, &mut properties).unwrap();
    verify(
        source,
        &output,
        properties.get::<Vec<usize>>("initial_layout").unwrap(),
        properties.get::<Vec<usize>>("final_layout").unwrap(),
    );
    (output, properties)
}

#[test]
fn diagonal_fanout_can_execute_nearby_interactions_before_distant_ones() {
    let mut circuit = Circuit::new(5, 0);
    for (target, angle) in [(4, 0.23), (1, -0.41), (3, 0.67), (2, -0.89)] {
        circuit.add_op(gate(GateType::CU1, &[0, target], &[angle]));
    }
    let (_, strict) = route(&circuit, 5, false);
    let (_, commuting) = route(&circuit, 5, true);
    let strict_swaps = *strict.get::<usize>("swaps_inserted").unwrap();
    let commuting_swaps = *commuting.get::<usize>("swaps_inserted").unwrap();
    assert!(
        commuting_swaps < strict_swaps,
        "commuting {commuting_swaps}, strict {strict_swaps}"
    );
}

#[test]
fn controlled_phase_reordering_respects_local_fences_and_barriers() {
    for theta in [0.0, 1e-15, -1e-15, 0.71, std::f64::consts::TAU] {
        for reversed in [false, true] {
            let mut circuit = Circuit::new(4, 0);
            for (a, b, angle) in [(0, 3, 0.31), (0, 1, -0.73), (2, 0, 1e-12)] {
                let wires = if reversed { [b, a] } else { [a, b] };
                circuit.add_op(gate(GateType::CU1, &wires, &[angle]));
            }
            circuit.add_op(gate(GateType::U, &[0], &[theta, 0.27, -0.19]));
            circuit.add_op(gate(GateType::CU1, &[2, 0], &[-0.53]));
            circuit.add_op(gate(GateType::CU1, &[0, 3], &[0.97]));
            circuit.add_op(Operation::Barrier {
                qubits: vec![0, 1, 2, 3],
            });
            circuit.add_op(gate(GateType::U, &[3], &[0.43, -0.17, 0.29]));
            circuit.add_op(gate(GateType::CU1, &[3, 1], &[-1.13]));
            route(&circuit, 5, true);
        }
    }
}

#[test]
fn mixed_controlled_phase_epochs_preserve_the_public_pipeline_isometry() {
    let mut circuit = Circuit::new(5, 0);
    for (target, angle) in [(4, 0.29), (1, -0.43), (3, 0.61), (2, -0.83)] {
        circuit.add_op(gate(GateType::CU1, &[0, target], &[angle]));
    }
    circuit.add_op(gate(GateType::U, &[0], &[0.79, 0.17, -0.37]));
    for target in [3, 1, 4, 2] {
        circuit.add_op(gate(GateType::CU1, &[target, 0], &[-0.31 * target as f64]));
    }
    circuit.add_op(gate(GateType::CX, &[1, 4], &[]));
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(6))
        .target_basis(["u", "cx"])
        .build();
    let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
    verify(
        &circuit,
        &output,
        report.initial_layout.as_ref().unwrap(),
        report.final_layout.as_ref().unwrap(),
    );
}
