//! Search-policy regressions verify the full clean-ancilla isometry, including
//! amplitudes outside the logical output subspace and one shared global phase.
use num_complex::Complex64;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::simulator::circuit_to_unitary;
use q_rust::transpiler::optimization::{
    ApplyOutputPermutationPass, ExtractOutputPermutationPass, TrailingSwapElisionPass,
};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::BeamSabrePass;
use q_rust::transpiler::{transpile_with_report, RoutingSearchConfig, TranspilerConfig};

fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    }
}

fn three_cycle() -> Circuit {
    let mut circuit = Circuit::new(3, 0);
    circuit.add_op(gate(GateType::H, &[0], &[]));
    circuit.add_op(gate(GateType::RX, &[1], &[0.37]));
    circuit.add_op(gate(GateType::CX, &[0, 2], &[]));
    // This permutation is not self-inverse: testing only QFT reversal would
    // miss applying the output permutation in the wrong direction.
    circuit.add_op(gate(GateType::SWAP, &[0, 1], &[]));
    circuit.add_op(gate(GateType::SWAP, &[1, 2], &[]));
    circuit
}

fn assert_full_isometry(
    original: &Circuit,
    output: &Circuit,
    initial: &[usize],
    final_layout: &[usize],
) {
    let source = circuit_to_unitary(original);
    let actual = circuit_to_unitary(output);
    let embed = |index: usize, layout: &[usize]| -> usize {
        layout
            .iter()
            .enumerate()
            .map(|(l, &p)| ((index >> l) & 1) << p)
            .sum()
    };
    let mut phase = None;
    for col in 0..source.ncols() {
        let physical_col = embed(col, initial);
        for physical_row in 0..actual.nrows() {
            let expected = (0..source.nrows())
                .find(|&row| embed(row, final_layout) == physical_row)
                .map(|row| source[(row, col)])
                .unwrap_or(Complex64::new(0.0, 0.0));
            let observed = actual[(physical_row, physical_col)];
            if phase.is_none() && expected.norm() > 1e-8 {
                phase = Some(observed / expected);
            }
            let factor = phase.unwrap_or(Complex64::new(1.0, 0.0));
            assert!(
                (observed - factor * expected).norm() < 1e-8,
                "isometry mismatch input={col} physical_output={physical_row}"
            );
        }
    }
    assert!((phase.unwrap().norm() - 1.0).abs() < 1e-8);
}

#[test]
fn logical_three_cycle_composes_with_nontrivial_layout_and_ancilla() {
    let source = three_cycle();
    let mut properties = PropertySet::new();
    let extracted = ExtractOutputPermutationPass.run(&source, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("logical_output_permutation"),
        Some(&vec![1, 2, 0])
    );
    properties.insert("sabre_initial_layout", vec![3usize, 1, 0]);
    let router = BeamSabrePass {
        backend: Backend::linear(4),
        bidir_iterations: 1,
        ..Default::default()
    };
    let routed = router.try_run(&extracted, &mut properties).unwrap();
    let lowered = TrailingSwapElisionPass.run(&routed, &mut properties);
    let output = ApplyOutputPermutationPass.run(&lowered, &mut properties);
    assert_full_isometry(
        &source,
        &output,
        properties.get::<Vec<usize>>("initial_layout").unwrap(),
        properties.get::<Vec<usize>>("final_layout").unwrap(),
    );
}

#[test]
fn all_search_ablation_policies_preserve_arbitrary_inputs() {
    let source = three_cycle();
    for lookahead_size in [0, 20] {
        for rank_progress in [false, true] {
            for elide_output_permutations in [false, true] {
                let options = RoutingSearchConfig {
                    layout_trials: Some(3),
                    layout_iterations: Some(2),
                    lookahead_size,
                    rank_progress,
                    elide_output_permutations,
                    orient_swaps: true,
                    commute_diagonal: true,
                };
                let config = TranspilerConfig::builder()
                    .optimization_level(3)
                    .backend(Backend::linear(4))
                    .target_basis(["u", "cx"])
                    .routing_search(options)
                    .build();
                let (output, report) = transpile_with_report(&source, Some(config)).unwrap();
                assert_full_isometry(
                    &source,
                    &output,
                    report.initial_layout.as_ref().unwrap(),
                    report.final_layout.as_ref().unwrap(),
                );
            }
        }
    }
}

#[test]
fn pure_output_permutation_requires_no_physical_gates() {
    let mut source = Circuit::new(4, 0);
    source.add_op(gate(GateType::SWAP, &[0, 3], &[]));
    source.add_op(gate(GateType::SWAP, &[1, 3], &[]));
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(5))
        .build();
    let (output, report) = transpile_with_report(&source, Some(config)).unwrap();
    assert!(output.operations.is_empty());
    assert_full_isometry(
        &source,
        &output,
        report.initial_layout.as_ref().unwrap(),
        report.final_layout.as_ref().unwrap(),
    );
}

#[test]
fn extraction_respects_barriers_and_terminal_classical_write_order() {
    let mut source = three_cycle();
    source.num_cbits = 2;
    source.operations.insert(
        4,
        Operation::Barrier {
            qubits: vec![0, 1, 2],
        },
    );
    source.add_op(Operation::Measure { qubit: 1, cbit: 0 });
    source.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    source.add_op(Operation::Measure { qubit: 2, cbit: 1 });
    let mut properties = PropertySet::new();
    let extracted = ExtractOutputPermutationPass.run(&source, &mut properties);
    // Direct pass has no permission to move a permutation across readout.
    assert_eq!(extracted.operations, source.operations);
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(4))
        .build();
    let (mut output, report) = transpile_with_report(&source, Some(config)).unwrap();
    let final_layout = report.final_layout.as_ref().unwrap();
    assert_eq!(
        &output.operations[output.operations.len() - 3..],
        &[
            Operation::Measure {
                qubit: final_layout[1],
                cbit: 0
            },
            Operation::Measure {
                qubit: final_layout[0],
                cbit: 0
            },
            Operation::Measure {
                qubit: final_layout[2],
                cbit: 1
            },
        ]
    );
    assert_eq!(
        output
            .operations
            .iter()
            .filter(|op| op.is_barrier())
            .count(),
        1
    );
    source.operations.truncate(source.operations.len() - 3);
    output.operations.truncate(output.operations.len() - 3);
    assert_full_isometry(
        &source,
        &output,
        report.initial_layout.as_ref().unwrap(),
        final_layout,
    );
}

#[test]
fn diagnostics_do_not_change_output_and_are_opt_in() {
    let source = three_cycle();
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .backend(Backend::linear(4))
        .build();
    let (plain, plain_report) = transpile_with_report(&source, Some(config.clone())).unwrap();
    assert!(plain_report.pass_diagnostics.is_empty());
    let mut profiled = config;
    profiled.profile_passes = true;
    let (observed, report) = transpile_with_report(&source, Some(profiled)).unwrap();
    assert_eq!(plain.operations, observed.operations);
    assert_eq!(plain_report.initial_layout, report.initial_layout);
    assert_eq!(plain_report.final_layout, report.final_layout);
    assert!(report
        .pass_diagnostics
        .iter()
        .any(|p| p.output.stage == "SabreLayoutPass"));
}

#[test]
fn swap_orientation_exposes_both_directed_cx_cancellations() {
    for pair in [[0, 1], [1, 0]] {
        for before in [false, true] {
            let mut source = Circuit::new(2, 0);
            let cx = gate(GateType::CX, &pair, &[]);
            let swap = gate(GateType::SWAP, &[0, 1], &[]);
            if before {
                source.add_op(cx);
                source.add_op(swap);
            } else {
                source.add_op(swap);
                source.add_op(cx);
            }
            // Keep the final SWAP semantic (not a free output permutation).
            source.add_op(gate(GateType::H, &[0], &[]));
            source.add_op(gate(GateType::H, &[1], &[]));
            let config = TranspilerConfig::builder()
                .optimization_level(3)
                .backend(Backend::all_to_all(2))
                .target_basis(["u", "cx"])
                .routing_search(RoutingSearchConfig {
                    elide_output_permutations: false,
                    ..Default::default()
                })
                .build();
            let (output, report) = transpile_with_report(&source, Some(config)).unwrap();
            let cx_count = output
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
                .count();
            assert_eq!(cx_count, 2, "direction={pair:?} before={before}");
            assert_full_isometry(
                &source,
                &output,
                report.initial_layout.as_ref().unwrap(),
                report.final_layout.as_ref().unwrap(),
            );
        }
    }
}

#[test]
fn swap_orientation_does_not_skip_noncommuting_single_qubit_gates() {
    use q_rust::transpiler::optimization::OrientSwapDecompositionPass;
    let mut circuit = Circuit::new(2, 0);
    circuit.add_op(gate(GateType::CX, &[1, 0], &[]));
    circuit.add_op(gate(GateType::RZ, &[0], &[0.37]));
    circuit.add_op(gate(GateType::SWAP, &[0, 1], &[]));
    circuit.add_op(gate(GateType::RX, &[1], &[0.29]));
    circuit.add_op(gate(GateType::CX, &[1, 0], &[]));
    let output = OrientSwapDecompositionPass.run(&circuit, &mut PropertySet::new());
    assert_eq!(output.operations, circuit.operations);
}

#[test]
fn permutation_metadata_is_cleared_on_reuse_and_consumed_once() {
    let mut properties = PropertySet::new();
    ExtractOutputPermutationPass.run(&three_cycle(), &mut properties);
    assert!(properties.contains("logical_output_permutation"));
    let mut next = Circuit::new(2, 0);
    next.add_op(gate(GateType::H, &[0], &[]));
    ExtractOutputPermutationPass.run(&next, &mut properties);
    assert!(!properties.contains("logical_output_permutation"));
    properties.insert("final_layout", vec![0usize, 1]);
    ApplyOutputPermutationPass.run(&next, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&vec![0, 1])
    );
    properties.insert("logical_output_permutation", vec![1usize, 0]);
    ApplyOutputPermutationPass.run(&next, &mut properties);
    ApplyOutputPermutationPass.run(&next, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&vec![1, 0])
    );
}
