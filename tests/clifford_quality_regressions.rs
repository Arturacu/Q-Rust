//! Independent arbitrary-input checks for final Clifford candidate selection.
//!
//! Literal matrices and complete clean-ancilla isometries keep relative phases,
//! layouts, spare-wire leakage and output quality independent of compiler code.
//!
//! The two QASMBench-derived fixtures below are Copyright © 2020, Battelle
//! Memorial Institute. Their license and notice are retained in
//! experiments/external_sources/qasmbench/{LICENSE,NOTICE}.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::clifford_simplification::{
    CliffordSimplificationPass, CliffordSimplificationStats,
};
use q_rust::transpiler::numerical_policy::{NumericalBudget, CIRCUIT_ERROR_ALLOWANCE};
use q_rust::transpiler::optimization::{ExtractOutputPermutationPass, GateCrystallizationPass};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::{transpile, transpile_with_report, TranspilationReport, TranspilerConfig};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

// QASMBench qec_en_n5, normalized unitary prefix. Retain the non-Clifford T
// preparation: the comparison is of arbitrary-input actions, not |0>-only states.
const QEC_BODY: &str = "h q[2]; t q[2]; h q[2]; h q[0]; h q[1]; h q[2];
    cx q[1],q[2]; cx q[0],q[2]; h q[0]; h q[1]; h q[3]; cx q[3],q[2];
    h q[2]; h q[3]; cx q[3],q[2]; cx q[0],q[2]; cx q[1],q[2];
    h q[2]; h q[4]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2];
    cx q[1],q[2]; cx q[3],q[2];";

// The normalized QASMBench error_correctiond3_n5 unitary prefix, from source
// revision 357b942396d5c2b7cbc1c229c585a6ef5ccaebac. Its twelve encoded
// permutations include one interleaved window; keep the complete gate order.
const ERROR_CORRECTION_BODY: &str = "h q[0]; h q[1]; id q[2]; h q[3]; h q[4];
cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2];
cx q[4],q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2];
sdg q[4]; cx q[4],q[2]; h q[2]; cx q[4],q[2]; h q[2];
cx q[0],q[2]; h q[0]; h q[2]; cx q[0],q[2]; h q[0]; h q[2]; cx q[0],q[2];
cx q[3],q[2]; cx q[0],q[2]; h q[0]; h q[2]; cx q[0],q[2]; h q[0]; h q[2]; cx q[0],q[2];
cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2];
cx q[3],q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2]; h q[1]; h q[2]; cx q[1],q[2];
cx q[3],q[2]; cx q[0],q[2]; h q[3]; h q[4];
cx q[3],q[2]; h q[2]; h q[3]; cx q[3],q[2]; h q[2]; h q[3]; cx q[3],q[2];
cx q[0],q[2]; cx q[3],q[2]; h q[2]; h q[3]; cx q[3],q[2]; h q[2]; h q[3]; cx q[3],q[2];
cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2];
cx q[0],q[2]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2];
cx q[1],q[2]; h q[2]; cx q[1],q[2]; h q[1];
cx q[4],q[2]; cx q[0],q[1]; h q[2]; h q[4]; h q[1]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2];
cx q[1],q[2]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2]; h q[2]; h q[4]; cx q[4],q[2];
sdg q[1]; h q[2]; h q[3]; h q[4];";

fn source(width: usize, body: &str) -> Circuit {
    parse_qasm(&format!(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[{width}]; {body}"
    ))
    .unwrap()
}

fn phase(theta: f64) -> C {
    C::from_polar(1.0, theta)
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
        match name {
            GateType::CX => {
                for low in 0..state.len() {
                    if low & (1 << qubits[0]) != 0 && low & (1 << qubits[1]) == 0 {
                        state.swap(low, low | (1 << qubits[1]));
                    }
                }
            }
            GateType::CZ => {
                for (basis, amp) in state.iter_mut().enumerate() {
                    if qubits.iter().all(|q| basis & (1 << q) != 0) {
                        *amp = -*amp;
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
                let matrix = match name {
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
                    _ => panic!("unsupported literal-oracle gate {name:?}"),
                };
                let target = *qubits.last().unwrap();
                for low in 0..state.len() {
                    if low & (1 << target) != 0
                        || (*name == GateType::CY && low & (1 << qubits[0]) == 0)
                    {
                        continue;
                    }
                    let high = low | (1 << target);
                    let (a, b) = (state[low], state[high]);
                    state[low] = matrix[0][0] * a + matrix[0][1] * b;
                    state[high] = matrix[1][0] * a + matrix[1][1] * b;
                }
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

// Independently recognize the local Clifford action by conjugating the three
// literal Pauli matrices. This does not use the compiler's gate classifier,
// angle snapping, or any of its stored Clifford matrices.
fn local_gate_preserves_paulis(operation: &Operation) -> bool {
    let Operation::Gate {
        name,
        params,
        qubits,
    } = operation
    else {
        return false;
    };
    if qubits.len() != 1 {
        return false;
    }
    let zero = C::new(0.0, 0.0);
    let one = C::new(1.0, 0.0);
    let i = C::new(0.0, 1.0);
    let mut local = Circuit::new(1, 0);
    local.add_op(Operation::Gate {
        name: name.clone(),
        params: params.clone(),
        qubits: vec![0],
    });
    let columns = [
        evolve(&local, vec![one, zero]),
        evolve(&local, vec![zero, one]),
    ];
    let u = [
        [columns[0][0], columns[1][0]],
        [columns[0][1], columns[1][1]],
    ];
    let paulis = [
        [[zero, one], [one, zero]],
        [[zero, -i], [i, zero]],
        [[one, zero], [zero, -one]],
    ];
    paulis.iter().all(|pauli| {
        let mut image = [[zero; 2]; 2];
        for (a, row) in image.iter_mut().enumerate() {
            for (b, entry) in row.iter_mut().enumerate() {
                for (c, pauli_row) in pauli.iter().enumerate() {
                    for (d, p) in pauli_row.iter().enumerate() {
                        *entry += u[a][c] * p * u[b][d].conj();
                    }
                }
            }
        }
        paulis.iter().any(|target| {
            [-1.0, 1.0].into_iter().any(|sign| {
                (0..2).all(|a| (0..2).all(|b| (image[a][b] - sign * target[a][b]).norm() < 1e-12))
            })
        })
    })
}

fn assert_isometry(input: &Circuit, output: &Circuit, report: &TranspilationReport) {
    let n = input.num_qubits;
    let identity: Vec<_> = (0..n).collect();
    let initial = report.initial_layout.as_ref().unwrap_or(&identity);
    let final_map = report.final_layout.as_ref().unwrap_or(&identity);
    for map in [initial, final_map] {
        assert!(map.len() >= n);
        for (index, wire) in map[..n].iter().enumerate() {
            assert!(*wire < output.num_qubits);
            assert!(!map[..index].contains(wire));
        }
    }
    let mut actual = Vec::new();
    let mut expected = Vec::new();
    for column in 0..1 << n {
        let mut logical_input = vec![C::new(0.0, 0.0); 1 << n];
        logical_input[column] = C::new(1.0, 0.0);
        let logical = evolve(input, logical_input);
        let mut physical_input = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        physical_input[embed(column, &initial[..n])] = C::new(1.0, 0.0);
        actual.extend(evolve(output, physical_input));
        let mut wanted = vec![C::new(0.0, 0.0); 1 << output.num_qubits];
        for (row, amplitude) in logical.into_iter().enumerate() {
            wanted[embed(row, &final_map[..n])] = amplitude;
        }
        expected.extend(wanted);
    }
    let overlap: C = expected
        .iter()
        .zip(&actual)
        .map(|(e, a)| e.conj() * a)
        .sum();
    assert!(overlap.norm() > 1e-12);
    let global = overlap / overlap.norm();
    let error = expected
        .iter()
        .zip(&actual)
        .map(|(e, a)| (*a - global * e).norm_sqr())
        .sum::<f64>()
        .sqrt();
    assert!(error < 2e-10, "full-isometry residual {error:e}");
}

fn cost(circuit: &Circuit) -> [usize; 3] {
    let mut count = 0;
    let mut total = vec![0; circuit.num_qubits];
    let mut entangling = total.clone();
    for operation in &circuit.operations {
        let wires = operation.qubits();
        let barrier = matches!(operation, Operation::Barrier { .. });
        let is_two = matches!(operation, Operation::Gate { qubits, .. } if qubits.len() == 2);
        let layer = wires.iter().map(|&q| total[q]).max().unwrap_or(0) + usize::from(!barrier);
        for &wire in wires {
            total[wire] = layer;
        }
        if is_two || barrier {
            let layer =
                wires.iter().map(|&q| entangling[q]).max().unwrap_or(0) + usize::from(is_two);
            for &wire in wires {
                entangling[wire] = layer;
            }
        }
        count += usize::from(is_two);
    }
    [
        count,
        *total.iter().max().unwrap(),
        *entangling.iter().max().unwrap(),
    ]
}

fn configuration(
    backend: Backend,
    enabled: bool,
    retention: bool,
    profile: bool,
) -> TranspilerConfig {
    TranspilerConfig::builder()
        .optimization_level(3)
        .backend(backend)
        .target_basis(["u", "cx"])
        .encoded_swap_recognition(true)
        .clifford_simplification(enabled)
        .retain_routing_candidates(retention)
        .profile_passes(profile)
        .build()
}

fn assert_no_regression(candidate: &Circuit, baseline: &Circuit) {
    let (a, b) = (cost(candidate), cost(baseline));
    assert!(
        a.iter().zip(b).all(|(a, b)| *a <= b),
        "candidate {a:?} regresses baseline {b:?}"
    );
}

fn assert_retention_costs(selected: &Circuit, baseline: &Circuit, report: &TranspilationReport) {
    let stats = report.clifford_retention_stats.as_ref().unwrap();
    assert!(!stats.retained_compilation_failed);
    assert!(stats.attempted_simplification.rewrites > 0);
    assert_eq!(stats.retained_output_cost, Some(cost(baseline)));
    assert_eq!(stats.selected_output_cost, cost(selected));
    if stats.selected_retained_input {
        assert_eq!(selected, baseline);
        assert!(report.clifford_simplification_stats.is_none());
    }
}

#[test]
fn qec_encoder_line_retains_final_entangling_depth() {
    let input = source(5, QEC_BODY);
    let (baseline, baseline_report) = transpile_with_report(
        &input,
        Some(configuration(Backend::linear(5), false, true, false)),
    )
    .unwrap();
    assert_eq!(cost(&baseline), [20, 21, 16]);
    let config = configuration(Backend::linear(5), true, true, false);
    let (selected, report) = transpile_with_report(&input, Some(config.clone())).unwrap();
    assert_no_regression(&selected, &baseline);
    assert_retention_costs(&selected, &baseline, &report);
    assert_isometry(&input, &baseline, &baseline_report);
    assert_isometry(&input, &selected, &report);
    assert_eq!(selected, transpile(&input, Some(config)).unwrap());
    assert!(report.pass_diagnostics.is_empty());
}

#[test]
fn error_correction_keeps_nine_cx_and_improves_the_routed_circuit() {
    let input = source(5, ERROR_CORRECTION_BODY);
    assert_eq!(cost(&input)[0], 49);
    for (backend, upper_cost) in [
        (Backend::all_to_all(5), [9, 11, 6]),
        (Backend::linear(5), [17, 19, 12]),
    ] {
        let (baseline, _) = transpile_with_report(
            &input,
            Some(configuration(backend.clone(), false, true, false)),
        )
        .unwrap();
        let (output, report) =
            transpile_with_report(&input, Some(configuration(backend, true, true, false))).unwrap();
        assert_no_regression(&output, &baseline);
        let result = cost(&output);
        assert!(
            result
                .iter()
                .zip(upper_cost)
                .all(|(got, bound)| *got <= bound),
            "error-correction quality regressed: {result:?}"
        );
        assert_retention_costs(&output, &baseline, &report);
        assert_isometry(&input, &output, &report);
    }
}

#[test]
fn quality_comparison_preserves_directed_spares_and_readout_writes() {
    let unitary = source(5, QEC_BODY);
    let mut input = unitary.clone();
    input.num_cbits = 2;
    for (qubit, cbit) in [(4, 0), (1, 0), (3, 1)] {
        input.add_op(Operation::Measure { qubit, cbit });
    }
    for directed in [false, true] {
        let edges: Vec<_> = if directed {
            (0..5).map(|q| (q, q + 1)).collect()
        } else {
            (0..5).flat_map(|q| [(q, q + 1), (q + 1, q)]).collect()
        };
        let mut backend = Backend::new("six-wire line", 6);
        backend.set_coupling_map(edges.clone());
        let (baseline, _) = transpile_with_report(
            &input,
            Some(configuration(backend.clone(), false, true, false)),
        )
        .unwrap();
        let config = configuration(backend, true, true, false);
        let (selected, report) = transpile_with_report(&input, Some(config.clone())).unwrap();
        assert_no_regression(&selected, &baseline);
        assert_retention_costs(&selected, &baseline, &report);
        let map = report.final_layout.as_ref().unwrap();
        let observed: Vec<_> = selected
            .operations
            .iter()
            .filter_map(|op| match op {
                Operation::Measure { qubit, cbit } => Some((*qubit, *cbit)),
                _ => None,
            })
            .collect();
        assert_eq!(observed, vec![(map[4], 0), (map[1], 0), (map[3], 1)]);
        let mut output_unitary = selected.clone();
        output_unitary
            .operations
            .retain(|op| !matches!(op, Operation::Measure { .. }));
        for op in &output_unitary.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                assert!(matches!(name, GateType::U | GateType::CX));
                if *name == GateType::CX {
                    assert!(edges.contains(&(qubits[0], qubits[1])));
                }
            }
        }
        assert_isometry(&unitary, &output_unitary, &report);
        assert_eq!(selected, transpile(&input, Some(config)).unwrap());
    }
}

#[test]
fn candidate_ledger_keeps_the_shared_prefix_and_selected_branch() {
    let prefix = "u(pi/2,0,pi) q[0]; rz(1e-9) q[1]; barrier q;";
    let input = source(5, &format!("{prefix} {QEC_BODY}"));
    let mut prefix_properties = PropertySet::new();
    GateCrystallizationPass::default().run(&input, &mut prefix_properties);
    let common = prefix_properties
        .get::<NumericalBudget>("numerical_budget")
        .unwrap();
    assert!(common.accepted > 0 && common.charged_error > 0.0);
    let (baseline, baseline_report) = transpile_with_report(
        &input,
        Some(configuration(Backend::linear(5), false, true, false)),
    )
    .unwrap();
    let (selected, report) = transpile_with_report(
        &input,
        Some(configuration(Backend::linear(5), true, true, false)),
    )
    .unwrap();
    assert_no_regression(&selected, &baseline);
    assert_retention_costs(&selected, &baseline, &report);
    assert!(report.numerical_budget.charged_error >= common.charged_error);
    assert!(report.numerical_budget.charged_error <= CIRCUIT_ERROR_ALLOWANCE);
    assert_isometry(&input, &selected, &report);
    if report
        .clifford_retention_stats
        .as_ref()
        .unwrap()
        .selected_retained_input
    {
        assert_eq!(report.initial_layout, baseline_report.initial_layout);
        assert_eq!(report.final_layout, baseline_report.final_layout);
        assert_eq!(
            report.numerical_budget, baseline_report.numerical_budget,
            "discarded branches must not debit the selected ledger"
        );
    }
}

#[test]
fn retention_opt_out_and_profiling_keep_their_contracts() {
    let input = source(5, QEC_BODY);
    let mut profile = configuration(Backend::linear(5), true, true, true);
    let (profiled, mut profiled_report) =
        transpile_with_report(&input, Some(profile.clone())).unwrap();
    assert!(profiled_report
        .pass_diagnostics
        .iter()
        .any(|p| p.output.stage.starts_with("pre-clifford-input::")));
    profile.profile_passes = false;
    let (ordinary, ordinary_report) = transpile_with_report(&input, Some(profile)).unwrap();
    assert_eq!(ordinary, profiled);
    profiled_report.pass_diagnostics.clear();
    assert_eq!(ordinary_report, profiled_report);
    for enabled in [false, true] {
        let config = configuration(Backend::linear(5), enabled, false, false);
        let (output, report) = transpile_with_report(&input, Some(config.clone())).unwrap();
        assert!(report.clifford_retention_stats.is_none());
        assert_isometry(&input, &output, &report);
        assert_eq!(output, transpile(&input, Some(config)).unwrap());
    }
}

#[test]
fn cz_target_counts_entanglers_and_does_not_assume_local_u_cx_fusion() {
    let input = source(5, QEC_BODY);
    for backend in [Backend::linear(5), Backend::all_to_all(5)] {
        let baseline_config = TranspilerConfig::builder()
            .optimization_level(3)
            .backend(backend.clone())
            .target_basis(["u", "cz"])
            .clifford_simplification(false)
            .build();
        let config = TranspilerConfig::builder()
            .optimization_level(3)
            .backend(backend)
            .target_basis(["u", "cz"])
            .build();
        let (baseline, _) = transpile_with_report(&input, Some(baseline_config)).unwrap();
        let (selected, report) = transpile_with_report(&input, Some(config)).unwrap();
        assert_no_regression(&selected, &baseline);
        assert_retention_costs(&selected, &baseline, &report);
        assert!(cost(&selected)[0] > 0);
        for op in &selected.operations {
            if let Operation::Gate { name, .. } = op {
                assert!(matches!(name, GateType::U | GateType::CZ));
            }
        }
        assert_isometry(&input, &selected, &report);
    }
}

#[test]
fn lower_levels_and_no_backend_do_not_enter_clifford_retention() {
    let input = source(5, QEC_BODY);
    for level in 0..3 {
        let mut enabled = configuration(Backend::linear(5), true, true, false);
        enabled.optimization_level = level;
        let mut disabled = enabled.clone();
        disabled.clifford_simplification = false;
        disabled.encoded_swap_recognition = false;
        let (a, ra) = transpile_with_report(&input, Some(enabled)).unwrap();
        let (b, rb) = transpile_with_report(&input, Some(disabled)).unwrap();
        assert_eq!(a, b);
        assert_eq!(ra, rb);
        assert!(ra.clifford_retention_stats.is_none());
        assert_isometry(&input, &a, &ra);
    }
    let config = TranspilerConfig::builder()
        .optimization_level(3)
        .target_basis(["u", "cx"])
        .build();
    let (output, report) = transpile_with_report(&input, Some(config)).unwrap();
    assert!(report.clifford_retention_stats.is_none());
    assert_isometry(&input, &output, &report);
}

#[test]
fn no_elision_and_named_gate_outputs_do_not_enter_fused_retention() {
    let input = source(5, QEC_BODY);
    let mut no_elision = configuration(Backend::linear(5), true, true, false);
    no_elision.routing_search.elide_output_permutations = false;
    let named = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(false)
        .backend(Backend::linear(5))
        .build();
    for config in [no_elision, named] {
        let (output, report) = transpile_with_report(&input, Some(config.clone())).unwrap();
        assert!(report.clifford_retention_stats.is_none());
        assert_isometry(&input, &output, &report);
        assert_eq!(output, transpile(&input, Some(config)).unwrap());
    }
}

#[test]
fn relabelled_qec_cases_compare_final_cost_and_complete_isometry() {
    let original = source(5, QEC_BODY);
    for shift in 0..5 {
        let mut input = original.clone();
        for operation in &mut input.operations {
            if let Operation::Gate { qubits, .. } = operation {
                for q in qubits {
                    *q = (*q + shift) % 5;
                }
            }
        }
        for complete in [false, true] {
            let backend = if complete {
                Backend::all_to_all(5)
            } else {
                Backend::linear(5)
            };
            let (baseline, _) = transpile_with_report(
                &input,
                Some(configuration(backend.clone(), false, true, false)),
            )
            .unwrap();
            let (selected, report) =
                transpile_with_report(&input, Some(configuration(backend, true, true, false)))
                    .unwrap();
            assert_no_regression(&selected, &baseline);
            assert_isometry(&input, &selected, &report);
        }
    }
}

fn all_four_wire_permutations() -> Vec<[usize; 4]> {
    let mut permutations = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let p = [a, b, c, d];
                    if p.iter().enumerate().all(|(i, q)| !p[..i].contains(q)) {
                        permutations.push(p);
                    }
                }
            }
        }
    }
    assert_eq!(permutations.len(), 24);
    permutations
}

fn meeting_witness(p: [usize; 4], blocked: bool) -> Circuit {
    let [a, b, c, d] = p;
    let first = format!("cx q[{c}],q[{b}];");
    let second = format!("cx q[{b}],q[{d}];");
    let middle = if blocked {
        format!("{second} {first}")
    } else {
        format!("{first} {second}")
    };
    source(4, &format!("cx q[{a}],q[{b}]; {middle} cx q[{b}],q[{a}];"))
}

#[test]
fn two_sided_meeting_reaches_an_interior_cut_with_exact_output_maps() {
    for permutation in all_four_wire_permutations() {
        let input = meeting_witness(permutation, false);
        let mut properties = PropertySet::new();
        let output = CliffordSimplificationPass {
            allow_swap_replacements: true,
        }
        .run(&input, &mut properties);
        assert_isometry(&input, &output, &TranspilationReport::default());
        let folded = ExtractOutputPermutationPass.run(&output, &mut properties);
        let map = properties
            .get::<Vec<usize>>("logical_output_permutation")
            .unwrap()
            .clone();
        let report = TranspilationReport {
            final_layout: Some(map),
            ..Default::default()
        };
        assert_eq!(cost(&folded)[0], 3, "interior meeting on {permutation:?}");
        assert_isometry(&input, &folded, &report);
        let fixed_wires =
            CliffordSimplificationPass::default().run(&input, &mut PropertySet::new());
        assert_eq!(cost(&fixed_wires)[0], 4);
        assert_isometry(&input, &fixed_wires, &TranspilationReport::default());
    }
}

#[test]
fn disjoint_reachable_intervals_and_barriers_block_the_meeting() {
    for permutation in all_four_wire_permutations() {
        let blocked = meeting_witness(permutation, true);
        let mut fenced = meeting_witness(permutation, false);
        // The public pass treats a barrier as a global optimization fence even
        // when its operand is disjoint from an endpoint's selected wire pair.
        fenced.operations.insert(
            1,
            Operation::Barrier {
                qubits: vec![permutation[2]],
            },
        );
        let mut near_blocker = meeting_witness(permutation, false);
        near_blocker.operations.insert(
            2,
            Operation::Gate {
                name: GateType::RZ,
                qubits: vec![permutation[1]],
                params: vec![1e-9],
            },
        );
        for input in [blocked, fenced, near_blocker] {
            let output = CliffordSimplificationPass {
                allow_swap_replacements: true,
            }
            .run(&input, &mut PropertySet::new());
            assert_eq!(output, input);
            assert_isometry(&input, &output, &TranspilationReport::default());
        }
    }
}

#[test]
fn interior_meeting_preserves_signed_frames_and_small_non_clifford_phases() {
    for p in all_four_wire_permutations() {
        for phase_gate in ["s", "sdg", "y", "rz(1e-9)"] {
            let mut input = meeting_witness(p, false);
            let prefix = source(4, &format!("{phase_gate} q[{}]; y q[{}];", p[0], p[1]));
            let suffix = source(4, &format!("x q[{}]; {phase_gate} q[{}];", p[1], p[0]));
            input.operations.splice(0..0, prefix.operations);
            input.operations.extend(suffix.operations);
            let mut properties = PropertySet::new();
            let output = CliffordSimplificationPass {
                allow_swap_replacements: true,
            }
            .run(&input, &mut properties);
            assert_isometry(&input, &output, &TranspilationReport::default());
            let folded = ExtractOutputPermutationPass.run(&output, &mut properties);
            let report = TranspilationReport {
                final_layout: properties
                    .get::<Vec<usize>>("logical_output_permutation")
                    .cloned(),
                ..Default::default()
            };
            assert_isometry(&input, &folded, &report);
        }
    }
}

#[test]
fn backward_local_transport_preserves_adjoint_order_and_pauli_signs() {
    // Local words adjacent to the right endpoint are crossed in reverse order,
    // with adjoints. Non-real S/Sdg/SX matrices make the wrong convention visible.
    let mut accepted = 0;
    for left in ["id", "h", "s", "sdg", "sx", "x", "y", "z"] {
        for right in ["id", "h", "s", "sdg", "sx", "x", "y", "z"] {
            for flipped in [false, true] {
                let [a, b, c, d] = if flipped { [3, 1, 2, 0] } else { [0, 1, 2, 3] };
                let input = source(
                    4,
                    &format!(
                        "cx q[{a}],q[{b}]; {left} q[{a}]; cx q[{c}],q[{b}];
                     cx q[{b}],q[{d}]; {right} q[{a}]; s q[{b}]; cx q[{b}],q[{a}];"
                    ),
                );
                let mut properties = PropertySet::new();
                let output = CliffordSimplificationPass {
                    allow_swap_replacements: true,
                }
                .run(&input, &mut properties);
                accepted += properties
                    .get::<CliffordSimplificationStats>("clifford_simplification_stats")
                    .unwrap()
                    .interior_rewrites;
                assert_isometry(&input, &output, &TranspilationReport::default());
                let folded = ExtractOutputPermutationPass.run(&output, &mut properties);
                let report = TranspilationReport {
                    final_layout: properties
                        .get::<Vec<usize>>("logical_output_permutation")
                        .cloned(),
                    ..Default::default()
                };
                assert_isometry(&input, &folded, &report);
            }
        }
    }
    assert!(
        accepted > 0,
        "the adjoint oracle must exercise accepted rewrites"
    );
}

#[test]
fn expanded_routing_swap_exposes_a_cross_block_cancellation() {
    for a in 0..3 {
        for b in 0..3 {
            if a == b {
                continue;
            }
            let c = 3 - a - b;
            let input = source(
                3,
                &format!("swap q[{a}],q[{b}]; cx q[{a}],q[{c}]; cx q[{a}],q[{b}];"),
            );
            for block_synthesis in [false, true] {
                let mut config = configuration(Backend::all_to_all(3), true, true, false);
                config.routing_search.elide_output_permutations = false;
                config.two_qubit_block_synthesis = block_synthesis;
                let mut disabled = config.clone();
                disabled.clifford_simplification = false;
                let (before, before_report) =
                    transpile_with_report(&input, Some(disabled)).unwrap();
                let (after, report) = transpile_with_report(&input, Some(config)).unwrap();
                // Expanding SWAP reveals A B A M A; A and M commute, so the
                // latter A pair disappears. A two-wire block cannot cross M.
                assert_eq!(cost(&before), [5, 5, 5]);
                assert_eq!(cost(&after), [3, 3, 3]);
                assert!(report.clifford_retention_stats.is_none());
                assert!(
                    report
                        .clifford_simplification_stats
                        .as_ref()
                        .unwrap()
                        .rewrites
                        > 0
                );
                assert_isometry(&input, &before, &before_report);
                assert_isometry(&input, &after, &report);
            }
        }
    }
}

#[test]
fn report_accumulates_optimization_and_post_expansion_clifford_rewrites() {
    let input = source(
        3,
        "cx q[0],q[1]; z q[1]; cx q[0],q[1]; barrier q;
        swap q[0],q[1]; cx q[0],q[2]; cx q[0],q[1];",
    );
    let mut config = configuration(Backend::all_to_all(3), true, true, false);
    config.routing_search.elide_output_permutations = false;
    let mut only_prefix = config.clone();
    only_prefix.post_routing_cleanup = false;
    let (prefix_output, prefix_report) = transpile_with_report(&input, Some(only_prefix)).unwrap();
    let (full_output, full_report) = transpile_with_report(&input, Some(config.clone())).unwrap();
    let prefix_stats = prefix_report
        .clifford_simplification_stats
        .as_ref()
        .unwrap();
    let full_stats = full_report.clifford_simplification_stats.as_ref().unwrap();
    assert_eq!((prefix_stats.rewrites, prefix_stats.cx_removed), (1, 2));
    assert_eq!((full_stats.rewrites, full_stats.cx_removed), (2, 4));
    assert_eq!(cost(&prefix_output)[0], 5);
    assert_eq!(cost(&full_output)[0], 3);
    assert_isometry(&input, &prefix_output, &prefix_report);
    assert_isometry(&input, &full_output, &full_report);
    assert_eq!(full_output, transpile(&input, Some(config)).unwrap());
}

#[test]
fn numerical_block_synthesis_keeps_a_gatewise_clifford_representation() {
    // This exact Clifford has a generic two-CX decomposition whose local
    // factors are not Cliffords. Rejecting that representation must retain the
    // input action without charging the numerical ledger for a rejected edit.
    for (a, b) in [(0, 1), (1, 0)] {
        let input = source(
            2,
            &format!(
                "cx q[{a}],q[{b}]; cx q[{b}],q[{a}]; h q[{b}]; s q[{a}];
                 cx q[{b}],q[{a}]; sdg q[{a}]; s q[{b}];"
            ),
        );
        for elide in [false, true] {
            let mut config = configuration(Backend::all_to_all(2), false, true, false);
            config.encoded_swap_recognition = false;
            config.routing_search.elide_output_permutations = elide;
            let mut disabled = config.clone();
            disabled.two_qubit_block_synthesis = false;
            let (baseline, baseline_report) =
                transpile_with_report(&input, Some(disabled)).unwrap();
            let (output, report) = transpile_with_report(&input, Some(config.clone())).unwrap();
            let stats = report.two_qubit_block_stats.as_ref().unwrap();
            assert_eq!(stats.rejected_non_clifford_representation, 1);
            assert_eq!(stats.accepted, 0);
            assert_eq!(stats.cx_removed, 0);
            assert_eq!(stats.residual_sum, 0.0);
            assert_eq!(stats.max_residual, 0.0);
            assert_eq!(cost(&output)[0], 3);
            assert_eq!(output, baseline);
            assert_eq!(report.initial_layout, baseline_report.initial_layout);
            assert_eq!(report.final_layout, baseline_report.final_layout);
            assert_eq!(report.numerical_budget, baseline_report.numerical_budget);
            assert_isometry(&input, &output, &report);
            for operation in &output.operations {
                match operation {
                    Operation::Gate {
                        name: GateType::CX, ..
                    } => (),
                    Operation::Gate { qubits, .. } if qubits.len() == 1 => {
                        assert!(local_gate_preserves_paulis(operation), "{operation:?}");
                    }
                    _ => panic!("unexpected output operation {operation:?}"),
                }
            }
            assert_eq!(output, transpile(&input, Some(config.clone())).unwrap());
            config.profile_passes = true;
            let (profiled, mut profiled_report) =
                transpile_with_report(&input, Some(config)).unwrap();
            assert_eq!(output, profiled);
            assert!(!profiled_report.pass_diagnostics.is_empty());
            profiled_report.pass_diagnostics.clear();
            assert_eq!(report, profiled_report);
        }
    }
}
