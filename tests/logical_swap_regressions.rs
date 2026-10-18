//! Logical SWAP absorption against independent analytic state evolution. The oracle below
//! does not use Q-Rust's gate definitions, decomposition, simulator or verifier.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::optimization::{ApplyOutputPermutationPass, ExtractOutputPermutationPass};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::{transpile_with_report, TranspilerConfig};

fn phase(theta: f64) -> C {
    C::from_polar(1.0, theta)
}

fn source(body: &str, width: usize) -> Circuit {
    parse_qasm(&format!(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[{width}]; {body}"
    ))
    .expect("qelib1 cu1 must parse")
}

fn evolve(circuit: &Circuit, mut state: Vec<C>) -> Vec<C> {
    for op in &circuit.operations {
        if matches!(op, Operation::Barrier { .. }) {
            continue;
        }
        let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        else {
            panic!("unexpected oracle operation {op:?}")
        };
        match name.to_qasm_name() {
            "swap" => {
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) == 0 && basis & (1 << qubits[1]) != 0 {
                        state.swap(basis, basis ^ (1 << qubits[0]) ^ (1 << qubits[1]));
                    }
                }
            }
            "cu1" => {
                // The only affected basis state has both operands equal to 1.
                // In particular, |01> does not acquire CRZ's control phase.
                for (basis, amplitude) in state.iter_mut().enumerate() {
                    if qubits.iter().all(|q| basis & (1 << q) != 0) {
                        *amplitude *= phase(params[0]);
                    }
                }
            }
            "cx" => {
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) != 0 && basis & (1 << qubits[1]) == 0 {
                        state.swap(basis, basis | (1 << qubits[1]));
                    }
                }
            }
            "ccx" => {
                for basis in 0..state.len() {
                    if basis & (1 << qubits[0]) != 0
                        && basis & (1 << qubits[1]) != 0
                        && basis & (1 << qubits[2]) == 0
                    {
                        state.swap(basis, basis | (1 << qubits[2]));
                    }
                }
            }
            "u" => {
                let (s, c) = (params[0] / 2.0).sin_cos();
                let matrix = [
                    [C::new(c, 0.0), -phase(params[2]) * s],
                    [phase(params[1]) * s, phase(params[1] + params[2]) * c],
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
            other => panic!("unsupported independent oracle gate {other}"),
        }
    }
    state
}

fn embed(basis: usize, layout: &[usize]) -> usize {
    layout
        .iter()
        .enumerate()
        .map(|(logical, physical)| ((basis >> logical) & 1) << physical)
        .sum()
}

fn full_isometry_error(
    expected: &Circuit,
    actual: &Circuit,
    initial: &[usize],
    final_map: &[usize],
) -> f64 {
    let logical_size = 1 << expected.num_qubits;
    let physical_size = 1 << actual.num_qubits;
    let mut columns = Vec::new();
    for basis in 0..logical_size {
        let mut input = vec![C::new(0.0, 0.0); logical_size];
        input[basis] = C::new(1.0, 0.0);
        let expected_column = evolve(expected, input);
        let mut physical_input = vec![C::new(0.0, 0.0); physical_size];
        physical_input[embed(basis, initial)] = C::new(1.0, 0.0);
        let actual_column = evolve(actual, physical_input);
        let mut physical_expected = vec![C::new(0.0, 0.0); physical_size];
        for (row, amplitude) in expected_column.into_iter().enumerate() {
            physical_expected[embed(row, final_map)] = amplitude;
        }
        columns.push((physical_expected, actual_column));
    }
    // One phase for all inputs, never a separate phase for each basis column.
    let overlap: C = columns
        .iter()
        .flat_map(|(a, b)| a.iter().zip(b).map(|(a, b)| a.conj() * b))
        .sum();
    let shared_phase = if overlap.norm() > 0.0 {
        overlap / overlap.norm()
    } else {
        C::new(1.0, 0.0)
    };
    columns
        .iter()
        .map(|(a, b)| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (*b - shared_phase * a).norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max)
}

#[test]
fn interior_three_cycle_relabels_noncommuting_suffix_and_barrier() {
    let original = source(
        "u(0.7,-0.2,0.4) q[0]; swap q[0],q[1];
        u(-0.8,0.3,0.6) q[0]; cx q[0],q[2]; swap q[1],q[2];
        barrier q[0],q[2]; cu1(-0.43) q[2],q[0]; u(0.1,0.2,-0.7) q[1];",
        3,
    );
    let mut properties = PropertySet::new();
    let output = ExtractOutputPermutationPass.run(&original, &mut properties);
    let permutation = properties
        .get::<Vec<usize>>("logical_output_permutation")
        .unwrap();
    assert_eq!(permutation, &vec![1, 2, 0]); // Deliberately not its own inverse.
    assert_eq!(output.operations.len(), original.operations.len() - 2);
    assert!(output
        .operations
        .iter()
        .any(|op| matches!(op,Operation::Barrier{qubits} if qubits==&vec![1,0])));
    assert!(full_isometry_error(&original, &output, &[0, 1, 2], permutation) < 1e-12);
    properties.insert("final_layout", vec![4usize, 0, 2]);
    ApplyOutputPermutationPass.run(&output, &mut properties);
    assert_eq!(
        properties.get::<Vec<usize>>("final_layout"),
        Some(&vec![0, 2, 4])
    );
    assert!(properties
        .get::<Vec<usize>>("logical_output_permutation")
        .is_none());
}

#[test]
fn generated_swap_networks_preserve_every_input_with_one_shared_phase() {
    for width in 3..=6 {
        for seed in 0..6usize {
            let mut body = String::new();
            for i in 0..15usize {
                let a = (i * 7 + seed) % width;
                let b = (a + 1 + (i + seed) % (width - 1)) % width;
                let theta = if i % 3 == 0 {
                    1e-9
                } else {
                    -0.37 + (i as f64) * 0.07
                };
                body += &format!("swap q[{a}],q[{b}]; u({theta},0.2,-0.1) q[{a}];");
                if i % 2 == 0 {
                    body += &format!("cu1({theta}) q[{b}],q[{a}];");
                } else {
                    body += &format!("cx q[{a}],q[{b}];");
                }
                if i % 4 == 0 {
                    body += &format!("barrier q[{a}];");
                }
            }
            let original = source(&body, width);
            let mut properties = PropertySet::new();
            let output = ExtractOutputPermutationPass.run(&original, &mut properties);
            let layout = (0..width).collect::<Vec<_>>();
            let final_map = properties
                .get::<Vec<usize>>("logical_output_permutation")
                .unwrap();
            assert!(full_isometry_error(&original, &output, &layout, final_map) < 1e-12);
        }
    }
}

#[test]
fn routed_interior_permutations_preserve_inputs_and_clean_idle_wires() {
    let original = source(
        "u(0.7,-0.2,0.4) q[0]; swap q[0],q[2];
        u(-0.8,0.3,0.6) q[0]; cx q[0],q[1]; swap q[1],q[2];
        barrier q[2]; cu1(-0.43) q[2],q[0]; u(0.1,0.2,-0.7) q[1];",
        3,
    );
    let mut directed = Backend::new("directed_line", 4);
    directed.set_coupling_map([(0, 1), (1, 2), (2, 3)]);
    for backend in [
        Backend::all_to_all(3),
        Backend::linear(3),
        Backend::linear(4),
        directed,
    ] {
        for level in 0..=3 {
            let config = TranspilerConfig::builder()
                .backend(backend.clone())
                .optimization_level(level)
                .target_basis(["u", "cx"])
                .build();
            let (output, report) = transpile_with_report(&original, Some(config)).unwrap();
            for operation in &output.operations {
                if let Operation::Gate {
                    name: GateType::CX,
                    qubits,
                    ..
                } = operation
                {
                    assert!(backend.has_directed_edge(qubits[0], qubits[1]));
                }
            }
            let error = full_isometry_error(
                &original,
                &output,
                &report.initial_layout.unwrap()[..3],
                &report.final_layout.unwrap()[..3],
            );
            assert!(error < 1e-8, "level {level}, residual {error:e}");
        }
    }
}

#[test]
fn internal_wire_permutation_does_not_require_cx_for_local_suffixes() {
    let original = source(
        "swap q[0],q[1]; u(0.7,0.2,-0.1) q[0];
        swap q[1],q[2]; u(-0.4,0.6,0.3) q[1]; u(0.2,0.9,-0.5) q[2];",
        3,
    );
    let config = TranspilerConfig::builder()
        .backend(Backend::all_to_all(3))
        .optimization_level(3)
        .target_basis(["u", "cx"])
        .build();
    let (output, report) = transpile_with_report(&original, Some(config)).unwrap();
    assert!(output
        .operations
        .iter()
        .all(|op| matches!(op,Operation::Gate{qubits,..} if qubits.len()==1)));
    assert!(
        full_isometry_error(
            &original,
            &output,
            &report.initial_layout.unwrap(),
            &report.final_layout.unwrap()
        ) < 1e-8
    );
}

#[test]
fn terminal_readout_uses_composed_mapping_and_preserves_write_order() {
    let circuit = parse_qasm(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[3]; creg c[2];
        swap q[0],q[1]; u(0.3,0.1,-0.2) q[0]; swap q[1],q[2]; cx q[0],q[2];
        measure q[0] -> c[1]; measure q[2] -> c[1]; measure q[1] -> c[0];",
    )
    .unwrap();
    let config = TranspilerConfig::builder()
        .backend(Backend::linear(4))
        .target_basis(["u", "cx"])
        .optimization_level(3)
        .build();
    let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
    let final_map = report.final_layout.unwrap();
    let measured = output
        .operations
        .iter()
        .filter_map(|op| {
            if let Operation::Measure { qubit, cbit } = op {
                Some((*qubit, *cbit))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    assert_eq!(
        measured,
        vec![(final_map[0], 1), (final_map[2], 1), (final_map[1], 0)]
    );
}

#[test]
fn custom_gate_operands_are_relabelled_and_direct_nonunitary_calls_are_conservative() {
    let circuit = parse_qasm(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[3];
        gate local(a) x { u(a,0.2,-0.4) x; }
        swap q[0],q[2]; local(-0.31) q[0]; local(0.7) q[2];",
    )
    .unwrap();
    let expected = source(
        "swap q[0],q[2]; u(-0.31,0.2,-0.4) q[0]; u(0.7,0.2,-0.4) q[2];",
        3,
    );
    let config = TranspilerConfig::builder()
        .backend(Backend::linear(4))
        .target_basis(["u", "cx"])
        .optimization_level(3)
        .build();
    let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
    assert!(
        full_isometry_error(
            &expected,
            &output,
            &report.initial_layout.unwrap()[..3],
            &report.final_layout.unwrap()[..3]
        ) < 1e-8
    );
    for suffix in [
        Operation::Reset { qubit: 1 },
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 9],
            params: vec![],
        },
    ] {
        let mut invalid = source("swap q[0],q[1]; u(0.3,0.1,-0.2) q[0];", 2);
        invalid.num_cbits = 1;
        invalid.operations.push(suffix);
        let mut properties = PropertySet::new();
        properties.insert("logical_output_permutation", vec![1usize, 0]);
        let output = ExtractOutputPermutationPass.run(&invalid, &mut properties);
        assert_eq!(output.operations, invalid.operations);
        assert!(properties
            .get::<Vec<usize>>("logical_output_permutation")
            .is_none());
    }
}

#[test]
fn internal_swaps_before_high_arity_and_custom_two_qubit_gates_preserve_all_inputs() {
    let circuit = parse_qasm(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[4];
        gate bridge(a) c,t { u(a,0.2,-0.1) c; cx c,t; cu1(a) t,c; }
        u(0.7,0.3,-0.2) q[0]; swap q[0],q[3]; ccx q[0],q[1],q[2];
        bridge(-0.31) q[3],q[1]; swap q[1],q[2];
        ccx q[2],q[3],q[0]; bridge(0.43) q[0],q[2];",
    )
    .unwrap();
    let expected = source(
        "u(0.7,0.3,-0.2) q[0]; swap q[0],q[3]; ccx q[0],q[1],q[2];
        u(-0.31,0.2,-0.1) q[3]; cx q[3],q[1]; cu1(-0.31) q[1],q[3];
        swap q[1],q[2]; ccx q[2],q[3],q[0];
        u(0.43,0.2,-0.1) q[0]; cx q[0],q[2]; cu1(0.43) q[2],q[0];",
        4,
    );
    let backend = Backend::linear(5);
    for level in 0..=3 {
        let config = TranspilerConfig::builder()
            .backend(backend.clone())
            .target_basis(["u", "cx"])
            .optimization_level(level)
            .build();
        let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
        for operation in &output.operations {
            if let Operation::Gate {
                name: GateType::CX,
                qubits,
                ..
            } = operation
            {
                assert!(backend.has_directed_edge(qubits[0], qubits[1]));
            }
        }
        let error = full_isometry_error(
            &expected,
            &output,
            &report.initial_layout.unwrap()[..4],
            &report.final_layout.unwrap()[..4],
        );
        assert!(error < 1e-8, "level {level}, residual {error:e}");
    }
}
