//! CU1 checks against the analytic controlled-phase action. The oracle below
//! does not use Q-Rust's gate definitions, decomposition, simulator or verifier.

use num_complex::Complex64 as C;
use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateDefinition, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::simulator::try_circuit_to_unitary;
use q_rust::transpiler::{transpile, transpile_with_report, TranspilerConfig};
use std::f64::consts::PI;

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
        let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        else {
            panic!("unexpected oracle operation {op:?}")
        };
        match name.to_qasm_name() {
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
fn qelib1_cu1_parses_symbolic_angles_and_round_trips() {
    let circuit = source("cu1(-pi/3) q[1],q[0]; cu1(2*pi) q[0],q[1];", 2);
    assert_eq!(circuit.operations.len(), 2);
    let Operation::Gate {
        name,
        qubits,
        params,
    } = &circuit.operations[0]
    else {
        panic!("expected gate")
    };
    assert_eq!(name.to_qasm_name(), "cu1");
    assert_eq!(name.static_qasm_name(), Some("cu1"));
    assert_eq!(name.num_qubits(), 2);
    assert_eq!(qubits, &[1, 0]);
    assert_eq!(params, &[-PI / 3.0]);
    let serialized = circuit.to_qasm(None);
    assert!(serialized.contains("cu1("));
    assert_eq!(
        parse_qasm(&serialized).unwrap().operations,
        circuit.operations
    );
}

#[test]
fn cu1_unitary_and_decomposition_match_the_analytic_diagonal() {
    let gate = "cu1".parse::<GateType>().unwrap();
    for theta in [0.0, 1e-9, -1e-9, 0.73, -0.73, PI, -PI, 2.0 * PI, 1e16] {
        let unitary = gate.unitary(&[theta]);
        assert_eq!(unitary.shape(), (4, 4));
        for row in 0..4 {
            for col in 0..4 {
                let expected = if row != col {
                    C::new(0.0, 0.0)
                } else if row == 3 {
                    phase(theta)
                } else {
                    C::new(1.0, 0.0)
                };
                assert!((unitary[(row, col)] - expected).norm() < 1e-12);
            }
        }
        for wires in [[0, 2], [2, 0]] {
            let original = source(&format!("cu1({theta}) q[{}],q[{}];", wires[0], wires[1]), 3);
            let mut decomposed = Circuit::new(3, 0);
            decomposed.operations = gate.decompose(&wires, &[theta]).expect("basis expansion");
            let error = full_isometry_error(&original, &decomposed, &[0, 1, 2], &[0, 1, 2]);
            assert!(error < 1e-12, "angle {theta}, residual {error:e}");
        }
    }
}

#[test]
fn cu1_preserves_the_control_phase_that_a_crz_alias_would_lose() {
    let theta = 0.73;
    let circuit = source(&format!("cu1({theta}) q[0],q[1];"), 2);
    let actual = try_circuit_to_unitary(&circuit).unwrap();
    assert!((actual[(1, 1)] - C::new(1.0, 0.0)).norm() < 1e-12);
    assert!((actual[(3, 3)] - phase(theta)).norm() < 1e-12);
    let crz = try_circuit_to_unitary(&source(&format!("crz({theta}) q[0],q[1];"), 2)).unwrap();
    // Shared |00> entry fixes the common phase. This ratio distinguishes CU1
    // from CRZ even when every computational-basis probability is identical.
    assert!((actual[(1, 1)] / actual[(0, 0)] - crz[(1, 1)] / crz[(0, 0)]).norm() > 0.3);
}

#[test]
fn cu1_register_broadcast_and_custom_symbolic_body_are_supported() {
    let circuit = parse_qasm(
        "OPENQASM 2.0; include \"qelib1.inc\";
         gate wrapped(a) control,target { cu1(-a/2) target,control; }
         qreg a[2]; qreg b[2];
         cu1(pi/7) a,b; wrapped(pi/3) a[1],b[0];",
    )
    .unwrap();
    let expected = source(
        "cu1(pi/7) q[0],q[2]; cu1(pi/7) q[1],q[3]; cu1(-pi/6) q[2],q[1];",
        4,
    );
    for level in 0..=3 {
        let output = transpile(
            &circuit,
            Some(
                TranspilerConfig::builder()
                    .optimization_level(level)
                    .target_basis(["u", "cx"])
                    .build(),
            ),
        )
        .unwrap();
        assert!(full_isometry_error(&expected, &output, &[0, 1, 2, 3], &[0, 1, 2, 3]) < 1e-8);
    }
}

#[test]
fn cu1_routing_preserves_all_inputs_and_idle_physical_wires() {
    for wires in [[0, 2], [2, 0]] {
        let circuit = source(
            &format!(
                "u(0.3,0.7,-0.2) q[0]; u(-0.8,0.1,0.5) q[2];
                 cu1(-pi/3) q[{}],q[{}]; cx q[1],q[0];
                 u(0.6,-0.4,0.8) q[2]; cu1(pi/7) q[2],q[1];
                 cu1(2*pi) q[0],q[2]; u(0.9,0.2,-0.1) q[1];",
                wires[0], wires[1]
            ),
            3,
        );
        for level in 0..=3 {
            for backend in [
                Backend::all_to_all(3),
                Backend::linear(3),
                Backend::linear(4),
            ] {
                let config = TranspilerConfig::builder()
                    .optimization_level(level)
                    .backend(backend.clone())
                    .target_basis(["u", "cx"])
                    .build();
                let (output, report) = transpile_with_report(&circuit, Some(config)).unwrap();
                for op in &output.operations {
                    let Operation::Gate { name, qubits, .. } = op else {
                        panic!("unexpected operation {op:?}")
                    };
                    assert!(matches!(name, GateType::U | GateType::CX));
                    if name == &GateType::CX {
                        assert!(backend.has_directed_edge(qubits[0], qubits[1]));
                    }
                }
                let error = full_isometry_error(
                    &circuit,
                    &output,
                    &report.initial_layout.unwrap()[..3],
                    &report.final_layout.unwrap()[..3],
                );
                assert!(
                    error < 1e-8,
                    "level {level}, {}, residual {error:e}",
                    backend.name
                );
            }
        }
    }
}

#[test]
fn cu1_full_turn_is_identity_but_half_turn_retains_relative_phase() {
    for theta in [2.0 * PI, -2.0 * PI, PI, -PI] {
        let circuit = source(&format!("cu1({theta}) q[0],q[1];"), 2);
        for level in 0..=3 {
            let output = transpile(
                &circuit,
                Some(
                    TranspilerConfig::builder()
                        .optimization_level(level)
                        .target_basis(["u", "cx"])
                        .build(),
                ),
            )
            .unwrap();
            assert!(full_isometry_error(&circuit, &output, &[0, 1], &[0, 1]) < 1e-8);
        }
    }
}

#[test]
fn invalid_cu1_shapes_and_nonfinite_angles_are_rejected_before_execution() {
    for body in [
        "cu1 q[0],q[1];",
        "cu1(0.1,0.2) q[0],q[1];",
        "cu1(0.1) q[0];",
        "cu1(0.1) q[0],q[1],q[2];",
        "cu1(0.1) q[1],q[1];",
    ] {
        // Built-in gate arity is checked at the public execution boundary,
        // just as for existing RX/CRZ gates, rather than during parsing.
        let circuit = source(body, 3);
        assert!(transpile(&circuit, None).is_err(), "accepted {body}");
        assert!(
            try_circuit_to_unitary(&circuit).is_err(),
            "simulated {body}"
        );
    }
    for theta in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut circuit = source("cu1(0.1) q[0],q[1];", 2);
        let Operation::Gate { params, .. } = &mut circuit.operations[0] else {
            unreachable!()
        };
        params[0] = theta;
        assert!(transpile(&circuit, None).is_err());
        assert!(try_circuit_to_unitary(&circuit).is_err());
    }
}
