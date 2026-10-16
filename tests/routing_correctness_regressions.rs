//! Regressions for routing failure propagation, dynamic-operation loss, and
//! hardware legality after basis lowering.

use q_rust::backend::Backend;
use q_rust::ir::{Circuit, ClassicalCondition, GateType, Operation};
use q_rust::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::BeamSabrePass;
use q_rust::transpiler::{transpile, transpile_with_report, TranspilerConfig};

fn gate(name: GateType, qubits: &[usize]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: vec![],
    }
}

fn config(backend: Backend) -> TranspilerConfig {
    TranspilerConfig::builder()
        .optimization_level(0)
        .decompose_basis(false)
        .backend(backend)
        .build()
}

#[test]
fn conditional_is_rejected_before_it_can_be_dropped() {
    let mut c = Circuit::new(3, 1);
    c.add_op(Operation::Conditional {
        condition: ClassicalCondition {
            creg: "c".into(),
            value: 0,
        },
        op: Box::new(gate(GateType::X, &[0])),
    });
    c.add_op(gate(GateType::CX, &[0, 1]));
    for level in 0..=3 {
        for backend in [None, Some(Backend::linear(3)), Some(Backend::all_to_all(3))] {
            let mut cfg = TranspilerConfig::builder()
                .optimization_level(level)
                .build();
            cfg.backend = backend;
            assert!(transpile(&c, Some(cfg.clone())).is_err());
            assert!(transpile_with_report(&c, Some(cfg)).is_err());
        }
    }
    let pass = BeamSabrePass {
        backend: Backend::linear(3),
        ..Default::default()
    };
    assert!(pass.try_run(&c, &mut PropertySet::new()).is_err());
}

#[test]
fn mid_circuit_measurement_is_rejected_instead_of_reordered() {
    let mut c = Circuit::new(3, 1);
    c.add_op(gate(GateType::X, &[0]));
    c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    c.add_op(gate(GateType::CX, &[0, 1]));
    let cfg = config(Backend::linear(3));
    assert!(transpile(&c, Some(cfg.clone())).is_err());
    assert!(transpile_with_report(&c, Some(cfg)).is_err());
}

#[test]
fn reset_is_rejected_before_optimization_or_routing() {
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::X, &[0]));
    c.add_op(Operation::Reset { qubit: 0 });
    c.add_op(gate(GateType::CX, &[0, 1]));
    assert!(transpile(&c, Some(config(Backend::linear(3)))).is_err());
    assert!(transpile(&c, None).is_err());
}

#[test]
fn disconnected_failure_reaches_both_public_entry_points() {
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::CX, &[0, 2]));
    for decompose in [false, true] {
        for level in 0..=3 {
            let mut cfg = config(Backend::new("disconnected", 3));
            cfg.decompose_basis = decompose;
            cfg.optimization_level = level;
            assert!(transpile(&c, Some(cfg.clone())).is_err());
            assert!(transpile_with_report(&c, Some(cfg)).is_err());
        }
    }
}

#[test]
fn insufficient_width_is_rejected_even_without_two_qubit_gates() {
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::H, &[2]));
    let cfg = config(Backend::linear(2));
    assert!(transpile(&c, Some(cfg.clone())).is_err());
    assert!(transpile_with_report(&c, Some(cfg)).is_err());
}

#[test]
fn terminal_measurements_follow_all_gates_and_keep_classical_mapping() {
    let mut c = Circuit::new(3, 3);
    c.add_op(gate(GateType::X, &[0]));
    c.add_op(gate(GateType::CX, &[0, 2]));
    c.add_op(gate(GateType::H, &[1]));
    for (qubit, cbit) in [(0, 2), (2, 0), (1, 1)] {
        c.add_op(Operation::Measure { qubit, cbit });
    }
    let (mut out, report) = transpile_with_report(&c, Some(config(Backend::linear(4)))).unwrap();
    let final_layout = report.final_layout.unwrap();
    let expected: Vec<_> = [(0, 2), (2, 0), (1, 1)]
        .iter()
        .map(|&(q, cbit)| Operation::Measure {
            qubit: final_layout[q],
            cbit,
        })
        .collect();
    assert_eq!(
        &out.operations[out.operations.len() - 3..],
        expected.as_slice()
    );
    c.operations.truncate(c.operations.len() - 3);
    out.operations.truncate(out.operations.len() - 3);
    let actual = extract_logical_unitary(
        &circuit_to_unitary(&out),
        3,
        &report.initial_layout.unwrap(),
        &final_layout,
    );
    assert!(unitary_fidelity(&circuit_to_unitary(&c), &actual) > 1.0 - 1e-10);
}

#[test]
fn directed_swap_lowering_finishes_legal_and_in_basis() {
    let bases = [
        vec!["u", "cx"],
        vec!["rz", "rx", "cx"],
        vec!["rz", "sx", "cx"],
        vec!["rz", "h", "cx"],
        vec!["u", "cz"],
    ];
    for basis in bases {
        for level in 0..=3 {
            let mut backend = Backend::new("oneway", 2);
            backend.set_coupling_map([(0, 1)]);
            for name in &basis {
                backend.add_basis_gate(*name);
            }
            let mut c = Circuit::new(2, 0);
            c.add_op(gate(GateType::SWAP, &[0, 1]));
            c.add_op(gate(GateType::H, &[0]));
            let cfg = TranspilerConfig::builder()
                .optimization_level(level)
                .backend(backend.clone())
                .build();
            let (out, report) = transpile_with_report(&c, Some(cfg)).unwrap();
            for op in &out.operations {
                if let Operation::Gate { name, qubits, .. } = op {
                    assert!(
                        basis.contains(&name.to_qasm_name()),
                        "{op:?} escaped {basis:?}"
                    );
                    if name == &GateType::CX {
                        assert!(backend.has_directed_edge(qubits[0], qubits[1]), "{op:?}");
                    }
                }
            }
            let actual = extract_logical_unitary(
                &circuit_to_unitary(&out),
                2,
                &report.initial_layout.unwrap(),
                &report.final_layout.unwrap(),
            );
            assert!(unitary_fidelity(&circuit_to_unitary(&c), &actual) > 1.0 - 1e-9);
        }
    }
}

#[test]
fn directed_sparse_routing_can_traverse_reverse_edges() {
    let mut backend = Backend::new("inward", 3);
    backend.set_coupling_map([(1, 0), (1, 2)]);
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::CX, &[0, 2]));
    let cfg = TranspilerConfig::builder()
        .optimization_level(0)
        .backend(backend.clone())
        .target_basis(["u", "cx"])
        .build();
    let (out, report) = transpile_with_report(&c, Some(cfg)).expect("undirected route exists");
    for op in &out.operations {
        if let Operation::Gate {
            name: GateType::CX,
            qubits,
            ..
        } = op
        {
            assert!(backend.has_directed_edge(qubits[0], qubits[1]));
        }
    }
    let actual = extract_logical_unitary(
        &circuit_to_unitary(&out),
        3,
        &report.initial_layout.unwrap(),
        &report.final_layout.unwrap(),
    );
    assert!(unitary_fidelity(&circuit_to_unitary(&c), &actual) > 1.0 - 1e-10);
}

#[test]
fn disabled_decomposition_cannot_return_gates_outside_backend_basis() {
    let mut backend = Backend::all_to_all(2);
    backend.add_basis_gate("u");
    backend.add_basis_gate("cx");
    let mut c = Circuit::new(2, 0);
    c.add_op(gate(GateType::H, &[0]));
    assert!(transpile(&c, Some(config(backend))).is_err());
}

#[test]
fn unresolved_high_arity_custom_gate_is_not_truncated_by_router() {
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::Custom("missing".into()), &[0, 1, 2]));
    assert!(transpile(&c, Some(config(Backend::linear(3)))).is_err());
}

#[test]
fn barriers_preserve_operation_order_and_layout_between_segments() {
    let mut c = Circuit::new(3, 1);
    c.add_op(gate(GateType::X, &[0]));
    c.add_op(Operation::Barrier { qubits: vec![0] });
    c.add_op(gate(GateType::CX, &[0, 2]));
    c.add_op(Operation::Barrier { qubits: vec![0, 2] });
    c.add_op(gate(GateType::H, &[1]));
    c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    let pass = BeamSabrePass {
        backend: Backend::linear(4),
        ..Default::default()
    };
    let mut properties = PropertySet::new();
    // A supplied placement also exercises the one-qubit-only first segment.
    properties.insert("sabre_initial_layout", vec![3usize, 1, 0]);
    let mut out = pass.try_run(&c, &mut properties).unwrap();
    assert_eq!(out.operations[0], gate(GateType::X, &[3]));
    assert_eq!(out.operations[1], Operation::Barrier { qubits: vec![3] });
    assert_eq!(
        out.operations.iter().filter(|op| op.is_barrier()).count(),
        2
    );
    let initial = properties.get::<Vec<usize>>("initial_layout").unwrap();
    let final_layout = properties.get::<Vec<usize>>("final_layout").unwrap();
    assert_eq!(
        out.operations.last(),
        Some(&Operation::Measure {
            qubit: final_layout[0],
            cbit: 0
        })
    );
    c.operations.pop();
    out.operations.pop();
    let actual = extract_logical_unitary(&circuit_to_unitary(&out), 3, initial, final_layout);
    assert!(unitary_fidelity(&circuit_to_unitary(&c), &actual) > 1.0 - 1e-10);
}

#[test]
fn requested_layout_is_respected_by_fast_paths() {
    for backend in [Backend::linear(4), Backend::all_to_all(4)] {
        let mut c = Circuit::new(2, 1);
        c.add_op(gate(GateType::X, &[0]));
        c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        let pass = BeamSabrePass {
            backend,
            ..Default::default()
        };
        let mut properties = PropertySet::new();
        properties.insert("sabre_initial_layout", vec![3usize, 2]);
        let out = pass.try_run(&c, &mut properties).unwrap();
        assert_eq!(out.num_qubits, 4);
        assert_eq!(out.operations[0], gate(GateType::X, &[3]));
        assert_eq!(out.operations[1], Operation::Measure { qubit: 3, cbit: 0 });
        assert_eq!(
            properties.get::<Vec<usize>>("initial_layout"),
            Some(&vec![3, 2])
        );
        assert_eq!(
            properties.get::<Vec<usize>>("final_layout"),
            Some(&vec![3, 2])
        );
    }
}

#[test]
fn malformed_layout_is_rejected_without_indexing_panic() {
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::CX, &[0, 2]));
    let pass = BeamSabrePass {
        backend: Backend::linear(3),
        ..Default::default()
    };
    for mapping in [vec![0usize], vec![0, 0, 1], vec![0, 1, 3]] {
        let mut properties = PropertySet::new();
        properties.insert("sabre_initial_layout", mapping);
        assert!(pass.try_run(&c, &mut properties).is_err());
    }
}

#[test]
fn malformed_ir_is_rejected_before_passes_can_panic_or_default_parameters() {
    let bad = [
        gate(GateType::CX, &[0]),
        gate(GateType::CX, &[0, 0]),
        gate(GateType::H, &[3]),
        gate(GateType::RX, &[0]),
        Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![f64::NAN],
        },
        Operation::Measure { qubit: 0, cbit: 1 },
    ];
    for op in bad {
        let mut c = Circuit::new(3, 1);
        c.add_op(op);
        assert!(transpile(&c, Some(config(Backend::linear(3)))).is_err());
    }
}

#[test]
fn checked_pass_manager_stops_at_error() {
    use q_rust::transpiler::pass::{Pass, PassManager};
    struct MustNotRun;
    impl Pass for MustNotRun {
        fn name(&self) -> &str {
            "MustNotRun"
        }
        fn run(&self, _: &Circuit, _: &mut PropertySet) -> Circuit {
            panic!("a pass ran after routing failed");
        }
    }
    let mut c = Circuit::new(3, 0);
    c.add_op(gate(GateType::CX, &[0, 2]));
    let mut pm = PassManager::new();
    pm.add_pass(Box::new(BeamSabrePass {
        backend: Backend::new("disconnected", 3),
        ..Default::default()
    }));
    pm.add_pass(Box::new(MustNotRun));
    assert!(pm.try_run(&c).is_err());
}

#[test]
fn decomposition_failure_reaches_caller_without_target_basis() {
    let mut c = Circuit::new(1, 0);
    c.add_op(gate(GateType::Custom("missing".into()), &[0]));
    assert!(transpile(&c, None).is_err());
    assert!(transpile_with_report(&c, None).is_err());
}

#[test]
fn recursive_custom_gate_returns_error_instead_of_overflowing_stack() {
    use q_rust::ir::ast::ParsedStatement;
    let mut c = Circuit::new(1, 0);
    c.custom_gates.register(
        "recursive".into(),
        vec![],
        vec!["a".into()],
        vec![ParsedStatement::Gate(
            "recursive".into(),
            vec![("a".into(), None)],
            vec![],
        )],
    );
    c.add_op(gate(GateType::Custom("recursive".into()), &[0]));
    assert!(transpile(&c, None).is_err());
}

#[test]
fn optimization_keeps_terminal_readout_order_and_records_elided_swaps() {
    use q_rust::verify::{verify_equivalence, verify_equivalence_with_layout};
    let mut c = Circuit::new(3, 2);
    c.add_op(gate(GateType::H, &[0]));
    c.add_op(gate(GateType::CX, &[0, 1]));
    c.add_op(gate(GateType::H, &[2]));
    c.add_op(gate(GateType::SWAP, &[1, 2]));
    // Repeated writes to one classical bit deliberately pin their order.
    c.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    c.add_op(Operation::Measure { qubit: 1, cbit: 0 });
    c.add_op(Operation::Barrier {
        qubits: vec![0, 1, 2],
    });
    c.add_op(Operation::Measure { qubit: 2, cbit: 1 });
    for level in 0..=3 {
        for backend in [None, Some(Backend::all_to_all(3)), Some(Backend::linear(4))] {
            let mut cfg = TranspilerConfig::builder()
                .optimization_level(level)
                .build();
            cfg.backend = backend;
            let (mut out, report) = transpile_with_report(&c, Some(cfg)).unwrap();
            let first_measure = out
                .operations
                .iter()
                .position(|op| matches!(op, Operation::Measure { .. }))
                .unwrap();
            assert!(out.operations[first_measure..]
                .iter()
                .all(|op| matches!(op, Operation::Measure { .. } | Operation::Barrier { .. })));
            let mapping = report.final_layout.clone().unwrap_or_else(|| vec![0, 1, 2]);
            assert_eq!(
                &out.operations[first_measure..],
                &[
                    Operation::Measure {
                        qubit: mapping[0],
                        cbit: 0
                    },
                    Operation::Measure {
                        qubit: mapping[1],
                        cbit: 0
                    },
                    Operation::Barrier {
                        qubits: mapping.clone()
                    },
                    Operation::Measure {
                        qubit: mapping[2],
                        cbit: 1
                    },
                ]
            );
            // Readout mapping/order was checked above. Verify the unitary
            // prefix separately because repeated classical writes are outside
            // the high-level verifier's deliberately narrower domain.
            out.operations.truncate(first_measure);
            let mut source_prefix = c.clone();
            source_prefix.operations.truncate(4);
            let verdict = match (&report.initial_layout, &report.final_layout) {
                (Some(initial), Some(final_layout)) => {
                    verify_equivalence_with_layout(&source_prefix, &out, initial, final_layout)
                }
                _ => verify_equivalence(&source_prefix, &out),
            }
            .unwrap();
            assert!(verdict.is_equivalent(), "{verdict:?}");
        }
    }
}
