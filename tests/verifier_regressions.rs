//! Adversarial checks with analytic expectations, independent of decomposition.
use nalgebra::{DMatrix, DVector};
use num_complex::Complex64 as C;
use q_rust::ir::{Circuit, ClassicalCondition, GateType, Operation};
use q_rust::simulator::{
    equivalence_by_sampling, equivalence_by_sampling_with_layout, evolve_state,
    extract_logical_unitary, try_circuit_to_unitary, try_extract_logical_unitary, unitary_fidelity,
};
use q_rust::verify::{
    verify_equivalence, verify_equivalence_configured, verify_equivalence_with_layout,
    verify_equivalence_with_layout_configured, Verdict, VerificationConfig,
};

fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    }
}

fn zero(n: usize) -> DVector<C> {
    let mut state = DVector::zeros(1 << n);
    state[0] = C::new(1.0, 0.0);
    state
}

#[test]
fn dirty_separable_ancilla_retains_leakage_penalty() {
    let logical = Circuit::new(1, 0);
    let mut routed = Circuit::new(2, 0);
    routed.add_op(gate(GateType::H, &[1], &[]));
    // I |psi> tensor |+>: overlap with I |psi> tensor |0> is 1/sqrt(2).
    let fidelity =
        equivalence_by_sampling_with_layout(&logical, &routed, &[0], &[0], 4, 123).unwrap();
    assert!((fidelity - 0.5).abs() < 1e-12, "got {fidelity}");
    assert!(
        matches!(verify_equivalence_with_layout(&logical, &routed, &[0], &[0]).unwrap(),
        Verdict::NotEquivalent { fidelity, .. } if (fidelity - 0.5).abs() < 1e-12)
    );
}

#[test]
fn dirty_ancilla_is_rejected_on_public_fifteen_qubit_sampling_path() {
    let logical = Circuit::new(1, 0);
    let mut routed = Circuit::new(15, 0);
    routed.add_op(gate(GateType::H, &[14], &[]));
    assert!(
        matches!(verify_equivalence_with_layout(&logical, &routed, &[0], &[0]).unwrap(),
        Verdict::NotEquivalent { fidelity, method: "layout-aware sampling" }
        if (fidelity - 0.5).abs() < 1e-12)
    );
}

#[test]
fn dirty_entangled_ancilla_and_complete_leakage_are_rejected() {
    let logical = Circuit::new(1, 0);
    for operation in [
        gate(GateType::CX, &[0, 1], &[]),
        gate(GateType::X, &[1], &[]),
    ] {
        let mut routed = Circuit::new(2, 0);
        routed.add_op(operation);
        assert!(
            equivalence_by_sampling_with_layout(&logical, &routed, &[0], &[0], 8, 123).unwrap()
                < 0.99
        );
        assert!(
            !verify_equivalence_with_layout(&logical, &routed, &[0], &[0])
                .unwrap()
                .is_equivalent()
        );
    }
}

#[test]
fn clean_padded_permuted_circuits_pass_exact_and_sampling() {
    let mut logical = Circuit::new(2, 0);
    logical.add_op(gate(GateType::H, &[0], &[]));
    logical.add_op(gate(GateType::CX, &[0, 1], &[]));
    for width in [3, 15] {
        let mut routed = Circuit::new(width, 0);
        routed.add_op(gate(GateType::H, &[2], &[]));
        routed.add_op(gate(GateType::CX, &[2, 0], &[]));
        // A final SWAP changes the physical output permutation.
        routed.add_op(gate(GateType::SWAP, &[0, 1], &[]));
        assert!(
            verify_equivalence_with_layout(&logical, &routed, &[2, 0], &[2, 1])
                .unwrap()
                .is_equivalent()
        );
    }
}

#[test]
fn unitary_simulators_reject_nonunitary_operations() {
    let ops = [
        Operation::Reset { qubit: 0 },
        Operation::Measure { qubit: 0, cbit: 0 },
        Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 0,
            },
            op: Box::new(gate(GateType::X, &[0], &[])),
        },
    ];
    for op in ops {
        let mut circuit = Circuit::new(1, 1);
        circuit.add_op(op);
        assert!(
            try_circuit_to_unitary(&circuit).is_err(),
            "{:?}",
            circuit.operations
        );
        assert!(evolve_state(&circuit, &zero(1)).is_err());
        assert!(equivalence_by_sampling(&circuit, &circuit, 1, 7).is_err());
        assert!(equivalence_by_sampling_with_layout(&circuit, &circuit, &[0], &[0], 1, 7).is_err());
    }
}

#[test]
fn public_verifiers_reject_reset_conditionals_and_midcircuit_measurements() {
    for width in [1, 15] {
        let identity = Circuit::new(width, 1);
        let mut invalid = vec![];
        let mut reset = identity.clone();
        reset.add_op(Operation::Reset { qubit: 0 });
        invalid.push(reset);
        let mut conditional = identity.clone();
        conditional.add_op(Operation::Conditional {
            condition: ClassicalCondition {
                creg: "c".into(),
                value: 0,
            },
            op: Box::new(gate(GateType::X, &[0], &[])),
        });
        invalid.push(conditional);
        let mut mid_measure = identity.clone();
        mid_measure.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        mid_measure.add_op(gate(GateType::X, &[0], &[]));
        invalid.push(mid_measure);
        let layout: Vec<_> = (0..width).collect();
        for circuit in invalid {
            assert!(
                verify_equivalence(&identity, &circuit).is_err(),
                "{:?}",
                circuit.operations
            );
            assert!(verify_equivalence_with_layout(&identity, &circuit, &layout, &layout).is_err());
        }
    }
}

#[test]
fn matching_terminal_readout_is_supported_and_mapping_mismatches_fail() {
    let mut logical = Circuit::new(2, 2);
    logical.add_op(gate(GateType::H, &[0], &[]));
    logical.add_op(Operation::Measure { qubit: 0, cbit: 1 });
    logical.add_op(Operation::Measure { qubit: 1, cbit: 0 });
    assert!(verify_equivalence(&logical, &logical)
        .unwrap()
        .is_equivalent());
    let mut routed = Circuit::new(3, 2);
    routed.add_op(gate(GateType::H, &[2], &[]));
    routed.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    routed.add_op(Operation::Measure { qubit: 2, cbit: 1 });
    assert!(
        verify_equivalence_with_layout(&logical, &routed, &[2, 0], &[2, 0])
            .unwrap()
            .is_equivalent()
    );
    routed.operations[1] = Operation::Measure { qubit: 0, cbit: 1 };
    routed.operations[2] = Operation::Measure { qubit: 2, cbit: 0 };
    assert!(verify_equivalence_with_layout(&logical, &routed, &[2, 0], &[2, 0]).is_err());
    let unitary = Circuit::new(2, 2);
    assert!(verify_equivalence(&unitary, &logical).is_err());
}

#[test]
fn terminal_readout_is_supported_on_public_sampling_paths() {
    let mut circuit = Circuit::new(15, 1);
    circuit.add_op(gate(GateType::H, &[14], &[]));
    circuit.add_op(Operation::Measure { qubit: 14, cbit: 0 });
    circuit.add_op(Operation::Barrier { qubits: vec![14] });
    assert!(matches!(
        verify_equivalence(&circuit, &circuit).unwrap(),
        Verdict::StatisticallyEquivalent { .. }
    ));
    let layout: Vec<_> = (0..15).collect();
    assert!(matches!(
        verify_equivalence_with_layout(&circuit, &circuit, &layout, &layout).unwrap(),
        Verdict::StatisticallyEquivalent { .. }
    ));
}

#[test]
fn invalid_readout_and_repeated_classical_writes_fail_closed() {
    for operations in [
        vec![Operation::Measure { qubit: 2, cbit: 0 }],
        vec![Operation::Measure { qubit: 0, cbit: 2 }],
        vec![
            Operation::Measure { qubit: 0, cbit: 0 },
            Operation::Measure { qubit: 1, cbit: 0 },
        ],
    ] {
        let mut circuit = Circuit::new(2, 2);
        circuit.operations = operations;
        assert!(verify_equivalence(&circuit, &circuit).is_err());
        assert!(verify_equivalence_with_layout(&circuit, &circuit, &[0, 1], &[0, 1]).is_err());
    }
}

#[test]
fn low_level_sampling_validates_layout_prefixes_and_injectivity() {
    let logical = Circuit::new(2, 0);
    let routed = Circuit::new(3, 0);
    for (initial, final_l) in [
        (vec![0, 0], vec![0, 1]),
        (vec![0, 1], vec![1, 1]),
        // The old chain(...).take(2*n) skipped final_l[1] for an extended initial layout.
        (vec![0, 1, 2], vec![0, 99]),
    ] {
        assert!(
            equivalence_by_sampling_with_layout(&logical, &routed, &initial, &final_l, 1, 7)
                .is_err()
        );
    }
    let too_narrow = Circuit::new(1, 0);
    assert!(
        equivalence_by_sampling_with_layout(&logical, &too_narrow, &[0, 0], &[0, 0], 1, 7).is_err()
    );
}

#[test]
fn extended_layout_metadata_does_not_turn_ancillas_into_logical_qubits() {
    let mut routed = Circuit::new(2, 0);
    routed.add_op(gate(GateType::H, &[1], &[]));
    let full = try_circuit_to_unitary(&routed).unwrap();
    let logical = extract_logical_unitary(&full, 1, &[0, 1], &[0, 1]);
    let expected = DMatrix::<C>::identity(2, 2) * C::new(1.0 / 2.0_f64.sqrt(), 0.0);
    assert!((logical - expected).norm() < 1e-12);
}

#[test]
fn fallible_extraction_rejects_invalid_matrix_dimensions_and_layouts() {
    for matrix in [
        DMatrix::<C>::zeros(0, 0),
        DMatrix::zeros(3, 3),
        DMatrix::zeros(2, 4),
    ] {
        assert!(try_extract_logical_unitary(&matrix, 1, &[0], &[0]).is_err());
    }
    let matrix = DMatrix::<C>::identity(4, 4);
    assert!(try_extract_logical_unitary(&matrix, 2, &[0, 0], &[0, 1]).is_err());
    assert!(try_extract_logical_unitary(&matrix, 3, &[0, 1, 2], &[0, 1, 2]).is_err());
}

#[test]
fn malformed_or_nonfinite_gates_fail_closed_in_both_simulators() {
    for op in [
        gate(GateType::RX, &[0], &[]),
        gate(GateType::U, &[0], &[0.1, 0.2]),
        gate(GateType::H, &[0, 1], &[]),
        gate(GateType::CX, &[0, 0], &[]),
        gate(GateType::RX, &[0], &[f64::NAN]),
        gate(GateType::RY, &[0], &[f64::INFINITY]),
    ] {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(op);
        assert!(
            try_circuit_to_unitary(&circuit).is_err(),
            "{:?}",
            circuit.operations
        );
        assert!(evolve_state(&circuit, &zero(2)).is_err());
        assert!(equivalence_by_sampling(&circuit, &circuit, 1, 7).is_err());
    }
}

#[test]
fn nan_sampling_never_returns_perfect_fidelity() {
    let mut circuit = Circuit::new(1, 0);
    circuit.add_op(gate(GateType::RX, &[0], &[f64::NAN]));
    assert!(equivalence_by_sampling(&circuit, &circuit, 1, 7).is_err());
    assert!(equivalence_by_sampling_with_layout(&circuit, &circuit, &[0], &[0], 1, 7).is_err());
}

#[test]
fn direct_sampling_checks_width_before_allocating_or_shifting() {
    let oversized = Circuit::new(usize::BITS as usize, 0);
    assert!(equivalence_by_sampling(&oversized, &oversized, 1, 7).is_err());
}

#[test]
fn resource_policy_uses_sampling_without_bypassing_semantic_checks() {
    let config = VerificationConfig {
        max_exact_qubits: 1,
        max_sampling_qubits: 3,
        samples: 2,
        ..VerificationConfig::default()
    };
    let mut circuit = Circuit::new(2, 1);
    circuit.add_op(gate(GateType::H, &[0], &[]));
    circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });
    assert!(matches!(
        verify_equivalence_configured(&circuit, &circuit, &config).unwrap(),
        Verdict::StatisticallyEquivalent { samples: 2, .. }
    ));
    assert!(matches!(
        verify_equivalence_with_layout_configured(&circuit, &circuit, &[0, 1], &[0, 1], &config)
            .unwrap(),
        Verdict::StatisticallyEquivalent { samples: 2, .. }
    ));

    let logical = Circuit::new(1, 0);
    let mut dirty = Circuit::new(2, 0);
    dirty.add_op(gate(GateType::H, &[1], &[]));
    assert!(
        matches!(verify_equivalence_with_layout_configured(&logical, &dirty, &[0], &[0], &config).unwrap(),
        Verdict::NotEquivalent { fidelity, .. } if (fidelity - 0.5).abs() < 1e-12)
    );
    circuit.add_op(Operation::Reset { qubit: 0 });
    assert!(verify_equivalence_configured(&circuit, &circuit, &config).is_err());
    assert!(verify_equivalence_with_layout_configured(
        &circuit,
        &circuit,
        &[0, 1],
        &[0, 1],
        &config
    )
    .is_err());
}

#[test]
fn resource_policy_rejects_invalid_limits_and_reports_inconclusive_above_budget() {
    let circuit = Circuit::new(2, 0);
    for config in [
        VerificationConfig {
            max_exact_qubits: 15,
            ..VerificationConfig::default()
        },
        VerificationConfig {
            max_sampling_qubits: 23,
            ..VerificationConfig::default()
        },
        VerificationConfig {
            samples: 0,
            ..VerificationConfig::default()
        },
    ] {
        assert!(verify_equivalence_configured(&circuit, &circuit, &config).is_err());
        assert!(verify_equivalence_with_layout_configured(
            &circuit,
            &circuit,
            &[0, 1],
            &[0, 1],
            &config
        )
        .is_err());
    }
    let config = VerificationConfig {
        max_exact_qubits: 0,
        max_sampling_qubits: 1,
        ..VerificationConfig::default()
    };
    assert!(matches!(
        verify_equivalence_configured(&circuit, &circuit, &config).unwrap(),
        Verdict::Unverifiable { .. }
    ));
    assert!(matches!(
        verify_equivalence_with_layout_configured(&circuit, &circuit, &[0, 1], &[0, 1], &config)
            .unwrap(),
        Verdict::Unverifiable { .. }
    ));
    let mut unsupported = circuit.clone();
    unsupported.add_op(Operation::Reset { qubit: 0 });
    // Semantic validation still runs when simulation would exceed the budget.
    assert!(verify_equivalence_configured(&unsupported, &unsupported, &config).is_err());
}

#[test]
fn rotations_and_controlled_rotations_match_independent_analytic_matrices() {
    let theta: f64 = 0.713;
    let c = C::new((theta / 2.0).cos(), 0.0);
    let s = (theta / 2.0).sin();
    let rotation_x = DMatrix::from_row_slice(2, 2, &[c, C::new(0.0, -s), C::new(0.0, -s), c]);
    let rotation_y = DMatrix::from_row_slice(2, 2, &[c, C::new(-s, 0.0), C::new(s, 0.0), c]);
    let rotation_z = DMatrix::from_diagonal(&DVector::from_vec(vec![
        C::from_polar(1.0, -theta / 2.0),
        C::from_polar(1.0, theta / 2.0),
    ]));
    for (single, controlled, analytic) in [
        (GateType::RX, GateType::CRX, rotation_x),
        (GateType::RY, GateType::CRY, rotation_y),
        (GateType::RZ, GateType::CRZ, rotation_z),
    ] {
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(gate(single, &[0], &[theta]));
        assert!(
            unitary_fidelity(&try_circuit_to_unitary(&circuit).unwrap(), &analytic) > 1.0 - 1e-12
        );
        let mut expected = DMatrix::<C>::identity(4, 4);
        for r in 0..2 {
            for col in 0..2 {
                expected[(1 + 2 * r, 1 + 2 * col)] = analytic[(r, col)];
            }
        }
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(gate(controlled, &[0, 1], &[theta]));
        assert!(
            unitary_fidelity(&try_circuit_to_unitary(&circuit).unwrap(), &expected) > 1.0 - 1e-12
        );
    }
}
