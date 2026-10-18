use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{
    numerical_policy::CIRCUIT_ERROR_ALLOWANCE,
    optimization::{
        CommutationCancellationPass, GateCrystallizationPass, InverseCancellationPass,
        ParameterSimplificationPass, RotationMergePass,
    },
    pass::Pass,
    property_set::PropertySet,
    transpile, transpile_with_report, TranspilerConfig,
};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
    Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    }
}

#[test]
fn default_simplification_preserves_nonzero_rotations_and_tiny_u_phases() {
    let mut source = Circuit::new(1, 0);
    for angle in [9e-10, -9e-10, 1e-16, -1e-16] {
        for name in [GateType::RX, GateType::RY, GateType::RZ] {
            source.add_op(gate(name, &[0], &[angle]));
        }
        source.add_op(gate(GateType::U, &[0], &[0.0, 0.0, angle]));
    }
    let output = ParameterSimplificationPass::default().run(&source, &mut PropertySet::new());
    assert_eq!(output.operations, source.operations);
}

#[test]
fn an_explicit_candidate_epsilon_cannot_bypass_residual_acceptance() {
    let mut source = Circuit::new(1, 0);
    for angle in [9e-10, -9e-10] {
        source.add_op(gate(GateType::RX, &[0], &[angle]));
        source.add_op(gate(GateType::RY, &[0], &[PI + angle]));
    }
    assert_eq!(
        ParameterSimplificationPass { epsilon: 0.1 }
            .run(&source, &mut PropertySet::new())
            .operations,
        source.operations
    );
    assert_eq!(
        GateCrystallizationPass { epsilon: 0.1 }
            .run(&source, &mut PropertySet::new())
            .operations,
        source.operations
    );
}

#[test]
fn near_special_angles_are_not_named_exact_gates() {
    let mut source = Circuit::new(1, 0);
    for angle in [PI, -PI, FRAC_PI_2, -FRAC_PI_2, FRAC_PI_4, -FRAC_PI_4] {
        for offset in [-9e-10, 9e-10] {
            source.add_op(gate(GateType::RZ, &[0], &[angle + offset]));
        }
    }
    source.add_op(gate(GateType::U, &[0], &[FRAC_PI_2, 9e-10, PI]));
    let output = GateCrystallizationPass::default().run(&source, &mut PropertySet::new());
    assert_eq!(output.operations, source.operations);
}

#[test]
fn almost_commuting_u_gates_do_not_authorize_cx_cancellation() {
    for (wire, params) in [
        (0, [9e-10, 0.0, 0.0]),
        (1, [0.7, -FRAC_PI_2 + 9e-10, FRAC_PI_2]),
    ] {
        let mut source = Circuit::new(2, 0);
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        source.add_op(gate(GateType::U, &[wire], &params));
        source.add_op(gate(GateType::CX, &[0, 1], &[]));
        let output = CommutationCancellationPass.run(&source, &mut PropertySet::new());
        assert_eq!(output.operations, source.operations);
    }
}

#[test]
fn almost_inverse_rotations_survive_for_both_period_conventions() {
    for (name, wires, period) in [
        (GateType::RX, vec![0], 2.0 * PI),
        (GateType::CRX, vec![0, 1], 4.0 * PI),
    ] {
        for sum in [9e-10, -9e-10, period + 9e-10, period - 9e-10] {
            let mut source = Circuit::new(2, 0);
            source.add_op(gate(name.clone(), &wires, &[0.3]));
            source.add_op(gate(name.clone(), &wires, &[sum - 0.3]));
            let output = InverseCancellationPass.run(&source, &mut PropertySet::new());
            assert_eq!(output.operations, source.operations);
        }
    }
}

#[test]
fn large_argument_reduction_and_overflowing_merges_fail_closed() {
    let mut source = Circuit::new(1, 0);
    source.add_op(gate(GateType::RX, &[0], &[1e16]));
    let output = ParameterSimplificationPass::default().run(&source, &mut PropertySet::new());
    // f64 remainder by a rounded 2pi has a substantial phase error at 1e16.
    assert_eq!(output.operations, source.operations);
    source.operations.clear();
    source.add_op(gate(GateType::RX, &[0], &[f64::MAX]));
    source.add_op(gate(GateType::RX, &[0], &[f64::MAX]));
    assert_eq!(
        RotationMergePass
            .run(&source, &mut PropertySet::new())
            .operations,
        source.operations
    );
}

#[test]
fn qasm_serialization_roundtrips_f64_parameters_including_subnormals() {
    let mut source = Circuit::new(1, 0);
    for value in [
        f64::from_bits(1),
        -f64::from_bits(1),
        4.9e-11,
        -4.9e-11,
        0.12345678901234567,
        PI + 9e-10,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        source.add_op(gate(GateType::RX, &[0], &[value]));
    }
    let output = parse_qasm(&source.to_qasm(None)).unwrap();
    for (before, after) in source.operations.iter().zip(&output.operations) {
        let (Operation::Gate { params: a, .. }, Operation::Gate { params: b, .. }) =
            (before, after)
        else {
            panic!("gate changed")
        };
        assert_eq!(a[0].to_bits(), b[0].to_bits());
    }
}

#[test]
fn shared_budget_survives_pipeline_split_and_plain_report_paths_agree() {
    let mut source = Circuit::new(1, 0);
    // The first segment creates numerical work before the stage split;
    // the second creates U fusion after it, with a barrier separating them.
    for _ in 0..1000 {
        source.add_op(gate(GateType::RX, &[0], &[9e-10]));
    }
    source.add_op(Operation::Barrier { qubits: vec![0] });
    for _ in 0..10 {
        source.add_op(gate(GateType::U, &[0], &[0.1, 0.2, 0.3]));
    }
    let config = TranspilerConfig::builder().optimization_level(3).build();
    let plain = transpile(&source, Some(config.clone())).unwrap();
    let (output, report) = transpile_with_report(&source, Some(config)).unwrap();
    assert_eq!(plain.operations, output.operations);
    assert!(report.numerical_budget.accepted >= 1008);
    assert!(report.numerical_budget.charged_error <= CIRCUIT_ERROR_ALLOWANCE);
    assert!(report.numerical_budget.charged_error >= report.numerical_budget.residual_sum);
}
