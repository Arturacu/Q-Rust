//! Fixture-driven transpiler correctness suite.

use q_rust::ir::{GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

fn config_default() -> TranspilerConfig {
    TranspilerConfig::default()
}

fn config_no_decompose() -> TranspilerConfig {
    TranspilerConfig::builder()
        .decompose_basis(false)
        .optimization_level(0)
        .build()
}

fn config_decompose_only() -> TranspilerConfig {
    TranspilerConfig::builder()
        .decompose_basis(true)
        .optimization_level(0)
        .build()
}

const FIDELITY_THRESHOLD: f64 = 1.0 - 1e-8;

fn run_and_check(qasm: &str, config: TranspilerConfig, fixture_name: &str) {
    let circuit = parse_qasm(qasm).unwrap_or_else(|e| panic!("[{fixture_name}] parse: {e}"));
    let original_qubits = circuit.num_qubits;
    let original_cbits = circuit.num_cbits;
    let decompose = config.decompose_basis;
    let u_orig = circuit_to_unitary(&circuit);

    let result = transpile(&circuit, Some(config)).expect("transpile failed");

    assert_eq!(
        result.num_qubits, original_qubits,
        "[{fixture_name}] num_qubits"
    );
    assert_eq!(
        result.num_cbits, original_cbits,
        "[{fixture_name}] num_cbits"
    );

    if decompose {
        for (i, op) in result.operations.iter().enumerate() {
            if let Operation::Gate { name, .. } = op {
                assert!(
                    matches!(name, GateType::U | GateType::CX),
                    "[{fixture_name}] non-basis gate at #{i}: {name:?}",
                );
            }
        }
    }

    let om = circuit
        .operations
        .iter()
        .filter(|op| matches!(op, Operation::Measure { .. }))
        .count();
    let rm = result
        .operations
        .iter()
        .filter(|op| matches!(op, Operation::Measure { .. }))
        .count();
    assert_eq!(rm, om, "[{fixture_name}] measurement count");

    let u_res = circuit_to_unitary(&result);
    let fid = unitary_fidelity(&u_orig, &u_res);
    assert!(
        fid > FIDELITY_THRESHOLD,
        "[{fixture_name}] fidelity {fid} below threshold"
    );
}

macro_rules! fixture_test {
    ($name:ident, $path:expr, $config_fn:ident) => {
        #[test]
        fn $name() {
            let qasm = include_str!($path);
            run_and_check(qasm, $config_fn(), stringify!($name));
        }
    };
}

fixture_test!(
    test_bell_state_default,
    "fixtures/bell_state.qasm",
    config_default
);
fixture_test!(
    test_bell_state_no_decompose,
    "fixtures/bell_state.qasm",
    config_no_decompose
);
fixture_test!(
    test_bell_state_decompose_only,
    "fixtures/bell_state.qasm",
    config_decompose_only
);

fixture_test!(test_ghz_3_default, "fixtures/ghz_3.qasm", config_default);
fixture_test!(
    test_ghz_3_no_decompose,
    "fixtures/ghz_3.qasm",
    config_no_decompose
);
fixture_test!(
    test_ghz_3_decompose_only,
    "fixtures/ghz_3.qasm",
    config_decompose_only
);

fixture_test!(
    test_rotations_default,
    "fixtures/rotations.qasm",
    config_default
);
fixture_test!(
    test_rotations_no_decompose,
    "fixtures/rotations.qasm",
    config_no_decompose
);
fixture_test!(
    test_rotations_decompose_only,
    "fixtures/rotations.qasm",
    config_decompose_only
);

fixture_test!(
    test_multi_gate_default,
    "fixtures/multi_gate.qasm",
    config_default
);
fixture_test!(
    test_multi_gate_no_decompose,
    "fixtures/multi_gate.qasm",
    config_no_decompose
);
fixture_test!(
    test_multi_gate_decompose_only,
    "fixtures/multi_gate.qasm",
    config_decompose_only
);

fixture_test!(
    test_toffoli_default,
    "fixtures/toffoli.qasm",
    config_default
);
fixture_test!(
    test_toffoli_no_decompose,
    "fixtures/toffoli.qasm",
    config_no_decompose
);
fixture_test!(
    test_toffoli_decompose_only,
    "fixtures/toffoli.qasm",
    config_decompose_only
);

fixture_test!(
    test_custom_gate_default,
    "fixtures/custom_gate.qasm",
    config_default
);
fixture_test!(
    test_custom_gate_no_decompose,
    "fixtures/custom_gate.qasm",
    config_no_decompose
);
fixture_test!(
    test_custom_gate_decompose_only,
    "fixtures/custom_gate.qasm",
    config_decompose_only
);

fixture_test!(
    test_identity_only_default,
    "fixtures/identity_only.qasm",
    config_default
);
fixture_test!(
    test_identity_only_no_decompose,
    "fixtures/identity_only.qasm",
    config_no_decompose
);
fixture_test!(
    test_identity_only_decompose_only,
    "fixtures/identity_only.qasm",
    config_decompose_only
);

fixture_test!(
    test_controlled_gates_default,
    "fixtures/controlled_gates.qasm",
    config_default
);
fixture_test!(
    test_controlled_gates_no_decompose,
    "fixtures/controlled_gates.qasm",
    config_no_decompose
);
fixture_test!(
    test_controlled_gates_decompose_only,
    "fixtures/controlled_gates.qasm",
    config_decompose_only
);

fixture_test!(
    test_csx_gate_default,
    "fixtures/csx_gate.qasm",
    config_default
);
fixture_test!(
    test_csx_gate_no_decompose,
    "fixtures/csx_gate.qasm",
    config_no_decompose
);
fixture_test!(
    test_csx_gate_decompose_only,
    "fixtures/csx_gate.qasm",
    config_decompose_only
);