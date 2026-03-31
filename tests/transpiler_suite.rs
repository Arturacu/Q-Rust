//! Fixture-driven transpiler test suite.
//!
//! This test file discovers QASM circuits from `tests/fixtures/` and runs each
//! through the transpiler with multiple configurations, asserting universal
//! invariants. To add a new test circuit, drop a `.qasm` file in `tests/fixtures/`.
//! To add a new transpiler mode, add a config variant and invoke the macro.

use q_rust::ir::{GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
use q_rust::transpiler::{transpile, TranspilerConfig};

// ---------------------------------------------------------------------------
// Transpiler config presets
// ---------------------------------------------------------------------------

/// Default config: decompose on, opt level 1.
fn config_default() -> TranspilerConfig {
    TranspilerConfig::default()
}

/// No decomposition, no optimization.
fn config_no_decompose() -> TranspilerConfig {
    TranspilerConfig {
        decompose_basis: false,
        optimization_level: 0,
    }
}

/// Decompose to basis, but no optimization passes.
fn config_decompose_only() -> TranspilerConfig {
    TranspilerConfig {
        decompose_basis: true,
        optimization_level: 0,
    }
}

// ---------------------------------------------------------------------------
// Universal invariant checks
// ---------------------------------------------------------------------------

/// Minimum fidelity threshold. A perfect transpilation yields 1.0.
const FIDELITY_THRESHOLD: f64 = 1.0 - 1e-8;

/// Run the transpiler on a QASM string with a given config and check invariants.
fn run_and_check(qasm: &str, config: TranspilerConfig, fixture_name: &str) {
    // 1. Parse must succeed
    let circuit = parse_qasm(qasm).unwrap_or_else(|e| {
        panic!("[{}] Failed to parse fixture: {}", fixture_name, e);
    });

    let original_qubits = circuit.num_qubits;
    let original_cbits = circuit.num_cbits;
    let decompose = config.decompose_basis;

    // Compute original unitary (before transpilation)
    let original_unitary = circuit_to_unitary(&circuit);

    // 2. Transpile must not panic
    let result = transpile(&circuit, Some(config));

    // 3. Qubit/cbit counts must be preserved
    assert_eq!(
        result.num_qubits, original_qubits,
        "[{}] num_qubits changed: {} -> {}",
        fixture_name, original_qubits, result.num_qubits
    );
    assert_eq!(
        result.num_cbits, original_cbits,
        "[{}] num_cbits changed: {} -> {}",
        fixture_name, original_cbits, result.num_cbits
    );

    // 4. When decompose is on, all gates must be U or CX
    if decompose {
        for (i, op) in result.operations.iter().enumerate() {
            if let Operation::Gate { name, .. } = op {
                match name {
                    GateType::U | GateType::CX => {} // basis gates — OK
                    other => panic!(
                        "[{}] Non-basis gate at position {}: {:?}",
                        fixture_name, i, other
                    ),
                }
            }
        }
    }

    // 5. Measurement count must be preserved
    let original_measures = circuit
        .operations
        .iter()
        .filter(|op| matches!(op, Operation::Measure { .. }))
        .count();
    let result_measures = result
        .operations
        .iter()
        .filter(|op| matches!(op, Operation::Measure { .. }))
        .count();
    assert_eq!(
        result_measures, original_measures,
        "[{}] Measurement count changed: {} -> {}",
        fixture_name, original_measures, result_measures
    );

    // 6. Unitary equivalence: the transpiled circuit must produce the same
    //    unitary as the original (up to global phase).
    let result_unitary = circuit_to_unitary(&result);
    let fidelity = unitary_fidelity(&original_unitary, &result_unitary);
    assert!(
        fidelity > FIDELITY_THRESHOLD,
        "[{}] Unitary fidelity too low: {:.10} (threshold: {:.10})",
        fixture_name,
        fidelity,
        FIDELITY_THRESHOLD
    );
}

// ---------------------------------------------------------------------------
// Test generation macro
// ---------------------------------------------------------------------------

/// Generates a `#[test]` for every (fixture, config) combination.
///
/// Usage:
///   fixture_test!(test_name, "fixtures/file.qasm", config_fn);
macro_rules! fixture_test {
    ($name:ident, $path:expr, $config_fn:ident) => {
        #[test]
        fn $name() {
            let qasm = include_str!($path);
            let config = $config_fn();
            run_and_check(qasm, config, stringify!($name));
        }
    };
}

// ---------------------------------------------------------------------------
// Bell state
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// GHZ 3-qubit
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Rotations
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Multi gate
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Toffoli
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Custom gate
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Identity only
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Controlled gates (CZ, CY, CH, CRX, CRY, CRZ, RXX, RYY, RZZ)
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// CSX gate
// ---------------------------------------------------------------------------
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
