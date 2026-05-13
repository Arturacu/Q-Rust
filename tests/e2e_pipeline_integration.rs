//! [E2E-NEW-FEATURE] End-to-end integration test for thesis Section 4.3.
//!
//! Mirrors the thesis' GHZ walkthrough: read QASM file → parse →
//! transpile → emit QASM → re-parse.
//!
//! Note on verification: post-routing, qubit indices in the transpiled
//! circuit are *physical*, not logical, so a direct `verify_equivalence`
//! between the source and the transpiled circuit only matches when SABRE
//! picks the identity layout. We therefore verify equivalence on the
//! pre-routing pipeline (no backend) and treat the routed pipeline as a
//! structural-correctness check (parses, round-trips, matches op count).

use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, transpile_with_report, TranspilerConfig};
use q_rust::verify::{verify_equivalence, Verdict};

const GHZ_QASM: &str = r#"
OPENQASM 2.0;
qreg q[3];
creg c[3];
h q[0];
cx q[0], q[1];
cx q[1], q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
"#;

#[test]
fn test_e2e_full_pipeline_ghz_no_backend_verifiable() {
    // Stage 1: Parsing.
    let circuit = parse_qasm(GHZ_QASM).expect("parse");
    assert_eq!(circuit.num_qubits, 3);
    assert!(circuit.operations.len() >= 6);

    // Stages 2–4: Transpilation WITHOUT a backend (no routing → indices
    // remain logical → safe to verify equivalence directly).
    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .build();
    let transpiled = transpile(&circuit, Some(cfg)).expect("transpile");

    // Strip measurements (non-unitary) before verifying.
    let mut c1 = circuit.clone();
    c1.operations
        .retain(|op| !matches!(op, q_rust::ir::Operation::Measure { .. }));
    let mut c2 = transpiled.clone();
    c2.operations
        .retain(|op| !matches!(op, q_rust::ir::Operation::Measure { .. }));

    let verdict = verify_equivalence(&c1, &c2).expect("verify");
    assert!(
        verdict.is_equivalent(),
        "GHZ pipeline broke equivalence: {}",
        verdict.describe()
    );
    assert!(matches!(verdict, Verdict::ExactlyEquivalent { .. }));

    // Output: emit QASM and re-parse, confirming round-trip.
    let emitted = transpiled.to_qasm(None);
    assert!(emitted.starts_with("OPENQASM 2.0;"));
    assert!(emitted.contains("qreg q[3];"));
    let reparsed = parse_qasm(&emitted).expect("re-parse emitted QASM");
    assert_eq!(reparsed.num_qubits, transpiled.num_qubits);
    assert_eq!(reparsed.operations.len(), transpiled.operations.len());
}

#[test]
fn test_e2e_full_pipeline_ghz_with_backend_structural() {
    // Same as above but with a backend — verifies the routed pipeline
    // produces a parseable, well-formed circuit. We do *not* verify
    // unitary equivalence here because qubit indices are physical.
    let circuit = parse_qasm(GHZ_QASM).expect("parse");
    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .backend(Backend::linear(3))
        .build();
    let transpiled = transpile(&circuit, Some(cfg)).expect("transpile");

    assert!(!transpiled.operations.is_empty());
    let emitted = transpiled.to_qasm(None);
    assert!(emitted.contains("OPENQASM 2.0;"));
    let reparsed = parse_qasm(&emitted).expect("re-parse");
    assert_eq!(reparsed.operations.len(), transpiled.operations.len());
}

#[test]
fn test_e2e_with_report_thesis_fig_4_1_annotations() {
    let circuit = parse_qasm(GHZ_QASM).expect("parse");
    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .decompose_basis(true)
        .backend(Backend::linear(3))
        .build();

    let (out, report) = transpile_with_report(&circuit, Some(cfg)).expect("transpile");
    assert_eq!(report.stages.len(), 3);
    assert_eq!(report.stages[0].stage, "1. parsed");
    assert_eq!(report.stages[1].stage, "2. optimized");
    assert!(report.stages[2].stage.starts_with("3."));

    let lines = report.format_lines();
    assert_eq!(lines.len(), 4);
    assert!(!out.operations.is_empty());
}

#[test]
fn test_e2e_qasm_to_qasm_via_builtin_backend() {
    let circuit = parse_qasm(GHZ_QASM).expect("parse");
    let backend = Backend::ibm_quito();
    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .decompose_basis(true)
        .backend(backend)
        .build();
    let transpiled = transpile(&circuit, Some(cfg)).expect("transpile");
    let emitted = transpiled.to_qasm(None);
    assert!(emitted.contains("OPENQASM 2.0;"));
    let _ = parse_qasm(&emitted).expect("re-parse");
}