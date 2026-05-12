//! End-to-end verification suite for well-known quantum algorithms.
//!
//! Each test builds a canonical circuit for a textbook algorithm, then runs
//! the complete Q-Rust pipeline (optimize → layout → route → basis translate)
//! against a specific hardware topology and target gate set.  Mathematical
//! correctness is verified by comparing the unitary of the original circuit
//! against the logical unitary extracted from the compiled output.
//!
//! ## Verification strategy
//!
//! Two complementary approaches are used:
//!
//! 1. **Basis-translation fidelity** (`verify_basis_only`): Uses an all-to-all
//!    topology with no routing so the qubit mapping is always the identity.
//!    Directly compares `circuit_to_unitary(original)` vs
//!    `circuit_to_unitary(compiled)`.  Validates that optimization passes and
//!    basis translation preserve the unitary exactly.
//!
//! 2. **Routed fidelity** (`verify_routed`): Runs the router through
//!    `BeamSabrePass` directly to obtain the `initial_layout` / `final_layout`
//!    from `PropertySet`, then uses `extract_logical_unitary` to account for
//!    qubit permutations introduced by routing.

use q_rust::backend::{Backend, BackendConfig};
use q_rust::error::QRustError;
use q_rust::ir::Circuit;
use q_rust::parser::parse_qasm;
use q_rust::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::BeamSabrePass;
use q_rust::transpiler::{transpile, TranspilerConfig};
use std::collections::HashSet;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn load_backend(name: &str) -> Backend {
    let path = format!("tests/fixtures/{name}");
    let json = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing fixture: {path}"));
    let cfg: BackendConfig = serde_json::from_str(&json).unwrap();
    Backend::from_config(cfg)
}

fn basis(gates: &[&str]) -> HashSet<String> {
    gates.iter().map(|s| s.to_string()).collect()
}

fn load_circuit(fixture: &str) -> Circuit {
    let path = format!("tests/fixtures/{fixture}");
    let qasm = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing fixture: {path}"));
    parse_qasm(&qasm).unwrap_or_else(|e| panic!("parse failed for {fixture}: {e:?}"))
}

/// Verify that optimize + basis-translation preserves the unitary exactly.
/// Uses an all-to-all backend so routing never permutes qubits.
fn verify_basis_only(label: &str, circuit: &Circuit, target_basis: HashSet<String>) {
    let u_orig = circuit_to_unitary(circuit);
    let b = Backend::all_to_all(circuit.num_qubits);

    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .decompose_basis(true)
        .backend(b)
        .target_basis(target_basis)
        .build();

    let compiled =
        transpile(circuit, Some(cfg)).unwrap_or_else(|e| panic!("[{label}] transpile failed: {e}"));

    let u_out = circuit_to_unitary(&compiled);
    let fid = unitary_fidelity(&u_orig, &u_out);
    assert!(
        (fid - 1.0).abs() < 1e-6,
        "[{label}] basis-translation fidelity = {fid:.10} (expected ≈ 1.0)"
    );
}

/// Verify fidelity after routing to a real hardware topology.
/// Uses the layout maps from PropertySet to extract the logical unitary.
fn verify_routed(label: &str, circuit: &Circuit, backend: &Backend, target_basis: HashSet<String>) {
    let u_orig = circuit_to_unitary(circuit);

    // 1. Run routing pass directly to capture layout metadata.
    let router = BeamSabrePass {
        backend: backend.clone(),
        beam_width: 4,
        branch_factor: 3,
        bidir_iterations: 2,
    };
    let mut props = PropertySet::new();
    let routed = router.run(circuit, &mut props);

    let initial: Vec<usize> = props
        .get::<Vec<usize>>("initial_layout")
        .cloned()
        .unwrap_or_else(|| (0..circuit.num_qubits).collect());
    let final_l: Vec<usize> = props
        .get::<Vec<usize>>("final_layout")
        .cloned()
        .unwrap_or_else(|| (0..circuit.num_qubits).collect());

    // 2. Apply basis translation on top of the routed circuit.
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .target_basis(target_basis)
        .build();
    let translated = transpile(&routed, Some(cfg))
        .unwrap_or_else(|e| panic!("[{label}] basis translation failed: {e}"));

    // 3. Pad to backend size and extract logical unitary.
    let mut padded = translated.clone();
    padded.num_qubits = backend.num_qubits;
    let u_phys = circuit_to_unitary(&padded);
    let u_log = extract_logical_unitary(&u_phys, circuit.num_qubits, &initial, &final_l);

    let fid = unitary_fidelity(&u_orig, &u_log);
    assert!(
        (fid - 1.0).abs() < 1e-6,
        "[{label}] routed fidelity = {fid:.10} (expected ≈ 1.0)"
    );
}

// ===========================================================================
// Bell State (2 qubits — H + CX, creates maximally entangled pair)
// ===========================================================================

#[test]
fn test_e2e_bell_basis_rz_rx_cx() {
    verify_basis_only(
        "bell/{rz,rx,cx}",
        &load_circuit("bell_state.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_bell_basis_u_cx() {
    verify_basis_only(
        "bell/{u,cx}",
        &load_circuit("bell_state.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_bell_routed_quito_rz_rx_cx() {
    verify_routed(
        "bell/ibm_quito/{rz,rx,cx}",
        &load_circuit("bell_state.qasm"),
        &load_backend("ibm_quito_5q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_bell_routed_linear3_u_cx() {
    verify_routed(
        "bell/linear-3/{u,cx}",
        &load_circuit("bell_state.qasm"),
        &Backend::linear(3),
        basis(&["u", "cx"]),
    );
}

// ===========================================================================
// GHZ-3 (3 qubits — H + 2×CX, creates 3-qubit entanglement)
// ===========================================================================

#[test]
fn test_e2e_ghz3_basis_rz_rx_cx() {
    verify_basis_only(
        "ghz3/{rz,rx,cx}",
        &load_circuit("ghz_3.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_ghz3_basis_u_cx() {
    verify_basis_only(
        "ghz3/{u,cx}",
        &load_circuit("ghz_3.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_ghz3_routed_nairobi_rz_rx_cx() {
    verify_routed(
        "ghz3/ibm_nairobi/{rz,rx,cx}",
        &load_circuit("ghz_3.qasm"),
        &load_backend("ibm_nairobi_7q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_ghz3_routed_linear5_u_cx() {
    verify_routed(
        "ghz3/linear-5/{u,cx}",
        &load_circuit("ghz_3.qasm"),
        &Backend::linear(5),
        basis(&["u", "cx"]),
    );
}

// ===========================================================================
// Grover 2-qubit (oracle + diffusion, marks state |11⟩)
// ===========================================================================

#[test]
fn test_e2e_grover2_basis_rz_rx_cx() {
    verify_basis_only(
        "grover2/{rz,rx,cx}",
        &load_circuit("grover_2.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_grover2_basis_u_cx() {
    verify_basis_only(
        "grover2/{u,cx}",
        &load_circuit("grover_2.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_grover2_routed_quito_rz_rx_cx() {
    verify_routed(
        "grover2/ibm_quito/{rz,rx,cx}",
        &load_circuit("grover_2.qasm"),
        &load_backend("ibm_quito_5q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

// ===========================================================================
// QFT-3 (Quantum Fourier Transform on 3 qubits)
// ===========================================================================

#[test]
fn test_e2e_qft3_basis_rz_rx_cx() {
    verify_basis_only(
        "qft3/{rz,rx,cx}",
        &load_circuit("qft_3.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_qft3_basis_u_cx() {
    verify_basis_only(
        "qft3/{u,cx}",
        &load_circuit("qft_3.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_qft3_routed_nairobi_rz_rx_cx() {
    verify_routed(
        "qft3/ibm_nairobi/{rz,rx,cx}",
        &load_circuit("qft_3.qasm"),
        &load_backend("ibm_nairobi_7q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

// ===========================================================================
// Deutsch-Jozsa (2 query + 1 ancilla, balanced oracle f(x)=x₀⊕x₁)
// ===========================================================================

#[test]
fn test_e2e_dj_basis_rz_rx_cx() {
    verify_basis_only(
        "dj/{rz,rx,cx}",
        &load_circuit("deutsch_jozsa.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_dj_basis_u_cx() {
    verify_basis_only(
        "dj/{u,cx}",
        &load_circuit("deutsch_jozsa.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_dj_routed_linear5_rz_rx_cx() {
    verify_routed(
        "dj/linear-5/{rz,rx,cx}",
        &load_circuit("deutsch_jozsa.qasm"),
        &Backend::linear(5),
        basis(&["rz", "rx", "cx"]),
    );
}

// ===========================================================================
// Bernstein-Vazirani (3 query + 1 ancilla, hidden string s=101)
// ===========================================================================

#[test]
fn test_e2e_bv_basis_rz_rx_cx() {
    verify_basis_only(
        "bv/{rz,rx,cx}",
        &load_circuit("bernstein_vazirani.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_bv_basis_u_cx() {
    verify_basis_only(
        "bv/{u,cx}",
        &load_circuit("bernstein_vazirani.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_bv_routed_nairobi_rz_rx_cx() {
    verify_routed(
        "bv/ibm_nairobi/{rz,rx,cx}",
        &load_circuit("bernstein_vazirani.qasm"),
        &load_backend("ibm_nairobi_7q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

// ===========================================================================
// VQE ansatz (4q variational circuit with parametric rotations)
// ===========================================================================

#[test]
fn test_e2e_vqe_basis_rz_rx_cx() {
    verify_basis_only(
        "vqe/{rz,rx,cx}",
        &load_circuit("vqe_ansatz.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_vqe_basis_u_cx() {
    verify_basis_only(
        "vqe/{u,cx}",
        &load_circuit("vqe_ansatz.qasm"),
        basis(&["u", "cx"]),
    );
}

// ===========================================================================
// Quantum Phase Estimation (toy 2q: estimates eigenphase of T gate)
// ===========================================================================

#[test]
fn test_e2e_qpe_basis_rz_rx_cx() {
    verify_basis_only(
        "qpe/{rz,rx,cx}",
        &load_circuit("qpe_2.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_qpe_routed_quito_rz_rx_cx() {
    verify_routed(
        "qpe/ibm_quito/{rz,rx,cx}",
        &load_circuit("qpe_2.qasm"),
        &load_backend("ibm_quito_5q.json"),
        basis(&["rz", "rx", "cx"]),
    );
}

// ===========================================================================
// Deep Clifford (H, S, CX, Sdg cascade — 3 qubits)
// ===========================================================================

#[test]
fn test_e2e_clifford_basis_rz_rx_cx() {
    verify_basis_only(
        "clifford/{rz,rx,cx}",
        &load_circuit("deep_clifford.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_clifford_basis_u_cx() {
    verify_basis_only(
        "clifford/{u,cx}",
        &load_circuit("deep_clifford.qasm"),
        basis(&["u", "cx"]),
    );
}

#[test]
fn test_e2e_clifford_routed_nairobi_u_cx() {
    verify_routed(
        "clifford/ibm_nairobi/{u,cx}",
        &load_circuit("deep_clifford.qasm"),
        &load_backend("ibm_nairobi_7q.json"),
        basis(&["u", "cx"]),
    );
}

// ===========================================================================
// QFT-5 (5-qubit Quantum Fourier Transform — deep parametric circuit)
// ===========================================================================

#[test]
fn test_e2e_qft5_basis_rz_rx_cx() {
    verify_basis_only(
        "qft5/{rz,rx,cx}",
        &load_circuit("qft_5.qasm"),
        basis(&["rz", "rx", "cx"]),
    );
}

#[test]
fn test_e2e_qft5_basis_u_cx() {
    verify_basis_only(
        "qft5/{u,cx}",
        &load_circuit("qft_5.qasm"),
        basis(&["u", "cx"]),
    );
}

// ===========================================================================
// Gate-set validation: universality checks
// ===========================================================================

#[test]
fn test_gateset_rejects_clifford_only() {
    let c = load_circuit("bell_state.qasm");
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .target_basis(["h", "s", "cx"].iter().map(|s| s.to_string()))
        .build();
    let result = transpile(&c, Some(cfg));
    assert!(
        matches!(result, Err(QRustError::NonUniversalBasisGateSet { .. })),
        "expected NonUniversalBasisGateSet, got: {result:?}"
    );
}

#[test]
fn test_gateset_rejects_no_entangler() {
    let c = load_circuit("bell_state.qasm");
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .target_basis(["rz", "rx", "ry"].iter().map(|s| s.to_string()))
        .build();
    let result = transpile(&c, Some(cfg));
    assert!(
        matches!(result, Err(QRustError::NonUniversalBasisGateSet { .. })),
        "expected NonUniversalBasisGateSet, got: {result:?}"
    );
}

#[test]
fn test_gateset_accepts_solovay_kitaev_basis() {
    // {H, T, CX}: T is non-Clifford → densely covers SU(2) → universal.
    let c = load_circuit("bell_state.qasm");
    let b = Backend::all_to_all(2);
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .backend(b)
        .target_basis(["h", "t", "cx"].iter().map(|s| s.to_string()))
        .build();
    assert!(transpile(&c, Some(cfg)).is_ok());
}

#[test]
fn test_gateset_accepts_rz_cx() {
    let c = load_circuit("bell_state.qasm");
    let b = Backend::all_to_all(2);
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(true)
        .backend(b)
        .target_basis(["rz", "cx"].iter().map(|s| s.to_string()))
        .build();
    assert!(transpile(&c, Some(cfg)).is_ok());
}

// ===========================================================================
// Permissive mode: no basis specified → circuit passes through unchanged
// ===========================================================================

#[test]
fn test_permissive_mode_no_basis() {
    let c = load_circuit("bell_state.qasm");
    let cfg = TranspilerConfig::builder()
        .optimization_level(1)
        .decompose_basis(false)
        .build();
    let result = transpile(&c, Some(cfg));
    assert!(
        result.is_ok(),
        "permissive mode should not error: {result:?}"
    );
    assert!(!result.unwrap().operations.is_empty());
}
