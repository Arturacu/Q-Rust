//! Target-basis validation and translation.
//!
//! This module provides:
//! - [`validate_universality`], a conservative recognizer for common universal
//!   gate-set families;
//! - [`TargetBasisPass`], an early, best-effort named-gate rewrite pass; and
//! - [`BasisClosurePass`], the final lowering step that guarantees the emitted
//!   circuit contains only gates from a supported target basis.
//!
//! Mathematical universality and implemented synthesis support are deliberately
//! separate checks. For example, `{H, T, CX}` is universal, but Q-Rust does not
//! yet implement approximate Clifford+T synthesis for arbitrary rotations, so
//! it is not accepted by [`BasisClosurePass`]. The exact output families
//! currently supported are `{U, CX/CZ}`, `{RZ, RX, CX/CZ}`,
//! `{RZ, SX, CX/CZ}`, and `{RZ, H, CX/CZ}` (additional named gates may
//! also be present).

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashSet;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Universality validation
// ---------------------------------------------------------------------------

/// Entanglers recognized by the conservative universality preflight.
const ENTANGLING_GATES: &[&str] = &[
    "cx", "cz", "cy", "ch", "csx", "ecr", "iswap", "dcx", "rzz", "rxx", "ryy", "crx", "cry", "crz",
];

fn canonical_gate_name(name: &str) -> String {
    match name.to_ascii_lowercase().as_str() {
        "cnot" => "cx".into(),
        "u3" => "u".into(),
        "u1" | "p" => "rz".into(),
        other => other.into(),
    }
}

fn canonical_basis(basis: &HashSet<String>) -> HashSet<String> {
    basis.iter().map(|name| canonical_gate_name(name)).collect()
}

fn sorted_basis(basis: &HashSet<String>) -> Vec<String> {
    let mut names: Vec<_> = basis.iter().cloned().collect();
    names.sort();
    names
}

/// Conservatively validates that `basis` matches a known universal family.
///
/// This is not an exhaustive decision procedure. It recognizes an entangler
/// together with one of the following single-qubit resources:
///
/// - an arbitrary `U`/`U3` gate;
/// - rotations about two distinct axes;
/// - arbitrary `RZ` plus `H` or `SX`; or
/// - the discrete Clifford+T pair `H` and `T`/`Tdg`.
///
/// A single continuous axis is insufficient: `{RZ, CX}`, for example, cannot
/// prepare an arbitrary one-qubit state. A Clifford-only set such as
/// `{H, S, CX}` is also rejected.
///
/// # Examples
/// ```rust
/// use q_rust::transpiler::target_basis::validate_universality;
/// use std::collections::HashSet;
///
/// let basis: HashSet<String> = ["rz", "rx", "cx"].iter().map(|s| s.to_string()).collect();
/// assert!(validate_universality(&basis).is_ok());
///
/// let clifford_only: HashSet<String> = ["h", "s", "cx"].iter().map(|s| s.to_string()).collect();
/// assert!(validate_universality(&clifford_only).is_err());
/// ```
pub fn validate_universality(basis: &HashSet<String>) -> Result<()> {
    let lower = canonical_basis(basis);
    let basis_names = sorted_basis(&lower);

    let has_entangler = ENTANGLING_GATES.iter().any(|&g| lower.contains(g));
    let rotation_axes = ["rx", "ry", "rz"]
        .iter()
        .filter(|&&axis| lower.contains(axis))
        .count();
    let has_universal_1q = lower.contains("u")
        || rotation_axes >= 2
        || (lower.contains("rz") && (lower.contains("h") || lower.contains("sx")))
        || (lower.contains("rx") && lower.contains("h"))
        || (lower.contains("h") && (lower.contains("t") || lower.contains("tdg")));

    if !has_entangler {
        return Err(QRustError::NonUniversalBasisGateSet {
            reason: format!(
                "no entangling 2-qubit gate found in basis {basis_names:?}; \
                 need at least one of {ENTANGLING_GATES:?}"
            ),
        });
    }

    if !has_universal_1q {
        return Err(QRustError::NonUniversalBasisGateSet {
            reason: format!(
                "single-qubit gates in basis {basis_names:?} do not match a recognized \
                 universal family; provide U, two independent rotation axes, \
                 RZ with H/SX, or H with T/Tdg"
            ),
        });
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Equivalence library (data-driven)
// ---------------------------------------------------------------------------

/// A single rewrite rule: expand gate `from` into a sequence of operations.
/// Parameters in `ops` use symbolic indices into the original gate's param
/// list; qubits use symbolic indices into the original gate's qubit list.
struct RewriteRule {
    /// Lower-case name of the gate to rewrite.
    from: &'static str,
    /// Sequence of (gate_name, qubit_indices, param_spec) triples.
    /// `param_spec` encodes how to compute each parameter:
    ///   `ParamSpec::Const(v)` — a fixed constant
    ///   `ParamSpec::Passthrough(i)` — take params[i] from the incoming gate
    ops: Vec<(&'static str, Vec<usize>, Vec<ParamSpec>)>,
}

#[derive(Clone)]
enum ParamSpec {
    Const(f64),
    Passthrough(usize),
    /// params[i] + constant offset.
    Sum(usize, f64),
}

impl ParamSpec {
    fn eval(&self, params: &[f64]) -> f64 {
        match self {
            ParamSpec::Const(v) => *v,
            ParamSpec::Passthrough(i) => params[*i],
            ParamSpec::Sum(i, off) => params[*i] + off,
        }
    }
}

/// Builds the equivalence library appropriate for a given target basis.
/// Each rule is only included when ALL its output gate names appear in the
/// basis — ensuring we only add rules that terminate.
fn build_equivalence_library(basis: &HashSet<String>) -> Vec<RewriteRule> {
    let lower = canonical_basis(basis);

    let has = |g: &str| lower.contains(g);

    let mut rules = Vec::new();

    // H = Rz(π/2) · Rx(π/2) · Rz(π/2)  [exact, no global phase]
    // Convention: Rz(φ) = diag(1, e^{iφ}), Rx(θ) = [[cos,⋅-i·sin],[-i·sin,cos]].
    // Verified: Rz(π/2)·Rx(π/2)·Rz(π/2) = (1/√2)[[1,1],[1,-1]] = H.
    if has("rz") && has("rx") {
        rules.push(RewriteRule {
            from: "h",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)]),
                ("rx", vec![0], vec![ParamSpec::Const(PI / 2.0)]),
                ("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)]),
            ],
        });
        // X = Rx(π)
        rules.push(RewriteRule {
            from: "x",
            ops: vec![("rx", vec![0], vec![ParamSpec::Const(PI)])],
        });
        // Y = Rx(π) · Rz(π)  (up to global phase)
        rules.push(RewriteRule {
            from: "y",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Const(PI)]),
                ("rx", vec![0], vec![ParamSpec::Const(PI)]),
            ],
        });
        // Z = Rz(π)
        rules.push(RewriteRule {
            from: "z",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI)])],
        });
        // S = Rz(π/2)
        rules.push(RewriteRule {
            from: "s",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)])],
        });
        // Sdg = Rz(-π/2)
        rules.push(RewriteRule {
            from: "sdg",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(-PI / 2.0)])],
        });
        // T = Rz(π/4)
        rules.push(RewriteRule {
            from: "t",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI / 4.0)])],
        });
        // Tdg = Rz(-π/4)
        rules.push(RewriteRule {
            from: "tdg",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(-PI / 4.0)])],
        });
        // In matrix order, RY(θ) = Rz(π/2) · Rx(θ) · Rz(-π/2).
        // Operations below are stored in circuit order (rightmost first).
        rules.push(RewriteRule {
            from: "ry",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Const(-PI / 2.0)]),
                ("rx", vec![0], vec![ParamSpec::Passthrough(0)]),
                ("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)]),
            ],
        });
        // U(θ,φ,λ) ZYZ form: Rz(φ) · Ry(θ) · Rz(λ)
        // Substituting Ry(θ) = Rz(π/2) · Rx(θ) · Rz(-π/2) gives:
        //   = Rz(φ + π/2) · Rx(θ) · Rz(λ - π/2)
        // Operations are stored in circuit order (rightmost matrix first).
        // params[0]=θ, params[1]=φ, params[2]=λ
        rules.push(RewriteRule {
            from: "u",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Sum(2, -PI / 2.0)]), // Rz(λ - π/2) first
                ("rx", vec![0], vec![ParamSpec::Passthrough(0)]),    // Rx(θ)
                ("rz", vec![0], vec![ParamSpec::Sum(1, PI / 2.0)]),  // Rz(φ + π/2) last
            ],
        });
    }

    // If basis has rz only (no rx), use ZYZ: RY replaceable by {rx,rz}
    // The canonical ZYZ form is used for U gate expansion above when rx is absent.
    if has("rz") && !has("rx") {
        // Z = Rz(π)
        rules.push(RewriteRule {
            from: "z",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI)])],
        });
        // S = Rz(π/2)
        rules.push(RewriteRule {
            from: "s",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)])],
        });
        // Sdg = Rz(-π/2)
        rules.push(RewriteRule {
            from: "sdg",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(-PI / 2.0)])],
        });
        // T = Rz(π/4)
        rules.push(RewriteRule {
            from: "t",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(PI / 4.0)])],
        });
        // Tdg = Rz(-π/4)
        rules.push(RewriteRule {
            from: "tdg",
            ops: vec![("rz", vec![0], vec![ParamSpec::Const(-PI / 4.0)])],
        });
        // RX(θ) = Rz(-π/2) · Ry(θ) · Rz(π/2)  — note ry must also be in basis
        // We skip this when rx is unavailable and ry is unavailable.
    }

    // CZ <-> CX (with H sandwich on target)
    if has("cx") && !has("cz") {
        rules.push(RewriteRule {
            from: "cz",
            ops: vec![
                ("h", vec![1], vec![]),
                ("cx", vec![0, 1], vec![]),
                ("h", vec![1], vec![]),
            ],
        });
    }
    if has("cz") && !has("cx") {
        rules.push(RewriteRule {
            from: "cx",
            ops: vec![
                ("h", vec![1], vec![]),
                ("cz", vec![0, 1], vec![]),
                ("h", vec![1], vec![]),
            ],
        });
    }

    // SWAP = CX(0,1) · CX(1,0) · CX(0,1)
    if has("cx") {
        rules.push(RewriteRule {
            from: "swap",
            ops: vec![
                ("cx", vec![0, 1], vec![]),
                ("cx", vec![1, 0], vec![]),
                ("cx", vec![0, 1], vec![]),
            ],
        });
        // CRZ: ParamSpec does not support θ/2 (division).
        // We let CRZ fall through to BasisDecompositionPass, which calls
        // GateType::CRZ.decompose() — a correct Rz(θ/2)·CX·Rz(-θ/2)·CX decomposition.
    }

    rules
}

// ---------------------------------------------------------------------------
// TargetBasisPass
// ---------------------------------------------------------------------------

/// Applies early, exact named-gate rewrites selected by `basis`.
///
/// This pass intentionally does **not** promise basis closure: unmatched gates
/// are retained for analytic decomposition/KAK, and some rewrites introduce
/// temporary gates such as `H`. Run [`BasisClosurePass`] after decomposition to
/// obtain and validate a circuit containing only target-basis gates.
#[derive(Debug, Clone)]
pub struct TargetBasisPass {
    pub basis: HashSet<String>,
}

impl TargetBasisPass {
    /// Creates a new `TargetBasisPass` and validates universality.
    pub fn new(basis: HashSet<String>) -> Result<Self> {
        validate_universality(&basis)?;
        Ok(Self {
            basis: canonical_basis(&basis),
        })
    }
}

impl Pass for TargetBasisPass {
    fn name(&self) -> &str {
        "TargetBasisPass"
    }

    fn run(&self, circuit: &Circuit, _props: &mut PropertySet) -> Circuit {
        let lower = canonical_basis(&self.basis);
        let lib = build_equivalence_library(&self.basis);

        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } => {
                    let gate_name = name.to_qasm_name().to_lowercase();
                    if lower.contains(&gate_name) {
                        // Already in target basis — pass through unchanged.
                        out.add_op(op.clone());
                    } else if let Some(rule) = lib.iter().find(|r| r.from == gate_name) {
                        // Rewrite using the equivalence rule.
                        for (sub_name_str, sub_q_idx, sub_params_spec) in &rule.ops {
                            let sub_q: Vec<usize> = sub_q_idx.iter().map(|&i| qubits[i]).collect();
                            let sub_p: Vec<f64> =
                                sub_params_spec.iter().map(|ps| ps.eval(params)).collect();
                            let sub_name: GateType = sub_name_str
                                .parse()
                                .unwrap_or_else(|_| GateType::Custom(sub_name_str.to_string()));
                            out.add_op(Operation::Gate {
                                name: sub_name,
                                qubits: sub_q,
                                params: sub_p,
                            });
                        }
                    } else {
                        // No early rule: a later decomposition/synthesis pass
                        // handles the residual gate before basis closure.
                        out.add_op(op.clone());
                    }
                }
                other => out.add_op(other.clone()),
            }
        }

        out
    }
}

// ---------------------------------------------------------------------------
// Final basis closure
// ---------------------------------------------------------------------------

/// Validates that Q-Rust has an exact final lowering for `basis`.
///
/// Universality alone is not enough: `{H, T, CX}` is universal only with an
/// approximation algorithm, which is outside the current implementation.
pub fn validate_exact_translation_support(basis: &HashSet<String>) -> Result<()> {
    validate_universality(basis)?;
    let basis = canonical_basis(basis);

    if !basis.contains("cx") && !basis.contains("cz") {
        return Err(QRustError::UntranslatableGate {
            gate: "cx".into(),
            basis: sorted_basis(&basis),
        });
    }

    let has_exact_1q = basis.contains("u")
        || (basis.contains("rz") && basis.contains("rx"))
        || (basis.contains("rz") && basis.contains("sx"))
        || (basis.contains("rz") && basis.contains("h"));
    if !has_exact_1q {
        return Err(QRustError::UntranslatableGate {
            gate: "u".into(),
            basis: sorted_basis(&basis),
        });
    }

    Ok(())
}

fn gate(name: GateType, qubits: Vec<usize>, params: Vec<f64>) -> Operation {
    Operation::Gate {
        name,
        qubits,
        params,
    }
}

fn lower_u_to_basis(q: usize, params: &[f64], basis: &HashSet<String>) -> Vec<Operation> {
    let theta = params.first().copied().unwrap_or(0.0);
    let phi = params.get(1).copied().unwrap_or(0.0);
    let lambda = params.get(2).copied().unwrap_or(0.0);

    if basis.contains("u") {
        return vec![gate(GateType::U, vec![q], vec![theta, phi, lambda])];
    }

    if basis.contains("rz") && basis.contains("rx") {
        return vec![
            gate(GateType::RZ, vec![q], vec![lambda - PI / 2.0]),
            gate(GateType::RX, vec![q], vec![theta]),
            gate(GateType::RZ, vec![q], vec![phi + PI / 2.0]),
        ];
    }

    if basis.contains("rz") && basis.contains("sx") {
        return vec![
            gate(GateType::RZ, vec![q], vec![lambda]),
            gate(GateType::SX, vec![q], vec![]),
            gate(GateType::RZ, vec![q], vec![theta + PI]),
            gate(GateType::SX, vec![q], vec![]),
            gate(GateType::RZ, vec![q], vec![phi + PI]),
        ];
    }

    // H RZ(theta) H = RX(theta), up to global phase.
    vec![
        gate(GateType::RZ, vec![q], vec![lambda - PI / 2.0]),
        gate(GateType::H, vec![q], vec![]),
        gate(GateType::RZ, vec![q], vec![theta]),
        gate(GateType::H, vec![q], vec![]),
        gate(GateType::RZ, vec![q], vec![phi + PI / 2.0]),
    ]
}

fn close_operation(op: &Operation, basis: &HashSet<String>) -> Vec<Operation> {
    match op {
        Operation::Gate {
            name: GateType::U,
            qubits,
            params,
        } if qubits.len() == 1 => lower_u_to_basis(qubits[0], params, basis),
        Operation::Gate {
            name: GateType::CX,
            qubits,
            ..
        } if qubits.len() == 2 && !basis.contains("cx") && basis.contains("cz") => {
            let h = [PI / 2.0, 0.0, PI];
            let mut out = lower_u_to_basis(qubits[1], &h, basis);
            out.push(gate(GateType::CZ, qubits.clone(), vec![]));
            out.extend(lower_u_to_basis(qubits[1], &h, basis));
            out
        }
        Operation::Conditional { condition, op } => close_operation(op, basis)
            .into_iter()
            .map(|inner| Operation::Conditional {
                condition: condition.clone(),
                op: Box::new(inner),
            })
            .collect(),
        other => vec![other.clone()],
    }
}

/// Final exact lowering from the internal `{U, CX}` representation to a
/// supported target basis.
#[derive(Debug, Clone)]
pub struct BasisClosurePass {
    pub basis: HashSet<String>,
}

impl BasisClosurePass {
    pub fn new(basis: HashSet<String>) -> Result<Self> {
        validate_exact_translation_support(&basis)?;
        Ok(Self {
            basis: canonical_basis(&basis),
        })
    }
}

impl Pass for BasisClosurePass {
    fn name(&self) -> &str {
        "BasisClosurePass"
    }

    fn run(&self, circuit: &Circuit, _props: &mut PropertySet) -> Circuit {
        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        for op in &circuit.operations {
            for lowered in close_operation(op, &self.basis) {
                out.add_op(lowered);
            }
        }
        out
    }
}

/// Checks that every emitted gate is a member of `basis`.
pub fn validate_circuit_basis(circuit: &Circuit, basis: &HashSet<String>) -> Result<()> {
    let basis = canonical_basis(basis);

    fn check_op(op: &Operation, basis: &HashSet<String>) -> Option<String> {
        match op {
            Operation::Gate { name, .. } => {
                let name = canonical_gate_name(name.to_qasm_name());
                (!basis.contains(&name)).then_some(name)
            }
            Operation::Conditional { op, .. } => check_op(op, basis),
            _ => None,
        }
    }

    if let Some(gate) = circuit
        .operations
        .iter()
        .find_map(|op| check_op(op, &basis))
    {
        return Err(QRustError::UntranslatableGate {
            gate,
            basis: sorted_basis(&basis),
        });
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn basis(gates: &[&str]) -> HashSet<String> {
        gates.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_universality_rejects_single_rotation_axis() {
        assert!(validate_universality(&basis(&["rz", "cx"])).is_err());
    }

    #[test]
    fn test_universality_accepts_two_rotation_axes() {
        assert!(validate_universality(&basis(&["rz", "rx", "cz"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_ibm_style_basis() {
        assert!(validate_universality(&basis(&["rz", "sx", "x", "cx"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_rz_h_basis() {
        assert!(validate_universality(&basis(&["rz", "h", "cz"])).is_ok());
        assert!(validate_exact_translation_support(&basis(&["rz", "h", "cz"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_u_cx() {
        assert!(validate_universality(&basis(&["u", "cx"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_h_t_cx() {
        // Mathematical universality is distinct from implemented exact
        // Clifford+T approximation support.
        assert!(validate_universality(&basis(&["h", "t", "cx"])).is_ok());
        assert!(matches!(
            validate_exact_translation_support(&basis(&["h", "t", "cx"])),
            Err(QRustError::UntranslatableGate { .. })
        ));
    }

    #[test]
    fn test_universality_rejects_clifford_only() {
        // {H, S, CX} is Clifford-only and NOT universal.
        let err = validate_universality(&basis(&["h", "s", "cx"])).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("non-Clifford") || msg.contains("universal"),
            "{msg}"
        );
    }

    #[test]
    fn test_universality_rejects_no_entangler() {
        let err = validate_universality(&basis(&["rz", "rx", "h"])).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("entangling") || msg.contains("2-qubit"),
            "{msg}"
        );
    }

    #[test]
    fn test_target_basis_h_to_rz_rx() {
        let pass = TargetBasisPass::new(basis(&["rz", "rx", "cx"])).unwrap();
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        // H must expand to RZ+RX sequence
        assert!(out.operations.len() >= 2);
        for op in &out.operations {
            if let Operation::Gate { name, .. } = op {
                let n = name.to_qasm_name().to_lowercase();
                assert!(
                    n == "rz" || n == "rx",
                    "unexpected gate {n} after H expansion"
                );
            }
        }
    }

    #[test]
    fn test_target_basis_cz_to_cx_h() {
        let pass = TargetBasisPass::new(basis(&["rz", "rx", "cx"])).unwrap();
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CZ,
            qubits: vec![0, 1],
            params: vec![],
        });
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        // CZ → H · CX · H
        let names: Vec<_> = out
            .operations
            .iter()
            .filter_map(|op| {
                if let Operation::Gate { name, .. } = op {
                    Some(name.to_qasm_name().to_lowercase())
                } else {
                    None
                }
            })
            .collect();
        assert!(
            names.contains(&"cx".to_string()),
            "expected CX in {names:?}"
        );
        assert!(names.contains(&"h".to_string()), "expected H in {names:?}");
    }

    #[test]
    fn test_basis_closure_rz_rx_cx_is_closed_and_equivalent() {
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};

        let mut c = Circuit::new(2, 0);
        c.add_op(gate(GateType::U, vec![0], vec![0.7, -0.2, 1.1]));
        c.add_op(gate(GateType::CX, vec![0, 1], vec![]));

        let pass = BasisClosurePass::new(basis(&["rz", "rx", "cx"])).unwrap();
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        validate_circuit_basis(&out, &pass.basis).unwrap();
        let fid = unitary_fidelity(&circuit_to_unitary(&c), &circuit_to_unitary(&out));
        assert!((fid - 1.0).abs() < 1e-9, "fidelity = {fid}");
    }

    #[test]
    fn test_basis_closure_rz_sx_cx_is_closed_and_equivalent() {
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};

        let mut c = Circuit::new(2, 0);
        c.add_op(gate(GateType::U, vec![0], vec![0.7, -0.2, 1.1]));
        c.add_op(gate(GateType::CX, vec![0, 1], vec![]));

        let pass = BasisClosurePass::new(basis(&["rz", "sx", "x", "cx"])).unwrap();
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        validate_circuit_basis(&out, &pass.basis).unwrap();
        let fid = unitary_fidelity(&circuit_to_unitary(&c), &circuit_to_unitary(&out));
        assert!((fid - 1.0).abs() < 1e-9, "fidelity = {fid}");
    }

    #[test]
    fn test_basis_closure_rz_h_cx_is_closed_and_equivalent() {
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};

        let mut c = Circuit::new(2, 0);
        c.add_op(gate(GateType::U, vec![0], vec![0.7, -0.2, 1.1]));
        c.add_op(gate(GateType::CX, vec![0, 1], vec![]));

        let pass = BasisClosurePass::new(basis(&["rz", "h", "cx"])).unwrap();
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        validate_circuit_basis(&out, &pass.basis).unwrap();
        let fid = unitary_fidelity(&circuit_to_unitary(&c), &circuit_to_unitary(&out));
        assert!((fid - 1.0).abs() < 1e-9, "fidelity = {fid}");
    }

    #[test]
    fn test_basis_closure_can_rebase_cx_to_cz() {
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};

        let mut c = Circuit::new(2, 0);
        c.add_op(gate(GateType::CX, vec![0, 1], vec![]));
        let pass = BasisClosurePass::new(basis(&["u", "cz"])).unwrap();
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        validate_circuit_basis(&out, &pass.basis).unwrap();
        let fid = unitary_fidelity(&circuit_to_unitary(&c), &circuit_to_unitary(&out));
        assert!((fid - 1.0).abs() < 1e-9, "fidelity = {fid}");
    }

    #[test]
    fn test_target_basis_passthrough_when_in_basis() {
        let pass = TargetBasisPass::new(basis(&["rz", "rx", "cx", "h"])).unwrap();
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let mut props = PropertySet::new();
        let out = pass.run(&c, &mut props);
        // H is already in basis → exactly 1 H gate, unchanged.
        assert_eq!(out.operations.len(), 1);
        assert!(matches!(
            &out.operations[0],
            Operation::Gate {
                name: GateType::H,
                ..
            }
        ));
    }

    #[test]
    fn test_no_basis_specified_error_variant_exists() {
        // Smoke-test that the error type compiles and formats.
        let e = QRustError::NoBasisSpecified;
        assert!(!e.to_string().is_empty());
    }
}
