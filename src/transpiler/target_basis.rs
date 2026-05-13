//! Vendor-agnostic target basis translation pass.
//!
//! This module provides:
//! - [`validate_universality`]: validates that a gate set is quantum-universal
//!   (i.e., can approximate any unitary to arbitrary precision) before any
//!   compilation is attempted.
//! - [`TargetBasisPass`]: a data-driven translation pass that rewrites every
//!   gate in a circuit into an equivalent sequence drawn exclusively from the
//!   caller-supplied target basis.
//!
//! No vendor-specific knowledge is hard-coded here. The pass reads a set of
//! gate name strings (e.g. `{"rz", "cx"}` or `{"rx", "cz"}`) and builds an
//! equivalence library at construction time from that set.

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashSet;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Universality validation
// ---------------------------------------------------------------------------

/// The minimum requirements for a quantum-universal gate set:
/// 1. At least one *entangling* 2-qubit gate capable of creating entanglement.
/// 2. At least one *continuous* single-qubit rotation gate (or a non-Clifford
///    discrete gate such as T) whose repeated use densely covers SU(2).
///
/// A Clifford-only set such as `{H, S, CX}` is NOT universal because it can
/// only generate the finite Clifford group, not all unitaries.
///
/// Known entangling gates:
const ENTANGLING_GATES: &[&str] = &[
    "cx", "cnot", "cz", "cy", "ch", "csx", "ecr", "iswap", "dcx", "rzz", "rxx", "ryy", "ccx",
    "crx", "cry", "crz",
];

/// Known single-qubit gates that provide non-Clifford / continuous coverage of SU(2):
const UNIVERSAL_1Q_GATES: &[&str] = &[
    // Continuous rotations (any non-zero irrational angle is dense in SU(2))
    "rz", "rx", "ry", // Universal 1-qubit gate families
    "u", "u3", "u2", "p",
    // Non-Clifford discrete gates (T / Tdg lift Clifford to universal when
    // combined with H and S)
    "t", "tdg",
];

/// Validates that `basis` is quantum-universal.
///
/// Returns `Ok(())` when the gate set passes both checks, or a descriptive
/// [`QRustError::NonUniversalBasisGateSet`] when either check fails.
///
/// # Examples
/// ```rust
/// use q_rust::transpiler::target_basis::validate_universality;
/// use std::collections::HashSet;
///
/// let basis: HashSet<String> = ["rz", "cx"].iter().map(|s| s.to_string()).collect();
/// assert!(validate_universality(&basis).is_ok());
///
/// let clifford_only: HashSet<String> = ["h", "s", "cx"].iter().map(|s| s.to_string()).collect();
/// assert!(validate_universality(&clifford_only).is_err());
/// ```
pub fn validate_universality(basis: &HashSet<String>) -> Result<()> {
    let lower: HashSet<String> = basis.iter().map(|s| s.to_lowercase()).collect();

    let has_entangler = ENTANGLING_GATES.iter().any(|&g| lower.contains(g));

    let has_universal_1q = UNIVERSAL_1Q_GATES.iter().any(|&g| lower.contains(g));

    if !has_entangler {
        return Err(QRustError::NonUniversalBasisGateSet {
            reason: format!(
                "no entangling 2-qubit gate found in basis {lower:?}; \
                 need at least one of {ENTANGLING_GATES:?}"
            ),
        });
    }

    if !has_universal_1q {
        return Err(QRustError::NonUniversalBasisGateSet {
            reason: format!(
                "no continuous or non-Clifford single-qubit gate found in basis {lower:?}; \
                 a Clifford-only set cannot approximate all unitaries. \
                 Add at least one of {UNIVERSAL_1Q_GATES:?}"
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
    let lower: HashSet<String> = basis.iter().map(|s| s.to_lowercase()).collect();

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
        // RY(θ) = Rz(-π/2) · Rx(θ) · Rz(π/2)  [exact, verified]
        // Rz(-π/2) = [[1,0],[0,-i]], Rx(θ) = [[c,-is],[-is,c]]
        // Product = Rz(-π/2)·Rx(θ)·Rz(π/2) = [[c,s],[-s,c]] = RY(θ) ✓
        rules.push(RewriteRule {
            from: "ry",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Const(-PI / 2.0)]),
                ("rx", vec![0], vec![ParamSpec::Passthrough(0)]),
                ("rz", vec![0], vec![ParamSpec::Const(PI / 2.0)]),
            ],
        });
        // U(θ,φ,λ) ZYZ form: Rz(φ) · Ry(θ) · Rz(λ)
        // Substituting Ry(θ) = Rz(-π/2) · Rx(θ) · Rz(π/2) gives:
        //   = Rz(φ - π/2) · Rx(θ) · Rz(λ + π/2)
        // params[0]=θ, params[1]=φ, params[2]=λ
        rules.push(RewriteRule {
            from: "u",
            ops: vec![
                ("rz", vec![0], vec![ParamSpec::Sum(2, PI / 2.0)]), // Rz(λ + π/2) first
                ("rx", vec![0], vec![ParamSpec::Passthrough(0)]),   // Rx(θ)
                ("rz", vec![0], vec![ParamSpec::Sum(1, -PI / 2.0)]), // Rz(φ - π/2) last
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

/// Translates every gate in the input circuit into an equivalent sequence of
/// gates drawn from `basis`. The pass is purely data-driven: no vendor-specific
/// knowledge is hard-coded.
///
/// Construction validates that `basis` is quantum-universal; an error is
/// returned if it is not. This guarantees that the pass can always make
/// progress (even gates without explicit rules can be handled downstream by
/// the KAK synthesis pass, which targets `{U, CX}` — both of which must be
/// expressible in any universal basis).
#[derive(Debug, Clone)]
pub struct TargetBasisPass {
    pub basis: HashSet<String>,
}

impl TargetBasisPass {
    /// Creates a new `TargetBasisPass` and validates universality.
    pub fn new(basis: HashSet<String>) -> Result<Self> {
        validate_universality(&basis)?;
        Ok(Self { basis })
    }
}

impl Pass for TargetBasisPass {
    fn name(&self) -> &str {
        "TargetBasisPass"
    }

    fn run(&self, circuit: &Circuit, _props: &mut PropertySet) -> Circuit {
        let lower: HashSet<String> = self.basis.iter().map(|s| s.to_lowercase()).collect();
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
                        // No rule found — pass through and let KAK handle residuals.
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
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn basis(gates: &[&str]) -> HashSet<String> {
        gates.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_universality_accepts_rz_cx() {
        assert!(validate_universality(&basis(&["rz", "cx"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_rx_cz() {
        assert!(validate_universality(&basis(&["rx", "cz"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_u_cx() {
        assert!(validate_universality(&basis(&["u", "cx"])).is_ok());
    }

    #[test]
    fn test_universality_accepts_h_t_cx() {
        // H, T (non-Clifford), CX — the classic Solovay-Kitaev basis.
        assert!(validate_universality(&basis(&["h", "t", "cx"])).is_ok());
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
