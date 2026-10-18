//! Circuit-equivalence verification.
//!
//! Exact checks compare full unitaries with the phase-invariant normalized
//! Hilbert--Schmidt (process) fidelity
//! `|tr(U†V)|² / d²`. Statistical checks compare output-state fidelities for
//! seeded Haar-random inputs; they provide finite-sample evidence, not a proof.
//!
//! [`verify_equivalence`] assumes both circuits use the same qubit register and
//! ordering. Routed circuits may use a larger physical register and a changed
//! output permutation; use [`verify_equivalence_with_layout`] for those.
//!
//! Matching terminal measurements are supported by comparing their classical
//! destinations (and the final layout), then checking their unitary prefixes.
//! Reset, classical control, mid-circuit measurements, and different readout
//! maps return an error. This is not a general quantum-channel equivalence test.
//! "Exact" denotes a full floating-point operator comparison, not a symbolic
//! or formal proof. Gate-matrix implementations still need independent tests.
//!
//! The automatic direct-verification policy is:
//!
//! 1. **n ≤ 14** — exact full-unitary process fidelity.
//! 2. **14 < n ≤ 22** — statistical: Haar-random sampling via
//!    [`equivalence_by_sampling`] (8 samples by default).
//! 3. **n > 22** — out of reach for both methods; returns
//!    [`Verdict::Unverifiable`].

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, Operation};
use crate::simulator::{
    equivalence_by_sampling, equivalence_by_sampling_with_layout, try_circuit_to_unitary,
    try_extract_logical_unitary, unitary_fidelity, validate_layout, MAX_QUBITS,
    MAX_STATE_VECTOR_QUBITS,
};
use std::collections::BTreeMap;

/// Threshold above which the exact unitary check is too memory-intensive.
pub const EXACT_VERIFY_LIMIT: usize = MAX_QUBITS;

/// Default number of Haar-random samples for statistical equivalence.
pub const DEFAULT_SAMPLE_COUNT: usize = 8;

/// Default fidelity threshold for the exact path.
pub const DEFAULT_EXACT_TOLERANCE: f64 = 1e-9;

/// Default fidelity threshold for the statistical path.
pub const DEFAULT_SAMPLING_TOLERANCE: f64 = 1e-6;

/// Conservative physical-width limit for layout-aware state-vector checks.
///
/// The simulator supports wider states, but verification keeps headroom for
/// the logical state, routed state, and scratch allocations held together.
pub const MAX_LAYOUT_VERIFY_QUBITS: usize = 20;

/// Resource policy for the high-level verifier. Lower the dense limit for
/// benchmark or memory-constrained runs instead of bypassing semantic checks.
/// The defaults preserve the original verification policy. These are upper
/// guards, not promises that every permitted dense problem fits in memory.
#[derive(Debug, Clone, PartialEq)]
pub struct VerificationConfig {
    /// Maximum width for full-operator comparison (at most 14).
    pub max_exact_qubits: usize,
    /// Maximum width for sampling (at most 22). Layout-aware verification
    /// additionally caps the physical width at [`MAX_LAYOUT_VERIFY_QUBITS`].
    pub max_sampling_qubits: usize,
    pub samples: usize,
    pub seed: u64,
    pub exact_tolerance: f64,
    pub sampling_tolerance: f64,
}

impl Default for VerificationConfig {
    fn default() -> Self {
        Self {
            max_exact_qubits: EXACT_VERIFY_LIMIT,
            max_sampling_qubits: MAX_STATE_VECTOR_QUBITS - 2,
            samples: DEFAULT_SAMPLE_COUNT,
            seed: 0x00C0_FFEE_DEAD_BEEF_u64,
            exact_tolerance: DEFAULT_EXACT_TOLERANCE,
            sampling_tolerance: DEFAULT_SAMPLING_TOLERANCE,
        }
    }
}

impl VerificationConfig {
    fn validate(&self) -> Result<()> {
        validate_verification_parameters(
            self.samples,
            self.exact_tolerance,
            self.sampling_tolerance,
        )?;
        if self.max_exact_qubits > EXACT_VERIFY_LIMIT
            || self.max_sampling_qubits > MAX_STATE_VECTOR_QUBITS - 2
        {
            return Err(QRustError::Simulation(format!(
                "verification limits exceed supported maxima (exact {EXACT_VERIFY_LIMIT}, sampling {})",
                MAX_STATE_VECTOR_QUBITS - 2,
            )));
        }
        Ok(())
    }
}

/// Outcome of an equivalence check.
#[derive(Debug, Clone, PartialEq)]
pub enum Verdict {
    /// Exact unitary fidelity ≥ `1 - tol`.
    ExactlyEquivalent { fidelity: f64 },
    /// Min sampled fidelity ≥ `1 - tol` over `samples` Haar-random states.
    StatisticallyEquivalent { min_fidelity: f64, samples: usize },
    /// Fidelity below threshold.
    NotEquivalent { fidelity: f64, method: &'static str },
    /// Circuit too large for any verification method.
    Unverifiable { reason: String },
}

impl Verdict {
    pub fn is_equivalent(&self) -> bool {
        matches!(
            self,
            Verdict::ExactlyEquivalent { .. } | Verdict::StatisticallyEquivalent { .. }
        )
    }

    pub fn describe(&self) -> String {
        match self {
            Verdict::ExactlyEquivalent { fidelity } => {
                format!("exactly equivalent (fidelity = {fidelity:.10})")
            }
            Verdict::StatisticallyEquivalent {
                min_fidelity,
                samples,
            } => format!(
                "statistically equivalent over {samples} Haar samples \
                 (min fidelity = {min_fidelity:.10})"
            ),
            Verdict::NotEquivalent { fidelity, method } => {
                format!("NOT equivalent ({method}, fidelity = {fidelity:.10})")
            }
            Verdict::Unverifiable { reason } => format!("unverifiable: {reason}"),
        }
    }
}

/// Verifies circuits on the same qubit register and in the same qubit order.
///
/// Routed output with a non-trivial layout must instead use
/// [`verify_equivalence_with_layout`].
pub fn verify_equivalence(c1: &Circuit, c2: &Circuit) -> Result<Verdict> {
    verify_equivalence_with(
        c1,
        c2,
        DEFAULT_SAMPLE_COUNT,
        0x00C0_FFEE_DEAD_BEEF_u64,
        DEFAULT_EXACT_TOLERANCE,
        DEFAULT_SAMPLING_TOLERANCE,
    )
}

/// Configurable variant of [`verify_equivalence`].
pub fn verify_equivalence_with(
    c1: &Circuit,
    c2: &Circuit,
    samples: usize,
    seed: u64,
    exact_tol: f64,
    sampling_tol: f64,
) -> Result<Verdict> {
    verify_equivalence_configured(
        c1,
        c2,
        &VerificationConfig {
            samples,
            seed,
            exact_tolerance: exact_tol,
            sampling_tolerance: sampling_tol,
            ..VerificationConfig::default()
        },
    )
}

/// Verify with explicit resource limits, retaining all readout and operation
/// validation performed by [`verify_equivalence`].
pub fn verify_equivalence_configured(
    c1: &Circuit,
    c2: &Circuit,
    config: &VerificationConfig,
) -> Result<Verdict> {
    config.validate()?;
    if c1.num_qubits != c2.num_qubits {
        return Err(QRustError::Simulation(format!(
            "verify_equivalence: qubit count mismatch ({} vs {})",
            c1.num_qubits, c2.num_qubits
        )));
    }
    let n = c1.num_qubits;
    let (c1, readout1) = unitary_prefix_and_readout(c1, None)?;
    let (c2, readout2) = unitary_prefix_and_readout(c2, None)?;
    validate_matching_readout(&c1, &c2, &readout1, &readout2)?;

    if n <= config.max_exact_qubits {
        let u1 = try_circuit_to_unitary(&c1)?;
        let u2 = try_circuit_to_unitary(&c2)?;
        let fid = unitary_fidelity(&u1, &u2);
        if fid >= 1.0 - config.exact_tolerance {
            return Ok(Verdict::ExactlyEquivalent { fidelity: fid });
        }
        return Ok(Verdict::NotEquivalent {
            fidelity: fid,
            method: "exact",
        });
    }

    // Leave a 2-qubit safety margin under MAX_STATE_VECTOR_QUBITS=24:
    // 24-qubit state-vector evolution costs ~256 MiB per copy, and we
    // need two copies plus a Haar sample.
    if n <= config.max_sampling_qubits {
        let min_fid = equivalence_by_sampling(&c1, &c2, config.samples, config.seed)?;
        if min_fid >= 1.0 - config.sampling_tolerance {
            return Ok(Verdict::StatisticallyEquivalent {
                min_fidelity: min_fid,
                samples: config.samples,
            });
        }
        return Ok(Verdict::NotEquivalent {
            fidelity: min_fid,
            method: "sampling",
        });
    }

    Ok(Verdict::Unverifiable {
        reason: format!(
            "{n} qubits exceeds both configured exact ({}) and sampling ({}) limits",
            config.max_exact_qubits, config.max_sampling_qubits,
        ),
    })
}

fn validate_verification_parameters(
    samples: usize,
    exact_tol: f64,
    sampling_tol: f64,
) -> Result<()> {
    if samples == 0 {
        return Err(QRustError::Simulation(
            "equivalence verification requires at least one statistical sample".into(),
        ));
    }
    for (name, tolerance) in [("exact", exact_tol), ("sampling", sampling_tol)] {
        if !tolerance.is_finite() || !(0.0..1.0).contains(&tolerance) {
            return Err(QRustError::Simulation(format!(
                "{name} verification tolerance must be finite and in [0, 1)"
            )));
        }
    }
    Ok(())
}

fn unitary_prefix_and_readout(
    circuit: &Circuit,
    logical_output_locations: Option<&[usize]>,
) -> Result<(Circuit, BTreeMap<usize, usize>)> {
    crate::transpiler::validate_circuit_structure(circuit)
        .map_err(|e| QRustError::Simulation(format!("invalid circuit: {e}")))?;
    let mut prefix = circuit.clone();
    prefix.operations.clear();
    let mut measured = false;
    let mut readout = BTreeMap::new();
    for op in &circuit.operations {
        match op {
            Operation::Gate { name: crate::ir::GateType::Barrier, .. } | Operation::Barrier { .. } => {
                prefix.operations.push(op.clone());
            }
            Operation::Gate { .. } if !measured => prefix.operations.push(op.clone()),
            Operation::Gate { .. } => return Err(QRustError::Simulation(
                "equivalence verification supports only terminal measurements; found a gate after measurement".into(),
            )),
            Operation::Measure { qubit, cbit } => {
                measured = true;
                let logical_qubit = match logical_output_locations {
                    Some(layout) => layout.iter().position(|p| p == qubit).ok_or_else(|| {
                        QRustError::Simulation("terminal measurement addresses a non-logical ancilla".into())
                    })?,
                    None => *qubit,
                };
                if readout.insert(*cbit, logical_qubit).is_some() {
                    return Err(QRustError::Simulation(
                        "equivalence verification does not support repeated writes to a classical bit".into(),
                    ));
                }
            }
            Operation::Reset { .. } | Operation::Conditional { .. } => return Err(QRustError::Simulation(
                "unitary equivalence verification does not support reset or conditional operations".into(),
            )),
        }
    }
    Ok((prefix, readout))
}

fn validate_matching_readout(
    c1: &Circuit,
    c2: &Circuit,
    readout1: &BTreeMap<usize, usize>,
    readout2: &BTreeMap<usize, usize>,
) -> Result<()> {
    if readout1 != readout2
        || ((!readout1.is_empty() || !readout2.is_empty()) && c1.num_cbits != c2.num_cbits)
    {
        return Err(QRustError::Simulation(
            "terminal measurement maps or classical widths differ; general measurement-channel equivalence is not supported".into(),
        ));
    }
    Ok(())
}

/// Verifies a routed circuit while accounting for initial/final layouts and
/// any idle physical qubits added by the backend.
///
/// `initial_layout[i]` and `final_layout[i]` are the physical locations of
/// logical qubit `i` at circuit input and output. Exact verification is used
/// through 14 physical qubits; widths 15--20 use layout-aware Haar sampling.
pub fn verify_equivalence_with_layout(
    logical: &Circuit,
    routed: &Circuit,
    initial_layout: &[usize],
    final_layout: &[usize],
) -> Result<Verdict> {
    verify_equivalence_with_layout_and(
        logical,
        routed,
        initial_layout,
        final_layout,
        DEFAULT_SAMPLE_COUNT,
        0x00C0_FFEE_DEAD_BEEF_u64,
        DEFAULT_EXACT_TOLERANCE,
        DEFAULT_SAMPLING_TOLERANCE,
    )
}

/// Configurable variant of [`verify_equivalence_with_layout`].
#[allow(clippy::too_many_arguments)]
pub fn verify_equivalence_with_layout_and(
    logical: &Circuit,
    routed: &Circuit,
    initial_layout: &[usize],
    final_layout: &[usize],
    samples: usize,
    seed: u64,
    exact_tol: f64,
    sampling_tol: f64,
) -> Result<Verdict> {
    verify_equivalence_with_layout_configured(
        logical,
        routed,
        initial_layout,
        final_layout,
        &VerificationConfig {
            samples,
            seed,
            exact_tolerance: exact_tol,
            sampling_tolerance: sampling_tol,
            ..VerificationConfig::default()
        },
    )
}

/// Layout-aware verification with configurable resource limits. The dense
/// limit is applied to the physical width; sampling uses the lower of the
/// configured limit and [`MAX_LAYOUT_VERIFY_QUBITS`].
pub fn verify_equivalence_with_layout_configured(
    logical: &Circuit,
    routed: &Circuit,
    initial_layout: &[usize],
    final_layout: &[usize],
    config: &VerificationConfig,
) -> Result<Verdict> {
    config.validate()?;
    let n_logical = logical.num_qubits;
    let n_physical = routed.num_qubits;
    validate_layout(n_logical, n_physical, initial_layout, final_layout)?;
    let (logical, logical_readout) = unitary_prefix_and_readout(logical, None)?;
    let (routed, routed_readout) =
        unitary_prefix_and_readout(routed, Some(&final_layout[..n_logical]))?;
    validate_matching_readout(&logical, &routed, &logical_readout, &routed_readout)?;

    if n_physical <= config.max_exact_qubits {
        let u_logical = try_circuit_to_unitary(&logical)?;
        let u_routed = try_circuit_to_unitary(&routed)?;
        let extracted =
            try_extract_logical_unitary(&u_routed, n_logical, initial_layout, final_layout)?;
        let fidelity = unitary_fidelity(&u_logical, &extracted);
        if fidelity >= 1.0 - config.exact_tolerance {
            return Ok(Verdict::ExactlyEquivalent { fidelity });
        }
        return Ok(Verdict::NotEquivalent {
            fidelity,
            method: "layout-aware exact",
        });
    }

    let sampling_limit = config.max_sampling_qubits.min(MAX_LAYOUT_VERIFY_QUBITS);
    if n_physical <= sampling_limit {
        let min_fidelity = equivalence_by_sampling_with_layout(
            &logical,
            &routed,
            initial_layout,
            final_layout,
            config.samples,
            config.seed,
        )?;
        if min_fidelity >= 1.0 - config.sampling_tolerance {
            return Ok(Verdict::StatisticallyEquivalent {
                min_fidelity,
                samples: config.samples,
            });
        }
        return Ok(Verdict::NotEquivalent {
            fidelity: min_fidelity,
            method: "layout-aware sampling",
        });
    }

    Ok(Verdict::Unverifiable {
        reason: format!(
            "{n_physical} physical qubits exceeds configured exact ({}) and layout-aware sampling ({sampling_limit}) limits",
            config.max_exact_qubits,
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{GateType, Operation};

    fn bell_pair() -> Circuit {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c
    }

    #[test]
    fn test_verify_exact_equivalent() {
        let c = bell_pair();
        let v = verify_equivalence(&c, &c).unwrap();
        assert!(matches!(v, Verdict::ExactlyEquivalent { .. }));
        assert!(v.is_equivalent());
    }

    #[test]
    fn test_verify_exact_not_equivalent() {
        let c1 = bell_pair();
        let mut c2 = Circuit::new(2, 0);
        c2.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        let v = verify_equivalence(&c1, &c2).unwrap();
        assert!(matches!(v, Verdict::NotEquivalent { .. }));
        assert!(!v.is_equivalent());
    }

    #[test]
    fn test_verify_statistical_path_for_large_circuits() {
        let mut c = Circuit::new(18, 0);
        for i in 0..18 {
            c.add_op(Operation::Gate {
                name: GateType::H,
                qubits: vec![i],
                params: vec![],
            });
        }
        for i in 0..17 {
            c.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![i, i + 1],
                params: vec![],
            });
        }
        let v = verify_equivalence(&c, &c).unwrap();
        assert!(
            matches!(v, Verdict::StatisticallyEquivalent { .. }),
            "got {v:?}"
        );
        assert!(v.is_equivalent());
    }

    #[test]
    fn test_verify_qubit_mismatch_errors() {
        let c1 = Circuit::new(2, 0);
        let c2 = Circuit::new(3, 0);
        assert!(verify_equivalence(&c1, &c2).is_err());
    }

    #[test]
    fn test_verify_rejects_invalid_parameters() {
        let c = bell_pair();
        assert!(verify_equivalence_with(&c, &c, 0, 7, 1e-9, 1e-6).is_err());
        assert!(verify_equivalence_with(&c, &c, 8, 7, f64::NAN, 1e-6).is_err());
        assert!(verify_equivalence_with(&c, &c, 8, 7, 1e-9, 1.0).is_err());
    }

    #[test]
    fn test_verify_layout_aware_permutation() {
        let logical = bell_pair();
        let mut routed = Circuit::new(2, 0);
        routed.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });
        routed.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 0],
            params: vec![],
        });

        let verdict = verify_equivalence_with_layout(&logical, &routed, &[1, 0], &[1, 0]).unwrap();
        assert!(matches!(verdict, Verdict::ExactlyEquivalent { .. }));
    }

    #[test]
    fn test_verify_layout_aware_rejects_non_injective_layout() {
        let logical = bell_pair();
        assert!(verify_equivalence_with_layout(&logical, &logical, &[0, 0], &[0, 1]).is_err());
    }

    #[test]
    fn test_verdict_describe_formats() {
        let v = Verdict::ExactlyEquivalent { fidelity: 1.0 };
        assert!(v.describe().contains("equivalent"));
    }
}
