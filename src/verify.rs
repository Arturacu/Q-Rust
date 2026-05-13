//! [E2E-NEW-FEATURE] Top-level circuit-equivalence verification.
//!
//! Thesis §1.2 Goal 3: "every pass and every benchmark circuit is
//! verified by reconstructing the induced unitary through an internal
//! simulator and comparing it to the input under operator-norm fidelity."
//!
//! For circuits beyond [`crate::simulator::MAX_QUBITS`] (=14), exact
//! unitary materialization is intractable. This module provides a
//! single [`verify_equivalence`] entry point that:
//!
//! 1. **n ≤ 14** — exact: full unitary fidelity via [`circuit_to_unitary`].
//! 2. **14 < n ≤ 22** — statistical: Haar-random sampling via
//!    [`equivalence_by_sampling`] (8 samples by default).
//! 3. **n > 22** — out of reach for both methods; returns
//!    [`Verdict::Unverifiable`].
//!
//! References: thesis §3.2 (Burgholzer & Wille 2021); Mele 2024 on
//! Haar-measure concentration.

use crate::error::{QRustError, Result};
use crate::ir::Circuit;
use crate::simulator::{
    circuit_to_unitary, equivalence_by_sampling, unitary_fidelity, MAX_QUBITS,
    MAX_STATE_VECTOR_QUBITS,
};

/// Threshold above which the exact unitary check is too memory-intensive.
pub const EXACT_VERIFY_LIMIT: usize = MAX_QUBITS;

/// Default number of Haar-random samples for statistical equivalence.
pub const DEFAULT_SAMPLE_COUNT: usize = 8;

/// Default fidelity threshold for the exact path.
pub const DEFAULT_EXACT_TOLERANCE: f64 = 1e-9;

/// Default fidelity threshold for the statistical path.
pub const DEFAULT_SAMPLING_TOLERANCE: f64 = 1e-6;

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

/// Verifies that two circuits implement equivalent unitaries (up to global
/// phase), automatically selecting exact or statistical methods based on
/// circuit size.
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
    if c1.num_qubits != c2.num_qubits {
        return Err(QRustError::Simulation(format!(
            "verify_equivalence: qubit count mismatch ({} vs {})",
            c1.num_qubits, c2.num_qubits
        )));
    }
    let n = c1.num_qubits;

    if n <= EXACT_VERIFY_LIMIT {
        let u1 = circuit_to_unitary(c1);
        let u2 = circuit_to_unitary(c2);
        let fid = unitary_fidelity(&u1, &u2);
        if (fid - 1.0).abs() < exact_tol {
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
    if n <= MAX_STATE_VECTOR_QUBITS - 2 {
        let min_fid = equivalence_by_sampling(c1, c2, samples, seed)?;
        if (min_fid - 1.0).abs() < sampling_tol {
            return Ok(Verdict::StatisticallyEquivalent {
                min_fidelity: min_fid,
                samples,
            });
        }
        return Ok(Verdict::NotEquivalent {
            fidelity: min_fid,
            method: "sampling",
        });
    }

    Ok(Verdict::Unverifiable {
        reason: format!(
            "{n} qubits exceeds both exact ({EXACT_VERIFY_LIMIT}) and sampling ({}) limits",
            MAX_STATE_VECTOR_QUBITS - 2
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
    fn test_verdict_describe_formats() {
        let v = Verdict::ExactlyEquivalent { fidelity: 1.0 };
        assert!(v.describe().contains("equivalent"));
    }
}
