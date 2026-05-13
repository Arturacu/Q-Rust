//! QSearch: A*-style approximate synthesis (stub).
//!
//! This is a placeholder for a future implementation that searches a
//! structure-and-parameter space using A* with a distance-to-target
//! heuristic.
//!
//! ## Loop 5 review §Finding 2 — silent-stub correctness bug
//!
//! Previously this stub returned `Some(Circuit::new(n, 0))` — an empty
//! identity circuit — for any well-shaped 2^n × 2^n input. That silently
//! corrupted any caller that assumed `Some(_)` implied a valid synthesis,
//! because the returned circuit is *not* equivalent to the input unitary
//! (except in the trivial identity case).
//!
//! Per the [`Synthesizer`] trait contract — "Returns `None` when synthesis
//! is not possible" — this stub now returns `None`. Callers that need
//! actual 1-qubit synthesis should use
//! [`crate::transpiler::synthesis::zyz::ZyzSynthesizer`] (exact, ZYZ Euler
//! decomposition) or
//! [`crate::transpiler::synthesis::qsd::NelderMead1qSynthesizer`]
//! (gradient-free numerical search).
//!
//! Reference: Davis et al. 2020, *Towards Optimal Topology Aware Quantum
//! Circuit Synthesis*, IEEE QCE 2020 — describes the genuine A*-over-
//! gate-structure algorithm we intend to implement here.

use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// QSearch A*-style synthesizer — **stub**.
///
/// Returns `None` until the actual A* search over gate-structure space is
/// implemented. For 1-qubit synthesis, prefer
/// [`crate::transpiler::synthesis::zyz::ZyzSynthesizer`]; for 2-qubit
/// synthesis, prefer [`crate::transpiler::synthesis::kak::KakSynthesizer`].
#[derive(Debug, Clone, Copy, Default)]
pub struct QSearchSynthesizer;

impl Synthesizer for QSearchSynthesizer {
    fn synthesize(&self, _unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        // Loop 5 §Finding 2: return None rather than an empty (identity)
        // Circuit so callers fail loudly rather than silently producing
        // a no-op.
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Loop 5 §Finding 2 — regression: the stub must NOT return a
    /// silently-wrong identity circuit for arbitrary 2^n × 2^n inputs.
    #[test]
    fn test_qsearch_stub_returns_none_for_well_shaped_input() {
        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
        assert!(QSearchSynthesizer.synthesize(&id2, &[]).is_none());
        let id4 = DMatrix::<Complex<f64>>::identity(4, 4);
        assert!(QSearchSynthesizer.synthesize(&id4, &[]).is_none());
    }

    #[test]
    fn test_qsearch_stub_returns_none_for_ill_shaped_input() {
        let nonsquare = DMatrix::<Complex<f64>>::zeros(2, 4);
        assert!(QSearchSynthesizer.synthesize(&nonsquare, &[]).is_none());
        let non_pow2 = DMatrix::<Complex<f64>>::identity(3, 3);
        assert!(QSearchSynthesizer.synthesize(&non_pow2, &[]).is_none());
    }
}