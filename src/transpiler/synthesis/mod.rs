//! Unitary synthesis algorithms.
//!
//! The [`Synthesizer`] trait is implemented by a family of algorithms:
//!
//! - [`zyz::ZyzSynthesizer`]   — exact single-qubit (Euler-angle) decomposition.
//! - [`kak::KakSynthesizer`]   — exact two-qubit (Cartan) decomposition.
//! - [`qsd::QsdSynthesizer`]   — Shannon decomposition for N≥3 (stub).
//! - [`qsearch::QSearchSynthesizer`] — A*-style approximate synthesis (stub).
//! - [`numerical::NumericalSynthesizer`] — gradient-based synthesis (stub).
//!
//! [`GlobalSynthesizer`] is a dispatcher that routes by matrix size.

pub mod kak;
pub mod numerical;
pub mod qsd;
pub mod qsearch;
pub mod zyz;

use crate::ir::{Circuit, GateType};
use nalgebra::DMatrix;
use num_complex::Complex;

/// Synthesizes a circuit that implements a given unitary.
pub trait Synthesizer {
    /// Returns a circuit equivalent (up to global phase) to `unitary`, using
    /// at least the gates in `basis`. Returns `None` when synthesis is not
    /// possible (wrong size, numerical failure, unsupported N).
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit>;
}

/// Dispatch by matrix size: 2×2 → ZYZ, 4×4 → KAK, otherwise `None`.
///
/// Does **not** panic on unsupported sizes; use [`qsd::QsdSynthesizer`] (still
/// a stub) or the numerical synthesizers for larger matrices.
#[derive(Debug, Clone, Copy)]
pub struct GlobalSynthesizer;

impl Synthesizer for GlobalSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != unitary.ncols() {
            return None;
        }
        match unitary.nrows() {
            2 => zyz::ZyzSynthesizer.synthesize(unitary, basis),
            4 => kak::KakSynthesizer.synthesize(unitary, basis),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use num_complex::Complex;

    #[test]
    fn test_global_synthesizer_routes_1q() {
        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
        let circ = GlobalSynthesizer.synthesize(&id2, &[]).unwrap();
        assert_eq!(circ.operations.len(), 1);
    }

    #[test]
    fn test_global_synthesizer_routes_2q() {
        let id4 = DMatrix::<Complex<f64>>::identity(4, 4);
        let circ = GlobalSynthesizer.synthesize(&id4, &[]).unwrap();
        assert_eq!(circ.operations.len(), 4);
    }

    #[test]
    fn test_global_synthesizer_returns_none_for_3q() {
        let id8 = DMatrix::<Complex<f64>>::identity(8, 8);
        assert!(GlobalSynthesizer.synthesize(&id8, &[]).is_none());
    }
}
