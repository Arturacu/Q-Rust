//! Gradient-based numerical synthesis (stub).
//!
//! Placeholder for a GRAPE/BFGS-style solver that fits a parameterized
//! ansatz to a target unitary. Unsupported synthesis returns `None`; an empty
//! circuit would silently misrepresent every target other than identity.

use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Numerical synthesizer (stub).
#[derive(Debug, Clone, Copy)]
pub struct NumericalSynthesizer;

impl Synthesizer for NumericalSynthesizer {
    fn synthesize(&self, _unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unsupported_numerical_synthesis_never_claims_identity() {
        let zero = Complex::new(0.0, 0.0);
        let one = Complex::new(1.0, 0.0);
        let x = DMatrix::from_row_slice(2, 2, &[zero, one, one, zero]);
        assert!(NumericalSynthesizer.synthesize(&x, &[]).is_none());
        for size in [2, 4, 8] {
            assert!(NumericalSynthesizer
                .synthesize(&DMatrix::identity(size, size), &[])
                .is_none());
        }
    }
}
