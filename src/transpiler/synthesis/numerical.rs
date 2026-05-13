//! Gradient-based numerical synthesis (stub).
//!
//! Placeholder for a GRAPE/BFGS-style solver that fits a parameterized
//! ansatz to a target unitary. Currently returns an empty circuit.

use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Numerical synthesizer (stub).
#[derive(Debug, Clone, Copy)]
pub struct NumericalSynthesizer;

impl Synthesizer for NumericalSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != unitary.ncols() || !unitary.nrows().is_power_of_two() {
            return None;
        }
        let n = (unitary.nrows() as f64).log2() as usize;
        Some(Circuit::new(n, 0))
    }
}
