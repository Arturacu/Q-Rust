//! QSearch: A*-style approximate synthesis (stub).
//!
//! This is a placeholder for a future implementation that searches a
//! structure-and-parameter space using A* with a distance-to-target
//! heuristic. Returning a trivially-empty circuit for now.

use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// QSearch synthesizer (stub).
#[derive(Debug, Clone, Copy)]
pub struct QSearchSynthesizer;

impl Synthesizer for QSearchSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != unitary.ncols() || !unitary.nrows().is_power_of_two() {
            return None;
        }
        let n = (unitary.nrows() as f64).log2() as usize;
        Some(Circuit::new(n, 0))
    }
}
