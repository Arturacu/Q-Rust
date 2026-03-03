pub mod kak;
pub mod numerical;
pub mod qsd;
pub mod qsearch;
pub mod zyz;

use crate::ir::{Circuit, GateType};
use nalgebra::DMatrix;
use num_complex::Complex;

/// A trait for unitary synthesis algorithms.
///
/// Implementors of this trait can synthesize a quantum circuit that approximates
/// a given unitary matrix using a specific set of basis gates.
pub trait Synthesizer {
    /// Synthesizes a unitary matrix into a Circuit using the allowed basis gates.
    ///
    /// # Arguments
    ///
    /// * `unitary` - The target unitary matrix to synthesize.
    /// * `basis` - The allowed set of basis gates to use in the output circuit.
    ///
    /// # Returns
    ///
    /// * `Option<Circuit>` - The synthesized circuit, or None if synthesis failed.
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit>;
}
