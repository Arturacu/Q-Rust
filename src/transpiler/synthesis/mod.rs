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

/// A global synthesis router that dispatches arbitrary unitaries to the correct exact algorithm.
/// Enforces a hard physical ceiling restricting exact synthesis to N <= 2 qubits.
pub struct GlobalSynthesizer;

impl Synthesizer for GlobalSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit> {
        let rows = unitary.nrows();
        let cols = unitary.ncols();

        if rows != cols {
            return None;
        }

        match rows {
            2 => zyz::ZyzSynthesizer.synthesize(unitary, basis),
            4 => kak::KakSynthesizer.synthesize(unitary, basis),
            _ => unimplemented!(
                "Exact matrix synthesis for N >= 3 qubits (matrix size {}x{}) is intentionally unsupported. \
                 Exponential O(4^N) scaling renders exact unrolling physically impractical. \
                 Please rely on rule-based decomposition or implement Phase 3 approximate compilers.",
                rows, cols
            ),
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
        let synth = GlobalSynthesizer;
        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
        let circ = synth.synthesize(&id2, &[]).unwrap();
        // ZYZ produces exactly 1 gate (U gate)
        assert_eq!(circ.operations.len(), 1);
    }

    #[test]
    fn test_global_synthesizer_routes_2q() {
        let synth = GlobalSynthesizer;
        let id4 = DMatrix::<Complex<f64>>::identity(4, 4);
        let circ = synth.synthesize(&id4, &[]).unwrap();
        // KAK naturally generates 21 operations
        assert_eq!(circ.operations.len(), 21);
    }

    #[test]
    #[should_panic(expected = "Exact matrix synthesis for N >= 3 qubits")]
    fn test_global_synthesizer_rejects_3q() {
        let synth = GlobalSynthesizer;
        let id8 = DMatrix::<Complex<f64>>::identity(8, 8);
        let _ = synth.synthesize(&id8, &[]);
    }
}
