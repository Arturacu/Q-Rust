use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Synthesizer using Quantum Shannon Decomposition for N-qubit unitaries.
pub struct QsdSynthesizer;

impl Synthesizer for QsdSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        let n = (unitary.nrows() as f64).log2() as usize;
        let circuit = Circuit::new(n, 0);

        // TODO: Implement recursive QSD
        // Base case: n=1 (ZYZ) or n=2 (KAK)

        Some(circuit)
    }
}
