use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Synthesizer using Numerical Optimization (GRAPE/BFGS).
pub struct NumericalSynthesizer;

impl Synthesizer for NumericalSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        let n = (unitary.nrows() as f64).log2() as usize;
        let circuit = Circuit::new(n, 0);

        // TODO: Implement Numerical Optimization
        // 1. Define Ansatz (parameterized circuit)
        // 2. Define Cost Function (1 - Fidelity)
        // 3. Optimize parameters

        Some(circuit)
    }
}
