use crate::ir::{Circuit, GateType};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Synthesizer using QSearch (A* Search) algorithm.
pub struct QSearchSynthesizer;

impl Synthesizer for QSearchSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        let n = (unitary.nrows() as f64).log2() as usize;
        let circuit = Circuit::new(n, 0);

        // TODO: Implement A* search
        // PriorityQueue<CircuitState>
        // Cost = distance(current_unitary, target_unitary) + depth_penalty

        Some(circuit)
    }
}
