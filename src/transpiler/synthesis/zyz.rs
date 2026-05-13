//! Exact single-qubit (Euler ZYZ) decomposition.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

/// A 2×2 unitary.
pub type Unitary2x2 = [[Complex<f64>; 2]; 2];

/// Decomposes a single-qubit unitary into Euler angles and a global phase.
///
/// Returns `(θ, φ, λ, γ)` such that `U = e^{iγ} U(θ, φ, λ)` where
/// `U(θ, φ, λ)` is the standard OpenQASM unitary gate.
pub fn zyz_decomposition(u: Unitary2x2) -> (f64, f64, f64, f64) {
    let v00 = u[0][0];
    let v01 = u[0][1];
    let v10 = u[1][0];
    let v11 = u[1][1];

    let mag_v00 = v00.norm().min(1.0);
    let theta = 2.0 * mag_v00.acos();
    let (gamma, phi, lambda);

    if theta.abs() < 1e-10 {
        gamma = v00.arg();
        lambda = 0.0;
        phi = v11.arg() - gamma;
    } else if (theta - PI).abs() < 1e-10 {
        lambda = 0.0;
        gamma = v01.arg() - PI;
        phi = v10.arg() - gamma;
    } else {
        gamma = v00.arg();
        phi = v10.arg() - gamma;
        lambda = v01.arg() - gamma - PI;
    }

    (
        theta,
        normalize_angle(phi),
        normalize_angle(lambda),
        normalize_angle(gamma),
    )
}

fn normalize_angle(angle: f64) -> f64 {
    let mut a = angle;
    while a <= -PI {
        a += 2.0 * PI;
    }
    while a > PI {
        a -= 2.0 * PI;
    }
    a
}

/// Converts Euler angles back to a 2×2 unitary (ignoring global phase).
pub fn u_to_matrix(theta: f64, phi: f64, lambda: f64) -> Unitary2x2 {
    let cos = (theta / 2.0).cos();
    let sin = (theta / 2.0).sin();
    let e_phi = Complex::from_polar(1.0, phi);
    let e_lambda = Complex::from_polar(1.0, lambda);
    let e_phi_lambda = Complex::from_polar(1.0, phi + lambda);
    [
        [Complex::new(cos, 0.0), -e_lambda * sin],
        [e_phi * sin, e_phi_lambda * cos],
    ]
}

/// Exact single-qubit (ZYZ) synthesizer.
#[derive(Debug, Clone, Copy)]
pub struct ZyzSynthesizer;

impl Synthesizer for ZyzSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != 2 || unitary.ncols() != 2 {
            return None;
        }
        let u_arr = [
            [unitary[(0, 0)], unitary[(0, 1)]],
            [unitary[(1, 0)], unitary[(1, 1)]],
        ];
        let (theta, phi, lambda, _) = zyz_decomposition(u_arr);
        let mut circuit = Circuit::new(1, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![theta, phi, lambda],
        });
        Some(circuit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex<f64> {
        Complex::new(re, im)
    }

    #[test]
    fn test_decompose_identity() {
        let u = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        let (theta, phi, lambda, gamma) = zyz_decomposition(u);
        assert!(theta.abs() < 1e-10);
        assert!((phi + lambda).abs() < 1e-10);
        assert!(gamma.abs() < 1e-10);
    }

    #[test]
    fn test_decompose_x() {
        let u = [[c(0.0, 0.0), c(1.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]];
        let (theta, _phi, _lambda, _gamma) = zyz_decomposition(u);
        assert!((theta - PI).abs() < 1e-10);
    }

    #[test]
    fn test_decompose_hadamard() {
        let s2 = 1.0 / 2.0_f64.sqrt();
        let u = [[c(s2, 0.0), c(s2, 0.0)], [c(s2, 0.0), c(-s2, 0.0)]];
        let (theta, _phi, _lambda, _gamma) = zyz_decomposition(u);
        assert!((theta - PI / 2.0).abs() < 1e-10);
    }
}