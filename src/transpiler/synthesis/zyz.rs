//! Exact single-qubit (Euler ZYZ) decomposition.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

/// A 2×2 unitary.
pub type Unitary2x2 = [[Complex<f64>; 2]; 2];

// These are admission/reconstruction guards for one synthesis call, not an
// allowance to discard small angles. A caller replacing many subcircuits
// must also account for their combined residual.
const INPUT_UNITARITY_TOLERANCE: f64 = 1e-12;
const OUTPUT_RESIDUAL_TOLERANCE: f64 = 1e-12;

pub(super) fn valid_unitary_1q(matrix: &DMatrix<Complex<f64>>) -> bool {
    if matrix.nrows() != 2
        || matrix.ncols() != 2
        || matrix
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return false;
    }
    let residual = (matrix.adjoint() * matrix - DMatrix::identity(2, 2)).norm();
    residual.is_finite() && residual <= INPUT_UNITARITY_TOLERANCE
}

/// Direct phase-aligned Frobenius residual. Unlike `1 - |trace| / 2`,
/// subtracting the matrices retains first-order small-angle discrepancies.
pub(super) fn residual_1q(target: &DMatrix<Complex<f64>>, actual: Unitary2x2) -> f64 {
    let overlap: Complex<f64> = (0..2)
        .flat_map(|row| (0..2).map(move |col| actual[row][col].conj() * target[(row, col)]))
        .sum();
    let magnitude = overlap.norm();
    if !magnitude.is_finite() || magnitude == 0.0 {
        return f64::INFINITY;
    }
    let phase = overlap / magnitude;
    (0..2)
        .flat_map(|row| {
            (0..2).map(move |col| (target[(row, col)] - actual[row][col] * phase).norm_sqr())
        })
        .sum::<f64>()
        .sqrt()
}

/// Decomposes a single-qubit unitary into Euler angles and a global phase.
///
/// Returns `(θ, φ, λ, γ)` such that `U = e^{iγ} U(θ, φ, λ)` where
/// `U(θ, φ, λ)` is the standard OpenQASM unitary gate.
pub fn zyz_decomposition(u: Unitary2x2) -> (f64, f64, f64, f64) {
    let v00 = u[0][0];
    let v01 = u[0][1];
    let v10 = u[1][0];
    let v11 = u[1][1];

    // atan2 preserves tiny off-diagonal rotations when |v00| rounds to 1.
    // acos(|v00|) would lose these rotations and amplify rounding near zero.
    let mag_v00 = v00.norm();
    let mag_v10 = v10.norm();
    let theta = 2.0 * mag_v10.atan2(mag_v00);
    let (gamma, phi, lambda);

    if mag_v10 < 4.0 * f64::EPSILON {
        gamma = v00.arg();
        lambda = 0.0;
        phi = v11.arg() - gamma;
    } else if mag_v00 < 4.0 * f64::EPSILON {
        lambda = 0.0;
        gamma = v01.arg() - PI;
        phi = v10.arg() - gamma;
    } else {
        gamma = v00.arg();
        phi = v10.arg() - gamma;
        // Near a diagonal unitary, independent phases of two tiny noisy
        // off-diagonal entries do not reliably determine phi+lambda. Fix
        // the large diagonal entries directly; near the other endpoint,
        // use the large off-diagonal entries instead. For exact unitaries
        // these formulas are identical modulo 2π.
        lambda = if mag_v00 >= mag_v10 {
            v11.arg() - v10.arg()
        } else {
            (-v01).arg() - gamma
        };
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

/// Single-qubit (ZYZ) synthesizer with finite-input, unitarity and direct
/// reconstruction guards. Numerical failures return `None`.
#[derive(Debug, Clone, Copy)]
pub struct ZyzSynthesizer;

impl Synthesizer for ZyzSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if !valid_unitary_1q(unitary) {
            return None;
        }
        let u_arr = [
            [unitary[(0, 0)], unitary[(0, 1)]],
            [unitary[(1, 0)], unitary[(1, 1)]],
        ];
        let (theta, phi, lambda, _) = zyz_decomposition(u_arr);
        let residual = residual_1q(unitary, u_to_matrix(theta, phi, lambda));
        if !residual.is_finite() || residual > OUTPUT_RESIDUAL_TOLERANCE {
            return None;
        }
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
    fn synthesis_rejects_nonunitary_and_nonfinite_inputs() {
        for target in [
            DMatrix::zeros(2, 2),
            DMatrix::from_element(2, 2, c(1.0, 0.0)),
            DMatrix::identity(2, 2) * c(1.0 + 1e-8, 0.0),
            DMatrix::identity(3, 3),
            DMatrix::zeros(2, 3),
        ] {
            assert!(ZyzSynthesizer.synthesize(&target, &[]).is_none());
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for imaginary in [false, true] {
                let mut target = DMatrix::identity(2, 2);
                target[(0, 1)] = if imaginary { c(0.0, bad) } else { c(bad, 0.0) };
                assert!(ZyzSynthesizer.synthesize(&target, &[]).is_none());
            }
        }
    }

    #[test]
    fn guarded_synthesis_retains_small_rotations_and_global_phase() {
        for angle in [9e-10, -4e-11, 1e-13, PI - 1e-11] {
            // Analytic RX target; no named-gate lowering or synthesis oracle.
            let (sine, cosine) = (angle / 2.0).sin_cos();
            let target = DMatrix::from_row_slice(
                2,
                2,
                &[c(cosine, 0.0), c(0.0, -sine), c(0.0, -sine), c(cosine, 0.0)],
            ) * Complex::from_polar(1.0, -0.47);
            let output = ZyzSynthesizer.synthesize(&target, &[]).unwrap();
            let actual = crate::simulator::circuit_to_unitary(&output);
            let overlap = (actual.adjoint() * &target).trace();
            let error = (&target - actual * (overlap / overlap.norm())).norm();
            assert!(error < 2e-15, "angle={angle}, residual={error}");
        }
    }

    #[test]
    fn reconstructs_tiny_rotations_and_near_pi_with_shared_phase() {
        for theta in [1e-13, -1e-11, 1e-9, 1e-7, PI - 1e-11, PI + 1e-9] {
            let phase = Complex::from_polar(1.0, 0.71);
            let mut target = u_to_matrix(theta, -0.83, 1.27);
            for row in &mut target {
                for entry in row {
                    *entry *= phase;
                }
            }
            let (t, p, l, g) = zyz_decomposition(target);
            let actual = u_to_matrix(t, p, l);
            let global = Complex::from_polar(1.0, g);
            let residual = (0..2)
                .flat_map(|i| {
                    (0..2).map(move |j| (target[i][j] - actual[i][j] * global).norm_sqr())
                })
                .sum::<f64>()
                .sqrt();
            assert!(residual < 2e-15, "theta={theta}, residual={residual}");
        }
    }

    #[test]
    fn tiny_roundoff_entries_do_not_corrupt_large_diagonal_phases() {
        // A nearly diagonal local tensor factor found by the W-transfer
        // regression. Its off-diagonal relative error is large despite its
        // absolute error being at floating-point roundoff scale.
        let target = [
            [c(6.12e-17, 1.0 + 4e-16), c(-1.138e-15, -8.33e-17)],
            [c(1.138e-15, 8.33e-17), c(4.98e-17, -1.0 + 3e-16)],
        ];
        let (theta, phi, lambda, gamma) = zyz_decomposition(target);
        let actual = u_to_matrix(theta, phi, lambda);
        let global = Complex::from_polar(1.0, gamma);
        let residual = (0..2)
            .flat_map(|i| (0..2).map(move |j| (target[i][j] - actual[i][j] * global).norm_sqr()))
            .sum::<f64>()
            .sqrt();
        assert!(residual < 2e-15, "near-diagonal residual {residual}");
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
