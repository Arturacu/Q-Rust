use num_complex::Complex;
use std::f64::consts::PI;

/// Represents a 2x2 unitary matrix.
/// [[a, b], [c, d]]
pub type Unitary2x2 = [[Complex<f64>; 2]; 2];

/// Decomposes a single-qubit unitary matrix into Euler angles and global phase.
///
/// Returns (theta, phi, lambda, gamma) such that:
/// U = e^{i * gamma} * U(theta, phi, lambda)
///
/// Where U(theta, phi, lambda) is the standard OpenQASM unitary gate:
/// U(theta, phi, lambda) =
/// [ cos(theta/2)          -e^{i*lambda} * sin(theta/2) ]
/// [ e^{i*phi} * sin(theta/2)  e^{i*(phi+lambda)} * cos(theta/2) ]
pub fn zyz_decomposition(u: Unitary2x2) -> (f64, f64, f64, f64) {
    let v00 = u[0][0];
    let v01 = u[0][1];
    let v10 = u[1][0];
    let v11 = u[1][1];

    // theta = 2 * arccos(|v00|)
    // |v00| should be between 0 and 1.
    // Due to floating point errors, we clamp it.
    let mag_v00 = v00.norm().min(1.0);
    let theta = 2.0 * mag_v00.acos();

    let gamma;
    let phi;
    let lambda;

    // Check if theta is close to 0 (diagonal matrix)
    if theta.abs() < 1e-10 {
        // U is effectively diagonal:
        // [ e^{i*gamma}     0 ]
        // [ 0    e^{i*(gamma + phi + lambda)} ]

        // v00 = e^{i*gamma}
        gamma = v00.arg();

        // v11 = e^{i*(gamma + phi + lambda)}
        // phi + lambda = arg(v11) - gamma
        let total_phase = v11.arg() - gamma;

        // We can arbitrarily split total_phase between phi and lambda.
        // Standard convention: lambda = 0, phi = total_phase?
        // Or symmetric?
        // Let's set lambda = 0 for simplicity.
        lambda = 0.0;
        phi = total_phase;
    } else if (theta - PI).abs() < 1e-10 {
        // U is effectively anti-diagonal (X-like):
        // [ 0   -e^{i*(gamma + lambda)} ]
        // [ e^{i*(gamma + phi)}    0 ]

        // v01 = -e^{i*(gamma + lambda)}
        // v10 = e^{i*(gamma + phi)}

        // We can choose gamma such that one of them is simplified, or just use v10.
        // Let's define gamma based on v10 (arbitrary choice, but consistent).
        // Actually, we have 3 unknowns (gamma, phi, lambda) and 2 equations (args of v01, v10).
        // We have one degree of freedom (global phase vs relative phase).
        // But gamma is global phase of the *whole* decomposition.

        // Let's set gamma = arg(v10) - phi? No.

        // Let's look at the standard decomposition again.
        // v10 = e^{i*gamma} * e^{i*phi} * sin(pi/2) = e^{i*(gamma+phi)}
        // v01 = -e^{i*gamma} * e^{i*lambda} * sin(pi/2) = -e^{i*(gamma+lambda)} = e^{i*(gamma+lambda+pi)}

        // sum = arg(v10) + arg(v01) = 2*gamma + phi + lambda + pi
        // diff = arg(v10) - arg(v01) = phi - lambda - pi

        // We can set phi = 0 (or lambda = 0) if we want to fix the gauge?
        // No, usually we want to minimize parameters or match standard forms.
        // Let's set lambda = 0?
        // If lambda = 0:
        // arg(v01) = gamma + pi -> gamma = arg(v01) - pi
        // phi = arg(v10) - gamma

        lambda = 0.0;
        gamma = v01.arg() - PI;
        phi = v10.arg() - gamma;
    } else {
        // General case
        // v00 = e^{i*gamma} * cos(theta/2)
        // gamma = arg(v00) (since cos(theta/2) > 0 for 0 < theta < pi)
        gamma = v00.arg();

        // v10 = e^{i*(gamma+phi)} * sin(theta/2)
        // phi = arg(v10) - gamma
        phi = v10.arg() - gamma;

        // v01 = -e^{i*(gamma+lambda)} * sin(theta/2)
        // arg(v01) = gamma + lambda + pi
        // lambda = arg(v01) - gamma - pi
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

/// Converts Euler angles to a 2x2 unitary matrix.
///
/// U(theta, phi, lambda) =
/// [ cos(theta/2)          -e^{i*lambda} * sin(theta/2) ]
/// [ e^{i*phi} * sin(theta/2)  e^{i*(phi+lambda)} * cos(theta/2) ]
pub fn u_to_matrix(theta: f64, phi: f64, lambda: f64) -> Unitary2x2 {
    let cos = (theta / 2.0).cos();
    let sin = (theta / 2.0).sin();

    let e_phi = Complex::from_polar(1.0, phi);
    let e_lambda = Complex::from_polar(1.0, lambda);
    let e_phi_lambda = Complex::from_polar(1.0, phi + lambda);

    let u00 = Complex::new(cos, 0.0);
    let u01 = -e_lambda * sin;
    let u10 = e_phi * sin;
    let u11 = e_phi_lambda * cos;

    [[u00, u01], [u10, u11]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex;

    fn c(re: f64, im: f64) -> Complex<f64> {
        Complex::new(re, im)
    }

    #[test]
    fn test_decompose_identity() {
        let u = [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]];
        let (theta, phi, lambda, gamma) = zyz_decomposition(u);

        assert!(theta.abs() < 1e-10);
        // For identity, phi+lambda should be 0 (mod 2pi)
        // Our implementation sets lambda=0, so phi should be 0.
        assert!((phi + lambda).abs() < 1e-10);
        assert!(gamma.abs() < 1e-10);
    }

    #[test]
    fn test_decompose_x() {
        // X = [[0, 1], [1, 0]]
        let u = [[c(0.0, 0.0), c(1.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]];
        let (theta, phi, lambda, gamma) = zyz_decomposition(u);

        // X ~ U(pi, 0, pi) * e^{i*pi/2} ?
        // U(pi, 0, pi) = [[0, -e^{i*pi}*1], [1, 0]] = [[0, 1], [1, 0]]
        // So gamma should be 0.

        assert!((theta - PI).abs() < 1e-10);
        // Our implementation for anti-diagonal sets lambda=0
        // If lambda=0, gamma = arg(1) - pi = -pi
        // phi = arg(1) - (-pi) = pi
        // So we get U(pi, pi, 0) * e^{-i*pi}
        // U(pi, pi, 0) = [[0, -1], [e^{i*pi}, 0]] = [[0, -1], [-1, 0]]
        // * -1 = [[0, 1], [1, 0]] -> Correct!

        // Verify reconstruction
        let recon = reconstruct(theta, phi, lambda, gamma);
        assert_matrices_close(u, recon);
    }

    #[test]
    fn test_decompose_hadamard() {
        // H = 1/sqrt(2) * [[1, 1], [1, -1]]
        let s2 = 1.0 / 2.0f64.sqrt();
        let u = [[c(s2, 0.0), c(s2, 0.0)], [c(s2, 0.0), c(-s2, 0.0)]];
        let (theta, phi, lambda, gamma) = zyz_decomposition(u);

        // H = U(pi/2, 0, pi) * e^{i*pi/2}? No, H is real.
        // U(pi/2, 0, pi) = 1/s2 * [[1, -e^{i*pi}], [1, e^{i*pi}]] = 1/s2 * [[1, 1], [1, -1]]
        // So gamma should be 0.

        assert!((theta - PI / 2.0).abs() < 1e-10);

        let recon = reconstruct(theta, phi, lambda, gamma);
        assert_matrices_close(u, recon);
    }

    fn reconstruct(theta: f64, phi: f64, lambda: f64, gamma: f64) -> Unitary2x2 {
        let cos = (theta / 2.0).cos();
        let sin = (theta / 2.0).sin();

        let e_gamma = Complex::from_polar(1.0, gamma);
        let e_phi = Complex::from_polar(1.0, phi);
        let e_lambda = Complex::from_polar(1.0, lambda);
        let e_phi_lambda = Complex::from_polar(1.0, phi + lambda);

        let u00 = cos;
        let u01 = -e_lambda * sin;
        let u10 = e_phi * sin;
        let u11 = e_phi_lambda * cos;

        [
            [e_gamma * u00, e_gamma * u01],
            [e_gamma * u10, e_gamma * u11],
        ]
    }

    fn assert_matrices_close(a: Unitary2x2, b: Unitary2x2) {
        for i in 0..2 {
            for j in 0..2 {
                let diff = (a[i][j] - b[i][j]).norm();
                assert!(
                    diff < 1e-10,
                    "Mismatch at [{}][{}]: {:?} vs {:?}",
                    i,
                    j,
                    a[i][j],
                    b[i][j]
                );
            }
        }
    }
}
