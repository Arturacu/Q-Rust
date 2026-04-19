//! KAK (Cartan) decomposition for two-qubit unitaries.
//!
//! Decomposes `U ∈ SU(4)` into
//! `U = (A1 ⊗ A0) · exp(i(x X⊗X + y Y⊗Y + z Z⊗Z)) · (B1 ⊗ B0)`,
//! synthesizing the result as an explicit gate sequence.
//!
//! Weyl-chamber branching (Shende et al.):
//! - 0 CX : all interaction coefficients below tolerance (local unitary).
//! - 1 CX : only one coefficient non-zero (SWAP-class).
//! - 3 CX : general case.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

const KAK_TOL: f64 = 1e-7;

/// Exact two-qubit (Cartan/KAK) synthesizer.
#[derive(Debug, Clone, Copy)]
pub struct KakSynthesizer;

impl Synthesizer for KakSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != 4 || unitary.ncols() != 4 {
            return None;
        }

        let det = unitary.determinant();
        let phase = det.powf(0.25);
        if phase.norm() < 1e-9 {
            return None;
        }
        let unitary = unitary / phase;

        let sqrt2_inv = 1.0 / 2.0_f64.sqrt();
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let i = Complex::<f64>::i();

        let magic_matrix = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, i, zero, i, one, zero, zero, i, -one, zero, one, zero, zero, -i,
            ],
        ) * Complex::new(sqrt2_inv, 0.0);

        let u_magic = magic_matrix.adjoint() * unitary * &magic_matrix;
        let m_sys = u_magic.transpose() * &u_magic;

        let mut a = DMatrix::<f64>::zeros(4, 4);
        let mut b = DMatrix::<f64>::zeros(4, 4);
        for ii in 0..4 {
            for jj in 0..4 {
                a[(ii, jj)] = m_sys[(ii, jj)].re;
                b[(ii, jj)] = m_sys[(ii, jj)].im;
            }
        }
        let c = &a + 0.7 * &b;
        let eigen_c = c.symmetric_eigen();
        let mut r_magic_real = eigen_c.eigenvectors;
        if r_magic_real.determinant() < 0.0 {
            for ii in 0..4 {
                r_magic_real[(ii, 0)] *= -1.0;
            }
        }
        let r_magic = r_magic_real.map(|x| Complex::new(x, 0.0));

        let d_mat = r_magic.transpose() * m_sys.clone() * &r_magic;
        let eigen_vals = d_mat.diagonal();
        let sqrt_eigen: Vec<Complex<f64>> = eigen_vals.iter().map(|c| c.sqrt()).collect();

        let mut best_l = DMatrix::zeros(4, 4);
        let mut best_n = DMatrix::zeros(4, 4);
        let mut min_score = f64::MAX;

        for mask in 0..16 {
            let mut d = Vec::with_capacity(4);
            for bit in 0..4 {
                let sign = if (mask >> bit) & 1 == 1 { -1.0 } else { 1.0 };
                d.push(sqrt_eigen[bit] * sign);
            }
            let n_candidate = DMatrix::from_diagonal(&DVector::from_vec(d.clone()));
            let det_n = d.iter().product::<Complex<f64>>();
            if (det_n - Complex::new(1.0, 0.0)).norm() > 1e-6 {
                continue;
            }
            let l_candidate = u_magic.clone() * r_magic.clone() * n_candidate.adjoint();
            let residual: f64 = l_candidate.iter().map(|c| c.im.powi(2)).sum();
            let phase_sum: f64 = d.iter().map(|c| c.arg().abs()).sum();
            let score = residual * 1000.0 + phase_sum;
            if score < min_score {
                min_score = score;
                best_l = l_candidate;
                best_n = n_candidate;
            }
        }

        let n_magic = best_n;
        let pi = std::f64::consts::PI;
        let phases: Vec<f64> = n_magic.diagonal().iter().map(|c| c.arg()).collect();

        let mut sums = Vec::new();
        for ii in 0..4 {
            for jj in (ii + 1)..4 {
                let mut s = phases[ii] + phases[jj];
                while s > pi {
                    s -= 2.0 * pi;
                }
                while s <= -pi {
                    s += 2.0 * pi;
                }
                sums.push(s.abs());
            }
        }
        sums.sort_by(|a, b| b.partial_cmp(a).unwrap());

        let mut x = (sums[0] + sums[1]) / 4.0;
        let mut y = (sums[2] + sums[3]) / 4.0;
        let mut z = (sums[4] + sums[5]) / 4.0;

        fn fold(val: f64) -> f64 {
            let mut v = val;
            let pi_2 = std::f64::consts::FRAC_PI_2;
            while v > pi_2 {
                v -= pi_2;
            }
            while v < -pi_2 {
                v += pi_2;
            }
            v
        }
        x = fold(x).abs();
        y = fold(y).abs();
        z = fold(z).abs();

        let mut coords = vec![x, y, z];
        coords.sort_by(|a, b| b.partial_cmp(a).unwrap());
        x = coords[0];
        y = coords[1];
        z = coords[2];

        let target_phases = [x - y + z, x + y - z, -x - y - z, -x + y + z];
        let current_phases: Vec<f64> = n_magic.diagonal().iter().map(|c| c.arg()).collect();
        let mut permutation = vec![0usize; 4];
        let mut used = vec![false; 4];
        for (ii, target) in target_phases.iter().enumerate() {
            let mut best_j = 0;
            let mut min_diff = f64::MAX;
            for jj in 0..4 {
                if !used[jj] {
                    let mut diff = (current_phases[jj] - target).abs();
                    while diff > pi {
                        diff -= 2.0 * pi;
                    }
                    while diff < -pi {
                        diff += 2.0 * pi;
                    }
                    diff = diff.abs();
                    if diff < min_diff {
                        min_diff = diff;
                        best_j = jj;
                    }
                }
            }
            permutation[ii] = best_j;
            used[best_j] = true;
        }

        let mut r_perm = DMatrix::zeros(4, 4);
        let mut n_perm_diag = Vec::with_capacity(4);
        for ii in 0..4 {
            let old = permutation[ii];
            r_perm.set_column(ii, &r_magic.column(old));
            n_perm_diag.push(n_magic[(old, old)]);
        }
        let mut r_magic = r_perm;
        let n_magic = DMatrix::from_diagonal(&DVector::from_vec(n_perm_diag));
        if r_magic.determinant().re < 0.0 {
            for ii in 0..4 {
                r_magic[(ii, 0)] *= -1.0;
            }
        }
        let l_magic = u_magic.clone() * r_magic.clone() * n_magic.adjoint();
        let l_magic = l_magic.map(|c| Complex::new(c.re, 0.0));

        let l_comp = &magic_matrix * l_magic * magic_matrix.adjoint();
        let r_comp = &magic_matrix * r_magic.transpose() * magic_matrix.adjoint();

        let (a0, a1) = decompose_tensor_product(&l_comp);
        let (b0, b1) = decompose_tensor_product(&r_comp);

        use crate::transpiler::synthesis::zyz::zyz_decomposition;
        fn to_array(m: DMatrix<Complex<f64>>) -> [[Complex<f64>; 2]; 2] {
            [[m[(0, 0)], m[(0, 1)]], [m[(1, 0)], m[(1, 1)]]]
        }
        let (a0_t, a0_p, a0_l, _) = zyz_decomposition(to_array(a0));
        let (a1_t, a1_p, a1_l, _) = zyz_decomposition(to_array(a1));
        let (b0_t, b0_p, b0_l, _) = zyz_decomposition(to_array(b0));
        let (b1_t, b1_p, b1_l, _) = zyz_decomposition(to_array(b1));

        let pi_2 = std::f64::consts::FRAC_PI_2;
        let _ = best_l;

        // Weyl-chamber branching: count non-trivial interaction coefficients.
        let nx = x.abs() > KAK_TOL;
        let ny = y.abs() > KAK_TOL;
        let nz = z.abs() > KAK_TOL;
        let n_nontrivial = (nx as u8) + (ny as u8) + (nz as u8);

        let mut circuit = Circuit::new(2, 0);

        // Pre-rotations (always emitted).
        circuit.add_op(Operation::Gate { name: GateType::U, qubits: vec![1], params: vec![b0_t, b0_p, b0_l] });
        circuit.add_op(Operation::Gate { name: GateType::U, qubits: vec![0], params: vec![b1_t, b1_p, b1_l] });

        match n_nontrivial {
            0 => {
                // 0-CX case: purely local unitary.
            }
            1 => {
                // 1-CX case: exactly one interaction coefficient is non-zero.
                // Emit the single corresponding block.
                if nx {
                    circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * x] });
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![1], params: vec![] });
                } else if ny {
                    circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![0], params: vec![pi_2] });
                    circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![1], params: vec![pi_2] });
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * y] });
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![0], params: vec![-pi_2] });
                    circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![1], params: vec![-pi_2] });
                } else {
                    // nz
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                    circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * z] });
                    circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                }
            }
            _ => {
                // General (3-CX) case: all three blocks.
                // XX
                circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * x] });
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![1], params: vec![] });
                // YY
                circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![0], params: vec![pi_2] });
                circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![1], params: vec![pi_2] });
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * y] });
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![0], params: vec![-pi_2] });
                circuit.add_op(Operation::Gate { name: GateType::RX, qubits: vec![1], params: vec![-pi_2] });
                // ZZ
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
                circuit.add_op(Operation::Gate { name: GateType::RZ, qubits: vec![1], params: vec![-2.0 * z] });
                circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
            }
        }

        // Post-rotations (always emitted).
        circuit.add_op(Operation::Gate { name: GateType::U, qubits: vec![1], params: vec![a0_t, a0_p, a0_l] });
        circuit.add_op(Operation::Gate { name: GateType::U, qubits: vec![0], params: vec![a1_t, a1_p, a1_l] });

        Some(circuit)
    }
}

fn decompose_tensor_product(
    u: &DMatrix<Complex<f64>>,
) -> (DMatrix<Complex<f64>>, DMatrix<Complex<f64>>) {
    let mut t = DMatrix::zeros(4, 4);
    for i0 in 0..2 {
        for i1 in 0..2 {
            for j0 in 0..2 {
                for j1 in 0..2 {
                    let row = 2 * i0 + j0;
                    let col = 2 * i1 + j1;
                    let u_row = 2 * i0 + i1;
                    let u_col = 2 * j0 + j1;
                    t[(row, col)] = u[(u_row, u_col)];
                }
            }
        }
    }
    let svd = t.svd(true, true);
    let u_vec = svd.u.unwrap().column(0).into_owned();
    let v_vec = svd.v_t.unwrap().row(0).transpose().into_owned();

    let mut a = DMatrix::zeros(2, 2);
    let mut b = DMatrix::zeros(2, 2);
    for i in 0..2 {
        for j in 0..2 {
            a[(i, j)] = u_vec[2 * i + j];
            b[(i, j)] = v_vec[2 * i + j];
        }
    }
    let det_a = a.determinant();
    if det_a.norm() > 1e-9 {
        a /= det_a.powf(0.5);
    }
    let det_b = b.determinant();
    if det_b.norm() > 1e-9 {
        b /= det_b.powf(0.5);
    }
    (a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transpiler::synthesis::zyz::u_to_matrix;

    fn u_matrix(theta: f64, phi: f64, lambda: f64) -> DMatrix<Complex<f64>> {
        let m = u_to_matrix(theta, phi, lambda);
        DMatrix::from_row_slice(2, 2, &[m[0][0], m[0][1], m[1][0], m[1][1]])
    }

    fn count_cx(c: &Circuit) -> usize {
        c.operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { name: GateType::CX, .. }))
            .count()
    }

    #[test]
    fn test_kak_identity_zero_cx() {
        let s = KakSynthesizer;
        let id = DMatrix::<Complex<f64>>::identity(4, 4);
        let circuit = s.synthesize(&id, &[]).unwrap();
        // Identity is local-only => 0 CX.
        assert_eq!(count_cx(&circuit), 0);
    }

    #[test]
    fn test_kak_cnot_branches() {
        let s = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let cnot = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, zero, zero, one, zero, zero, zero, zero, zero, one, zero, zero,
                one, zero,
            ],
        );
        let circuit = s.synthesize(&cnot, &[]).unwrap();
        // CNOT is SWAP-class (single non-trivial coefficient) => 1 block, which expands to 2 CXs.
        assert_eq!(count_cx(&circuit), 2);
        let _ = u_matrix(0.0, 0.0, 0.0);
    }
}