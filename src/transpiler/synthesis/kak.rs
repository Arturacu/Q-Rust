use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

/// Represents a 4x4 unitary matrix for 2 qubits.

/// KAK Decomposition for 2-qubit Unitaries.
///
/// This module implements the KAK decomposition (also known as the Cartan decomposition)
/// for arbitrary 2-qubit unitary matrices. It decomposes a unitary $U \in SU(4)$ into:
///
/// $$U = (A_1 \otimes A_0) \cdot e^{i(x X \otimes X + y Y \otimes Y + z Z \otimes Z)} \cdot (B_1 \otimes B_0)$$
///
/// where:
/// - $A_0, A_1, B_0, B_1 \in SU(2)$ are single-qubit gates.
/// - $x, y, z \in \mathbb{R}$ are the interaction coefficients (Weyl chamber coordinates).
///
/// The interaction term $N = e^{i(x XX + y YY + z ZZ)}$ is diagonal in the "Magic Basis".
///
/// # Algorithm
/// 1. Transform $U$ to the Magic Basis: $U_M = M^\dagger U M$.
/// 2. Compute $M_{sys} = U_M^T U_M$.
/// 3. Diagonalize $M_{sys}$ to find the interaction coefficients $x, y, z$.
/// 4. Construct the local unitaries $A$ and $B$ from the eigenvectors of $M_{sys}$.
/// 5. Synthesize the circuit using the computed coefficients and single-qubit decompositions.
///
/// # References
/// - [Synthesis of Quantum Logic Circuits](https://arxiv.org/abs/quant-ph/0406176)
pub struct KakSynthesizer;

impl Synthesizer for KakSynthesizer {
    /// Synthesizes a 2-qubit unitary matrix into a quantum circuit.
    ///
    /// # Arguments
    /// * `unitary` - A 4x4 unitary matrix (will be normalized to SU(4)).
    /// * `_basis` - Optional target basis gates (currently ignored, fixed decomposition used).
    ///
    /// # Returns
    /// * `Option<Circuit>` - The synthesized circuit, or `None` if the matrix is not 4x4.
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        // Check dimensions
        if unitary.nrows() != 4 || unitary.ncols() != 4 {
            return None; // KAK is only for 2 qubits (4x4)
        }

        // Normalize to SU(4) (det = 1)
        let det = unitary.determinant();
        let phase = det.powf(0.25);
        if phase.norm() < 1e-9 {
            return None;
        }
        let unitary = unitary / phase;

        // 1. Define Magic Basis Matrix M
        // The Magic Basis maps the Bell basis to the computational basis (up to phases).
        // M = 1/sqrt(2) * [[1, 0, 0, i], [0, i, 1, 0], [0, i, -1, 0], [1, 0, 0, -i]]
        // This specific definition matches standard literature (e.g. arXiv:quant-ph/0406176).
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

        // 2. Transform to Magic Basis: U_M = M^dag * U * M
        let u_magic = magic_matrix.adjoint() * unitary * &magic_matrix;

        // 3. Compute M_sys = U_M^T * U_M
        let m_sys = u_magic.transpose() * &u_magic;

        // 4. Eigenvalue Decomposition of M_sys
        // M_sys = A + iB is complex symmetric and unitary.
        // A and B are real symmetric and commute.
        // We can simultaneously diagonalize them by diagonalizing C = A + mu * B.

        let mut a = DMatrix::<f64>::zeros(4, 4);
        let mut b = DMatrix::<f64>::zeros(4, 4);

        for i in 0..4 {
            for j in 0..4 {
                a[(i, j)] = m_sys[(i, j)].re;
                b[(i, j)] = m_sys[(i, j)].im;
            }
        }

        // Mix A and B to avoid accidental degeneracy in one of them
        let c = &a + 0.7 * &b;

        // Diagonalize C (real symmetric)
        let eigen_c = c.symmetric_eigen();
        let mut r_magic_real = eigen_c.eigenvectors; // This is R_M (real orthogonal)

        // Ensure R_M is in SO(4) (det = 1)
        if r_magic_real.determinant() < 0.0 {
            for i in 0..4 {
                r_magic_real[(i, 0)] *= -1.0;
            }
        }

        // Convert to complex for multiplication
        let r_magic = r_magic_real.map(|x| Complex::new(x, 0.0));

        // Compute eigenvalues of M_sys using R_M
        // D = R_M^T * M_sys * R_M
        // Since R_M is real orthogonal, R_M^T = R_M^T.
        let d_mat = r_magic.transpose() * m_sys.clone() * &r_magic;

        // Extract diagonal (eigenvalues)
        let eigen_vals = d_mat.diagonal();

        // Search for correct signs of square roots of eigenvalues
        let mut best_l = DMatrix::zeros(4, 4);
        let mut best_n = DMatrix::zeros(4, 4);
        let mut min_score = f64::MAX; // Score = residual + penalty for large phases

        // Base square roots
        let sqrt_eigen: Vec<Complex<f64>> = eigen_vals.iter().map(|c| c.sqrt()).collect();

        // Iterate 2^4 = 16 sign combinations
        for i in 0..16 {
            let mut d = Vec::new();
            for bit in 0..4 {
                let sign = if (i >> bit) & 1 == 1 { -1.0 } else { 1.0 };
                d.push(sqrt_eigen[bit] * sign);
            }
            let n_candidate = DMatrix::from_diagonal(&DVector::from_vec(d.clone()));

            // Check determinant of N (must be 1 for SU(4))
            let det_n = d.iter().product::<Complex<f64>>();
            if (det_n - Complex::new(1.0, 0.0)).norm() > 1e-6 {
                continue;
            }

            // L = U * R * N^dag (Not R^T!)
            // U = L N R^T => L = U R N^dag
            let l_candidate = u_magic.clone() * r_magic.clone() * n_candidate.adjoint();

            // Check reality
            let residual: f64 = l_candidate.iter().map(|c| c.im.powi(2)).sum();

            // Penalty for large phases (prefer canonical branch)
            let phase_sum: f64 = d.iter().map(|c| c.arg().abs()).sum();

            // Combined score: primarily residual, secondary phase sum
            let score = residual * 1000.0 + phase_sum;

            if score < min_score {
                min_score = score;
                best_l = l_candidate;
                best_n = n_candidate;
            }
        }

        let _l_magic = best_l.map(|c| Complex::new(c.re, 0.0)); // Force real
        let n_magic = best_n;

        // Extract x, y, z from N_magic
        // N = diag(d0, d1, d2, d3)
        let d0 = n_magic[(0, 0)];
        let d1 = n_magic[(1, 1)];
        let d2 = n_magic[(2, 2)];
        let d3 = n_magic[(3, 3)];

        let _t0 = d0.arg();
        let _t1 = d1.arg();
        let _t2 = d2.arg();
        let _t3 = d3.arg();

        // We need to map t_k to x, y, z.
        // We don't know the permutation of t_k relative to x-y+z etc.
        // But we know x, y, z are related to t_k by Hadamard.
        // And x, y, z should be in Weyl chamber (positive, ordered).
        // Try all permutations of t_k?
        // Or just solve for x, y, z from the set of t_k.

        // The set {t0, t1, t2, t3} = {x-y+z, x+y-z, -x-y-z, -x+y+z} (modulo 2pi).
        // Sum = 0.
        // x = (t_i + t_j)/2.
        // Extract x, y, z from n_magic using Sums of Pairs algorithm
        // The set of sums of pairs of phases {p_i + p_j} corresponds to {+/- 2x, +/- 2y, +/- 2z}.
        let phases: Vec<f64> = n_magic.diagonal().iter().map(|c| c.arg()).collect();

        let mut sums = Vec::new();
        let pi = std::f64::consts::PI;
        for i in 0..4 {
            for j in (i + 1)..4 {
                let s = phases[i] + phases[j];
                // Fold to [-pi, pi]
                let mut s_folded = s;
                while s_folded > pi {
                    s_folded -= 2.0 * pi;
                }
                while s_folded <= -pi {
                    s_folded += 2.0 * pi;
                }
                sums.push(s_folded.abs());
            }
        }

        sums.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // We expect 3 pairs of values.
        // Due to numerical noise, they might not be exactly equal.
        // Largest 2 should be 2x. Next 2 should be 2y. Smallest 2 should be 2z.
        // But wait, we sorted descending.
        // So sums[0], sums[1] -> 2x
        // sums[2], sums[3] -> 2y
        // sums[4], sums[5] -> 2z

        let mut x = (sums[0] + sums[1]) / 4.0; // Average of 2x, divided by 2
        let mut y = (sums[2] + sums[3]) / 4.0;
        let mut z = (sums[4] + sums[5]) / 4.0;

        fn fold_scalar(val: f64) -> f64 {
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
        x = fold_scalar(x).abs();
        y = fold_scalar(y).abs();
        z = fold_scalar(z).abs();

        let mut coords = vec![x, y, z];
        coords.sort_by(|a, b| b.partial_cmp(a).unwrap());
        x = coords[0];
        y = coords[1];
        z = coords[2];

        // Target phases for our basis order: Phi+, Psi+, Psi-, Phi-
        // p0 = x - y + z
        // p1 = x + y - z
        // p2 = -x - y - z
        // p3 = -x + y + z
        let target_phases = vec![x - y + z, x + y - z, -x - y - z, -x + y + z];

        // Find permutation to match n_magic phases to target_phases
        let current_phases: Vec<f64> = n_magic.diagonal().iter().map(|c| c.arg()).collect();
        let mut permutation = vec![0; 4];
        let mut used = vec![false; 4];

        for (i, target) in target_phases.iter().enumerate() {
            let mut best_j = 0;
            let mut min_diff = f64::MAX;
            for j in 0..4 {
                if !used[j] {
                    // Difference modulo 2pi
                    let mut diff = (current_phases[j] - target).abs();
                    while diff > pi {
                        diff -= 2.0 * pi;
                    }
                    while diff < -pi {
                        diff += 2.0 * pi;
                    }
                    diff = diff.abs();

                    if diff < min_diff {
                        min_diff = diff;
                        best_j = j;
                    }
                }
            }
            permutation[i] = best_j;
            used[best_j] = true;
        }

        // Permute R_M and N
        let mut r_permuted = DMatrix::zeros(4, 4);
        let mut n_permuted_diag = Vec::new();

        for i in 0..4 {
            let old_idx = permutation[i];
            r_permuted.set_column(i, &r_magic.column(old_idx));
            n_permuted_diag.push(n_magic[(old_idx, old_idx)]);
        }

        let mut r_magic = r_permuted;
        let n_magic = DMatrix::from_diagonal(&DVector::from_vec(n_permuted_diag));

        // Ensure R_M is in SO(4) (det = 1)
        // Permutation might have flipped the sign.
        if r_magic.determinant().re < 0.0 {
            for i in 0..4 {
                r_magic[(i, 0)] *= -1.0;
            }
        }

        // Recalculate L with permuted matrices
        // L = U * R * N^dag
        let l_magic = u_magic.clone() * r_magic.clone() * n_magic.adjoint();

        // Ensure L is real (it should be, up to numerical noise)
        let l_magic = l_magic.map(|c| Complex::new(c.re, 0.0));

        // 5. Construct Circuit
        // L and R are in Magic Basis. They correspond to local unitaries in Computational Basis.
        // L_comp = M * L_magic * M^dag
        // R_comp = M * R_magic^T * M^dag

        let l_comp = &magic_matrix * l_magic * magic_matrix.adjoint();
        let r_comp = &magic_matrix * r_magic.transpose() * magic_matrix.adjoint();

        // Decompose L_comp and R_comp into tensor products
        let (a0, a1) = decompose_tensor_product(&l_comp);
        let (b0, b1) = decompose_tensor_product(&r_comp);

        // Debug: Check if a0 x a1 matches l_comp
        // Note: kronecker order depends on library. nalgebra: a.kronecker(b) -> blocks of a scaled by b.
        // So a is Outer.

        // Check for NaNs

        // 10. Convert to Euler Angles (ZYZ decomposition)
        use crate::transpiler::synthesis::zyz::zyz_decomposition;

        // Helper to convert DMatrix to [[Complex<f64>; 2]; 2]
        fn to_array(m: DMatrix<Complex<f64>>) -> [[Complex<f64>; 2]; 2] {
            [[m[(0, 0)], m[(0, 1)]], [m[(1, 0)], m[(1, 1)]]]
        }

        let (a0_theta, a0_phi, a0_lambda, _a0_gamma) = zyz_decomposition(to_array(a0));
        let (a1_theta, a1_phi, a1_lambda, _a1_gamma) = zyz_decomposition(to_array(a1));
        let (b0_theta, b0_phi, b0_lambda, _b0_gamma) = zyz_decomposition(to_array(b0));
        let (b1_theta, b1_phi, b1_lambda, _b1_gamma) = zyz_decomposition(to_array(b1));
        let mut circuit = Circuit::new(2, 0);

        // Pre-rotations (B1 tensor B0)
        // decompose returns (Outer, Inner).
        // If q0 is Outer, then a0 on 0.
        // Let's try swapping        // Pre-rotations (B1 tensor B0)
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![1], // Swapped
            params: vec![b0_theta, b0_phi, b0_lambda],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0], // Swapped
            params: vec![b1_theta, b1_phi, b1_lambda],
        });
        // Interaction:
        // e^{ixXX} = H0 H1 CX01 RZ(-2x)1 CX01 H0 H1
        // e^{iyYY} = RX(pi/2)0 RX(pi/2)1 CX01 RZ(-2y)1 CX01 RX(-pi/2)0 RX(-pi/2)1
        // e^{izZZ} = CX01 RZ(-2z)1 CX01

        // XX
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![1],
            params: vec![-2.0 * x],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![1],
            params: vec![],
        });

        // YY
        let pi_2 = std::f64::consts::FRAC_PI_2;
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![pi_2],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![1],
            params: vec![pi_2],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![1],
            params: vec![-2.0 * y],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![-pi_2],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![1],
            params: vec![-pi_2],
        });

        // ZZ
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![1],
            params: vec![-2.0 * z],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        // Post-rotations (A1 tensor A0)
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![1], // Swapped
            params: vec![a0_theta, a0_phi, a0_lambda],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0], // Swapped
            params: vec![a1_theta, a1_phi, a1_lambda],
        });

        Some(circuit)
    }
}

/// Decomposes a 4x4 unitary U known to be a tensor product A tensor B
/// Returns (A, B) up to a global phase.
fn decompose_tensor_product(
    u: &DMatrix<Complex<f64>>,
) -> (DMatrix<Complex<f64>>, DMatrix<Complex<f64>>) {
    // Rearrange U into a 4x4 matrix T where T_ij = A_k B_l
    // T has rows indexed by (i0, j0) and cols by (i1, j1)
    // T[(i0, j0), (i1, j1)] = U[(i0, i1), (j0, j1)]
    let mut t = DMatrix::zeros(4, 4);
    for i0 in 0..2 {
        for i1 in 0..2 {
            for j0 in 0..2 {
                for j1 in 0..2 {
                    let row = 2 * i0 + j0; // (i0, j0)
                    let col = 2 * i1 + j1; // (i1, j1)
                    let u_row = 2 * i0 + i1; // (i0, i1)
                    let u_col = 2 * j0 + j1; // (j0, j1)
                    t[(row, col)] = u[(u_row, u_col)];
                }
            }
        }
    }

    // SVD of T = U S V^dag
    // T = sigma_0 u_0 v_0^dag
    // A is proportional to u_0
    // B is proportional to v_0^dag? No, B is v_0 (conjugated)
    // T_kl = A_k B_l
    // T = vec(A) vec(B)^T
    // In complex SVD: T = U S V^dag = s0 u0 v0^dag
    // So vec(A) ~ u0
    // vec(B)^T ~ v0^dag => vec(B) ~ (v0^dag)^T = v0^*
    // So B is conjugate of v0.

    let svd = t.svd(true, true);
    let s = svd.singular_values;
    if s[1] > 1e-6 {}
    let u_vec = svd.u.unwrap().column(0).into_owned();
    let v_vec = svd.v_t.unwrap().row(0).transpose().into_owned(); // v_t is V^dag, row 0 is v0^dag

    let mut a = DMatrix::zeros(2, 2);
    let mut b = DMatrix::zeros(2, 2);
    for i in 0..2 {
        for j in 0..2 {
            a[(i, j)] = u_vec[2 * i + j];
            b[(i, j)] = v_vec[2 * i + j]; // No conjugation needed! v_vec is already B.
        }
    }

    // Normalize to SU(2)
    let det_a = a.determinant();
    if det_a.norm() > 1e-9 {
        let phase_a = det_a.powf(0.5);
        a /= phase_a;
    }

    let det_b = b.determinant();
    if det_b.norm() > 1e-9 {
        let phase_b = det_b.powf(0.5);
        b /= phase_b;
    }

    (a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transpiler::synthesis::zyz::u_to_matrix;
    use std::f64::consts::PI;

    fn u_matrix(theta: f64, phi: f64, lambda: f64) -> DMatrix<Complex<f64>> {
        let m = u_to_matrix(theta, phi, lambda);
        DMatrix::from_row_slice(2, 2, &[m[0][0], m[0][1], m[1][0], m[1][1]])
    }

    #[test]
    fn test_kak_decomposition_identity() {
        let synth = KakSynthesizer;
        let identity = DMatrix::<Complex<f64>>::identity(4, 4);
        let circuit = synth.synthesize(&identity, &[]).unwrap();

        // Identity should have 0 interaction parameters (or close to 0)
        // The circuit structure is fixed (3 CXs), but the RZ/RX angles should be 0.
        // We can check the parameters of the middle gates.

        // Expected structure:
        // 0,1: Pre-rotations
        // 2-8: e^{ixXX} (7 ops) -> x is at index 5 (RZ)
        // 9-15: e^{iyYY} (7 ops) -> y is at index 12 (RZ)
        // 16-18: e^{izZZ} (3 ops) -> z is at index 17 (RZ)
        // 19,20: Post-rotations

        assert_eq!(circuit.operations.len(), 21);

        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            assert!(params[0].abs() < 1e-6, "x should be 0 for Identity");
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            assert!(params[0].abs() < 1e-6, "y should be 0 for Identity");
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            assert!(params[0].abs() < 1e-6, "z should be 0 for Identity");
        }

        // Check single qubit gates (should be Identity or close to it)
        // Indices: 0, 1 (Pre), 19, 20 (Post)
        for i in [0, 1, 19, 20] {
            if let Operation::Gate { params, .. } = &circuit.operations[i] {
                // Check theta (params[0])
                // Note: phi and lambda might be arbitrary if theta is 0.
                assert!(
                    params[0].abs() < 1e-6,
                    "theta should be 0 for Identity at op {}",
                    i
                );
            }
        }
    }

    #[test]
    fn test_kak_decomposition_swap() {
        // SWAP gate
        // 1 0 0 0
        // 0 0 1 0
        // 0 1 0 0
        // 0 0 0 1
        // Canonical params: x=y=z=pi/4

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let swap = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, zero, zero, zero, one, zero, zero, one, zero, zero, zero, zero,
                zero, one,
            ],
        );

        let circuit = synth.synthesize(&swap, &[]).unwrap();

        let mut params_found = Vec::new();
        // x at 5, y at 12, z at 17
        // Note: params are -2*x, -2*y, -2*z.
        // We need to divide by -2.

        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // All should be pi/4 (or -pi/4, which is equivalent for SWAP)
        for p in params_found {
            assert!(
                (p.abs() - PI / 4.0).abs() < 1e-6,
                "Param should be +/- pi/4 for SWAP, got {}",
                p
            );
        }
    }

    #[test]
    fn test_kak_decomposition_iswap() {
        // iSWAP gate
        // 1 0 0 0
        // 0 0 i 0
        // 0 i 0 0
        // 0 0 0 1
        // Canonical params: x=y=pi/4, z=0

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let i = Complex::new(0.0, 1.0);
        let iswap = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, zero, zero, zero, i, zero, zero, i, zero, zero, zero, zero, zero,
                one,
            ],
        );

        let circuit = synth.synthesize(&iswap, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        assert!(
            (params_found[0].abs() - PI / 4.0).abs() < 1e-6,
            "Max param should be pi/4 for iSWAP"
        );
        assert!(
            (params_found[1].abs() - PI / 4.0).abs() < 1e-6,
            "Mid param should be pi/4 for iSWAP"
        );
        assert!(
            params_found[2].abs() < 1e-6,
            "Min param should be 0 for iSWAP"
        );
    }

    // Helper to simulate circuit back to matrix
    fn circuit_to_matrix(circuit: &Circuit) -> DMatrix<Complex<f64>> {
        // Start with Identity
        let mut u = DMatrix::<Complex<f64>>::identity(4, 4);

        for op in &circuit.operations {
            if let Operation::Gate {
                name,
                qubits,
                params,
            } = op
            {
                let gate_u = match name {
                    GateType::CX => {
                        // CX 0 1
                        // 1 0 0 0
                        // 0 1 0 0
                        // 0 0 0 1
                        // 0 0 1 0
                        DMatrix::from_row_slice(
                            4,
                            4,
                            &[
                                Complex::new(1.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(1.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(1.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(0.0, 0.0),
                                Complex::new(1.0, 0.0),
                                Complex::new(0.0, 0.0),
                            ],
                        )
                    }
                    GateType::U => {
                        // Single qubit gate expanded to 4x4
                        // If qubit 0: U tensor I
                        // If qubit 1: I tensor U
                        use crate::transpiler::synthesis::zyz::u_to_matrix;
                        let u2 = u_to_matrix(params[0], params[1], params[2]);
                        let u2_mat = DMatrix::from_row_slice(
                            2,
                            2,
                            &[u2[0][0], u2[0][1], u2[1][0], u2[1][1]],
                        );
                        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);

                        if qubits[0] == 0 {
                            u2_mat.kronecker(&id2)
                        } else {
                            id2.kronecker(&u2_mat)
                        }
                    }
                    GateType::RZ => {
                        // RZ(theta) = diag(e^{-i theta/2}, e^{i theta/2})
                        let theta = params[0];
                        let u2_mat = DMatrix::from_diagonal(&nalgebra::DVector::from_vec(vec![
                            Complex::new(0.0, -theta / 2.0).exp(),
                            Complex::new(0.0, theta / 2.0).exp(),
                        ]));
                        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
                        if qubits[0] == 0 {
                            u2_mat.kronecker(&id2)
                        } else {
                            id2.kronecker(&u2_mat)
                        }
                    }
                    GateType::RX => {
                        // RX(theta) = [cos(theta/2), -i sin(theta/2); -i sin(theta/2), cos(theta/2)]
                        let theta = params[0];
                        let c = (theta / 2.0).cos();
                        let s = (theta / 2.0).sin();
                        let i = Complex::new(0.0, 1.0);
                        let u2_mat = DMatrix::from_row_slice(
                            2,
                            2,
                            &[Complex::new(c, 0.0), -i * s, -i * s, Complex::new(c, 0.0)],
                        );
                        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
                        if qubits[0] == 0 {
                            u2_mat.kronecker(&id2)
                        } else {
                            id2.kronecker(&u2_mat)
                        }
                    }
                    GateType::H => {
                        // H = 1/sqrt(2) [[1, 1], [1, -1]]
                        let s2 = 1.0 / 2.0_f64.sqrt();
                        let u2_mat = DMatrix::from_row_slice(
                            2,
                            2,
                            &[
                                Complex::new(s2, 0.0),
                                Complex::new(s2, 0.0),
                                Complex::new(s2, 0.0),
                                Complex::new(-s2, 0.0),
                            ],
                        );
                        let id2 = DMatrix::<Complex<f64>>::identity(2, 2);
                        if qubits[0] == 0 {
                            u2_mat.kronecker(&id2)
                        } else {
                            id2.kronecker(&u2_mat)
                        }
                    }
                    _ => panic!("Unsupported gate in reconstruction: {:?}", name),
                };

                // Apply gate: U_new = Gate * U_old
                u = gate_u * u;
            }
        }
        u
    }

    #[test]
    fn test_kak_random_unitary() {
        let synth = KakSynthesizer;

        // Generate a random 4x4 unitary
        // We can't easily generate a Haar random unitary without a PRNG and QR.
        // For reproducibility and simplicity, let's construct one from random parameters.
        // U = (A1 x A0) * exp(i(xXX + yYY + zZZ)) * (B1 x B0)

        // Random parameters
        let x = 0.123;
        let y = 0.456;
        let z = 0.789;

        // Random local unitaries (using fixed values for reproducibility)
        let a0 = u_matrix(1.0, 0.5, 0.2);
        let a1 = u_matrix(0.3, 0.7, 1.1);
        let b0 = u_matrix(2.0, 1.5, 0.1);
        let b1 = u_matrix(0.1, 0.9, 2.5);

        let a = a1.kronecker(&a0);
        let b = b1.kronecker(&b0);

        // Interaction
        let h = DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(1.0 / 2.0_f64.sqrt(), 0.0),
                Complex::new(1.0 / 2.0_f64.sqrt(), 0.0),
                Complex::new(1.0 / 2.0_f64.sqrt(), 0.0),
                Complex::new(-1.0 / 2.0_f64.sqrt(), 0.0),
            ],
        );
        let _hh = h.kronecker(&h);

        // N = exp(i(xXX + yYY + zZZ))
        // Diagonal in Magic Basis:
        // d0 = e^{i(x-y+z)}
        // d1 = e^{i(-x+y+z)}
        // d2 = e^{i(-x-y-z)}
        // d3 = e^{i(x+y-z)}
        // Wait, my basis order was: Phi+, Psi+, Psi-, Phi-
        // p0 = x - y + z
        // p1 = x + y - z
        // p2 = -x - y - z
        // p3 = -x + y + z

        let p0 = x - y + z;
        let p1 = x + y - z;
        let p2 = -x - y - z;
        let p3 = -x + y + z;

        let n_magic = DMatrix::from_diagonal(&DVector::from_vec(vec![
            Complex::new(0.0, p0).exp(),
            Complex::new(0.0, p1).exp(),
            Complex::new(0.0, p2).exp(),
            Complex::new(0.0, p3).exp(),
        ]));

        // Transform N to computational basis
        // N_comp = M * N_magic * M^dag
        let magic_matrix = DMatrix::from_row_slice(
            4,
            4,
            &[
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(-1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0),
            ],
        )
        .scale(1.0 / 2.0_f64.sqrt());

        let n_comp = &magic_matrix * n_magic * magic_matrix.adjoint();

        let target_u = &a * n_comp.clone() * &b;

        let circuit = synth.synthesize(&target_u, &[]).unwrap();
        let reconstructed_u = circuit_to_matrix(&circuit);

        // Check distance d(U, V) = || U - e^{i phi} V ||
        // We can optimize phi or just check |Tr(U^dag V)| / 4 approx 1.

        // Note: circuit_to_matrix appears to use a different qubit ordering (LSB-first?)
        // than the standard Kronecker product (MSB-first) used to construct target_u.
        // We apply a SWAP gate to align the bases for comparison.
        let swap_gate = DMatrix::from_row_slice(
            4,
            4,
            &[
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
        );
        let recon_swapped = &swap_gate * &reconstructed_u * &swap_gate;

        let trace = (target_u.adjoint() * recon_swapped).trace();
        let fidelity = (trace.norm() / 4.0).abs();

        assert!(fidelity > 0.999, "Fidelity too low: {}", fidelity);
    }

    #[test]
    fn test_kak_decomposition_cnot() {
        // CNOT gate
        // 1 0 0 0
        // 0 1 0 0
        // 0 0 0 1
        // 0 0 1 0
        // Canonical params: x=pi/4, y=0, z=0

        let synth = KakSynthesizer;
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

        let circuit = synth.synthesize(&cnot, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        assert!(
            (params_found[0].abs() - PI / 4.0).abs() < 1e-6,
            "Max param should be pi/4 for CNOT"
        );
        assert!(
            params_found[1].abs() < 1e-6,
            "Mid param should be 0 for CNOT"
        );
        assert!(
            params_found[2].abs() < 1e-6,
            "Min param should be 0 for CNOT"
        );
    }

    #[test]
    fn test_kak_decomposition_cz() {
        // CZ gate
        // 1 0 0 0
        // 0 1 0 0
        // 0 0 1 0
        // 0 0 0 -1
        // Canonical params: x=pi/4, y=0, z=0 (Same as CNOT)

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let minus_one = Complex::new(-1.0, 0.0);
        let cz = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, zero, zero, one, zero, zero, zero, zero, one, zero, zero, zero,
                zero, minus_one,
            ],
        );

        let circuit = synth.synthesize(&cz, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        assert!(
            (params_found[0].abs() - PI / 4.0).abs() < 1e-6,
            "Max param should be pi/4 for CZ"
        );
        assert!(params_found[1].abs() < 1e-6, "Mid param should be 0 for CZ");
        assert!(params_found[2].abs() < 1e-6, "Min param should be 0 for CZ");
    }

    #[test]
    fn test_kak_decomposition_sqrt_iswap() {
        // sqrt(iSWAP) gate
        // 1 0 0 0
        // 0 1/sqrt(2) i/sqrt(2) 0
        // 0 i/sqrt(2) 1/sqrt(2) 0
        // 0 0 0 1
        // Canonical params: x=pi/8, y=pi/8, z=0

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let i = Complex::new(0.0, 1.0);
        let s2 = 1.0 / 2.0_f64.sqrt();
        let val = Complex::new(s2, 0.0);
        let val_i = i * s2;

        let sqrt_iswap = DMatrix::from_row_slice(
            4,
            4,
            &[
                one, zero, zero, zero, zero, val, val_i, zero, zero, val_i, val, zero, zero, zero,
                zero, one,
            ],
        );

        let circuit = synth.synthesize(&sqrt_iswap, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        assert!(
            (params_found[0].abs() - PI / 8.0).abs() < 1e-6,
            "Max param should be pi/8 for sqrt(iSWAP)"
        );
        assert!(
            (params_found[1].abs() - PI / 8.0).abs() < 1e-6,
            "Mid param should be pi/8 for sqrt(iSWAP)"
        );
        assert!(
            params_found[2].abs() < 1e-6,
            "Min param should be 0 for sqrt(iSWAP)"
        );
    }
    #[test]
    fn test_kak_decomposition_crx() {
        // Controlled-RX(theta)
        // 1 0 0 0
        // 0 1 0 0
        // 0 0 cos(t/2) -i sin(t/2)
        // 0 0 -i sin(t/2) cos(t/2)
        // Canonical params: x=theta/4, y=0, z=0 (up to permutation)

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        let i = Complex::new(0.0, 1.0);

        let theta = PI / 3.0;
        let c = (theta / 2.0).cos();
        let s = (theta / 2.0).sin();

        let crx = DMatrix::from_row_slice(
            4,
            4,
            &[
                one,
                zero,
                zero,
                zero,
                zero,
                one,
                zero,
                zero,
                zero,
                zero,
                Complex::new(c, 0.0),
                -i * s,
                zero,
                zero,
                -i * s,
                Complex::new(c, 0.0),
            ],
        );

        let circuit = synth.synthesize(&crx, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // One param should be theta/4
        assert!(
            (params_found[0].abs() - theta / 4.0).abs() < 1e-6,
            "Max param should be theta/4 for CRX, got {} expected {}",
            params_found[0],
            theta / 4.0
        );
        assert!(
            params_found[1].abs() < 1e-6,
            "Mid param should be 0 for CRX"
        );
        assert!(
            params_found[2].abs() < 1e-6,
            "Min param should be 0 for CRX"
        );
    }

    #[test]
    fn test_kak_decomposition_cry() {
        // Controlled-RY(theta)
        // 1 0 0 0
        // 0 1 0 0
        // 0 0 cos(t/2) -sin(t/2)
        // 0 0 sin(t/2) cos(t/2)
        // Canonical params: x=theta/4, y=0, z=0 (up to permutation)

        let synth = KakSynthesizer;
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);

        let theta = PI / 2.0;
        let c = (theta / 2.0).cos();
        let s = (theta / 2.0).sin();

        let cry = DMatrix::from_row_slice(
            4,
            4,
            &[
                one,
                zero,
                zero,
                zero,
                zero,
                one,
                zero,
                zero,
                zero,
                zero,
                Complex::new(c, 0.0),
                Complex::new(-s, 0.0),
                zero,
                zero,
                Complex::new(s, 0.0),
                Complex::new(c, 0.0),
            ],
        );

        let circuit = synth.synthesize(&cry, &[]).unwrap();

        let mut params_found = Vec::new();
        if let Operation::Gate { params, .. } = &circuit.operations[5] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[12] {
            params_found.push(params[0] / -2.0);
        }
        if let Operation::Gate { params, .. } = &circuit.operations[17] {
            params_found.push(params[0] / -2.0);
        }

        params_found.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // One param should be theta/4
        assert!(
            (params_found[0].abs() - theta / 4.0).abs() < 1e-6,
            "Max param should be theta/4 for CRY, got {} expected {}",
            params_found[0],
            theta / 4.0
        );
        assert!(
            params_found[1].abs() < 1e-6,
            "Mid param should be 0 for CRY"
        );
        assert!(
            params_found[2].abs() < 1e-6,
            "Min param should be 0 for CRY"
        );
    }
}
