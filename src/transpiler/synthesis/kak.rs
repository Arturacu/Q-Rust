//! KAK (Cartan) decomposition for two-qubit unitaries.
//!
//! Decomposes `U ∈ SU(4)` into
//! `U = (A1 ⊗ A0) · exp(i(x X⊗X + y Y⊗Y + z Z⊗Z)) · (B1 ⊗ B0)`,
//! synthesizing the result as an explicit gate sequence.
//!
//! The non-local term is emitted axis by axis. Each non-zero Cartan
//! coefficient costs two CX gates, so this construction uses 0, 2, 4, or 6 CX.
//!
//! The optimal circuit counts from Shende et al. 2004 are 0/1/2/3 CX; reaching
//! those requires a combined synthesis of the full interaction term rather than
//! decomposing X, Y, Z axes independently. Optimal synthesis is left as future work.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

const CARTAN_ZERO_TOL: f64 = 1e-7;
const INPUT_UNITARITY_TOL: f64 = 1e-8;
const DECOMPOSITION_TOL: f64 = 1e-6;
const OUTPUT_FIDELITY_TOL: f64 = 1e-9;
const MAX_PHASE_LIFT: i32 = 2;

fn complex_matrix_is_finite(matrix: &DMatrix<Complex<f64>>) -> bool {
    matrix
        .iter()
        .all(|value| value.re.is_finite() && value.im.is_finite())
}

fn real_matrix_is_finite(matrix: &DMatrix<f64>) -> bool {
    matrix.iter().all(|value| value.is_finite())
}

fn is_unitary(matrix: &DMatrix<Complex<f64>>, tolerance: f64) -> bool {
    if matrix.nrows() != matrix.ncols() || !complex_matrix_is_finite(matrix) {
        return false;
    }
    let identity = DMatrix::<Complex<f64>>::identity(matrix.nrows(), matrix.ncols());
    let residual = (matrix.adjoint() * matrix - identity).norm();
    residual.is_finite() && residual <= tolerance * matrix.nrows() as f64
}

fn off_diagonal_norm(matrix: &DMatrix<f64>) -> f64 {
    let squared: f64 = (0..matrix.nrows())
        .flat_map(|row| (0..matrix.ncols()).map(move |col| (row, col)))
        .filter(|(row, col)| row != col)
        .map(|(row, col)| matrix[(row, col)].powi(2))
        .sum();
    squared.sqrt()
}

fn complex_off_diagonal_norm(matrix: &DMatrix<Complex<f64>>) -> f64 {
    let squared: f64 = (0..matrix.nrows())
        .flat_map(|row| (0..matrix.ncols()).map(move |col| (row, col)))
        .filter(|(row, col)| row != col)
        .map(|(row, col)| matrix[(row, col)].norm_sqr())
        .sum();
    squared.sqrt()
}

fn magic_basis() -> DMatrix<Complex<f64>> {
    let scale = Complex::new(1.0 / 2.0_f64.sqrt(), 0.0);
    let one = Complex::new(1.0, 0.0);
    let zero = Complex::new(0.0, 0.0);
    let i = Complex::<f64>::i();
    DMatrix::from_row_slice(
        4,
        4,
        &[
            one, zero, zero, i, zero, i, one, zero, zero, i, -one, zero, one, zero, zero, -i,
        ],
    ) * scale
}

/// Simultaneously diagonalise two commuting real-symmetric matrices.
///
/// Diagonalising only one matrix is insufficient when it has repeated
/// eigenvalues. Within each degenerate eigenspace of `a`, we therefore
/// diagonalise the restriction of `b`.
fn simultaneous_diagonalizer(a: &DMatrix<f64>, b: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    if a.shape() != b.shape()
        || a.nrows() != a.ncols()
        || !real_matrix_is_finite(a)
        || !real_matrix_is_finite(b)
    {
        return None;
    }

    let scale = 1.0 + a.norm() + b.norm();
    if (a - a.transpose()).norm() > INPUT_UNITARITY_TOL * scale
        || (b - b.transpose()).norm() > INPUT_UNITARITY_TOL * scale
        || (a * b - b * a).norm() > DECOMPOSITION_TOL * scale
    {
        return None;
    }

    let eig_a = a.clone().symmetric_eigen();
    let mut order: Vec<usize> = (0..a.nrows()).collect();
    order.sort_by(|left, right| eig_a.eigenvalues[*left].total_cmp(&eig_a.eigenvalues[*right]));

    let mut eigenvalues = Vec::with_capacity(a.nrows());
    let mut basis = DMatrix::<f64>::zeros(a.nrows(), a.ncols());
    for (new_index, old_index) in order.into_iter().enumerate() {
        eigenvalues.push(eig_a.eigenvalues[old_index]);
        basis.set_column(new_index, &eig_a.eigenvectors.column(old_index));
    }

    let mut start = 0;
    while start < eigenvalues.len() {
        let reference = eigenvalues[start];
        let mut end = start + 1;
        while end < eigenvalues.len()
            && (eigenvalues[end] - reference).abs() <= DECOMPOSITION_TOL * scale
        {
            end += 1;
        }

        if end - start > 1 {
            let subspace = basis.columns(start, end - start).into_owned();
            let projected_b = subspace.transpose() * b * &subspace;
            let eig_b = projected_b.symmetric_eigen();
            let rotated = subspace * eig_b.eigenvectors;
            for column in 0..(end - start) {
                basis.set_column(start + column, &rotated.column(column));
            }
        }
        start = end;
    }

    if basis.determinant() < 0.0 {
        for row in 0..basis.nrows() {
            basis[(row, 0)] *= -1.0;
        }
    }

    let identity = DMatrix::<f64>::identity(basis.nrows(), basis.ncols());
    let orthogonality = (&basis.transpose() * &basis - identity).norm();
    let diagonal_a = basis.transpose() * a * &basis;
    let diagonal_b = basis.transpose() * b * &basis;
    if orthogonality > DECOMPOSITION_TOL * scale
        || off_diagonal_norm(&diagonal_a) > DECOMPOSITION_TOL * scale
        || off_diagonal_norm(&diagonal_b) > DECOMPOSITION_TOL * scale
    {
        return None;
    }

    Some(basis)
}

fn select_cartan_root(
    u_magic: &DMatrix<Complex<f64>>,
    r_magic: &DMatrix<Complex<f64>>,
    diagonal_system: &DMatrix<Complex<f64>>,
) -> Option<DMatrix<Complex<f64>>> {
    if complex_off_diagonal_norm(diagonal_system) > DECOMPOSITION_TOL {
        return None;
    }

    let square_roots: Vec<Complex<f64>> = diagonal_system
        .diagonal()
        .iter()
        .map(|value| value.sqrt())
        .collect();
    if square_roots
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return None;
    }

    let mut best: Option<(f64, f64, DMatrix<Complex<f64>>)> = None;
    for mask in 0..16 {
        let diagonal: Vec<Complex<f64>> = square_roots
            .iter()
            .enumerate()
            .map(|(bit, value)| {
                if (mask >> bit) & 1 == 1 {
                    -*value
                } else {
                    *value
                }
            })
            .collect();
        let determinant = diagonal.iter().product::<Complex<f64>>();
        if (determinant - Complex::new(1.0, 0.0)).norm() > DECOMPOSITION_TOL {
            continue;
        }

        let candidate = DMatrix::from_diagonal(&DVector::from_vec(diagonal.clone()));
        let local = u_magic * r_magic * candidate.adjoint();
        let imaginary_residual = local
            .iter()
            .map(|value| value.im.powi(2))
            .sum::<f64>()
            .sqrt();
        let principal_phase_norm = diagonal.iter().map(|value| value.arg().abs()).sum::<f64>();
        if !imaginary_residual.is_finite() || !principal_phase_norm.is_finite() {
            continue;
        }

        let replace = match &best {
            None => true,
            Some((best_residual, best_phase_norm, _)) => {
                imaginary_residual < *best_residual - INPUT_UNITARITY_TOL
                    || ((imaginary_residual - *best_residual).abs() <= INPUT_UNITARITY_TOL
                        && principal_phase_norm < *best_phase_norm)
            }
        };
        if replace {
            best = Some((imaginary_residual, principal_phase_norm, candidate));
        }
    }

    let (residual, _, root) = best?;
    (residual <= DECOMPOSITION_TOL).then_some(root)
}

#[derive(Clone, Copy)]
struct CartanChoice {
    permutation: [usize; 4],
    coordinates: [f64; 3],
    nonzero_coordinates: usize,
    max_magnitude: f64,
    total_magnitude: f64,
}

fn cartan_choice_is_better(candidate: &CartanChoice, current: &CartanChoice) -> bool {
    candidate.nonzero_coordinates < current.nonzero_coordinates
        || (candidate.nonzero_coordinates == current.nonzero_coordinates
            && (candidate.max_magnitude < current.max_magnitude - INPUT_UNITARITY_TOL
                || ((candidate.max_magnitude - current.max_magnitude).abs()
                    <= INPUT_UNITARITY_TOL
                    && candidate.total_magnitude < current.total_magnitude)))
}

/// Choose phase lifts and a Bell-basis ordering without changing the
/// represented diagonal unitary. This avoids folding Cartan coordinates
/// unless the corresponding local factors are transformed as well.
fn select_cartan_coordinates(diagonal: &[Complex<f64>]) -> Option<CartanChoice> {
    if diagonal.len() != 4
        || diagonal
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return None;
    }

    let phases: Vec<f64> = diagonal.iter().map(|value| value.arg()).collect();
    let tau = std::f64::consts::TAU;
    let winding = (phases.iter().sum::<f64>() / tau).round() as i32;
    if (phases.iter().sum::<f64>() - f64::from(winding) * tau).abs() > DECOMPOSITION_TOL {
        return None;
    }

    let mut best: Option<CartanChoice> = None;
    for p0 in 0..4 {
        for p1 in 0..4 {
            if p1 == p0 {
                continue;
            }
            for p2 in 0..4 {
                if p2 == p0 || p2 == p1 {
                    continue;
                }
                for p3 in 0..4 {
                    if p3 == p0 || p3 == p1 || p3 == p2 {
                        continue;
                    }
                    let permutation = [p0, p1, p2, p3];
                    for k0 in -MAX_PHASE_LIFT..=MAX_PHASE_LIFT {
                        for k1 in -MAX_PHASE_LIFT..=MAX_PHASE_LIFT {
                            for k2 in -MAX_PHASE_LIFT..=MAX_PHASE_LIFT {
                                let k3 = -winding - k0 - k1 - k2;
                                if !(-MAX_PHASE_LIFT..=MAX_PHASE_LIFT).contains(&k3) {
                                    continue;
                                }
                                let lifts = [k0, k1, k2, k3];
                                let lifted: [f64; 4] = std::array::from_fn(|index| {
                                    phases[permutation[index]] + f64::from(lifts[index]) * tau
                                });
                                if lifted.iter().sum::<f64>().abs() > DECOMPOSITION_TOL {
                                    continue;
                                }

                                let x = (lifted[0] + lifted[1]) / 2.0;
                                let y = (lifted[1] + lifted[3]) / 2.0;
                                let z = (lifted[0] + lifted[3]) / 2.0;
                                let coordinates = [x, y, z];
                                let magnitudes = [x.abs(), y.abs(), z.abs()];
                                let candidate = CartanChoice {
                                    permutation,
                                    coordinates,
                                    nonzero_coordinates: magnitudes
                                        .iter()
                                        .filter(|value| **value > CARTAN_ZERO_TOL)
                                        .count(),
                                    max_magnitude: magnitudes
                                        .into_iter()
                                        .max_by(f64::total_cmp)
                                        .unwrap_or(0.0),
                                    total_magnitude: magnitudes.into_iter().sum(),
                                };
                                let replace = match best.as_ref() {
                                    None => true,
                                    Some(current) => cartan_choice_is_better(&candidate, current),
                                };
                                if replace {
                                    best = Some(candidate);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    best
}

fn real_orthogonal_component(matrix: &DMatrix<Complex<f64>>) -> Option<DMatrix<f64>> {
    if matrix.nrows() != matrix.ncols() || !complex_matrix_is_finite(matrix) {
        return None;
    }
    let imaginary_residual = matrix
        .iter()
        .map(|value| value.im.powi(2))
        .sum::<f64>()
        .sqrt();
    if imaginary_residual > DECOMPOSITION_TOL {
        return None;
    }
    let real = matrix.map(|value| value.re);
    let identity = DMatrix::<f64>::identity(real.nrows(), real.ncols());
    let orthogonality = (real.transpose() * &real - identity).norm();
    if orthogonality > DECOMPOSITION_TOL || real.determinant() < 0.0 {
        return None;
    }
    Some(real)
}

fn add_gate(circuit: &mut Circuit, name: GateType, qubits: &[usize], params: &[f64]) {
    circuit.add_op(Operation::Gate {
        name,
        qubits: qubits.to_vec(),
        params: params.to_vec(),
    });
}

fn emit_xx(circuit: &mut Circuit, angle: f64) {
    add_gate(circuit, GateType::H, &[0], &[]);
    add_gate(circuit, GateType::H, &[1], &[]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RZ, &[1], &[-2.0 * angle]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::H, &[0], &[]);
    add_gate(circuit, GateType::H, &[1], &[]);
}

fn emit_yy(circuit: &mut Circuit, angle: f64) {
    let half_pi = std::f64::consts::FRAC_PI_2;
    add_gate(circuit, GateType::RX, &[0], &[half_pi]);
    add_gate(circuit, GateType::RX, &[1], &[half_pi]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RZ, &[1], &[-2.0 * angle]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RX, &[0], &[-half_pi]);
    add_gate(circuit, GateType::RX, &[1], &[-half_pi]);
}

fn emit_zz(circuit: &mut Circuit, angle: f64) {
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RZ, &[1], &[-2.0 * angle]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
}

/// Exact two-qubit (Cartan/KAK) synthesizer.
///
/// Input unitarity, decomposition structure, tensor factors, and final process
/// fidelity are checked before a circuit is returned. Numerical failures yield
/// `None` through the [`Synthesizer`] contract.
#[derive(Debug, Clone, Copy)]
pub struct KakSynthesizer;

impl Synthesizer for KakSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, _basis: &[GateType]) -> Option<Circuit> {
        if unitary.nrows() != 4 || unitary.ncols() != 4 || !is_unitary(unitary, INPUT_UNITARITY_TOL)
        {
            return None;
        }

        let det = unitary.determinant();
        let phase = det.powf(0.25);
        if !phase.re.is_finite() || !phase.im.is_finite() || phase.norm() < INPUT_UNITARITY_TOL {
            return None;
        }
        let su4_unitary = unitary / phase;
        let magic_matrix = magic_basis();

        let u_magic = magic_matrix.adjoint() * su4_unitary * &magic_matrix;
        let m_sys = u_magic.transpose() * &u_magic;

        let mut a = DMatrix::<f64>::zeros(4, 4);
        let mut b = DMatrix::<f64>::zeros(4, 4);
        for ii in 0..4 {
            for jj in 0..4 {
                a[(ii, jj)] = m_sys[(ii, jj)].re;
                b[(ii, jj)] = m_sys[(ii, jj)].im;
            }
        }
        let r_magic_real = simultaneous_diagonalizer(&a, &b)?;
        let r_magic = r_magic_real.map(|x| Complex::new(x, 0.0));

        let diagonal_system = r_magic.transpose() * m_sys * &r_magic;
        let n_magic = select_cartan_root(&u_magic, &r_magic, &diagonal_system)?;
        let diagonal: Vec<Complex<f64>> = n_magic.diagonal().iter().copied().collect();
        let choice = select_cartan_coordinates(&diagonal)?;
        let [x, y, z] = choice.coordinates;

        let mut r_perm = DMatrix::zeros(4, 4);
        let mut n_perm_diag = Vec::with_capacity(4);
        for ii in 0..4 {
            let old = choice.permutation[ii];
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
        let l_magic = real_orthogonal_component(&l_magic)?.map(|value| Complex::new(value, 0.0));

        let l_comp = &magic_matrix * l_magic * magic_matrix.adjoint();
        let r_comp = &magic_matrix * r_magic.transpose() * magic_matrix.adjoint();

        let (a0, a1) = decompose_tensor_product(&l_comp)?;
        let (b0, b1) = decompose_tensor_product(&r_comp)?;

        use crate::transpiler::synthesis::zyz::zyz_decomposition;
        fn to_array(m: DMatrix<Complex<f64>>) -> [[Complex<f64>; 2]; 2] {
            [[m[(0, 0)], m[(0, 1)]], [m[(1, 0)], m[(1, 1)]]]
        }
        let (a0_t, a0_p, a0_l, _) = zyz_decomposition(to_array(a0));
        let (a1_t, a1_p, a1_l, _) = zyz_decomposition(to_array(a1));
        let (b0_t, b0_p, b0_l, _) = zyz_decomposition(to_array(b0));
        let (b1_t, b1_p, b1_l, _) = zyz_decomposition(to_array(b1));

        let nx = x.abs() > CARTAN_ZERO_TOL;
        let ny = y.abs() > CARTAN_ZERO_TOL;
        let nz = z.abs() > CARTAN_ZERO_TOL;

        let mut circuit = Circuit::new(2, 0);
        add_gate(&mut circuit, GateType::U, &[1], &[b0_t, b0_p, b0_l]);
        add_gate(&mut circuit, GateType::U, &[0], &[b1_t, b1_p, b1_l]);
        if nx {
            emit_xx(&mut circuit, x);
        }
        if ny {
            emit_yy(&mut circuit, y);
        }
        if nz {
            emit_zz(&mut circuit, z);
        }
        add_gate(&mut circuit, GateType::U, &[1], &[a0_t, a0_p, a0_l]);
        add_gate(&mut circuit, GateType::U, &[0], &[a1_t, a1_p, a1_l]);

        let synthesized = crate::simulator::try_circuit_to_unitary(&circuit).ok()?;
        let fidelity = crate::simulator::unitary_fidelity(unitary, &synthesized);
        (fidelity >= 1.0 - OUTPUT_FIDELITY_TOL).then_some(circuit)
    }
}

fn decompose_tensor_product(
    u: &DMatrix<Complex<f64>>,
) -> Option<(DMatrix<Complex<f64>>, DMatrix<Complex<f64>>)> {
    if u.nrows() != 4 || u.ncols() != 4 || !complex_matrix_is_finite(u) {
        return None;
    }

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
    let mut pivot = (0, 0);
    for row in 0..4 {
        for col in 0..4 {
            if t[(row, col)].norm_sqr() > t[pivot].norm_sqr() {
                pivot = (row, col);
            }
        }
    }
    let pivot_value = t[pivot];
    if pivot_value.norm() <= INPUT_UNITARITY_TOL {
        return None;
    }
    let u_vec = t.column(pivot.1).into_owned();
    let mut v_vec = DVector::zeros(4);
    for col in 0..4 {
        v_vec[col] = t[(pivot.0, col)] / pivot_value;
    }
    let rank_one = &u_vec * v_vec.transpose();
    if (&t - rank_one).norm() > DECOMPOSITION_TOL * (1.0 + t.norm()) {
        return None;
    }

    let mut a = DMatrix::zeros(2, 2);
    let mut b = DMatrix::zeros(2, 2);
    for i in 0..2 {
        for j in 0..2 {
            a[(i, j)] = u_vec[2 * i + j];
            b[(i, j)] = v_vec[2 * i + j];
        }
    }
    let det_a = a.determinant();
    if !det_a.re.is_finite() || !det_a.im.is_finite() || det_a.norm() <= INPUT_UNITARITY_TOL {
        return None;
    }
    a /= det_a.powf(0.5);
    let det_b = b.determinant();
    if !det_b.re.is_finite() || !det_b.im.is_finite() || det_b.norm() <= INPUT_UNITARITY_TOL {
        return None;
    }
    b /= det_b.powf(0.5);

    if !is_unitary(&a, DECOMPOSITION_TOL) || !is_unitary(&b, DECOMPOSITION_TOL) {
        return None;
    }
    let reconstructed = tensor_product_2x2(&a, &b);
    let trace = (u.adjoint() * reconstructed).trace();
    let fidelity = trace.norm_sqr() / 16.0;
    (fidelity >= 1.0 - DECOMPOSITION_TOL).then_some((a, b))
}

fn tensor_product_2x2(
    left: &DMatrix<Complex<f64>>,
    right: &DMatrix<Complex<f64>>,
) -> DMatrix<Complex<f64>> {
    let mut output = DMatrix::zeros(4, 4);
    for left_row in 0..2 {
        for right_row in 0..2 {
            for left_col in 0..2 {
                for right_col in 0..2 {
                    output[(2 * left_row + right_row, 2 * left_col + right_col)] =
                        left[(left_row, left_col)] * right[(right_row, right_col)];
                }
            }
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulator::{circuit_to_unitary, unitary_fidelity};

    fn count_cx(c: &Circuit) -> usize {
        c.operations
            .iter()
            .filter(|op| {
                matches!(
                    op,
                    Operation::Gate {
                        name: GateType::CX,
                        ..
                    }
                )
            })
            .count()
    }

    fn assert_synthesis_fidelity(target: &DMatrix<Complex<f64>>, circuit: &Circuit) {
        let actual = circuit_to_unitary(circuit);
        let fidelity = unitary_fidelity(target, &actual);
        assert!(
            fidelity >= 1.0 - 1e-9,
            "KAK synthesis fidelity {fidelity:.12} is below tolerance"
        );
    }

    fn add_u(circuit: &mut Circuit, qubit: usize, theta: f64, phi: f64, lambda: f64) {
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![qubit],
            params: vec![theta, phi, lambda],
        });
    }

    #[test]
    fn test_kak_identity_zero_cx() {
        let s = KakSynthesizer;
        let id = DMatrix::<Complex<f64>>::identity(4, 4);
        let circuit = s.synthesize(&id, &[]).unwrap();
        assert_eq!(count_cx(&circuit), 0);
        assert_synthesis_fidelity(&id, &circuit);
    }

    #[test]
    fn test_kak_local_product_zero_cx() {
        let mut source = Circuit::new(2, 0);
        add_u(&mut source, 0, 0.37, -0.21, 0.44);
        add_u(&mut source, 1, -0.18, 0.32, -0.29);
        let target = circuit_to_unitary(&source);

        let circuit = KakSynthesizer.synthesize(&target, &[]).unwrap();
        assert_eq!(count_cx(&circuit), 0);
        assert_synthesis_fidelity(&target, &circuit);
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
        assert_eq!(count_cx(&circuit), 2);
        assert_synthesis_fidelity(&cnot, &circuit);
    }

    #[test]
    fn test_kak_generic_unitary_fidelity() {
        let mut source = Circuit::new(2, 0);
        add_u(&mut source, 0, 0.31, -0.27, 0.19);
        add_u(&mut source, 1, -0.42, 0.11, 0.37);
        source.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        add_u(&mut source, 0, 0.23, 0.41, -0.17);
        add_u(&mut source, 1, -0.29, 0.13, 0.47);
        source.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 0],
            params: vec![],
        });
        add_u(&mut source, 0, -0.34, 0.21, 0.09);
        add_u(&mut source, 1, 0.28, -0.16, 0.33);
        source.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let target = circuit_to_unitary(&source);
        let circuit = KakSynthesizer.synthesize(&target, &[]).unwrap();
        assert_eq!(count_cx(&circuit), 6);
        assert_synthesis_fidelity(&target, &circuit);
    }

    #[test]
    fn test_kak_deterministic_generic_corpus() {
        fn next_angle(state: &mut u64) -> f64 {
            *state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let unit_interval = (*state >> 11) as f64 / ((1_u64 << 53) as f64);
            (2.0 * unit_interval - 1.0) * std::f64::consts::PI
        }

        let mut state = 0x4b41_4b2d_5152_5354;
        for case in 0..32 {
            let mut source = Circuit::new(2, 0);
            for layer in 0..4 {
                for qubit in 0..2 {
                    add_u(
                        &mut source,
                        qubit,
                        next_angle(&mut state),
                        next_angle(&mut state),
                        next_angle(&mut state),
                    );
                }
                if layer < 3 {
                    let qubits = if layer % 2 == 0 {
                        vec![0, 1]
                    } else {
                        vec![1, 0]
                    };
                    source.add_op(Operation::Gate {
                        name: GateType::CX,
                        qubits,
                        params: vec![],
                    });
                }
            }

            let target = circuit_to_unitary(&source);
            let circuit = KakSynthesizer
                .synthesize(&target, &[])
                .unwrap_or_else(|| panic!("KAK synthesis rejected deterministic case {case}"));
            assert!(count_cx(&circuit) <= 6);
            assert_synthesis_fidelity(&target, &circuit);
        }
    }

    #[test]
    fn test_kak_two_axis_interaction_uses_four_cx() {
        let mut source = Circuit::new(2, 0);
        source.add_op(Operation::Gate {
            name: GateType::RXX,
            qubits: vec![0, 1],
            params: vec![0.46],
        });
        source.add_op(Operation::Gate {
            name: GateType::RYY,
            qubits: vec![0, 1],
            params: vec![-0.34],
        });
        let target = circuit_to_unitary(&source);

        let circuit = KakSynthesizer.synthesize(&target, &[]).unwrap();
        assert_eq!(count_cx(&circuit), 4);
        assert_synthesis_fidelity(&target, &circuit);
    }

    #[test]
    fn test_simultaneous_diagonalizer_resolves_degenerate_subspace() {
        let a = DMatrix::from_diagonal(&DVector::from_vec(vec![1.0, 1.0, 2.0, 3.0]));
        let b = DMatrix::from_row_slice(
            4,
            4,
            &[
                0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 5.0,
            ],
        );

        let basis = simultaneous_diagonalizer(&a, &b).unwrap();
        assert!(off_diagonal_norm(&(basis.transpose() * a * &basis)) <= DECOMPOSITION_TOL);
        assert!(off_diagonal_norm(&(basis.transpose() * b * &basis)) <= DECOMPOSITION_TOL);
    }

    #[test]
    fn test_kak_rejects_invalid_input() {
        let wrong_shape = DMatrix::<Complex<f64>>::identity(2, 2);
        assert!(KakSynthesizer.synthesize(&wrong_shape, &[]).is_none());

        let mut non_finite = DMatrix::<Complex<f64>>::identity(4, 4);
        non_finite[(0, 0)] = Complex::new(f64::NAN, 0.0);
        assert!(KakSynthesizer.synthesize(&non_finite, &[]).is_none());

        let non_unitary = DMatrix::<Complex<f64>>::from_element(4, 4, Complex::new(1.0, 0.0));
        assert!(KakSynthesizer.synthesize(&non_unitary, &[]).is_none());
    }
}
