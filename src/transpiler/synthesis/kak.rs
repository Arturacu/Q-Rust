//! KAK (Cartan) decomposition for two-qubit unitaries.
//!
//! Decomposes `U ∈ SU(4)` into
//! `U = (A1 ⊗ A0) · exp(i(x X⊗X + y Y⊗Y + z Z⊗Z)) · (B1 ⊗ B0)`,
//! synthesizing the result as an explicit gate sequence.
//!
//! The interaction is synthesized jointly with at most three CX gates. Local,
//! CNOT-equivalent, and two-CX families use fewer entanglers. The constructions
//! are established Cartan synthesis, not a new decomposition theorem:
//! Vatan & Williams, PRA 69, 032315 (2004), arXiv:quant-ph/0308006;
//! Shende, Bullock & Markov, PRA 70, 012310 (2004), arXiv:quant-ph/0308045.
//! Every returned circuit passes a phase-aligned operator-residual check.

use crate::ir::gate_def::GateDefinition;
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

// Classification only. A reduced-count candidate must also pass the tighter
// residual bound; otherwise retry the untruncated three-CX construction.
const CARTAN_ZERO_TOL: f64 = 1e-12;
const REDUCED_RESIDUAL_TOL: f64 = 1e-12;
const INPUT_UNITARITY_TOL: f64 = 1e-11;
const DECOMPOSITION_TOL: f64 = 1e-11;
const OUTPUT_RESIDUAL_TOL: f64 = 1e-10;

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
/// eigenvalues. Within each degenerate eigenspace of a real projection,
/// diagonalise the restriction of `b` and check both original matrices.
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

    // A single projection can have an almost-colliding spectrum even when
    // the complex eigenvalues are well separated (cosines coincide near 0
    // and π, for example). Try deterministic real projections, validating
    // against both original matrices after resolving repeated eigenspaces.
    // Usually the first projection suffices; the alternatives avoid dividing
    // by an artificially tiny real eigenvalue gap.
    for weight in [std::f64::consts::SQRT_2, 0.0, 1.0, -1.0] {
        let projection = a + b * weight;
        let eig = projection.symmetric_eigen();
        let mut order: Vec<usize> = (0..a.nrows()).collect();
        order.sort_by(|left, right| eig.eigenvalues[*left].total_cmp(&eig.eigenvalues[*right]));

        let mut eigenvalues = Vec::with_capacity(a.nrows());
        let mut basis = DMatrix::<f64>::zeros(a.nrows(), a.ncols());
        for (new_index, old_index) in order.into_iter().enumerate() {
            eigenvalues.push(eig.eigenvalues[old_index]);
            basis.set_column(new_index, &eig.eigenvectors.column(old_index));
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
        if orthogonality <= DECOMPOSITION_TOL
            && off_diagonal_norm(&diagonal_a) <= DECOMPOSITION_TOL / 2.0
            && off_diagonal_norm(&diagonal_b) <= DECOMPOSITION_TOL / 2.0
        {
            return Some(basis);
        }
    }
    None
}

/// Principal square roots already give a Takagi factor, up to column signs.
/// Fix the determinant with one sign; enumerating all 16 choices cannot improve
/// the imaginary residual because changing a root just flips one real column.
fn select_cartan_root(diagonal_system: &DMatrix<Complex<f64>>) -> Option<DMatrix<Complex<f64>>> {
    if complex_off_diagonal_norm(diagonal_system) > DECOMPOSITION_TOL {
        return None;
    }
    let mut roots: Vec<Complex<f64>> = diagonal_system
        .diagonal()
        .iter()
        .map(|value| Complex::from_polar(1.0, value.arg() / 2.0))
        .collect();
    let determinant = roots.iter().product::<Complex<f64>>();
    if (determinant + Complex::new(1.0, 0.0)).norm() < DECOMPOSITION_TOL {
        roots[3] = -roots[3];
    } else if (determinant - Complex::new(1.0, 0.0)).norm() > DECOMPOSITION_TOL {
        return None;
    }
    Some(DMatrix::from_diagonal(&DVector::from_vec(roots)))
}

/// Recover one valid Cartan representative without exhaustive phase lifts.
/// For magic-basis phases d0..d3, x=(d0+d1)/2, y=(d1+d3)/2,
/// z=(d0+d3)/2 reproduces all four eigenvalues: the fourth phase differs
/// only by the determinant's integer multiple of 2π.
fn cartan_coordinates(diagonal: &DMatrix<Complex<f64>>) -> [f64; 3] {
    let phase = |i| diagonal[(i, i)].arg();
    [
        (phase(0) + phase(1)) / 2.0,
        (phase(1) + phase(3)) / 2.0,
        (phase(0) + phase(3)) / 2.0,
    ]
}

/// One shared global phase; no fidelity-to-amplitude conversion or
/// per-column phase alignment. The unnormalized Frobenius residual also
/// upper-bounds the error of every normalized input state.
fn phase_aligned_residual(target: &DMatrix<Complex<f64>>, actual: &DMatrix<Complex<f64>>) -> f64 {
    let overlap = (actual.adjoint() * target).trace();
    if !overlap.re.is_finite() || !overlap.im.is_finite() || overlap.norm() == 0.0 {
        return f64::INFINITY;
    }
    let phase = overlap / overlap.norm();
    (target - actual * phase).norm()
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

fn emit_three_cx(circuit: &mut Circuit, [x, y, z]: [f64; 3]) {
    // Joint Cartan network (Vatan–Williams three-entangler construction,
    // in an equivalent RX/RY/RZ local frame). Valid for arbitrary angles.
    let half_pi = std::f64::consts::FRAC_PI_2;
    add_gate(circuit, GateType::RX, &[0], &[half_pi]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RX, &[0], &[half_pi - 2.0 * x]);
    add_gate(circuit, GateType::RY, &[1], &[half_pi - 2.0 * y]);
    add_gate(circuit, GateType::CX, &[1, 0], &[]);
    add_gate(circuit, GateType::RX, &[1], &[-half_pi]);
    add_gate(circuit, GateType::RZ, &[1], &[half_pi - 2.0 * z]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
}

fn emit_two_cx(circuit: &mut Circuit, x: f64, y: f64) {
    // CX conjugates X(control) to XX and Y(target) to ZY. The outer
    // control RX rotations turn ZY into YY, yielding exp(i(xXX+yYY)).
    let half_pi = std::f64::consts::FRAC_PI_2;
    add_gate(circuit, GateType::RX, &[0], &[half_pi]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RX, &[0], &[-2.0 * x]);
    add_gate(circuit, GateType::RY, &[1], &[-2.0 * y]);
    add_gate(circuit, GateType::CX, &[0, 1], &[]);
    add_gate(circuit, GateType::RX, &[0], &[-half_pi]);
}

fn pair_frame(circuit: &mut Circuit, axis: usize, inverse: bool) {
    // Frame ZZ as XX (H) or YY (RX). Both Pauli signs cancel in P⊗P.
    for q in 0..2 {
        match axis {
            0 => add_gate(circuit, GateType::H, &[q], &[]),
            1 => add_gate(
                circuit,
                GateType::RX,
                &[q],
                &[if inverse {
                    -std::f64::consts::FRAC_PI_2
                } else {
                    std::f64::consts::FRAC_PI_2
                }],
            ),
            _ => {}
        }
    }
}

fn emit_reduced_cartan(circuit: &mut Circuit, coordinates: [f64; 3]) {
    let active: Vec<_> = (0..3)
        .filter(|&i| coordinates[i].abs() > CARTAN_ZERO_TOL)
        .collect();
    if active.is_empty() {
        return;
    }
    if active.len() == 1 {
        let axis = active[0];
        let angle = coordinates[axis];
        if (angle.abs() - std::f64::consts::FRAC_PI_4).abs() <= CARTAN_ZERO_TOL {
            // exp(±iπZZ/4) is CZ with local Z quarter turns, up to phase.
            pair_frame(circuit, axis, false);
            add_gate(circuit, GateType::H, &[1], &[]);
            add_gate(circuit, GateType::CX, &[0, 1], &[]);
            add_gate(circuit, GateType::H, &[1], &[]);
            for q in 0..2 {
                add_gate(
                    circuit,
                    GateType::RZ,
                    &[q],
                    &[-angle.signum() * std::f64::consts::FRAC_PI_2],
                );
            }
            pair_frame(circuit, axis, true);
            return;
        }
    }
    if active.len() <= 2 {
        // Move a zero coefficient to ZZ using a local frame, then use the
        // exact two-CX XX+YY network. No input-state assumptions are made.
        if coordinates[2].abs() <= CARTAN_ZERO_TOL {
            emit_two_cx(circuit, coordinates[0], coordinates[1]);
        } else if coordinates[0].abs() <= CARTAN_ZERO_TOL {
            pair_frame(circuit, 0, false);
            emit_two_cx(circuit, coordinates[2], coordinates[1]);
            pair_frame(circuit, 0, true);
        } else {
            pair_frame(circuit, 1, false);
            emit_two_cx(circuit, coordinates[0], coordinates[2]);
            pair_frame(circuit, 1, true);
        }
    } else {
        emit_three_cx(circuit, coordinates);
    }
}

/// Stable U parameters: atan2 retains tiny off-diagonal rotations that acos
/// of a rounded |u00| can lose. Singular endpoints are handled separately.
fn append_local(circuit: &mut Circuit, qubit: usize, matrix: &DMatrix<Complex<f64>>) {
    let (theta, phi, lambda, _) = super::zyz::zyz_decomposition([
        [matrix[(0, 0)], matrix[(0, 1)]],
        [matrix[(1, 0)], matrix[(1, 1)]],
    ]);
    add_gate(circuit, GateType::U, &[qubit], &[theta, phi, lambda]);
}

fn fuse_local_gates(circuit: &Circuit) -> Circuit {
    let mut output = Circuit::new(2, 0);
    let identity = DMatrix::<Complex<f64>>::identity(2, 2);
    let mut pending = [identity.clone(), identity.clone()];
    let flush = |q: usize, pending: &mut [DMatrix<Complex<f64>>; 2], output: &mut Circuit| {
        if phase_aligned_residual(&identity, &pending[q]) > 4.0 * f64::EPSILON {
            append_local(output, q, &pending[q]);
        }
        pending[q] = identity.clone();
    };
    for op in &circuit.operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            if qubits.len() == 1 {
                let q = qubits[0];
                pending[q] = name.unitary(params) * &pending[q];
            } else {
                flush(0, &mut pending, &mut output);
                flush(1, &mut pending, &mut output);
                output.add_op(op.clone());
            }
        }
    }
    flush(0, &mut pending, &mut output);
    flush(1, &mut pending, &mut output);
    output
}

/// Exact two-qubit (Cartan/KAK) synthesizer.
///
/// Input unitarity, decomposition structure, tensor factors, and final
/// phase-aligned operator residual are checked before returning U/CX output.
/// Reduced-count classifications are checked at a tighter residual tolerance
/// and retried without coordinate truncation if necessary. Numerical failures yield
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
        let n_magic = select_cartan_root(&diagonal_system)?;
        let raw_coordinates = cartan_coordinates(&n_magic);
        // exp(i(a+kπ/2)P⊗P) = i^k(P⊗P)^k exp(iaP⊗P).
        // Retain these local Pauli factors explicitly; coordinate folding
        // without them would silently change the represented unitary.
        let half_pi = std::f64::consts::FRAC_PI_2;
        let shifts: [i32; 3] = raw_coordinates.map(|a| (a / half_pi).round() as i32);
        let coordinates: [f64; 3] =
            std::array::from_fn(|i| raw_coordinates[i] - f64::from(shifts[i]) * half_pi);
        let l_magic = u_magic.clone() * r_magic.clone() * n_magic.adjoint();
        let l_magic = real_orthogonal_component(&l_magic)?.map(|value| Complex::new(value, 0.0));

        let l_comp = &magic_matrix * l_magic * magic_matrix.adjoint();
        let r_comp = &magic_matrix * r_magic.transpose() * magic_matrix.adjoint();

        let (a0, a1) = decompose_tensor_product(&l_comp)?;
        let (b0, b1) = decompose_tensor_product(&r_comp)?;

        for reduced in [true, false] {
            let mut circuit = Circuit::new(2, 0);
            append_local(&mut circuit, 1, &b0);
            append_local(&mut circuit, 0, &b1);
            if reduced {
                emit_reduced_cartan(&mut circuit, coordinates);
            } else {
                emit_three_cx(&mut circuit, coordinates);
            }
            for (axis, shift) in shifts.iter().enumerate() {
                if shift.rem_euclid(2) != 0 {
                    let pauli = [GateType::X, GateType::Y, GateType::Z][axis].clone();
                    for q in 0..2 {
                        add_gate(&mut circuit, pauli.clone(), &[q], &[]);
                    }
                }
            }
            append_local(&mut circuit, 1, &a0);
            append_local(&mut circuit, 0, &a1);
            let circuit = fuse_local_gates(&circuit);
            let synthesized = crate::simulator::try_circuit_to_unitary(&circuit).ok()?;
            let residual = phase_aligned_residual(unitary, &synthesized);
            let tolerance = if reduced {
                REDUCED_RESIDUAL_TOL
            } else {
                OUTPUT_RESIDUAL_TOL
            };
            if residual <= tolerance {
                return Some(circuit);
            }
        }
        None
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
    (phase_aligned_residual(u, &reconstructed) <= OUTPUT_RESIDUAL_TOL).then_some((a, b))
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
    use crate::simulator::circuit_to_unitary;

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
        let residual = phase_aligned_residual(target, &actual);
        assert!(
            residual <= OUTPUT_RESIDUAL_TOL,
            "KAK phase-aligned operator residual {residual:.3e} exceeds tolerance"
        );
        assert!(circuit.operations.iter().all(|op| matches!(
            op,
            Operation::Gate {
                name: GateType::U | GateType::CX,
                ..
            }
        )));
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
        assert_eq!(count_cx(&circuit), 1);
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
        assert_eq!(count_cx(&circuit), 3);
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
            assert!(count_cx(&circuit) <= 3);
            assert_synthesis_fidelity(&target, &circuit);
        }
    }

    #[test]
    fn test_kak_two_axis_interaction_uses_two_cx() {
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
        assert_eq!(count_cx(&circuit), 2);
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

    // Targets below are built directly from analytic matrices, independently
    // of both the emitted network and the simulator's gate decompositions.
    fn analytic_local(theta: f64, phi: f64, lambda: f64) -> DMatrix<Complex<f64>> {
        let (sin, cos) = (theta / 2.0).sin_cos();
        DMatrix::from_row_slice(
            2,
            2,
            &[
                Complex::new(cos, 0.0),
                -Complex::from_polar(sin, lambda),
                Complex::from_polar(sin, phi),
                Complex::from_polar(cos, phi + lambda),
            ],
        )
    }

    fn analytic_cartan(coordinates: [f64; 3]) -> DMatrix<Complex<f64>> {
        let identity = DMatrix::<Complex<f64>>::identity(4, 4);
        let mut output = identity.clone();
        for (axis, angle) in coordinates.into_iter().enumerate() {
            let mut pauli = DMatrix::<Complex<f64>>::zeros(4, 4);
            for col in 0..4 {
                let row = if axis == 2 { col } else { col ^ 3 };
                let sign = match axis {
                    0 => 1.0,
                    1 => {
                        if col == 0 || col == 3 {
                            -1.0
                        } else {
                            1.0
                        }
                    }
                    _ => {
                        if col == 0 || col == 3 {
                            1.0
                        } else {
                            -1.0
                        }
                    }
                };
                pauli[(row, col)] = Complex::new(sign, 0.0);
            }
            output = (&identity * Complex::new(angle.cos(), 0.0)
                + pauli * Complex::new(0.0, angle.sin()))
                * output;
        }
        output
    }

    fn dressed_cartan(coordinates: [f64; 3], phase: f64) -> DMatrix<Complex<f64>> {
        let before = tensor_product_2x2(
            &analytic_local(0.39, -0.67, 0.51),
            &analytic_local(-0.87, 0.23, 1.11),
        );
        let after = tensor_product_2x2(
            &analytic_local(1.03, 0.19, -0.43),
            &analytic_local(-0.33, -0.91, 0.71),
        );
        after * analytic_cartan(coordinates) * before * Complex::from_polar(1.0, phase)
    }

    #[test]
    fn analytic_reduced_families_and_noncanonical_angles() {
        let h = std::f64::consts::FRAC_PI_2;
        let q = std::f64::consts::FRAC_PI_4;
        let cases = [
            ([0.0, 0.0, 0.0], 0),
            ([h, -h, 2.0 * h], 0),
            ([q, 0.0, 0.0], 1),
            ([0.0, -q, 0.0], 1),
            ([0.0, 0.0, q + h], 1),
            ([0.34, 0.0, 0.0], 2),
            ([0.0, -0.23, 0.17], 2),
            ([0.27, 0.0, -0.31], 2),
            ([0.42 + h, -0.14 - h, h], 2),
            ([q, q, 0.0], 2), // iSWAP interaction class.
            ([q, q, q], 3),   // SWAP interaction class.
            ([0.22, -0.37, 0.16], 3),
            ([1.91, -2.73, 3.27], 3),
        ];
        for (case, (coordinates, expected_cx)) in cases.into_iter().enumerate() {
            for phase in [0.0, 0.29, -2.7, std::f64::consts::PI] {
                let target = dressed_cartan(coordinates, phase);
                let circuit = KakSynthesizer
                    .synthesize(&target, &[])
                    .unwrap_or_else(|| panic!("rejected analytic case {case}, phase {phase}"));
                assert_eq!(count_cx(&circuit), expected_cx, "analytic case {case}");
                assert_synthesis_fidelity(&target, &circuit);
            }
        }
    }

    #[test]
    fn degenerate_spectra_and_threshold_neighbors() {
        let q = std::f64::consts::FRAC_PI_4;
        let cases = [
            [0.23, 0.23, 0.23],
            [0.23, 0.23, 0.0],
            [0.23, -0.23, 0.23],
            [q, q, q - 1e-9],
            [q + 2e-11, 0.0, 0.0],
            [q - 2e-11, 0.0, 0.0],
            [1e-8, 2e-8, 3e-8],
            [1e-5, 2e-5, -3e-5],
            [0.27, -0.19, 2e-11],
            [0.27, -0.19, -2e-11],
            [0.27, -0.19, 8e-13],
            [0.27, -0.19, -8e-13],
            [0.0, 0.0, 8e-13],
            [0.0, 0.0, -8e-13],
        ];
        for (case, coordinates) in cases.into_iter().enumerate() {
            let target = dressed_cartan(coordinates, 0.67);
            let circuit = KakSynthesizer
                .synthesize(&target, &[])
                .unwrap_or_else(|| panic!("rejected near-degenerate case {case}: {coordinates:?}"));
            assert!(count_cx(&circuit) <= 3);
            assert_synthesis_fidelity(&target, &circuit);
            assert!(phase_aligned_residual(&target, &circuit_to_unitary(&circuit)) < 1e-12);
        }
    }

    fn uniform(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        // Open interval, suitable for the Box–Muller logarithm.
        ((*state >> 11) as f64 + 0.5) / ((1_u64 << 53) as f64)
    }

    fn gaussian(state: &mut u64) -> Complex<f64> {
        let radius = (-2.0 * uniform(state).ln()).sqrt();
        Complex::from_polar(radius, std::f64::consts::TAU * uniform(state))
    }

    fn haar_unitary(state: &mut u64) -> DMatrix<Complex<f64>> {
        // Independent complex Gaussian QR (positive real column norms).
        let mut output = DMatrix::<Complex<f64>>::zeros(4, 4);
        for col in 0..4 {
            let mut column = DVector::from_fn(4, |_, _| gaussian(state));
            for _ in 0..2 {
                // Reorthogonalize to keep the oracle well conditioned.
                for previous in 0..col {
                    let basis = output.column(previous);
                    let projection = basis.dotc(&column);
                    column -= basis * projection;
                }
            }
            column /= Complex::new(column.norm(), 0.0);
            output.set_column(col, &column);
        }
        output
    }

    #[test]
    fn seeded_haar_unitaries_have_three_cx_and_small_operator_residual() {
        let mut state = 0x4841_4152_5355_3401;
        for case in 0..256 {
            let target = haar_unitary(&mut state);
            assert!(is_unitary(&target, 1e-14));
            let circuit = KakSynthesizer
                .synthesize(&target, &[])
                .unwrap_or_else(|| panic!("rejected Haar case {case}"));
            assert_eq!(count_cx(&circuit), 3);
            assert_synthesis_fidelity(&target, &circuit);
        }
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
