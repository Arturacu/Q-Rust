use nalgebra::{ComplexField, DMatrix};
use num_complex::Complex;

fn main() {
    let sqrt2_inv = 1.0 / 2.0_f64.sqrt();
    let one = Complex::new(1.0, 0.0);
    let zero = Complex::new(0.0, 0.0);
    let i = Complex::new(0.0, 1.0);

    let magic_matrix = DMatrix::from_row_slice(
        4,
        4,
        &[
            one, zero, zero, i, zero, i, one, zero, zero, i, -one, zero, one, zero, zero, -i,
        ],
    ) * Complex::new(sqrt2_inv, 0.0);

    let swap = DMatrix::from_row_slice(
        4,
        4,
        &[
            one, zero, zero, zero, zero, zero, one, zero, zero, one, zero, zero, zero, zero, zero,
            one,
        ],
    );

    let u_magic = magic_matrix.adjoint() * &swap * &magic_matrix;
    println!("U_magic:");
    println!("{}", u_magic);

    let m_sys = u_magic.transpose() * &u_magic;
    println!("M_sys:");
    println!("{}", m_sys);

    let schur = m_sys.try_schur(1e-9, 0).unwrap();
    let eigen = schur.unpack().1.diagonal();
    println!("Eigenvalues:");
    println!("{}", eigen);

    let phases: Vec<f64> = eigen.iter().map(|c| c.arg()).collect();
    println!("Phases: {:?}", phases);

    let theta: Vec<f64> = phases.iter().map(|p| p / 2.0).collect();
    println!("Theta (psi): {:?}", theta);
}
