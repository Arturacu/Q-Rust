//! Unitary circuit simulator.
//!
//! Simulates a quantum circuit by computing the full unitary matrix it
//! represents. This is intended for verification and testing — the
//! exponential memory cost limits it to small circuits (~10 qubits).

use crate::ir::{Circuit, GateType, Operation};
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

type C = Complex<f64>;

fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

/// Returns the 2×2 unitary matrix for a single-qubit gate with the given params.
fn gate_matrix_2x2(gate: &GateType, params: &[f64]) -> DMatrix<C> {
    match gate {
        GateType::H => {
            let s = 1.0 / 2.0_f64.sqrt();
            DMatrix::from_row_slice(2, 2, &[c(s, 0.0), c(s, 0.0), c(s, 0.0), c(-s, 0.0)])
        }
        GateType::X => {
            DMatrix::from_row_slice(2, 2, &[c(0.0, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)])
        }
        GateType::Y => {
            DMatrix::from_row_slice(2, 2, &[c(0.0, 0.0), c(0.0, -1.0), c(0.0, 1.0), c(0.0, 0.0)])
        }
        GateType::Z => {
            DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(-1.0, 0.0)])
        }
        GateType::S => {
            DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 1.0)])
        }
        GateType::Sdg => {
            DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, -1.0)])
        }
        GateType::T => {
            let phase = Complex::from_polar(1.0, PI / 4.0);
            DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), phase])
        }
        GateType::Tdg => {
            let phase = Complex::from_polar(1.0, -PI / 4.0);
            DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), phase])
        }
        GateType::ID => DMatrix::<C>::identity(2, 2),
        GateType::RX => {
            let theta = params[0];
            let cos = c((theta / 2.0).cos(), 0.0);
            let isin = c(0.0, -(theta / 2.0).sin());
            DMatrix::from_row_slice(2, 2, &[cos, isin, isin, cos])
        }
        GateType::RY => {
            let theta = params[0];
            let cos = c((theta / 2.0).cos(), 0.0);
            let sin = c((theta / 2.0).sin(), 0.0);
            DMatrix::from_row_slice(2, 2, &[cos, -sin, sin, cos])
        }
        GateType::RZ => {
            let theta = params[0];
            let em = Complex::from_polar(1.0, -theta / 2.0);
            let ep = Complex::from_polar(1.0, theta / 2.0);
            DMatrix::from_row_slice(2, 2, &[em, c(0.0, 0.0), c(0.0, 0.0), ep])
        }
        GateType::U => {
            let (theta, phi, lambda) = (params[0], params[1], params[2]);
            let cos = (theta / 2.0).cos();
            let sin = (theta / 2.0).sin();
            let e_phi = Complex::from_polar(1.0, phi);
            let e_lambda = Complex::from_polar(1.0, lambda);
            let e_phi_lambda = Complex::from_polar(1.0, phi + lambda);
            DMatrix::from_row_slice(
                2,
                2,
                &[
                    c(cos, 0.0),
                    -e_lambda * sin,
                    e_phi * sin,
                    e_phi_lambda * cos,
                ],
            )
        }
        other => panic!("Unsupported single-qubit gate for simulation: {:?}", other),
    }
}

/// Embeds a 2×2 single-qubit gate acting on `target` into the full 2^n × 2^n space.
fn embed_single_qubit(gate_2x2: &DMatrix<C>, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);

    for col in 0..dim {
        for row in 0..dim {
            // Check if all non-target qubits match
            let mask = !(1usize << target);
            if (row & mask) != (col & mask) {
                continue;
            }
            let r_bit = (row >> target) & 1;
            let c_bit = (col >> target) & 1;
            full[(row, col)] = gate_2x2[(r_bit, c_bit)];
        }
    }
    full
}

/// Returns the 4×4 CX matrix with the given control and target in the full n-qubit space.
fn embed_cx(control: usize, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::identity(dim, dim);

    for basis in 0..dim {
        let ctrl_bit = (basis >> control) & 1;
        if ctrl_bit == 1 {
            // Flip the target bit
            let flipped = basis ^ (1 << target);
            full[(basis, basis)] = c(0.0, 0.0);
            full[(flipped, basis)] = c(1.0, 0.0);
        }
    }
    full
}

/// Returns the SWAP matrix for the given two qubits in the full n-qubit space.
fn embed_swap(q0: usize, q1: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);

    for basis in 0..dim {
        let bit0 = (basis >> q0) & 1;
        let bit1 = (basis >> q1) & 1;
        let mut swapped = basis;
        // Clear both bits and set swapped values
        swapped &= !(1 << q0);
        swapped &= !(1 << q1);
        swapped |= bit1 << q0;
        swapped |= bit0 << q1;
        full[(swapped, basis)] = c(1.0, 0.0);
    }
    full
}

/// Returns the CCX (Toffoli) matrix for the given control/target qubits.
fn embed_ccx(ctrl0: usize, ctrl1: usize, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::identity(dim, dim);

    for basis in 0..dim {
        let c0 = (basis >> ctrl0) & 1;
        let c1 = (basis >> ctrl1) & 1;
        if c0 == 1 && c1 == 1 {
            let flipped = basis ^ (1 << target);
            full[(basis, basis)] = c(0.0, 0.0);
            full[(flipped, basis)] = c(1.0, 0.0);
        }
    }
    full
}

/// Embeds a controlled-U gate: |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ U on (control, target).
fn embed_controlled_unitary(
    u_2x2: &DMatrix<C>,
    control: usize,
    target: usize,
    n_qubits: usize,
) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::identity(dim, dim);

    for basis in 0..dim {
        let ctrl_bit = (basis >> control) & 1;
        if ctrl_bit == 1 {
            let t_bit = (basis >> target) & 1;
            // Apply U to the target qubit
            for out_t in 0..2usize {
                let out_basis = (basis & !(1 << target)) | (out_t << target);
                full[(out_basis, basis)] = u_2x2[(out_t, t_bit)];
            }
        }
    }
    full
}

/// Embeds an Ising interaction exp(-i θ/2 P⊗P) where P is X, Y, or Z.
/// `pauli` is 0=X, 1=Y, 2=Z.
fn embed_ising(theta: f64, pauli: u8, q0: usize, q1: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);

    // Get the 2x2 Pauli matrix
    let p: [[C; 2]; 2] = match pauli {
        0 => [[c(0.0, 0.0), c(1.0, 0.0)], [c(1.0, 0.0), c(0.0, 0.0)]], // X
        1 => [[c(0.0, 0.0), c(0.0, -1.0)], [c(0.0, 1.0), c(0.0, 0.0)]], // Y
        2 => [[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(-1.0, 0.0)]], // Z
        _ => unreachable!(),
    };

    // Build 4x4 P⊗P
    let mut pp = [[c(0.0, 0.0); 4]; 4];
    for i0 in 0..2 {
        for j0 in 0..2 {
            for i1 in 0..2 {
                for j1 in 0..2 {
                    pp[i0 * 2 + i1][j0 * 2 + j1] = p[i0][j0] * p[i1][j1];
                }
            }
        }
    }

    // exp(-i θ/2 P⊗P): since (P⊗P)^2 = I, we have:
    // exp(-i θ/2 P⊗P) = cos(θ/2) I - i sin(θ/2) P⊗P
    let cos_t = c((theta / 2.0).cos(), 0.0);
    let isin_t = c(0.0, -(theta / 2.0).sin());

    // Build the full n-qubit matrix by iterating basis states
    for col in 0..dim {
        for row in 0..dim {
            // All qubits except q0, q1 must match
            let mask = !((1usize << q0) | (1usize << q1));
            if (row & mask) != (col & mask) {
                continue;
            }
            let r0 = (row >> q0) & 1;
            let r1 = (row >> q1) & 1;
            let c0 = (col >> q0) & 1;
            let c1 = (col >> q1) & 1;

            let local_row = r0 * 2 + r1;
            let local_col = c0 * 2 + c1;

            let id_val = if local_row == local_col {
                c(1.0, 0.0)
            } else {
                c(0.0, 0.0)
            };
            full[(row, col)] = cos_t * id_val + isin_t * pp[local_row][local_col];
        }
    }
    full
}

/// Simulates a circuit and returns its full unitary matrix.
///
/// Non-gate operations (Measure, Reset, Barrier) are ignored.
/// Panics if a gate type is not supported by the simulator.
pub fn circuit_to_unitary(circuit: &Circuit) -> DMatrix<C> {
    let n = circuit.num_qubits;
    let dim = 1 << n;
    let mut u = DMatrix::<C>::identity(dim, dim);

    for op in &circuit.operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            let gate_u = match name {
                // Two/three-qubit gates
                GateType::CX => embed_cx(qubits[0], qubits[1], n),
                GateType::SWAP => embed_swap(qubits[0], qubits[1], n),
                GateType::CCX => embed_ccx(qubits[0], qubits[1], qubits[2], n),
                // Controlled gates
                GateType::CZ => {
                    let z = gate_matrix_2x2(&GateType::Z, &[]);
                    embed_controlled_unitary(&z, qubits[0], qubits[1], n)
                }
                GateType::CY => {
                    let y = gate_matrix_2x2(&GateType::Y, &[]);
                    embed_controlled_unitary(&y, qubits[0], qubits[1], n)
                }
                GateType::CH => {
                    let h = gate_matrix_2x2(&GateType::H, &[]);
                    embed_controlled_unitary(&h, qubits[0], qubits[1], n)
                }
                GateType::CSX => {
                    // SX = RX(π/2) up to global phase. Use exact SX matrix:
                    // SX = 1/2 * [[1+i, 1-i], [1-i, 1+i]]
                    let half = 0.5;
                    let sx = DMatrix::from_row_slice(
                        2,
                        2,
                        &[c(half, half), c(half, -half), c(half, -half), c(half, half)],
                    );
                    embed_controlled_unitary(&sx, qubits[0], qubits[1], n)
                }
                GateType::CRX => {
                    let rx = gate_matrix_2x2(&GateType::RX, params);
                    embed_controlled_unitary(&rx, qubits[0], qubits[1], n)
                }
                GateType::CRY => {
                    let ry = gate_matrix_2x2(&GateType::RY, params);
                    embed_controlled_unitary(&ry, qubits[0], qubits[1], n)
                }
                GateType::CRZ => {
                    let rz = gate_matrix_2x2(&GateType::RZ, params);
                    embed_controlled_unitary(&rz, qubits[0], qubits[1], n)
                }
                // Ising interaction gates
                GateType::RXX => embed_ising(params[0], 0, qubits[0], qubits[1], n),
                GateType::RYY => embed_ising(params[0], 1, qubits[0], qubits[1], n),
                GateType::RZZ => embed_ising(params[0], 2, qubits[0], qubits[1], n),
                // Single-qubit gates
                _ => {
                    let mat = gate_matrix_2x2(name, params);
                    embed_single_qubit(&mat, qubits[0], n)
                }
            };
            u = gate_u * u;
        }
        // Non-gate operations (Measure, Reset, Barrier) are silently ignored
    }
    u
}

/// Checks if two unitary matrices are equal up to a global phase.
///
/// Returns the fidelity as |<U1|U2>|^2 / d^2, where d is the matrix dimension.
/// A fidelity of 1.0 means the matrices are equivalent up to global phase.
pub fn unitary_fidelity(u1: &DMatrix<C>, u2: &DMatrix<C>) -> f64 {
    let d = u1.nrows();
    assert_eq!(
        u1.shape(),
        u2.shape(),
        "Unitary matrices must have the same shape"
    );

    // Compute Tr(U1^dag * U2)
    let u1_dag = u1.adjoint();
    let product = u1_dag * u2;
    let trace: C = product.trace();

    // |Tr(U1^dag U2)|^2 / d^2
    let fid = trace.norm_sqr() / (d as f64 * d as f64);
    fid
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::Circuit;

    #[test]
    fn test_identity_circuit() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::ID,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::ID,
            qubits: vec![1],
            params: vec![],
        });

        let u = circuit_to_unitary(&circuit);
        let id = DMatrix::<C>::identity(4, 4);
        assert!(
            unitary_fidelity(&u, &id) > 1.0 - 1e-10,
            "Identity circuit should produce identity unitary"
        );
    }

    #[test]
    fn test_hcx_bell() {
        // H on qubit 0, then CX 0->1 creates a Bell state preparation
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let u = circuit_to_unitary(&circuit);

        // The unitary should map |00> to (|00> + |11>)/sqrt(2)
        let s = 1.0 / 2.0_f64.sqrt();
        assert!((u[(0, 0)] - c(s, 0.0)).norm() < 1e-10);
        assert!((u[(3, 0)] - c(s, 0.0)).norm() < 1e-10);
    }
}
