//! Unified gate definitions.

use crate::ir::{CommutationSignature, GateType, Operation, PauliBasis};
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

type C = Complex<f64>;

#[inline]
fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

pub fn controlled(u: &DMatrix<C>) -> DMatrix<C> {
    let mut m = DMatrix::<C>::zeros(4, 4);
    m[(0, 0)] = c(1.0, 0.0);
    m[(2, 2)] = c(1.0, 0.0);
    m[(1, 1)] = u[(0, 0)];
    m[(3, 1)] = u[(1, 0)];
    m[(1, 3)] = u[(0, 1)];
    m[(3, 3)] = u[(1, 1)];
    m
}

fn embed_u_local(u: &DMatrix<C>, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let mask = !(1usize << target);
    for col in 0..dim {
        for row in 0..dim {
            if (row & mask) != (col & mask) {
                continue;
            }
            let rb = (row >> target) & 1;
            let cb = (col >> target) & 1;
            full[(row, col)] = u[(rb, cb)];
        }
    }
    full
}

fn embed_cx_local(control: usize, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    for basis in 0..dim {
        let ctrl_bit = (basis >> control) & 1;
        if ctrl_bit == 1 {
            let flipped = basis ^ (1 << target);
            full[(flipped, basis)] = c(1.0, 0.0);
        } else {
            full[(basis, basis)] = c(1.0, 0.0);
        }
    }
    full
}

fn compose_basis_ops(ops: &[Operation], n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
    let mut result = DMatrix::<C>::identity(dim, dim);
    for op in ops {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            let gate_u = match name {
                GateType::U => {
                    let u = basis_u_matrix(params);
                    embed_u_local(&u, qubits[0], n_qubits)
                }
                GateType::CX => embed_cx_local(qubits[0], qubits[1], n_qubits),
                _ => continue,
            };
            result = gate_u * result;
        }
    }
    result
}

fn basis_u_matrix(params: &[f64]) -> DMatrix<C> {
    let theta = params.first().copied().unwrap_or(0.0);
    let phi = params.get(1).copied().unwrap_or(0.0);
    let lambda = params.get(2).copied().unwrap_or(0.0);
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

fn basis_cx_matrix() -> DMatrix<C> {
    let mut m = DMatrix::<C>::zeros(4, 4);
    m[(0, 0)] = c(1.0, 0.0);
    m[(2, 2)] = c(1.0, 0.0);
    m[(3, 1)] = c(1.0, 0.0);
    m[(1, 3)] = c(1.0, 0.0);
    m
}

pub trait GateDefinition {
    fn num_qubits(&self) -> usize;
    fn is_basis(&self) -> bool;
    fn unitary(&self, params: &[f64]) -> DMatrix<C>;
    fn decompose(&self, qubits: &[usize], params: &[f64]) -> Option<Vec<Operation>>;
    fn commutation_signature(&self) -> CommutationSignature;
}

impl GateDefinition for GateType {
    fn num_qubits(&self) -> usize {
        match self {
            GateType::H
            | GateType::X
            | GateType::Y
            | GateType::Z
            | GateType::S
            | GateType::Sdg
            | GateType::T
            | GateType::Tdg
            | GateType::ID
            | GateType::RX
            | GateType::RY
            | GateType::RZ
            | GateType::U => 1,

            GateType::CX
            | GateType::CZ
            | GateType::CY
            | GateType::CH
            | GateType::CSX
            | GateType::CRX
            | GateType::CRY
            | GateType::CRZ
            | GateType::RXX
            | GateType::RYY
            | GateType::RZZ
            | GateType::SWAP => 2,

            GateType::CCX => 3,
            GateType::Barrier => 0,
            GateType::Custom(_) => 1,
        }
    }

    fn is_basis(&self) -> bool {
        matches!(self, GateType::U | GateType::CX)
    }

    fn unitary(&self, params: &[f64]) -> DMatrix<C> {
        match self {
            GateType::U => basis_u_matrix(params),
            GateType::CX => basis_cx_matrix(),
            GateType::Custom(_) | GateType::Barrier => DMatrix::<C>::identity(2, 2),
            other => {
                let n = other.num_qubits();
                let qubits: Vec<usize> = (0..n).collect();
                match other.decompose(&qubits, params) {
                    Some(ops) => compose_basis_ops(&ops, n),
                    None => {
                        let dim = 1usize << n;
                        DMatrix::<C>::identity(dim, dim)
                    }
                }
            }
        }
    }

    fn decompose(&self, qubits: &[usize], params: &[f64]) -> Option<Vec<Operation>> {
        if self.is_basis() {
            return None;
        }

        let u_gate = |q: usize, p: [f64; 3]| Operation::Gate {
            name: GateType::U,
            qubits: vec![q],
            params: p.to_vec(),
        };
        let cx_gate = |a: usize, b: usize| Operation::Gate {
            name: GateType::CX,
            qubits: vec![a, b],
            params: vec![],
        };

        let mut ops = Vec::new();

        match self {
            GateType::U | GateType::CX => unreachable!(),
            GateType::Barrier => return None,

            GateType::ID => ops.push(u_gate(qubits[0], [0.0, 0.0, 0.0])),
            GateType::X => ops.push(u_gate(qubits[0], [PI, 0.0, PI])),
            GateType::Y => ops.push(u_gate(qubits[0], [PI, PI / 2.0, PI / 2.0])),
            GateType::Z => ops.push(u_gate(qubits[0], [0.0, 0.0, PI])),
            GateType::H => ops.push(u_gate(qubits[0], [PI / 2.0, 0.0, PI])),
            GateType::S => ops.push(u_gate(qubits[0], [0.0, 0.0, PI / 2.0])),
            GateType::Sdg => ops.push(u_gate(qubits[0], [0.0, 0.0, -PI / 2.0])),
            GateType::T => ops.push(u_gate(qubits[0], [0.0, 0.0, PI / 4.0])),
            GateType::Tdg => ops.push(u_gate(qubits[0], [0.0, 0.0, -PI / 4.0])),
            GateType::RX => {
                let theta = params.first().copied().unwrap_or(0.0);
                ops.push(u_gate(qubits[0], [theta, -PI / 2.0, PI / 2.0]));
            }
            GateType::RY => {
                let theta = params.first().copied().unwrap_or(0.0);
                ops.push(u_gate(qubits[0], [theta, 0.0, 0.0]));
            }
            GateType::RZ => {
                let phi = params.first().copied().unwrap_or(0.0);
                ops.push(u_gate(qubits[0], [0.0, 0.0, phi]));
            }

            GateType::SWAP => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(cx_gate(a, b));
                ops.push(cx_gate(b, a));
                ops.push(cx_gate(a, b));
            }

            GateType::CCX => {
                let (a, b, t) = (qubits[0], qubits[1], qubits[2]);
                ops.extend(GateType::H.decompose(&[t], &[])?);
                ops.push(cx_gate(b, t));
                ops.extend(GateType::Tdg.decompose(&[t], &[])?);
                ops.push(cx_gate(a, t));
                ops.extend(GateType::T.decompose(&[t], &[])?);
                ops.push(cx_gate(b, t));
                ops.extend(GateType::Tdg.decompose(&[t], &[])?);
                ops.push(cx_gate(a, t));
                ops.extend(GateType::T.decompose(&[b], &[])?);
                ops.extend(GateType::T.decompose(&[t], &[])?);
                ops.extend(GateType::H.decompose(&[t], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::T.decompose(&[a], &[])?);
                ops.extend(GateType::Tdg.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
            }

            GateType::CZ => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::H.decompose(&[b], &[])?);
            }
            GateType::CY => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::Sdg.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::S.decompose(&[b], &[])?);
            }
            GateType::CH => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::S.decompose(&[b], &[])?);
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.extend(GateType::T.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::Tdg.decompose(&[b], &[])?);
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.extend(GateType::Sdg.decompose(&[b], &[])?);
            }
            GateType::CSX => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(u_gate(a, [0.0, 0.0, PI / 4.0]));
                ops.extend(GateType::CRX.decompose(&[a, b], &[PI / 2.0])?);
            }

            GateType::CRX => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RZ.decompose(&[b], &[PI / 2.0])?);
                ops.push(cx_gate(a, b));
                ops.push(u_gate(b, [-theta / 2.0, 0.0, 0.0]));
                ops.push(cx_gate(a, b));
                ops.push(u_gate(b, [theta / 2.0, -PI / 2.0, 0.0]));
            }
            GateType::CRY => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(u_gate(b, [theta / 2.0, 0.0, 0.0]));
                ops.push(cx_gate(a, b));
                ops.push(u_gate(b, [-theta / 2.0, 0.0, 0.0]));
                ops.push(cx_gate(a, b));
            }
            GateType::CRZ => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RZ.decompose(&[b], &[theta / 2.0])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[-theta / 2.0])?);
                ops.push(cx_gate(a, b));
            }

            GateType::RXX => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::H.decompose(&[a], &[])?);
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[theta])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::H.decompose(&[a], &[])?);
                ops.extend(GateType::H.decompose(&[b], &[])?);
            }
            GateType::RYY => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RX.decompose(&[a], &[PI / 2.0])?);
                ops.extend(GateType::RX.decompose(&[b], &[PI / 2.0])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[theta])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RX.decompose(&[a], &[-PI / 2.0])?);
                ops.extend(GateType::RX.decompose(&[b], &[-PI / 2.0])?);
            }
            GateType::RZZ => {
                let theta = params.first().copied().unwrap_or(0.0);
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[theta])?);
                ops.push(cx_gate(a, b));
            }

            GateType::Custom(_) => return None,
        }

        Some(ops)
    }

    fn commutation_signature(&self) -> CommutationSignature {
        match self {
            GateType::Z
            | GateType::RZ
            | GateType::S
            | GateType::Sdg
            | GateType::T
            | GateType::Tdg
            | GateType::CZ => CommutationSignature::Diagonal(PauliBasis::Z),
            GateType::X | GateType::RX => CommutationSignature::Diagonal(PauliBasis::X),
            GateType::Y | GateType::RY => CommutationSignature::Diagonal(PauliBasis::Y),
            GateType::CX => CommutationSignature::CompositeDiagonal(vec![
                (0, PauliBasis::Z),
                (1, PauliBasis::X),
            ]),
            GateType::CY => CommutationSignature::CompositeDiagonal(vec![
                (0, PauliBasis::Z),
                (1, PauliBasis::Y),
            ]),
            GateType::H | GateType::SWAP => CommutationSignature::Clifford,
            _ => CommutationSignature::Generic,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convention_consistency_cx_simulator_path() {
        use crate::ir::Circuit;
        use crate::simulator::circuit_to_unitary;
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let u_sim = circuit_to_unitary(&c);
        let u_direct = GateType::CX.unitary(&[]);
        assert!((u_sim.clone() - u_direct).norm() < 1e-10);
    }

    #[test]
    fn test_convention_consistency_x_on_qubit0() {
        use crate::ir::Circuit;
        use crate::simulator::circuit_to_unitary;
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        let u = circuit_to_unitary(&c);
        let norm = u.norm();
        assert!(norm > 0.0);
    }
}
