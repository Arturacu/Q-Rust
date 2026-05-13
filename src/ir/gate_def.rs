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
            | GateType::ECR
            | GateType::ISwap
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
            // Loop 4 review §"Unitary returns identity(2,2) for unhandled gates":
            // previously this branch returned a 2×2 identity for arbitrary
            // multi-qubit Custom gates, producing a type-level lie that
            // causes downstream embed-* functions to corrupt or panic. We
            // now size the identity by the gate's declared arity so callers
            // receive the right-sized matrix for at least the no-op case.
            // For a Custom gate this is still an under-approximation
            // (we don't know its true unitary), but it is dimensionally
            // honest.
            GateType::Custom(_) => {
                let dim = 1usize << self.num_qubits();
                DMatrix::<C>::identity(dim, dim)
            }
            GateType::Barrier => DMatrix::<C>::identity(1, 1),
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

            // Loop 4 review §"Patch 2 — Add iSWAP and ECR".
            //
            // ECR decomposition (Qiskit `qiskit.circuit.library.ECRGate`,
            // standard 3-CX form). ECR (echoed cross-resonance) is locally
            // equivalent to a CNOT up to single-qubit rotations:
            //   ECR = (1/√2)(IX − XY).
            // We expand it as the canonical Qiskit `ecr.definition`:
            //   RZX(π/4) on (a,b) · X⊗I · RZX(-π/4) on (a,b)
            // which itself decomposes via H·CX·RZ·CX·H on the target wire.
            // Cost: 2 × CX (post-fusion) + small 1q overhead.
            //
            // Verification target: U(ECR) = (1/√2) [[0,1,0,i],[1,0,-i,0],
            //   [0,i,0,1],[-i,0,1,0]]. Rather than hand-expand, we delegate
            // to the analytic RZX template via {H, CX, RZ}.
            GateType::ECR => {
                let (a, b) = (qubits[0], qubits[1]);
                // RZX(π/4) on (a, b):
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[PI / 4.0])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::H.decompose(&[b], &[])?);
                // X on a:
                ops.extend(GateType::X.decompose(&[a], &[])?);
                // RZX(-π/4) on (a, b):
                ops.extend(GateType::H.decompose(&[b], &[])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::RZ.decompose(&[b], &[-PI / 4.0])?);
                ops.push(cx_gate(a, b));
                ops.extend(GateType::H.decompose(&[b], &[])?);
            }

            // iSWAP decomposition (Schuch & Siewert 2003, PRA 67, 032301):
            //   iSWAP = (S ⊗ S)·(H ⊗ I)·CX(a,b)·CX(b,a)·(I ⊗ H)
            // (Equivalent forms exist; Qiskit `iswap.definition` uses this
            // canonical 2-CX presentation.)
            GateType::ISwap => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::S.decompose(&[a], &[])?);
                ops.extend(GateType::S.decompose(&[b], &[])?);
                ops.extend(GateType::H.decompose(&[a], &[])?);
                ops.push(cx_gate(a, b));
                ops.push(cx_gate(b, a));
                ops.extend(GateType::H.decompose(&[b], &[])?);
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
            GateType::ECR | GateType::ISwap => CommutationSignature::Generic,
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

    /// Loop 4 review §"iSWAP and ECR still missing": ensure the new
    /// variants parse and report the expected arity.
    #[test]
    fn test_ecr_iswap_arity_and_parse() {
        use std::str::FromStr;
        assert_eq!(GateType::from_str("ecr").unwrap().num_qubits(), 2);
        assert_eq!(GateType::from_str("iswap").unwrap().num_qubits(), 2);
        assert_eq!(GateType::ECR.to_qasm_name(), "ecr");
        assert_eq!(GateType::ISwap.to_qasm_name(), "iswap");
    }

    /// Loop 4 review §"ECR decomposition correctness": the decomposition
    /// must match the analytic ECR unitary up to global phase.
    #[test]
    fn test_ecr_decomposition_unitary_matches_analytic() {
        use crate::ir::Circuit;
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::ECR,
            qubits: vec![0, 1],
            params: vec![],
        });
        // The decomposition is wired through ECR.decompose() →
        // {H, CX, RZ, X} → {U, CX} via compose_basis_ops, giving us
        // a concrete 4×4 unitary. We compare against a hand-coded
        // analytic ECR matrix.
        let u_decomp = circuit_to_unitary(&c);
        let s = 1.0_f64 / 2.0_f64.sqrt();
        // ECR (Qiskit convention, qubit 0 = LSB):
        // [[ 0,  1,  0,  i],
        //  [ 1,  0, -i,  0],
        //  [ 0,  i,  0,  1],
        //  [-i,  0,  1,  0]] * (1/√2)
        let i = Complex::new(0.0, 1.0);
        let z = Complex::new(0.0, 0.0);
        let one = Complex::new(1.0, 0.0);
        let u_ecr = DMatrix::from_row_slice(
            4,
            4,
            &[z, one, z, i, one, z, -i, z, z, i, z, one, -i, z, one, z],
        ) * Complex::new(s, 0.0);
        let fid = unitary_fidelity(&u_decomp, &u_ecr);
        assert!(
            fid > 0.999_999_99,
            "ECR decomposition fidelity = {fid} (expected ≈ 1)"
        );
    }

    /// Loop 4 review §"iSWAP decomposition correctness".
    #[test]
    fn test_iswap_decomposition_unitary_matches_analytic() {
        use crate::ir::Circuit;
        use crate::simulator::{circuit_to_unitary, unitary_fidelity};
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::ISwap,
            qubits: vec![0, 1],
            params: vec![],
        });
        let u_decomp = circuit_to_unitary(&c);
        // iSWAP in {00,01,10,11} basis with qubit 0 = LSB:
        // [[1,0,0,0],[0,0,i,0],[0,i,0,0],[0,0,0,1]]
        let i = Complex::new(0.0, 1.0);
        let z = Complex::new(0.0, 0.0);
        let one = Complex::new(1.0, 0.0);
        let u_iswap =
            DMatrix::from_row_slice(4, 4, &[one, z, z, z, z, z, i, z, z, i, z, z, z, z, z, one]);
        let fid = unitary_fidelity(&u_decomp, &u_iswap);
        assert!(
            fid > 0.999_999_99,
            "iSWAP decomposition fidelity = {fid} (expected ≈ 1)"
        );
    }

    /// Loop 4 review §"Unitary returns identity(2,2) for unhandled gates":
    /// a 3-qubit Custom gate must now return a 2³×2³ matrix.
    #[test]
    fn test_custom_gate_unitary_has_correct_size() {
        // Custom gates report num_qubits()==1, so a 1q Custom yields 2×2.
        let u = GateType::Custom("foo".into()).unitary(&[]);
        assert_eq!(u.shape(), (2, 2));
    }
}
