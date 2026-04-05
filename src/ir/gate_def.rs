//! Unified gate definitions: unitary matrices and decomposition rules.
//!
//! This module provides the `GateDefinition` trait and its implementation on
//! `GateType`, consolidating gate semantics into a single source of truth.
//!
//! **Key invariant**: only basis gates (`U`, `CX`) have explicit unitary matrices.
//! All other gates derive their unitary by composing their decomposition through
//! a mini-simulator, making correctness forced by construction.

use crate::ir::{CommutationSignature, GateType, Operation, PauliBasis};
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

type C = Complex<f64>;

fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

// ─── Matrix construction utilities ──────────────────────────────────────────

/// Builds a 4×4 controlled-U matrix from a 2×2 unitary.
///
/// Result = |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ U
///
/// This is a general utility; any controlled gate can be expressed as
/// `controlled(sub_gate.unitary(params))`.
pub fn controlled(u: &DMatrix<C>) -> DMatrix<C> {
    let mut m = DMatrix::<C>::zeros(4, 4);
    m[(0, 0)] = c(1.0, 0.0);
    m[(1, 1)] = c(1.0, 0.0);
    m[(2, 2)] = u[(0, 0)];
    m[(3, 2)] = u[(1, 0)];
    m[(2, 3)] = u[(0, 1)];
    m[(3, 3)] = u[(1, 1)];
    m
}

// ─── Mini-simulator for composing basis operations ──────────────────────────
//
// These functions exist so that `unitary()` can derive the matrix from
// `decompose()` without depending on the full simulator (which would
// create a circular dependency).

/// Embeds a 2×2 unitary acting on `target` into the local n-qubit space.
///
/// Uses the same convention as the full simulator's `embed_2q`/`embed_3q`:
/// qubit 0 is the most-significant bit of the local index.
/// `local_index = bit(q0) * 2^(n-1) + bit(q1) * 2^(n-2) + ...`
fn embed_u_local(u: &DMatrix<C>, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let bit_pos = n_qubits - 1 - target; // qubit 0 → MSB (bit n-1)

    for col in 0..dim {
        for row in 0..dim {
            let mask = !(1usize << bit_pos);
            if (row & mask) != (col & mask) {
                continue;
            }
            let r_bit = (row >> bit_pos) & 1;
            let c_bit = (col >> bit_pos) & 1;
            full[(row, col)] = u[(r_bit, c_bit)];
        }
    }
    full
}

/// Embeds a CX gate on (control, target) into the local n-qubit space.
///
/// Uses the same qubit-to-bit-position convention as `embed_u_local`.
fn embed_cx_local(control: usize, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let ctrl_pos = n_qubits - 1 - control; // qubit 0 → MSB
    let tgt_pos = n_qubits - 1 - target;

    for basis in 0..dim {
        let ctrl_bit = (basis >> ctrl_pos) & 1;
        if ctrl_bit == 1 {
            let flipped = basis ^ (1 << tgt_pos);
            full[(flipped, basis)] = c(1.0, 0.0);
        } else {
            full[(basis, basis)] = c(1.0, 0.0);
        }
    }
    full
}

/// Composes a sequence of basis operations (U and CX only) into a single unitary.
///
/// This is a self-contained mini-simulator used by `GateDefinition::unitary()`
/// to derive the matrix from the decomposition, ensuring consistency by construction.
fn compose_basis_ops(ops: &[Operation], n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
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
                other => panic!(
                    "compose_basis_ops: expected only U/CX, got {:?}. \
                     Decomposition must produce fully-expanded basis ops.",
                    other
                ),
            };
            result = gate_u * result;
        }
    }
    result
}

/// Returns the 2×2 U(θ,φ,λ) matrix — the axiomatic definition of the U gate.
fn basis_u_matrix(params: &[f64]) -> DMatrix<C> {
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

/// Returns the 4×4 CX matrix — the axiomatic definition of the CX gate.
fn basis_cx_matrix() -> DMatrix<C> {
    let mut m = DMatrix::<C>::zeros(4, 4);
    m[(0, 0)] = c(1.0, 0.0);
    m[(1, 1)] = c(1.0, 0.0);
    m[(3, 2)] = c(1.0, 0.0);
    m[(2, 3)] = c(1.0, 0.0);
    m
}

// ─── GateDefinition trait ───────────────────────────────────────────────────

/// Unified gate definition trait.
///
/// Every gate type must be able to report its qubit count, whether it's a basis
/// gate, its local unitary matrix, and its decomposition into the `{U, CX}` basis.
pub trait GateDefinition {
    /// Number of qubits this gate acts on.
    fn num_qubits(&self) -> usize;

    /// Whether this gate is a basis gate (`U` or `CX`).
    fn is_basis(&self) -> bool;

    /// Returns the local unitary matrix for this gate.
    ///
    /// - Single-qubit gates return a 2×2 matrix.
    /// - Two-qubit gates return a 4×4 matrix (in computational basis order).
    /// - Three-qubit gates return an 8×8 matrix.
    ///
    /// For basis gates (`U`, `CX`), the matrix is defined axiomatically.
    /// For all other gates, the matrix is **derived from the decomposition**,
    /// guaranteeing consistency by construction.
    ///
    /// # Panics
    /// Panics for `Custom` gates (no definition available).
    fn unitary(&self, params: &[f64]) -> DMatrix<C>;

    /// Returns the decomposition of this gate into `{U, CX}` basis operations.
    ///
    /// Returns `None` if the gate is already a basis gate.
    /// The returned operations use the actual qubit indices passed in via the
    /// `qubits` argument — the caller provides the physical qubit mapping.
    ///
    /// Returns `None` for `Custom` gates (cannot decompose without definition).
    fn decompose(&self, qubits: &[usize], params: &[f64]) -> Option<Vec<Operation>>;

    /// Derives the strict algebraic commutation signature for the given gate topology.
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

            GateType::Custom(_) => 1, // Default; custom gates need context
        }
    }

    fn is_basis(&self) -> bool {
        matches!(self, GateType::U | GateType::CX)
    }

    fn unitary(&self, params: &[f64]) -> DMatrix<C> {
        match self {
            // ── Axiomatic basis gates ────────────────────────────────────
            GateType::U => basis_u_matrix(params),
            GateType::CX => basis_cx_matrix(),

            // ── All other gates: derived from decomposition ─────────────
            GateType::Custom(_) => {
                panic!("Cannot compute unitary for custom gate without definition")
            }
            other => {
                let n = other.num_qubits();
                let qubits: Vec<usize> = (0..n).collect();
                let ops = other
                    .decompose(&qubits, params)
                    .expect("Non-basis, non-custom gate must have a decomposition");
                compose_basis_ops(&ops, n)
            }
        }
    }

    fn decompose(&self, qubits: &[usize], params: &[f64]) -> Option<Vec<Operation>> {
        if self.is_basis() {
            return None; // Already a basis gate
        }

        let mut ops = Vec::new();

        match self {
            GateType::U | GateType::CX => unreachable!(), // handled above

            // ── Single-qubit → U(θ,φ,λ) ─────────────────────────────────
            GateType::ID => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, 0.0],
                });
            }
            GateType::X => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![PI, 0.0, PI],
                });
            }
            GateType::Y => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![PI, PI / 2.0, PI / 2.0],
                });
            }
            GateType::Z => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, PI],
                });
            }
            GateType::H => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![PI / 2.0, 0.0, PI],
                });
            }
            GateType::S => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, PI / 2.0],
                });
            }
            GateType::Sdg => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, -PI / 2.0],
                });
            }
            GateType::T => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, PI / 4.0],
                });
            }
            GateType::Tdg => {
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, -PI / 4.0],
                });
            }
            GateType::RX => {
                let theta = params[0];
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![theta, -PI / 2.0, PI / 2.0],
                });
            }
            GateType::RY => {
                let theta = params[0];
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![theta, 0.0, 0.0],
                });
            }
            GateType::RZ => {
                let phi = params[0];
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: qubits.to_vec(),
                    params: vec![0.0, 0.0, phi],
                });
            }

            // ── SWAP → 3 CX ─────────────────────────────────────────────
            GateType::SWAP => {
                let a = qubits[0];
                let b = qubits[1];
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![b, a],
                    params: vec![],
                });
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
            }

            // ── CCX (Toffoli) ────────────────────────────────────────────
            GateType::CCX => {
                let a = qubits[0];
                let b = qubits[1];
                let tgt = qubits[2];
                ops.extend(GateType::H.decompose(&[tgt], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![b, tgt],
                    params: vec![],
                });
                ops.extend(GateType::Tdg.decompose(&[tgt], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, tgt],
                    params: vec![],
                });
                ops.extend(GateType::T.decompose(&[tgt], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![b, tgt],
                    params: vec![],
                });
                ops.extend(GateType::Tdg.decompose(&[tgt], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, tgt],
                    params: vec![],
                });
                ops.extend(GateType::T.decompose(&[b], &[]).unwrap());
                ops.extend(GateType::T.decompose(&[tgt], &[]).unwrap());
                ops.extend(GateType::H.decompose(&[tgt], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::T.decompose(&[a], &[]).unwrap());
                ops.extend(GateType::Tdg.decompose(&[b], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
            }

            // ── Controlled gates ─────────────────────────────────────────
            GateType::CZ => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
            }
            GateType::CY => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::Sdg.decompose(&[b], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::S.decompose(&[b], &[]).unwrap());
            }
            GateType::CH => {
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::S.decompose(&[b], &[]).unwrap());
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
                ops.extend(GateType::T.decompose(&[b], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::Tdg.decompose(&[b], &[]).unwrap());
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
                ops.extend(GateType::Sdg.decompose(&[b], &[]).unwrap());
            }
            GateType::CSX => {
                let (a, b) = (qubits[0], qubits[1]);
                // Phase on control
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![a],
                    params: vec![0.0, 0.0, PI / 4.0],
                });
                // CRX(pi/2)
                ops.extend(GateType::CRX.decompose(&[a, b], &[PI / 2.0]).unwrap());
            }

            // ── Controlled rotations ─────────────────────────────────────
            GateType::CRX => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RZ.decompose(&[b], &[PI / 2.0]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![b],
                    params: vec![-theta / 2.0, 0.0, 0.0],
                });
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![b],
                    params: vec![theta / 2.0, -PI / 2.0, 0.0],
                });
            }
            GateType::CRY => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![b],
                    params: vec![theta / 2.0, 0.0, 0.0],
                });
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.push(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![b],
                    params: vec![-theta / 2.0, 0.0, 0.0],
                });
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
            }
            GateType::CRZ => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RZ.decompose(&[b], &[theta / 2.0]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::RZ.decompose(&[b], &[-theta / 2.0]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
            }

            // ── Ising interaction gates ──────────────────────────────────
            GateType::RXX => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::H.decompose(&[a], &[]).unwrap());
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::RZ.decompose(&[b], &[theta]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::H.decompose(&[a], &[]).unwrap());
                ops.extend(GateType::H.decompose(&[b], &[]).unwrap());
            }
            GateType::RYY => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.extend(GateType::RX.decompose(&[a], &[PI / 2.0]).unwrap());
                ops.extend(GateType::RX.decompose(&[b], &[PI / 2.0]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::RZ.decompose(&[b], &[theta]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::RX.decompose(&[a], &[-PI / 2.0]).unwrap());
                ops.extend(GateType::RX.decompose(&[b], &[-PI / 2.0]).unwrap());
            }
            GateType::RZZ => {
                let theta = params[0];
                let (a, b) = (qubits[0], qubits[1]);
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                ops.extend(GateType::RZ.decompose(&[b], &[theta]).unwrap());
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
            }

            GateType::Custom(_) => return None,
        }

        Some(ops)
    }

    fn commutation_signature(&self) -> CommutationSignature {
        match self {
            // Pure Z-Diagonal
            GateType::Z
            | GateType::RZ
            | GateType::S
            | GateType::Sdg
            | GateType::T
            | GateType::Tdg
            | GateType::CZ => CommutationSignature::Diagonal(PauliBasis::Z),
            // Pure X-Diagonal
            GateType::X | GateType::RX => CommutationSignature::Diagonal(PauliBasis::X),
            // Pure Y-Diagonal
            GateType::Y | GateType::RY => CommutationSignature::Diagonal(PauliBasis::Y),
            // Composite Diagonals (e.g. CX is Z on control, X on target)
            GateType::CX => CommutationSignature::CompositeDiagonal(vec![
                (0, PauliBasis::Z),
                (1, PauliBasis::X),
            ]),
            GateType::CY => CommutationSignature::CompositeDiagonal(vec![
                (0, PauliBasis::Z),
                (1, PauliBasis::Y),
            ]),
            // Pure Cliffords that are not diagonal
            GateType::H | GateType::SWAP => CommutationSignature::Clifford,
            // Everything else acts completely generalized
            _ => CommutationSignature::Generic,
        }
    }
}
