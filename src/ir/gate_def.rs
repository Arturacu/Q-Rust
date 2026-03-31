//! Unified gate definitions: unitary matrices and decomposition rules.
//!
//! This module provides the `GateDefinition` trait and its implementation on
//! `GateType`, consolidating gate semantics into a single source of truth.
//! Both the simulator and the decomposition pass consume this trait, eliminating
//! the risk of drift between the mathematical and logical definitions of a gate.

use crate::ir::{GateType, Operation};
use nalgebra::DMatrix;
use num_complex::Complex;
use std::f64::consts::PI;

type C = Complex<f64>;

fn c(re: f64, im: f64) -> C {
    Complex::new(re, im)
}

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
    ///
    /// The caller (simulator) handles embedding into the full 2^n Hilbert space.
    ///
    /// # Panics
    /// Panics for `Custom` gates (no definition available).
    fn unitary(&self, params: &[f64]) -> DMatrix<C>;

    /// Returns the decomposition of this gate into `{U, CX}` basis operations.
    ///
    /// Returns `None` if the gate is already a basis gate.
    /// The returned operations use **placeholder qubit indices** (0, 1, 2, …)
    /// that must be remapped by the caller to the actual circuit qubits.
    ///
    /// Returns `None` for `Custom` gates (cannot decompose without definition).
    fn decompose(&self, qubits: &[usize], params: &[f64]) -> Option<Vec<Operation>>;
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
            // ── Single-qubit fixed gates ─────────────────────────────────
            GateType::H => {
                let s = 1.0 / 2.0_f64.sqrt();
                DMatrix::from_row_slice(2, 2, &[c(s, 0.0), c(s, 0.0), c(s, 0.0), c(-s, 0.0)])
            }
            GateType::X => {
                DMatrix::from_row_slice(2, 2, &[c(0.0, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)])
            }
            GateType::Y => DMatrix::from_row_slice(
                2,
                2,
                &[c(0.0, 0.0), c(0.0, -1.0), c(0.0, 1.0), c(0.0, 0.0)],
            ),
            GateType::Z => DMatrix::from_row_slice(
                2,
                2,
                &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(-1.0, 0.0)],
            ),
            GateType::S => {
                DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 1.0)])
            }
            GateType::Sdg => DMatrix::from_row_slice(
                2,
                2,
                &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, -1.0)],
            ),
            GateType::T => {
                let phase = Complex::from_polar(1.0, PI / 4.0);
                DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), phase])
            }
            GateType::Tdg => {
                let phase = Complex::from_polar(1.0, -PI / 4.0);
                DMatrix::from_row_slice(2, 2, &[c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), phase])
            }
            GateType::ID => DMatrix::<C>::identity(2, 2),

            // ── Single-qubit parametric gates ────────────────────────────
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

            // ── Two-qubit fixed gates (4×4) ──────────────────────────────
            GateType::CX => {
                // |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|11⟩, |11⟩→|10⟩
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(3, 2)] = c(1.0, 0.0);
                m[(2, 3)] = c(1.0, 0.0);
                m
            }
            GateType::CZ => {
                // diag(1, 1, 1, -1)
                let mut m = DMatrix::<C>::identity(4, 4);
                m[(3, 3)] = c(-1.0, 0.0);
                m
            }
            GateType::CY => {
                // |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ Y
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(3, 2)] = c(0.0, 1.0);
                m[(2, 3)] = c(0.0, -1.0);
                m
            }
            GateType::CH => {
                // |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ H
                let s = 1.0 / 2.0_f64.sqrt();
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(2, 2)] = c(s, 0.0);
                m[(3, 2)] = c(s, 0.0);
                m[(2, 3)] = c(s, 0.0);
                m[(3, 3)] = c(-s, 0.0);
                m
            }
            GateType::CSX => {
                // |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ SX
                // SX = 1/2 * [[1+i, 1-i], [1-i, 1+i]]
                let half = 0.5;
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(2, 2)] = c(half, half);
                m[(3, 2)] = c(half, -half);
                m[(2, 3)] = c(half, -half);
                m[(3, 3)] = c(half, half);
                m
            }
            GateType::SWAP => {
                // Permutation: |01⟩↔|10⟩
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(2, 1)] = c(1.0, 0.0);
                m[(1, 2)] = c(1.0, 0.0);
                m[(3, 3)] = c(1.0, 0.0);
                m
            }

            // ── Two-qubit parametric gates (4×4) ─────────────────────────
            GateType::CRX => {
                // |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ RX(θ)
                let rx = GateType::RX.unitary(params);
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(2, 2)] = rx[(0, 0)];
                m[(3, 2)] = rx[(1, 0)];
                m[(2, 3)] = rx[(0, 1)];
                m[(3, 3)] = rx[(1, 1)];
                m
            }
            GateType::CRY => {
                let ry = GateType::RY.unitary(params);
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(2, 2)] = ry[(0, 0)];
                m[(3, 2)] = ry[(1, 0)];
                m[(2, 3)] = ry[(0, 1)];
                m[(3, 3)] = ry[(1, 1)];
                m
            }
            GateType::CRZ => {
                let rz = GateType::RZ.unitary(params);
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = c(1.0, 0.0);
                m[(1, 1)] = c(1.0, 0.0);
                m[(2, 2)] = rz[(0, 0)];
                m[(3, 2)] = rz[(1, 0)];
                m[(2, 3)] = rz[(0, 1)];
                m[(3, 3)] = rz[(1, 1)];
                m
            }

            // ── Ising interaction gates (4×4) ────────────────────────────
            // exp(-iθ/2 P⊗P) = cos(θ/2) I₄ - i sin(θ/2) P⊗P
            GateType::RXX => {
                let theta = params[0];
                let cos = c((theta / 2.0).cos(), 0.0);
                let isin = c(0.0, -(theta / 2.0).sin());
                // X⊗X = [[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]]
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = cos;
                m[(1, 1)] = cos;
                m[(2, 2)] = cos;
                m[(3, 3)] = cos;
                m[(0, 3)] = isin;
                m[(1, 2)] = isin;
                m[(2, 1)] = isin;
                m[(3, 0)] = isin;
                m
            }
            GateType::RYY => {
                let theta = params[0];
                let cos = c((theta / 2.0).cos(), 0.0);
                let isin = c(0.0, -(theta / 2.0).sin());
                // Y⊗Y = [[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]]
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = cos;
                m[(1, 1)] = cos;
                m[(2, 2)] = cos;
                m[(3, 3)] = cos;
                m[(0, 3)] = -isin; // -i·sin · (-1) = i·sin
                m[(1, 2)] = isin;
                m[(2, 1)] = isin;
                m[(3, 0)] = -isin;
                m
            }
            GateType::RZZ => {
                let theta = params[0];
                // Z⊗Z = diag(1, -1, -1, 1)
                // exp(-iθ/2 Z⊗Z) = diag(e^{-iθ/2}, e^{iθ/2}, e^{iθ/2}, e^{-iθ/2})
                let em = Complex::from_polar(1.0, -theta / 2.0);
                let ep = Complex::from_polar(1.0, theta / 2.0);
                let mut m = DMatrix::<C>::zeros(4, 4);
                m[(0, 0)] = em;
                m[(1, 1)] = ep;
                m[(2, 2)] = ep;
                m[(3, 3)] = em;
                m
            }

            // ── Three-qubit gates (8×8) ──────────────────────────────────
            GateType::CCX => {
                // Toffoli: flip target iff both controls are |1⟩
                // Standard ordering: [ctrl0, ctrl1, target]
                let mut m = DMatrix::<C>::identity(8, 8);
                // Swap |110⟩ ↔ |111⟩ (indices 6 and 7)
                m[(6, 6)] = c(0.0, 0.0);
                m[(7, 7)] = c(0.0, 0.0);
                m[(6, 7)] = c(1.0, 0.0);
                m[(7, 6)] = c(1.0, 0.0);
                m
            }

            GateType::Custom(_) => {
                panic!("Cannot compute unitary for custom gate without definition")
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
                // h tgt
                ops.extend(GateType::H.decompose(&[tgt], &[]).unwrap());
                // cx b,tgt
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![b, tgt],
                    params: vec![],
                });
                // tdg tgt
                ops.extend(GateType::Tdg.decompose(&[tgt], &[]).unwrap());
                // cx a,tgt
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, tgt],
                    params: vec![],
                });
                // t tgt
                ops.extend(GateType::T.decompose(&[tgt], &[]).unwrap());
                // cx b,tgt
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![b, tgt],
                    params: vec![],
                });
                // tdg tgt
                ops.extend(GateType::Tdg.decompose(&[tgt], &[]).unwrap());
                // cx a,tgt
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, tgt],
                    params: vec![],
                });
                // t b
                ops.extend(GateType::T.decompose(&[b], &[]).unwrap());
                // t tgt
                ops.extend(GateType::T.decompose(&[tgt], &[]).unwrap());
                // h tgt
                ops.extend(GateType::H.decompose(&[tgt], &[]).unwrap());
                // cx a,b
                ops.push(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![a, b],
                    params: vec![],
                });
                // t a
                ops.extend(GateType::T.decompose(&[a], &[]).unwrap());
                // tdg b
                ops.extend(GateType::Tdg.decompose(&[b], &[]).unwrap());
                // cx a,b
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
}
