//! Unitary circuit simulator.
//!
//! Simulates a quantum circuit by computing the full unitary matrix it
//! represents. This is intended for verification and testing — the
//! exponential memory cost limits it to small circuits (~10 qubits).

use crate::ir::{Circuit, GateDefinition, Operation};
use nalgebra::DMatrix;
use num_complex::Complex;

type C = Complex<f64>;

// gate_matrix_2x2 has been moved to GateDefinition::unitary() in ir/gate_def.rs

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

/// Embeds a 4×4 two-qubit gate acting on `q0`, `q1` into the full 2^n space.
///
/// The local unitary uses qubit ordering [q0, q1] where q0 is the
/// most-significant local bit. For each basis state, we extract the
/// local 2-bit index, look up the 4×4 matrix element, and write it
/// into the full matrix.
fn embed_2q(gate_4x4: &DMatrix<C>, q0: usize, q1: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let mask = !((1usize << q0) | (1usize << q1));

    for col in 0..dim {
        for row in 0..dim {
            // Non-gate qubits must match
            if (row & mask) != (col & mask) {
                continue;
            }
            // Extract local qubit bits
            let r0 = (row >> q0) & 1;
            let r1 = (row >> q1) & 1;
            let c0 = (col >> q0) & 1;
            let c1 = (col >> q1) & 1;
            let local_row = r0 * 2 + r1;
            let local_col = c0 * 2 + c1;
            full[(row, col)] = gate_4x4[(local_row, local_col)];
        }
    }
    full
}

/// Embeds an 8×8 three-qubit gate acting on `q0`, `q1`, `q2` into the full 2^n space.
fn embed_3q(gate_8x8: &DMatrix<C>, q0: usize, q1: usize, q2: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1 << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let mask = !((1usize << q0) | (1usize << q1) | (1usize << q2));

    for col in 0..dim {
        for row in 0..dim {
            if (row & mask) != (col & mask) {
                continue;
            }
            let r0 = (row >> q0) & 1;
            let r1 = (row >> q1) & 1;
            let r2 = (row >> q2) & 1;
            let c0 = (col >> q0) & 1;
            let c1 = (col >> q1) & 1;
            let c2 = (col >> q2) & 1;
            let local_row = r0 * 4 + r1 * 2 + r2;
            let local_col = c0 * 4 + c1 * 2 + c2;
            full[(row, col)] = gate_8x8[(local_row, local_col)];
        }
    }
    full
}

/// Simulates a circuit and returns its full unitary matrix.
///
/// Uses `GateDefinition::unitary()` as the single source of truth for gate
/// matrices, then embeds them into the full 2^n Hilbert space.
///
/// Non-gate operations (Measure, Reset, Barrier) are ignored.
pub fn circuit_to_unitary(circuit: &Circuit) -> DMatrix<C> {
    let unrolled = crate::transpiler::decomposition::unroll_custom_gates(circuit);
    let n = unrolled.num_qubits;
    let dim = 1 << n;
    let mut u = DMatrix::<C>::identity(dim, dim);

    for op in &unrolled.operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            let local_u = name.unitary(params);

            let gate_u = match name.num_qubits() {
                1 => embed_single_qubit(&local_u, qubits[0], n),
                2 => embed_2q(&local_u, qubits[0], qubits[1], n),
                3 => embed_3q(&local_u, qubits[0], qubits[1], qubits[2], n),
                nq => panic!("Unsupported gate arity: {} qubits", nq),
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
    use crate::ir::{Circuit, GateType};

    fn c(re: f64, im: f64) -> C {
        Complex::new(re, im)
    }

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
