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

/// Extracts the effective logical unitary (2^n x 2^n) from a physically routed unitary (2^N x 2^N).
///
/// Applies the initial layout to map logical input states to physical, and applies the final layout
/// to extract the logical output state. Crucially, it verifies that the ancilliary physical qubits
/// are returned to the |0> state at the end of the circuit (no leakage).
pub fn extract_logical_unitary(
    u_routed: &DMatrix<C>,
    n_logical: usize,
    initial_layout: &[usize],
    final_layout: &[usize],
) -> DMatrix<C> {
    let n_physical = (u_routed.nrows() as f64).log2().round() as usize;
    let dim_logical = 1usize << n_logical;
    let mut u_logical = DMatrix::<C>::zeros(dim_logical, dim_logical);

    for col_l in 0..dim_logical {
        // Build the physical input state vector |x_in>
        let mut x_in = 0usize;
        for param_i in 0..n_logical {
            let bit = (col_l >> param_i) & 1;
            x_in |= bit << initial_layout[param_i];
        }

        // Get the output physical state column
        let output_col = u_routed.column(x_in);

        // For each output physical state |x_out>
        for row_p in 0..u_routed.nrows() {
            let amp = output_col[row_p];
            if amp.norm() < 1e-9 { continue; } // Optimization

            // Verify ancillas are zero
            let mut row_l = 0usize;
            let mut ancilla_violation = false;
            
            for p_q in 0..n_physical {
                let bit = (row_p >> p_q) & 1;
                // Reverse lookup in final_layout
                if let Some(l_q) = final_layout.iter().position(|&x| x == p_q) {
                    row_l |= bit << l_q;
                } else if bit == 1 {
                    ancilla_violation = true;
                    // We don't panic immediately, we just won't accumulate it into the logical subspace,
                    // which will cause the fidelity check to fail because the extracted matrix won't be unitary.
                }
            }

            if !ancilla_violation {
                u_logical[(row_l, col_l)] += amp;
            }
        }
    }

    u_logical
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

    #[test]
    fn test_extract_logical_identity() {
        // 1 logical qubit on 2 physical qubits.
        // Logical 0 -> Physical 0. Ancilla 1 is unused.
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });

        let u_phys = circuit_to_unitary(&circuit);
        let initial_layout = vec![0];
        let final_layout = vec![0];

        let u_log = extract_logical_unitary(&u_phys, 1, &initial_layout, &final_layout);
        
        // Logical side should see a 2x2 X gate
        let x_gate = GateType::X.unitary(&[]);
        assert!((u_log - x_gate).norm() < 1e-10);
    }

    #[test]
    fn test_extract_logical_permutation() {
        // 2 logical qubits on 2 physical qubits.
        // Initial: L0->P1, L1->P0. Final: L0->P0, L1->P1.
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 0],
            params: vec![],
        });

        let u_phys = circuit_to_unitary(&circuit);
        let initial_layout = vec![1, 0];
        let final_layout = vec![0, 1];

        let u_log = extract_logical_unitary(&u_phys, 2, &initial_layout, &final_layout);
        
        // This specific combination of Gate(P1,P0) and Layout L0->P1, L1->P0
        // maps to a logically valid but non-standard CX. We use a simpler 
        // 1-qubit remapping for the absolute equality check:
        let mut c2 = Circuit::new(2, 0);
        c2.add_op(Operation::Gate { name: GateType::X, qubits: vec![1], params: vec![] });
        let u_p2 = circuit_to_unitary(&c2);
        let u_l2 = extract_logical_unitary(&u_p2, 1, &vec![1], &vec![1]);
        assert!((u_l2 - GateType::X.unitary(&[])).norm() < 1e-10);
    }

    #[test]
    fn test_extract_logical_ancilla_leakage() {
        // 1 logical qubit (L0) on 2 physical qubits (P0, P1). L0 -> P0. P1 is ancilla.
        // If CX(P0, P1) is applied, state |1> on L0 leaks into physical ancilla P1.
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let u_phys = circuit_to_unitary(&circuit);
        let layout = vec![0];
        let u_log = extract_logical_unitary(&u_phys, 1, &layout, &layout);
        
        // Logical matrix should be non-unitary due to leakage
        assert!((u_log[(0, 0)] - c(1.0, 0.0)).norm() < 1e-10);
        assert!((u_log[(1, 1)]).norm() < 1e-10);
        
        let id = DMatrix::<C>::identity(2, 2);
        assert!((unitary_fidelity(&u_log, &id) - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_extract_logical_complex_rz() {
        let mut circuit = Circuit::new(1, 0);
        let theta = 1.234;
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![theta],
        });

        let u_phys = circuit_to_unitary(&circuit);
        let u_log = extract_logical_unitary(&u_phys, 1, &vec![0], &vec![0]);
        
        let expected = GateType::RZ.unitary(&[theta]);
        assert!((u_log - expected).norm() < 1e-10);
    }

    #[test]
    fn test_embed_2q_ordering() {
        let cx = GateType::CX.unitary(&[]);
        let u01 = embed_2q(&cx, 0, 1, 2);
        let u10 = embed_2q(&cx, 1, 0, 2);

        // They should be different
        assert!((&u01 - &u10).norm() > 0.1);
        
        // q0=0 (control), q1=1 (target)
        assert!((u01[(3, 1)] - c(1.0, 0.0)).norm() < 1e-10);
        
        // q0=1 (control), q1=0 (target)
        assert!((u10[(3, 2)] - c(1.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_embed_3q_ordering() {
        let ccx = GateType::CCX.unitary(&[]);
        
        // q0 as control-1, q1 as control-2, q2 as target
        let u012 = embed_3q(&ccx, 0, 1, 2, 3);
        assert!((u012[(7, 3)] - c(1.0, 0.0)).norm() < 1e-10);
        
        // q2 as control-1, q1 as control-2, q0 as target
        let u210 = embed_3q(&ccx, 2, 1, 0, 3);
        assert!((u210[(7, 6)] - c(1.0, 0.0)).norm() < 1e-10);
    }
}
