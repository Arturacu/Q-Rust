//! Unitary circuit simulator for verification.

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateDefinition, Operation};
use nalgebra::DMatrix;
use num_complex::Complex;

type C = Complex<f64>;

pub const MAX_QUBITS: usize = 14;

fn embed_single_qubit(gate_2x2: &DMatrix<C>, target: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let mask = !(1usize << target);
    for col in 0..dim {
        for row in 0..dim {
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

fn embed_2q(gate_4x4: &DMatrix<C>, q0: usize, q1: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
    let mut full = DMatrix::<C>::zeros(dim, dim);
    let mask = !((1usize << q0) | (1usize << q1));
    for col in 0..dim {
        for row in 0..dim {
            if (row & mask) != (col & mask) {
                continue;
            }
            let r0 = (row >> q0) & 1;
            let r1 = (row >> q1) & 1;
            let c0 = (col >> q0) & 1;
            let c1 = (col >> q1) & 1;
            let local_row = r0 + 2 * r1;
            let local_col = c0 + 2 * c1;
            full[(row, col)] = gate_4x4[(local_row, local_col)];
        }
    }
    full
}

fn embed_3q(gate_8x8: &DMatrix<C>, q0: usize, q1: usize, q2: usize, n_qubits: usize) -> DMatrix<C> {
    let dim = 1usize << n_qubits;
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
            let local_row = r0 + 2 * r1 + 4 * r2;
            let local_col = c0 + 2 * c1 + 4 * c2;
            full[(row, col)] = gate_8x8[(local_row, local_col)];
        }
    }
    full
}

pub fn try_circuit_to_unitary(circuit: &Circuit) -> Result<DMatrix<C>> {
    if circuit.num_qubits > MAX_QUBITS {
        return Err(QRustError::Simulation(format!(
            "circuit has {} qubits; simulator supports at most {}",
            circuit.num_qubits, MAX_QUBITS
        )));
    }
    let unrolled = crate::transpiler::decomposition::unroll_custom_gates(circuit);
    let n = unrolled.num_qubits;
    let dim = 1usize << n;
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
                0 => continue, // Barrier pseudo-gate; no unitary action.
                1 => {
                    if qubits.is_empty() {
                        return Err(QRustError::Simulation("1-qubit gate has no target".into()));
                    }
                    embed_single_qubit(&local_u, qubits[0], n)
                }
                2 => {
                    if qubits.len() < 2 {
                        return Err(QRustError::Simulation(
                            "2-qubit gate has insufficient qubits".into(),
                        ));
                    }
                    embed_2q(&local_u, qubits[0], qubits[1], n)
                }
                3 => {
                    if qubits.len() < 3 {
                        return Err(QRustError::Simulation(
                            "3-qubit gate has insufficient qubits".into(),
                        ));
                    }
                    embed_3q(&local_u, qubits[0], qubits[1], qubits[2], n)
                }
                nq => {
                    return Err(QRustError::Simulation(format!(
                        "unsupported gate arity: {}",
                        nq
                    )));
                }
            };
            u = gate_u * u;
        }
    }
    Ok(u)
}

pub fn circuit_to_unitary(circuit: &Circuit) -> DMatrix<C> {
    try_circuit_to_unitary(circuit).expect("circuit_to_unitary failed")
}

pub fn unitary_fidelity(u1: &DMatrix<C>, u2: &DMatrix<C>) -> f64 {
    let d = u1.nrows();
    assert_eq!(u1.shape(), u2.shape(), "unitary_fidelity: shape mismatch");
    let trace: C = (u1.adjoint() * u2).trace();
    trace.norm_sqr() / (d as f64 * d as f64)
}

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
        let mut x_in = 0usize;
        for param_i in 0..n_logical {
            let bit = (col_l >> param_i) & 1;
            x_in |= bit << initial_layout[param_i];
        }
        let output_col = u_routed.column(x_in);

        for row_p in 0..u_routed.nrows() {
            let amp = output_col[row_p];
            if amp.norm() < 1e-9 {
                continue;
            }
            let mut row_l = 0usize;
            let mut ancilla_violation = false;
            for p_q in 0..n_physical {
                let bit = (row_p >> p_q) & 1;
                if let Some(l_q) = final_layout.iter().position(|&x| x == p_q) {
                    row_l |= bit << l_q;
                } else if bit == 1 {
                    ancilla_violation = true;
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
        assert!(unitary_fidelity(&u, &id) > 1.0 - 1e-10);
    }

    #[test]
    fn test_hcx_bell() {
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
        let s = 1.0 / 2.0_f64.sqrt();
        assert!((u[(0, 0)] - c(s, 0.0)).norm() < 1e-10);
        assert!((u[(3, 0)] - c(s, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_max_qubits_guard() {
        let c = Circuit::new(MAX_QUBITS + 1, 0);
        assert!(try_circuit_to_unitary(&c).is_err());
    }
}
// /// Statevector simulator
