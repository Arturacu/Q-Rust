//! Unitary circuit simulator for verification.

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateDefinition, Operation};
use crate::SHIFT_GUARD;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

type C = Complex<f64>;

pub const MAX_QUBITS: usize = 14;

/// Maximum qubit count for state-vector evolution.
///
/// State-vector evolution is `O(g · 2^n)` per gate (vs `O(g · 4^n)` for the
/// full unitary materialization in [`circuit_to_unitary`]). We allow up to
/// 24 qubits (~ 256 MiB for the state vector) — well within reach of a
/// 2024-era workstation.
pub const MAX_STATE_VECTOR_QUBITS: usize = 24;

/// Maximum qubit count at which a 3-qubit gate may be applied via the
/// `embed_3q` fallback. The fallback materializes a full `2^n × 2^n`
/// matrix; at n=14 this would be 4 GiB. We cap at n=12 (256 MiB peak).
///
/// Loop 3 (Loop 2 review §"3-qubit gate ceiling"): the previous limit
/// `n >= MAX_QUBITS` (i.e. n=14 still allowed) silently OOMed during
/// `equivalence_by_sampling` on circuits with one Toffoli at ~14 qubits,
/// contradicting the documented `MAX_STATE_VECTOR_QUBITS=24` ceiling.
/// A native in-place 3q kernel is tracked as future work.
pub const MAX_3Q_EMBED_QUBITS: usize = 12;

fn embed_single_qubit(gate_2x2: &DMatrix<C>, target: usize, n_qubits: usize) -> DMatrix<C> {
    debug_assert!(n_qubits < SHIFT_GUARD);
    debug_assert!(target < n_qubits);
    if n_qubits >= SHIFT_GUARD {
        return DMatrix::<C>::identity(1, 1);
    }
    if target >= n_qubits {
        let dim = 1usize << n_qubits;
        return DMatrix::<C>::identity(dim, dim);
    }
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
    debug_assert!(n_qubits < SHIFT_GUARD);
    debug_assert!(q0 < n_qubits && q1 < n_qubits && q0 != q1);
    if n_qubits >= SHIFT_GUARD {
        return DMatrix::<C>::identity(1, 1);
    }
    if q0 >= n_qubits || q1 >= n_qubits || q0 == q1 {
        let dim = 1usize << n_qubits;
        return DMatrix::<C>::identity(dim, dim);
    }
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
    debug_assert!(n_qubits < SHIFT_GUARD);
    debug_assert!(q0 < n_qubits && q1 < n_qubits && q2 < n_qubits);
    debug_assert!(q0 != q1 && q1 != q2 && q0 != q2);
    if n_qubits >= SHIFT_GUARD {
        return DMatrix::<C>::identity(1, 1);
    }
    if q0 >= n_qubits || q1 >= n_qubits || q2 >= n_qubits || q0 == q1 || q1 == q2 || q0 == q2 {
        let dim = 1usize << n_qubits;
        return DMatrix::<C>::identity(dim, dim);
    }
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
    // Loop 5 §NI-1: prefer the fallible unroll path so missing
    // custom-gate registry entries propagate as Simulation errors
    // rather than silently producing a wrong-but-finite unitary.
    let unrolled = crate::transpiler::decomposition::try_unroll_custom_gates(circuit)
        .map_err(|e| QRustError::Simulation(format!("custom-gate unroll failed: {e}")))?;
    let n = unrolled.num_qubits;
    if n > MAX_QUBITS || n >= SHIFT_GUARD {
        return Err(QRustError::Simulation(format!(
            "unrolled circuit has {n} qubits; cannot simulate (max = {})",
            MAX_QUBITS.min(SHIFT_GUARD - 1)
        )));
    }
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
                0 => continue,
                1 => {
                    if qubits.is_empty() {
                        return Err(QRustError::Simulation("1-qubit gate has no target".into()));
                    }
                    if qubits[0] >= n {
                        return Err(QRustError::Simulation(format!(
                            "1-qubit gate target {} out of range (n={})",
                            qubits[0], n
                        )));
                    }
                    embed_single_qubit(&local_u, qubits[0], n)
                }
                2 => {
                    if qubits.len() < 2 {
                        return Err(QRustError::Simulation(
                            "2-qubit gate has insufficient qubits".into(),
                        ));
                    }
                    if qubits[0] >= n || qubits[1] >= n || qubits[0] == qubits[1] {
                        return Err(QRustError::Simulation(format!(
                            "invalid 2-qubit gate qubits {:?} (n={})",
                            &qubits[..2],
                            n
                        )));
                    }
                    embed_2q(&local_u, qubits[0], qubits[1], n)
                }
                3 => {
                    if qubits.len() < 3 {
                        return Err(QRustError::Simulation(
                            "3-qubit gate has insufficient qubits".into(),
                        ));
                    }
                    if qubits[0] >= n
                        || qubits[1] >= n
                        || qubits[2] >= n
                        || qubits[0] == qubits[1]
                        || qubits[1] == qubits[2]
                        || qubits[0] == qubits[2]
                    {
                        return Err(QRustError::Simulation(format!(
                            "invalid 3-qubit gate qubits {:?} (n={})",
                            &qubits[..3],
                            n
                        )));
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

#[track_caller]
pub fn circuit_to_unitary(circuit: &Circuit) -> DMatrix<C> {
    try_circuit_to_unitary(circuit)
        .unwrap_or_else(|e| panic!("circuit_to_unitary: simulation failed: {e}"))
}

pub fn unitary_fidelity(u1: &DMatrix<C>, u2: &DMatrix<C>) -> f64 {
    let d = u1.nrows();
    assert_eq!(u1.shape(), u2.shape(), "unitary_fidelity: shape mismatch");
    let trace: C = (u1.adjoint() * u2).trace();
    trace.norm_sqr() / (d as f64 * d as f64)
}

#[track_caller]
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

// ─── State-vector evolution (Loop 3 EXPLORATORY #4 Phase 1) ────────────────
//
// Unblocks verification at >14 qubits by avoiding full unitary materialization.
// Memory: O(2^n) for the state vector vs O(4^n) for the unitary.
// Time: O(g · 2^n) per gate vs O(g · 4^n).

/// Apply a 1-qubit gate to a state vector in-place.
/// Touches each amplitude exactly once: O(2^n).
fn apply_1q_gate(state: &mut DVector<C>, u: &DMatrix<C>, q: usize, n: usize) {
    debug_assert!(q < n);
    let dim = 1usize << n;
    let stride = 1usize << q;
    let u00 = u[(0, 0)];
    let u01 = u[(0, 1)];
    let u10 = u[(1, 0)];
    let u11 = u[(1, 1)];

    let mut i = 0;
    while i < dim {
        for j in 0..stride {
            let idx0 = i + j;
            let idx1 = idx0 + stride;
            let a = state[idx0];
            let b = state[idx1];
            state[idx0] = u00 * a + u01 * b;
            state[idx1] = u10 * a + u11 * b;
        }
        i += stride * 2;
    }
}

/// Apply a 2-qubit gate to a state vector in-place. `q0` is the LSB of
/// the local 4×4 matrix, `q1` is the MSB. O(2^n) per gate.
///
/// Loop 3 (review §"apply_2q_gate allocates O(2^n) bool buffer"): the
/// previous implementation allocated a `Vec<bool>` of length `2^n` per
/// gate invocation, leading to ~200 MiB of allocator churn over an
/// 18-qubit, 100-gate circuit. Replaced with a triple-nested stride
/// iteration (Stim-style; Gidney 2021, Quantum 5, 497) that visits each
/// 4-amplitude block exactly once with zero heap allocation.
fn apply_2q_gate(state: &mut DVector<C>, u: &DMatrix<C>, q0: usize, q1: usize, n: usize) {
    debug_assert!(q0 < n && q1 < n && q0 != q1);
    let dim = 1usize << n;
    let (lo, hi) = if q0 < q1 { (q0, q1) } else { (q1, q0) };
    let lo_stride = 1usize << lo;
    let hi_stride = 1usize << hi;

    // Cache the four matrix rows so the inner loop avoids repeated indexing.
    let u00 = u[(0, 0)];
    let u01 = u[(0, 1)];
    let u02 = u[(0, 2)];
    let u03 = u[(0, 3)];
    let u10 = u[(1, 0)];
    let u11 = u[(1, 1)];
    let u12 = u[(1, 2)];
    let u13 = u[(1, 3)];
    let u20 = u[(2, 0)];
    let u21 = u[(2, 1)];
    let u22 = u[(2, 2)];
    let u23 = u[(2, 3)];
    let u30 = u[(3, 0)];
    let u31 = u[(3, 1)];
    let u32 = u[(3, 2)];
    let u33 = u[(3, 3)];

    let bit_q0 = 1usize << q0;
    let bit_q1 = 1usize << q1;

    // Iterate over base indices that have bits `lo` and `hi` cleared.
    // The stride pattern: hi_block over [0, dim) in steps of 2*hi_stride;
    // within each, lo_block over [0, hi_stride) in steps of 2*lo_stride;
    // within each, tail over [0, lo_stride).
    let mut hi_block = 0usize;
    while hi_block < dim {
        let mut lo_block = 0usize;
        while lo_block < hi_stride {
            let block_base = hi_block + lo_block;
            for tail in 0..lo_stride {
                let base = block_base + tail;
                let i00 = base;
                let i01 = base | bit_q0;
                let i10 = base | bit_q1;
                let i11 = i01 | bit_q1;
                let a00 = state[i00];
                let a01 = state[i01];
                let a10 = state[i10];
                let a11 = state[i11];
                state[i00] = u00 * a00 + u01 * a01 + u02 * a10 + u03 * a11;
                state[i01] = u10 * a00 + u11 * a01 + u12 * a10 + u13 * a11;
                state[i10] = u20 * a00 + u21 * a01 + u22 * a10 + u23 * a11;
                state[i11] = u30 * a00 + u31 * a01 + u32 * a10 + u33 * a11;
            }
            lo_block += 2 * lo_stride;
        }
        hi_block += 2 * hi_stride;
    }
}

/// Apply a 3-qubit gate to a state vector. Falls back to embed_3q +
/// matrix-vector multiply.
///
/// Loop 3 (review §"apply_3q_gate ceiling"): the embed fallback
/// materializes a `2^n × 2^n` matrix. We cap n at [`MAX_3Q_EMBED_QUBITS`]
/// (=12, ~256 MiB peak) rather than the previous `MAX_QUBITS=14` which
/// would have peaked at 4 GiB. A native in-place 3q kernel is future work.
fn apply_3q_gate(
    state: &mut DVector<C>,
    u: &DMatrix<C>,
    q0: usize,
    q1: usize,
    q2: usize,
    n: usize,
) -> Result<()> {
    if n > MAX_3Q_EMBED_QUBITS {
        return Err(QRustError::Simulation(format!(
            "3-qubit gate at n={n} exceeds embed-fallback limit \
             ({MAX_3Q_EMBED_QUBITS}). Decompose to 1q+2q gates first; \
             a native 3q in-place kernel is tracked as future work."
        )));
    }
    let full = embed_3q(u, q0, q1, q2, n);
    *state = &full * &*state;
    Ok(())
}

fn apply_gate_to_state(
    state: &mut DVector<C>,
    name: &crate::ir::GateType,
    qubits: &[usize],
    params: &[f64],
    n: usize,
) -> Result<()> {
    let local_u = name.unitary(params);
    match name.num_qubits() {
        0 => Ok(()),
        1 => {
            if qubits.is_empty() || qubits[0] >= n {
                return Err(QRustError::Simulation(format!(
                    "invalid 1-qubit gate target (n={n}, qubits={qubits:?})"
                )));
            }
            apply_1q_gate(state, &local_u, qubits[0], n);
            Ok(())
        }
        2 => {
            if qubits.len() < 2 || qubits[0] >= n || qubits[1] >= n || qubits[0] == qubits[1] {
                return Err(QRustError::Simulation(format!(
                    "invalid 2-qubit gate (n={n}, qubits={qubits:?})"
                )));
            }
            apply_2q_gate(state, &local_u, qubits[0], qubits[1], n);
            Ok(())
        }
        3 => {
            if qubits.len() < 3
                || qubits[0] >= n
                || qubits[1] >= n
                || qubits[2] >= n
                || qubits[0] == qubits[1]
                || qubits[1] == qubits[2]
                || qubits[0] == qubits[2]
            {
                return Err(QRustError::Simulation(format!(
                    "invalid 3-qubit gate (n={n}, qubits={qubits:?})"
                )));
            }
            apply_3q_gate(state, &local_u, qubits[0], qubits[1], qubits[2], n)
        }
        nq => Err(QRustError::Simulation(format!(
            "unsupported gate arity: {nq}"
        ))),
    }
}

/// Evolves an initial state vector through a circuit.
///
/// `O(g · 2^n)` time, `O(2^n)` memory — vs `O(g · 4^n)` / `O(4^n)` for
/// `try_circuit_to_unitary`. Supports up to [`MAX_STATE_VECTOR_QUBITS`].
pub fn evolve_state(circuit: &Circuit, init: &DVector<C>) -> Result<DVector<C>> {
    let n = circuit.num_qubits;
    if n > MAX_STATE_VECTOR_QUBITS {
        return Err(QRustError::Simulation(format!(
            "evolve_state: {n} qubits exceeds practical limit ({MAX_STATE_VECTOR_QUBITS})"
        )));
    }
    if n >= SHIFT_GUARD {
        return Err(QRustError::Simulation(format!(
            "evolve_state: {n} qubits exceeds SHIFT_GUARD"
        )));
    }
    let expected_dim = 1usize << n;
    if init.len() != expected_dim {
        return Err(QRustError::Simulation(format!(
            "evolve_state: init has dimension {} (expected {})",
            init.len(),
            expected_dim
        )));
    }
    // Loop 5 §NI-1: same fail-loud-on-missing-custom-gate stance as
    // `try_circuit_to_unitary` above.
    let unrolled = crate::transpiler::decomposition::try_unroll_custom_gates(circuit)
        .map_err(|e| QRustError::Simulation(format!("custom-gate unroll failed: {e}")))?;
    let mut state = init.clone();
    for op in &unrolled.operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            apply_gate_to_state(&mut state, name, qubits, params, n)?;
        }
    }
    Ok(state)
}

/// Splitmix64 — simple, fast, good-quality PRNG with full 64-bit state.
struct SplitMix64(u64);
impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_add(0x9E3779B97F4A7C15))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn next_unit(&mut self) -> f64 {
        let bits = (self.next_u64() >> 11) as f64;
        bits * (1.0 / ((1u64 << 53) as f64))
    }
    fn next_gaussian_pair(&mut self) -> (f64, f64) {
        let u1 = (self.next_unit()).max(f64::MIN_POSITIVE);
        let u2 = self.next_unit();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
    }
}

/// Generates a Haar-random pure state via complex Gaussians + normalization.
fn haar_random_state(n: usize, rng: &mut SplitMix64) -> DVector<C> {
    let dim = 1usize << n;
    let mut psi = DVector::<C>::zeros(dim);
    let mut i = 0;
    while i < dim {
        let (g1, g2) = rng.next_gaussian_pair();
        let (g3, g4) = rng.next_gaussian_pair();
        psi[i] = C::new(g1, g2);
        if i + 1 < dim {
            psi[i + 1] = C::new(g3, g4);
        }
        i += 2;
    }
    let norm = psi.norm();
    if norm > 0.0 {
        psi /= C::new(norm, 0.0);
    }
    psi
}

/// Random-Haar input-state equivalence check.
///
/// Samples `k` Haar-random pure states, evolves each through both circuits
/// via `evolve_state`, and returns the *minimum* state fidelity observed
/// (`|⟨ψ₁|ψ₂⟩|²`). A return value of 1.0 (within tolerance) is strong
/// evidence — but not proof — that the circuits are equivalent up to a
/// global phase.
///
/// Cost: `O(k · g · 2^n)` time, `O(2^n)` memory. Practical to ~22 qubits.
///
/// **Statistical caveat**: this is a probabilistic check, not a proof.
/// For two unitaries differing by an `O(1)`-norm perturbation, a single
/// Haar sample's `1 - fidelity` is concentrated around its mean with
/// std-dev `O(1/√d)` where `d = 2^n` (Mele 2024, §III.A). So `k`
/// independent samples give a per-sample false-positive rate of roughly
/// `O(d^{-1/2}) = 2^{-n/2}`, and `k` samples accumulate (informally)
/// to `O(2^{-nk/2})` — e.g. `k=4, n=8` yields `< 2^{-16}`.
///
/// References:
/// - Burgholzer & Wille 2021 (IEEE TCAD 40), §III: motivation for
///   sampling-based equivalence at scale.
/// - Mele 2024 (Quantum 8, 1340), "Introduction to Haar measure tools":
///   concentration arguments for Haar-random sampling.
pub fn equivalence_by_sampling(c1: &Circuit, c2: &Circuit, k: usize, seed: u64) -> Result<f64> {
    if c1.num_qubits != c2.num_qubits {
        return Err(QRustError::Simulation(format!(
            "equivalence_by_sampling: qubit count mismatch ({} vs {})",
            c1.num_qubits, c2.num_qubits
        )));
    }
    if k == 0 {
        return Err(QRustError::Simulation(
            "equivalence_by_sampling: k must be >= 1".into(),
        ));
    }
    let mut rng = SplitMix64::new(seed);
    let mut min_fid = 1.0_f64;
    for _ in 0..k {
        let psi = haar_random_state(c1.num_qubits, &mut rng);
        let psi1 = evolve_state(c1, &psi)?;
        let psi2 = evolve_state(c2, &psi)?;
        let overlap = psi1.dotc(&psi2).norm_sqr();
        if overlap < min_fid {
            min_fid = overlap;
        }
    }
    Ok(min_fid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Circuit, GateType};

    /// Test-only seed; named to look intentional in test failure output.
    const QRUST_SEED: u64 = 0xCAFE_F00D_DEAD_BEEF;

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

    #[test]
    fn test_out_of_range_qubit_errors() {
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![5],
            params: vec![],
        });
        assert!(try_circuit_to_unitary(&circuit).is_err());
    }

    /// Loop 5 §NI-1: a circuit referencing a Custom gate that's not in
    /// the registry must produce a Simulation error rather than silently
    /// returning a wrong-but-finite unitary (previously this returned the
    /// 1-qubit identity because the `unroll_custom_gates` infallible
    /// helper swallowed the error).
    #[test]
    fn test_simulator_surfaces_missing_custom_gate() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: crate::ir::GateType::Custom("phantom_gate".into()),
            qubits: vec![0],
            params: vec![],
        });
        let res = try_circuit_to_unitary(&c);
        match res {
            Err(QRustError::Simulation(msg)) => {
                assert!(
                    msg.contains("phantom_gate") || msg.contains("custom-gate unroll"),
                    "expected error to identify the unroll failure, got: {msg}"
                );
            }
            other => panic!("expected Simulation error from missing Custom gate, got: {other:?}"),
        }
    }

    /// Loop 5 §NI-1: same contract for `evolve_state` (used by
    /// `equivalence_by_sampling` — see Loop 1 EP-4 Phase 1).
    #[test]
    fn test_evolve_state_surfaces_missing_custom_gate() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: crate::ir::GateType::Custom("phantom_gate".into()),
            qubits: vec![0],
            params: vec![],
        });
        let mut psi = DVector::<C>::zeros(2);
        psi[0] = C::new(1.0, 0.0);
        let res = evolve_state(&c, &psi);
        assert!(
            matches!(res, Err(QRustError::Simulation(_))),
            "expected Simulation error, got: {res:?}"
        );
    }

    // ─── State-vector evolution tests ──────────────────────────────────

    #[test]
    fn test_evolve_state_matches_unitary_on_basis_state() {
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
        let mut zero = DVector::<C>::zeros(4);
        zero[0] = C::new(1.0, 0.0);
        let final_state = evolve_state(&circuit, &zero).unwrap();
        let expected = u.column(0);
        for i in 0..4 {
            assert!(
                (final_state[i] - expected[i]).norm() < 1e-10,
                "mismatch at {i}: got {}, expected {}",
                final_state[i],
                expected[i]
            );
        }
    }

    #[test]
    fn test_evolve_state_preserves_norm() {
        let mut circuit = Circuit::new(3, 0);
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
        circuit.add_op(Operation::Gate {
            name: GateType::RY,
            qubits: vec![2],
            params: vec![0.7],
        });
        let mut rng = SplitMix64::new(42);
        let psi = haar_random_state(3, &mut rng);
        let phi = evolve_state(&circuit, &psi).unwrap();
        assert!((phi.norm() - 1.0).abs() < 1e-10, "norm = {}", phi.norm());
    }

    #[test]
    fn test_evolve_state_dimension_mismatch_errors() {
        let circuit = Circuit::new(3, 0);
        let bad_init = DVector::<C>::zeros(4); // expects 8
        assert!(evolve_state(&circuit, &bad_init).is_err());
    }

    #[test]
    fn test_equivalence_by_sampling_identical_circuits() {
        let mut c1 = Circuit::new(3, 0);
        c1.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c1.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c1.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 2],
            params: vec![],
        });
        let c2 = c1.clone();
        let fid = equivalence_by_sampling(&c1, &c2, 4, QRUST_SEED).unwrap();
        assert!(
            (fid - 1.0).abs() < 1e-10,
            "identical circuits should have fidelity 1, got {fid}"
        );
    }

    #[test]
    fn test_equivalence_by_sampling_distinct_circuits() {
        let mut c1 = Circuit::new(2, 0);
        c1.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let mut c2 = Circuit::new(2, 0);
        c2.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        let fid = equivalence_by_sampling(&c1, &c2, 8, 7).unwrap();
        assert!(
            fid < 0.99,
            "H ≠ X should give fidelity well below 1, got {fid}"
        );
    }

    #[test]
    fn test_equivalence_by_sampling_works_above_max_qubits() {
        let mut c = Circuit::new(18, 0);
        for i in 0..18 {
            c.add_op(Operation::Gate {
                name: GateType::H,
                qubits: vec![i],
                params: vec![],
            });
        }
        for i in 0..17 {
            c.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![i, i + 1],
                params: vec![],
            });
        }
        let fid = equivalence_by_sampling(&c, &c, 2, 12345).unwrap();
        assert!(
            (fid - 1.0).abs() < 1e-9,
            "self-equivalence at 18q failed: fid={fid}"
        );
    }

    #[test]
    fn test_evolve_state_above_limit_errors() {
        let circuit = Circuit::new(MAX_STATE_VECTOR_QUBITS + 1, 0);
        let small = DVector::<C>::zeros(2);
        assert!(evolve_state(&circuit, &small).is_err());
    }

    /// Loop 3 review §"apply_3q_gate ceiling": at n > MAX_3Q_EMBED_QUBITS,
    /// 3q gates must error out cleanly rather than OOM.
    #[test]
    fn test_apply_3q_gate_above_limit_errors() {
        // n = MAX_3Q_EMBED_QUBITS + 1 = 13. Build a circuit with one CCX.
        let n = MAX_3Q_EMBED_QUBITS + 1;
        assert!(n <= MAX_QUBITS, "test assumes n stays within MAX_QUBITS");
        let mut c = Circuit::new(n, 0);
        c.add_op(Operation::Gate {
            name: GateType::CCX,
            qubits: vec![0, 1, 2],
            params: vec![],
        });
        let dim = 1usize << n;
        let mut psi = DVector::<C>::zeros(dim);
        psi[0] = C::new(1.0, 0.0);
        let res = evolve_state(&c, &psi);
        assert!(
            matches!(res, Err(QRustError::Simulation(_))),
            "expected Simulation error for 3q gate above limit, got {res:?}"
        );
    }

    /// Loop 3 review §"apply_2q_gate stride iteration": correctness
    /// regression test. Apply a 2q gate at non-contiguous (q0, q1) and
    /// verify the result matches the embed-fallback path.
    #[test]
    fn test_apply_2q_gate_non_contiguous_matches_embed() {
        let mut c = Circuit::new(5, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 4],
            params: vec![],
        });
        let dim = 1usize << 5;
        // Initial state: random superposition (use a deterministic seed).
        let mut rng = SplitMix64::new(0xDEAD_BEEF);
        let psi = haar_random_state(5, &mut rng);
        // Path A: stride kernel.
        let phi_stride = evolve_state(&c, &psi).unwrap();
        // Path B: full unitary materialization.
        let u = circuit_to_unitary(&c);
        let phi_embed = &u * &psi;
        for i in 0..dim {
            assert!(
                (phi_stride[i] - phi_embed[i]).norm() < 1e-10,
                "stride/embed mismatch at i={i}: stride={}, embed={}",
                phi_stride[i],
                phi_embed[i]
            );
        }
    }
}
