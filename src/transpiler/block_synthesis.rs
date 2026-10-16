//! Conservative resynthesis of two-wire blocks in the lowered U/CX circuit.
//!
//! A block follows both wire dependencies and ends at the first interaction
//! with a third wire. Disjoint U/CX operations may be passed because they
//! commute with the entire block. Every other operation, including a barrier
//! on disjoint wires, is a global collection boundary. Layouts are untouched.
//! Unclaimed leading U runs on either wire move forward to the first CX;
//! the joint replacement is never moved before that CX's dependencies.
//!
//! Replacement requires fewer CX gates and an independently reconstructed,
//! phase-aligned matrix residual. The unnormalized Frobenius norm bounds the
//! operator norm; summing accepted residuals bounds composition error for
//! unitary blocks in exact arithmetic. This budget does not include later
//! passes, floating-point roundoff in the checker, or QASM serialization.

use super::decomposition::BasisDecompositionPass;
use super::optimization::GateFusionPass;
use super::pass::Pass;
use super::property_set::PropertySet;
use super::synthesis::kak::KakSynthesizer;
use super::synthesis::Synthesizer;
use crate::ir::{Circuit, GateType, Operation};
use nalgebra::DMatrix;
use num_complex::Complex;

const BLOCK_RESIDUAL_LIMIT: f64 = 1e-10;
const CIRCUIT_RESIDUAL_BUDGET: f64 = 1e-9;

/// Lightweight counters, written to `two_qubit_block_stats` after each run.
/// No clocks or dense matrices larger than four by four are used by this pass.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TwoQubitBlockStats {
    pub blocks_considered: usize,
    pub synthesis_attempted: usize,
    pub pure_cx_fast_paths: usize,
    pub accepted: usize,
    pub cx_removed: usize,
    pub rejected_numerical: usize,
    pub rejected_no_cx_improvement: usize,
    pub residual_sum: f64,
    pub max_residual: f64,
}

/// Consolidate dependency-safe U/CX blocks and accept strict CX reductions.
/// Numerical failures retain the original block without changing metadata.
#[derive(Debug, Clone, Copy)]
pub struct TwoQubitBlockSynthesisPass;

/// Only well-formed finite lowered operations can participate. Unsupported or
/// malformed operations are retained and act as boundaries, rather than being
/// interpreted using a shared named-gate decomposition oracle.
fn lowered_wires(op: &Operation, width: usize) -> Option<&[usize]> {
    match op {
        Operation::Gate {
            name,
            qubits,
            params,
        } if qubits.iter().all(|&q| q < width) && params.iter().all(|p| p.is_finite()) => {
            match name {
                GateType::U if qubits.len() == 1 && params.len() == 3 => Some(qubits),
                GateType::CX
                    if qubits.len() == 2 && qubits[0] != qubits[1] && params.is_empty() =>
                {
                    Some(qubits)
                }
                _ => None,
            }
        }
        _ => None,
    }
}

fn is_cx(op: &Operation) -> bool {
    matches!(
        op,
        Operation::Gate {
            name: GateType::CX,
            ..
        }
    )
}

/// Construct next-on-wire dependencies once, avoiding repeated scans through
/// unrelated gates. A non-U/CX operation clears every wire's dependency tail.
fn wire_successors(circuit: &Circuit) -> Vec<[Option<usize>; 2]> {
    let mut successors = vec![[None, None]; circuit.operations.len()];
    let mut next = vec![None; circuit.num_qubits];
    for (index, op) in circuit.operations.iter().enumerate().rev() {
        if let Some(wires) = lowered_wires(op, circuit.num_qubits) {
            for (operand, &wire) in wires.iter().enumerate() {
                successors[index][operand] = next[wire];
                next[wire] = Some(index);
            }
        } else {
            next.fill(None);
        }
    }
    successors
}

fn wire_predecessors(circuit: &Circuit) -> Vec<[Option<usize>; 2]> {
    let mut predecessors = vec![[None, None]; circuit.operations.len()];
    let mut previous = vec![None; circuit.num_qubits];
    for (index, op) in circuit.operations.iter().enumerate() {
        if let Some(wires) = lowered_wires(op, circuit.num_qubits) {
            for (operand, &wire) in wires.iter().enumerate() {
                predecessors[index][operand] = previous[wire];
                previous[wire] = Some(index);
            }
        } else {
            previous.fill(None);
        }
    }
    predecessors
}

/// Absorb only local runs that have not already belonged to another block.
/// Each U moves forward across operations disjoint from its own wire. The
/// replacement must still be emitted at `start`, not the earliest U: the
/// other block wire may have an intervening interaction with a third wire.
fn leading_u_indices(
    circuit: &Circuit,
    predecessors: &[[Option<usize>; 2]],
    start: usize,
    processed: &[bool],
) -> Vec<usize> {
    let mut indices = Vec::new();
    for mut previous in predecessors[start] {
        while let Some(index) = previous {
            if processed[index]
                || !matches!(
                    circuit.operations[index],
                    Operation::Gate {
                        name: GateType::U,
                        ..
                    }
                )
            {
                break;
            }
            indices.push(index);
            previous = predecessors[index][0];
        }
    }
    indices.sort_unstable();
    indices
}

fn collect_block(
    circuit: &Circuit,
    successors: &[[Option<usize>; 2]],
    start: usize,
    pair: [usize; 2],
) -> Vec<usize> {
    let mut indices = Vec::new();
    let mut pending = [Some(start), Some(start)];
    while let Some(index) = pending.iter().flatten().copied().min() {
        let Some(wires) = lowered_wires(&circuit.operations[index], circuit.num_qubits) else {
            break;
        };
        if wires.iter().any(|q| !pair.contains(q)) {
            break;
        }
        indices.push(index);
        for (operand, &wire) in wires.iter().enumerate() {
            let local = usize::from(wire == pair[1]);
            pending[local] = successors[index][operand];
        }
    }
    indices
}

/// Direct OpenQASM U/CX matrices in little-endian order, independent of the
/// compiler's gate-definition, decomposition and simulator implementations.
fn block_matrix(ops: &[Operation], pair: [usize; 2]) -> Option<DMatrix<Complex<f64>>> {
    let mut matrix = DMatrix::<Complex<f64>>::identity(4, 4);
    for op in ops {
        let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        else {
            return None;
        };
        if params.iter().any(|p| !p.is_finite()) || qubits.iter().any(|q| !pair.contains(q)) {
            return None;
        }
        match name {
            GateType::U if qubits.len() == 1 && params.len() == 3 => {
                let bit = usize::from(qubits[0] == pair[1]);
                let mask = 1 << bit;
                let (sin, cos) = (params[0] / 2.0).sin_cos();
                let one = Complex::new(cos, 0.0);
                let a01 = -Complex::from_polar(sin, params[2]);
                let a10 = Complex::from_polar(sin, params[1]);
                let a11 = Complex::from_polar(cos, params[1] + params[2]);
                for row in (0..4).filter(|r| r & mask == 0) {
                    for col in 0..4 {
                        let low = matrix[(row, col)];
                        let high = matrix[(row | mask, col)];
                        matrix[(row, col)] = one * low + a01 * high;
                        matrix[(row | mask, col)] = a10 * low + a11 * high;
                    }
                }
            }
            GateType::CX if qubits.len() == 2 && qubits[0] != qubits[1] && params.is_empty() => {
                let control = 1 << usize::from(qubits[0] == pair[1]);
                let target = 1 << usize::from(qubits[1] == pair[1]);
                for row in (0..4).filter(|r| r & control != 0 && r & target == 0) {
                    matrix.swap_rows(row, row | target);
                }
            }
            _ => return None,
        }
    }
    matrix
        .iter()
        .all(|z| z.re.is_finite() && z.im.is_finite())
        .then_some(matrix)
}

fn phase_aligned_residual(
    original: &DMatrix<Complex<f64>>,
    candidate: &DMatrix<Complex<f64>>,
) -> Option<f64> {
    let inner: Complex<f64> = original
        .iter()
        .zip(candidate.iter())
        .map(|(a, b)| a.conj() * b)
        .sum();
    let magnitude = inner.norm();
    if !magnitude.is_finite() || magnitude <= f64::EPSILON {
        return None;
    }
    let residual = (candidate - original * (inner / magnitude)).norm();
    residual.is_finite().then_some(residual)
}

fn cx_op(control: usize, target: usize) -> Operation {
    Operation::Gate {
        name: GateType::CX,
        qubits: vec![control, target],
        params: vec![],
    }
}

/// GL(2,2) has only six elements. Enumerating zero to three CX gates avoids
/// numerical synthesis of cancellation chains and irreducible SWAPs.
fn pure_cx_candidate(ops: &[Operation], pair: [usize; 2]) -> Vec<Operation> {
    fn permutation(ops: &[Operation], pair: [usize; 2]) -> [usize; 4] {
        let mut perm = [0, 1, 2, 3];
        for op in ops {
            let Operation::Gate { qubits, .. } = op else {
                unreachable!()
            };
            let control = 1 << usize::from(qubits[0] == pair[1]);
            let target = 1 << usize::from(qubits[1] == pair[1]);
            for value in &mut perm {
                if *value & control != 0 {
                    *value ^= target;
                }
            }
        }
        perm
    }
    let desired = permutation(ops, pair);
    for count in 0..=3 {
        for word in 0..(1 << count) {
            let candidate: Vec<_> = (0..count)
                .map(|index| {
                    let reverse = (word >> index) & 1;
                    cx_op(pair[reverse], pair[1 - reverse])
                })
                .collect();
            if permutation(&candidate, pair) == desired {
                return candidate;
            }
        }
    }
    // Unreachable for well-formed CX gates, but retaining the original is safe.
    ops.to_vec()
}

fn mapped_candidate(matrix: &DMatrix<Complex<f64>>, pair: [usize; 2]) -> Option<Vec<Operation>> {
    let synthesized = KakSynthesizer.synthesize(matrix, &[GateType::U, GateType::CX])?;
    if synthesized.num_qubits != 2 || synthesized.num_cbits != 0 {
        return None;
    }
    let mut properties = PropertySet::new();
    let lowered = BasisDecompositionPass.run(&synthesized, &mut properties);
    let fused = GateFusionPass.run(&lowered, &mut properties);
    fused
        .operations
        .into_iter()
        .map(|op| {
            lowered_wires(&op, 2)?;
            let Operation::Gate {
                name,
                qubits,
                params,
            } = op
            else {
                return None;
            };
            Some(Operation::Gate {
                name,
                qubits: qubits.into_iter().map(|q| pair[q]).collect(),
                params,
            })
        })
        .collect()
}

impl Pass for TwoQubitBlockSynthesisPass {
    fn name(&self) -> &str {
        "TwoQubitBlockSynthesisPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        let successors = wire_successors(circuit);
        let predecessors = wire_predecessors(circuit);
        let mut processed = vec![false; circuit.operations.len()];
        let mut removed = vec![false; circuit.operations.len()];
        let mut replacements: Vec<Option<Vec<Operation>>> = vec![None; circuit.operations.len()];
        let mut stats = TwoQubitBlockStats::default();
        for (index, op) in circuit.operations.iter().enumerate() {
            if processed[index] || !is_cx(op) {
                continue;
            }
            let Some(wires) = lowered_wires(op, circuit.num_qubits) else {
                continue;
            };
            let pair = [wires[0].min(wires[1]), wires[0].max(wires[1])];
            let body = collect_block(circuit, &successors, index, pair);
            for &i in &body {
                processed[i] = true;
            }
            let cx_before = body
                .iter()
                .filter(|&&i| is_cx(&circuit.operations[i]))
                .count();
            if cx_before < 2 || stats.residual_sum >= CIRCUIT_RESIDUAL_BUDGET {
                continue;
            }
            let mut indices = leading_u_indices(circuit, &predecessors, index, &processed);
            for &i in &indices {
                processed[i] = true;
            }
            let prefix_len = indices.len();
            indices.extend(body);
            stats.blocks_considered += 1;
            let original: Vec<_> = indices
                .iter()
                .map(|&i| circuit.operations[i].clone())
                .collect();
            let Some(matrix) = block_matrix(&original, pair) else {
                stats.rejected_numerical += 1;
                continue;
            };
            let candidate = if original[prefix_len..].iter().all(is_cx) {
                stats.pure_cx_fast_paths += 1;
                let mut candidate = original[..prefix_len].to_vec();
                candidate.extend(pure_cx_candidate(&original[prefix_len..], pair));
                Some(candidate)
            } else {
                stats.synthesis_attempted += 1;
                mapped_candidate(&matrix, pair)
            };
            let Some(candidate) = candidate else {
                stats.rejected_numerical += 1;
                continue;
            };
            let cx_after = candidate.iter().filter(|op| is_cx(op)).count();
            if cx_after >= cx_before {
                stats.rejected_no_cx_improvement += 1;
                continue;
            }
            let residual = block_matrix(&candidate, pair)
                .and_then(|out| phase_aligned_residual(&matrix, &out));
            let Some(residual) = residual.filter(|r| {
                *r <= BLOCK_RESIDUAL_LIMIT && stats.residual_sum + r <= CIRCUIT_RESIDUAL_BUDGET
            }) else {
                stats.rejected_numerical += 1;
                continue;
            };
            stats.accepted += 1;
            stats.cx_removed += cx_before - cx_after;
            stats.residual_sum += residual;
            stats.max_residual = stats.max_residual.max(residual);
            for i in indices {
                removed[i] = true;
            }
            replacements[index] = Some(candidate);
        }
        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        for (index, op) in circuit.operations.iter().enumerate() {
            if let Some(replacement) = replacements[index].take() {
                out.operations.extend(replacement);
            }
            if !removed[index] {
                out.operations.push(op.clone());
            }
        }
        property_set.insert("two_qubit_block_stats", stats);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_cx_group_fast_path_is_exact_and_avoids_kak() {
        for word in 0..64 {
            let mut circuit = Circuit::new(2, 0);
            for index in 0..6 {
                let reverse = (word >> index) & 1;
                circuit.add_op(cx_op(reverse, 1 - reverse));
            }
            let mut properties = PropertySet::new();
            let out = TwoQubitBlockSynthesisPass.run(&circuit, &mut properties);
            assert_eq!(
                block_matrix(&circuit.operations, [0, 1]),
                block_matrix(&out.operations, [0, 1])
            );
            assert!(out.operations.len() <= 3);
            assert_eq!(
                properties
                    .get::<TwoQubitBlockStats>("two_qubit_block_stats")
                    .unwrap()
                    .synthesis_attempted,
                0
            );
        }
    }

    #[test]
    fn block_collection_respects_third_wire_and_global_boundaries() {
        let mut circuit = Circuit::new(4, 0);
        circuit.operations = vec![
            cx_op(0, 1),
            cx_op(2, 3),
            cx_op(0, 1),
            cx_op(1, 2),
            cx_op(0, 1),
            Operation::Barrier { qubits: vec![3] },
            cx_op(0, 1),
        ];
        let successors = wire_successors(&circuit);
        assert_eq!(collect_block(&circuit, &successors, 0, [0, 1]), vec![0, 2]);
        assert_eq!(collect_block(&circuit, &successors, 4, [0, 1]), vec![4]);
        let out = TwoQubitBlockSynthesisPass.run(&circuit, &mut PropertySet::new());
        assert_eq!(
            out.operations,
            circuit.operations[1..]
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != 1)
                .map(|(_, op)| op.clone())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn residual_checks_phase_and_coherent_error_not_fidelity_only() {
        let identity = DMatrix::<Complex<f64>>::identity(4, 4);
        let phased = &identity * Complex::from_polar(1.0, 0.71);
        assert!(phase_aligned_residual(&identity, &phased).unwrap() < 1e-14);
        let mut wrong = identity.clone();
        wrong[(3, 3)] = Complex::from_polar(1.0, 1e-6);
        assert!(phase_aligned_residual(&identity, &wrong).unwrap() > BLOCK_RESIDUAL_LIMIT);
    }

    #[test]
    fn leading_runs_move_to_first_cx_without_crossing_other_wire_dependencies() {
        let u = |q| Operation::Gate {
            name: GateType::U,
            qubits: vec![q],
            params: vec![0.23, 0.41, -0.37],
        };
        let mut circuit = Circuit::new(3, 0);
        circuit.operations = vec![u(0), cx_op(1, 2), u(1), cx_op(0, 1), cx_op(0, 1)];
        let predecessors = wire_predecessors(&circuit);
        let mut processed = vec![false; circuit.operations.len()];
        assert_eq!(
            leading_u_indices(&circuit, &predecessors, 3, &processed),
            vec![0, 2]
        );
        processed[0] = true;
        assert_eq!(
            leading_u_indices(&circuit, &predecessors, 3, &processed),
            vec![2]
        );
        let mut properties = PropertySet::new();
        let out = TwoQubitBlockSynthesisPass.run(&circuit, &mut properties);
        // U(1) already belongs to the earlier one-CX block and stays in place.
        // The unclaimed U(0) moves to the original first CX at index 3.
        assert_eq!(out.operations, vec![cx_op(1, 2), u(1), u(0)]);
        assert_eq!(
            properties
                .get::<TwoQubitBlockStats>("two_qubit_block_stats")
                .unwrap()
                .synthesis_attempted,
            0
        );
        circuit
            .operations
            .insert(2, Operation::Barrier { qubits: vec![2] });
        assert_eq!(
            leading_u_indices(
                &circuit,
                &wire_predecessors(&circuit),
                4,
                &vec![false; circuit.operations.len()]
            ),
            vec![3]
        );
    }

    #[test]
    fn public_toggle_only_enables_level_three_lowered_cleanup() {
        use crate::transpiler::{transpile_with_report, TranspilerConfig};
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(cx_op(0, 1));
        let configurations = [
            (TranspilerConfig::default(), false),
            (
                TranspilerConfig::builder().optimization_level(2).build(),
                false,
            ),
            (
                TranspilerConfig::builder().optimization_level(3).build(),
                true,
            ),
            (
                TranspilerConfig::builder()
                    .optimization_level(3)
                    .two_qubit_block_synthesis(false)
                    .build(),
                false,
            ),
            (
                TranspilerConfig::builder()
                    .optimization_level(3)
                    .post_routing_cleanup(false)
                    .build(),
                false,
            ),
            (
                TranspilerConfig::builder()
                    .optimization_level(3)
                    .decompose_basis(false)
                    .build(),
                false,
            ),
        ];
        for (config, enabled) in configurations {
            let (_, report) = transpile_with_report(&circuit, Some(config)).unwrap();
            assert_eq!(report.two_qubit_block_stats.is_some(), enabled);
            assert!(report.pass_diagnostics.is_empty());
        }
    }
}
