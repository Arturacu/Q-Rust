//! Guarded cleanup of local gates introduced by two-qubit block synthesis.
//!
//! Only adjacent U gates on the same wire are combined. They may be separated
//! by operations on disjoint wires, but never by another operation on that
//! wire. Every non-U/CX operation is a global boundary. The scan and matrix
//! work are linear in the circuit size; no commutation approximation is used.

use super::pass::Pass;
use super::property_set::PropertySet;
use super::synthesis::zyz::zyz_decomposition;
use crate::ir::{Circuit, GateType, Operation};
use num_complex::Complex64;

const RUN_RESIDUAL_LIMIT: f64 = 1e-12;
const CIRCUIT_RESIDUAL_BUDGET: f64 = 1e-10;
const IDENTITY_RESIDUAL_LIMIT: f64 = 1e-14;

type Matrix = [[Complex64; 2]; 2];

/// Cheap counters stored under `post_block_fusion_stats`; no clocks are read.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PostBlockFusionStats {
    pub runs_attempted: usize,
    pub accepted: usize,
    pub u_gates_removed: usize,
    pub rejected_numerical: usize,
    pub rejected_budget: usize,
    pub residual_sum: f64,
    pub max_residual: f64,
}

/// Fuse lowered local runs with independently checked numerical error budgets.
/// Layout metadata and all unsupported operations are preserved.
#[derive(Debug, Clone, Copy)]
pub struct PostBlockFusionPass;

struct LocalRun {
    indices: Vec<usize>,
    matrix: Matrix,
}

fn finite(matrix: &Matrix) -> bool {
    matrix
        .iter()
        .flatten()
        .all(|z| z.re.is_finite() && z.im.is_finite())
}

/// Analytic OpenQASM U matrix, separate from GateDefinition and the shared
/// simulator. In particular, finite parameters need not have a finite sum.
fn analytic_u(params: &[f64]) -> Option<Matrix> {
    if params.len() != 3 || params.iter().any(|p| !p.is_finite()) {
        return None;
    }
    let (sine, cosine) = (params[0] / 2.0).sin_cos();
    let matrix = [
        [
            Complex64::new(cosine, 0.0),
            -Complex64::from_polar(sine, params[2]),
        ],
        [
            Complex64::from_polar(sine, params[1]),
            Complex64::from_polar(cosine, params[1] + params[2]),
        ],
    ];
    finite(&matrix).then_some(matrix)
}

fn multiply(left: &Matrix, right: &Matrix) -> Matrix {
    std::array::from_fn(|row| {
        std::array::from_fn(|col| left[row][0] * right[0][col] + left[row][1] * right[1][col])
    })
}

fn residual(original: &Matrix, candidate: &Matrix) -> Option<f64> {
    if !finite(original) || !finite(candidate) {
        return None;
    }
    let overlap: Complex64 = original
        .iter()
        .flatten()
        .zip(candidate.iter().flatten())
        .map(|(a, b)| a.conj() * b)
        .sum();
    let magnitude = overlap.norm();
    if !magnitude.is_finite() || magnitude <= f64::EPSILON {
        return None;
    }
    let phase = overlap / magnitude;
    let value = original
        .iter()
        .flatten()
        .zip(candidate.iter().flatten())
        .map(|(a, b)| (b - a * phase).norm_sqr())
        .sum::<f64>()
        .sqrt();
    value.is_finite().then_some(value)
}

fn flush(
    qubit: usize,
    pending: &mut [Option<LocalRun>],
    output: &mut [Option<Operation>],
    stats: &mut PostBlockFusionStats,
) {
    let Some(run) = pending[qubit].take() else {
        return;
    };
    // Avoid numerical work and representation changes for isolated U gates.
    if run.indices.len() < 2 {
        return;
    }
    stats.runs_attempted += 1;
    if stats.residual_sum >= CIRCUIT_RESIDUAL_BUDGET {
        stats.rejected_budget += 1;
        return;
    }
    if !finite(&run.matrix) {
        stats.rejected_numerical += 1;
        return;
    }
    let identity = [
        [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
    ];
    let identity_error = residual(&run.matrix, &identity);
    let (replacement, error) = if identity_error.is_some_and(|error| {
        error <= IDENTITY_RESIDUAL_LIMIT && stats.residual_sum + error <= CIRCUIT_RESIDUAL_BUDGET
    }) {
        (None, identity_error.unwrap())
    } else {
        let (theta, phi, lambda, _) = zyz_decomposition(run.matrix);
        let params = vec![theta, phi, lambda];
        let Some(error) = analytic_u(&params)
            .and_then(|candidate| residual(&run.matrix, &candidate))
            .filter(|&error| error <= RUN_RESIDUAL_LIMIT)
        else {
            stats.rejected_numerical += 1;
            return;
        };
        (
            Some(Operation::Gate {
                name: GateType::U,
                qubits: vec![qubit],
                params,
            }),
            error,
        )
    };
    if stats.residual_sum + error > CIRCUIT_RESIDUAL_BUDGET {
        stats.rejected_budget += 1;
        return;
    }

    // Materialize at the last U of this run. Earlier U gates move right only
    // across operations disjoint from this wire. No multiwire operation moves.
    let anchor = *run.indices.last().unwrap();
    stats.accepted += 1;
    stats.u_gates_removed += run.indices.len() - usize::from(replacement.is_some());
    stats.residual_sum += error;
    stats.max_residual = stats.max_residual.max(error);
    for index in run.indices {
        output[index] = None;
    }
    output[anchor] = replacement;
}

impl Pass for PostBlockFusionPass {
    fn name(&self) -> &str {
        "PostBlockFusionPass"
    }

    fn run(&self, circuit: &Circuit, properties: &mut PropertySet) -> Circuit {
        let mut pending: Vec<Option<LocalRun>> = (0..circuit.num_qubits).map(|_| None).collect();
        // Stale entries are harmless: each run adds one entry, so all global
        // boundary flushes together cost O(number of operations), not O(n*w).
        let mut active = Vec::new();
        let mut output: Vec<_> = circuit.operations.iter().cloned().map(Some).collect();
        let mut stats = PostBlockFusionStats::default();
        for (index, operation) in circuit.operations.iter().enumerate() {
            match operation {
                Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                } if qubits.len() == 1 && qubits[0] < circuit.num_qubits => {
                    if let Some(matrix) = analytic_u(params) {
                        let qubit = qubits[0];
                        if let Some(run) = &mut pending[qubit] {
                            run.indices.push(index);
                            run.matrix = multiply(&matrix, &run.matrix);
                        } else {
                            pending[qubit] = Some(LocalRun {
                                indices: vec![index],
                                matrix,
                            });
                            active.push(qubit);
                        }
                        continue;
                    }
                }
                Operation::Gate {
                    name: GateType::CX,
                    qubits,
                    params,
                } if qubits.len() == 2
                    && qubits[0] != qubits[1]
                    && qubits.iter().all(|&q| q < circuit.num_qubits)
                    && params.is_empty() =>
                {
                    for &qubit in qubits {
                        flush(qubit, &mut pending, &mut output, &mut stats);
                    }
                    continue;
                }
                _ => {}
            }
            // All non-lowered, malformed, nonfinite or nonunitary operations
            // remain in place and conservatively delimit every pending run.
            while let Some(qubit) = active.pop() {
                flush(qubit, &mut pending, &mut output, &mut stats);
            }
        }
        while let Some(qubit) = active.pop() {
            flush(qubit, &mut pending, &mut output, &mut stats);
        }
        let mut result = circuit.clone();
        result.operations = output.into_iter().flatten().collect();
        properties.insert("post_block_fusion_stats", stats);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn u(qubit: usize, theta: f64, phi: f64, lambda: f64) -> Operation {
        Operation::Gate {
            name: GateType::U,
            qubits: vec![qubit],
            params: vec![theta, phi, lambda],
        }
    }

    #[test]
    fn boundary_fusion_retains_other_wires_and_resets_stats() {
        let mut circuit = Circuit::new(3, 0);
        circuit.operations = vec![
            u(0, 0.2, 0.4, 0.7),
            Operation::Gate {
                name: GateType::CX,
                qubits: vec![1, 2],
                params: vec![],
            },
            u(0, -0.3, 0.1, 0.5),
        ];
        let mut properties = PropertySet::new();
        properties.insert("initial_layout", vec![2_usize, 1, 0]);
        let output = PostBlockFusionPass.run(&circuit, &mut properties);
        assert_eq!(output.operations.len(), 2);
        assert_eq!(output.operations[0], circuit.operations[1]);
        let stats = properties
            .get::<PostBlockFusionStats>("post_block_fusion_stats")
            .unwrap();
        assert_eq!(stats.runs_attempted, 1);
        assert_eq!(stats.accepted, 1);
        assert_eq!(stats.u_gates_removed, 1);
        assert!(stats.max_residual < 2e-15);
        assert_eq!(
            properties.get::<Vec<usize>>("initial_layout"),
            Some(&vec![2, 1, 0])
        );
        PostBlockFusionPass.run(&output, &mut properties);
        assert_eq!(
            properties.get::<PostBlockFusionStats>("post_block_fusion_stats"),
            Some(&PostBlockFusionStats::default())
        );
    }

    #[test]
    fn a_disjoint_barrier_and_overflowing_phase_are_global_boundaries() {
        for boundary in [
            Operation::Barrier { qubits: vec![1] },
            u(1, 0.2, f64::MAX, f64::MAX),
        ] {
            let mut circuit = Circuit::new(2, 0);
            circuit.operations = vec![u(0, 0.2, 0.0, 0.0), boundary, u(0, 0.3, 0.0, 0.0)];
            let output = PostBlockFusionPass.run(&circuit, &mut PropertySet::new());
            assert_eq!(output.operations, circuit.operations);
        }
    }

    #[test]
    fn scalar_identity_is_removed_but_small_rotation_is_kept() {
        let mut circuit = Circuit::new(1, 0);
        circuit.operations = vec![u(0, std::f64::consts::PI, 0.0, 0.0); 2];
        assert!(PostBlockFusionPass
            .run(&circuit, &mut PropertySet::new())
            .operations
            .is_empty());
        circuit.operations = vec![u(0, 4e-11, 0.0, 0.0); 100];
        let output = PostBlockFusionPass.run(&circuit, &mut PropertySet::new());
        assert_eq!(output.operations.len(), 1);
        let Operation::Gate { params, .. } = &output.operations[0] else {
            unreachable!()
        };
        assert!((params[0] - 4e-9).abs() < 1e-22);
    }

    #[test]
    fn exhausted_budget_retains_the_entire_original_run() {
        let operations = vec![u(0, 0.3, -0.2, 0.7), u(0, -0.6, 0.5, -0.1)];
        let mut output: Vec<_> = operations.iter().cloned().map(Some).collect();
        let mut pending = vec![Some(LocalRun {
            indices: vec![0, 1],
            matrix: multiply(
                &analytic_u(&[-0.6, 0.5, -0.1]).unwrap(),
                &analytic_u(&[0.3, -0.2, 0.7]).unwrap(),
            ),
        })];
        let mut stats = PostBlockFusionStats {
            residual_sum: CIRCUIT_RESIDUAL_BUDGET,
            ..PostBlockFusionStats::default()
        };
        flush(0, &mut pending, &mut output, &mut stats);
        assert_eq!(output.into_iter().flatten().collect::<Vec<_>>(), operations);
        assert_eq!(stats.rejected_budget, 1);
        assert_eq!(stats.accepted, 0);
    }
}
