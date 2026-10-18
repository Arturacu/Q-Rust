//! Structural certificates and final-output costs for retained routing candidates.

use crate::backend::Backend;
use crate::ir::{Circuit, GateType, Operation};

/// Certify that retaining a logical SWAP sequence needs no routing SWAPs in
/// the identity placement. This is a sufficient condition, not a search for
/// all possible zero-routing placements. Direction is repaired after lowering.
///
/// A terminal-only permutation is already handled by trailing-SWAP elision;
/// require an interior SWAP so a second lowering has a distinct opportunity.
pub(crate) fn has_connected_interior_swaps(circuit: &Circuit, backend: &Backend) -> bool {
    // Most inputs have no explicit logical SWAP. Avoid a second structural
    // validation and adjacency scan when there is no candidate to retain.
    if !circuit.operations.iter().any(|op| {
        matches!(
            op,
            Operation::Gate {
                name: GateType::SWAP,
                ..
            }
        )
    }) {
        return false;
    }
    if circuit.num_qubits > backend.num_qubits
        || super::validate_circuit_structure(circuit).is_err()
    {
        return false;
    }
    for op in &circuit.operations {
        match op {
            Operation::Gate { name, qubits, .. } => {
                if matches!(name, GateType::Custom(_))
                    || qubits.len() > 2
                    || (qubits.len() == 2 && !backend.is_adjacent(qubits[0], qubits[1]))
                {
                    return false;
                }
            }
            Operation::Barrier { .. } => {}
            _ => return false,
        }
    }
    let mut touched = vec![false; circuit.num_qubits];
    for op in circuit.operations.iter().rev() {
        match op {
            Operation::Gate {
                name: GateType::SWAP,
                qubits,
                ..
            } => {
                if touched[qubits[0]] || touched[qubits[1]] {
                    return true;
                }
                // An entirely trailing SWAP does not make earlier SWAPs
                // interior: the whole trailing permutation can be elided.
            }
            Operation::Gate { qubits, .. } => {
                for &q in qubits {
                    touched[q] = true;
                }
            }
            // Match the conservative scheduling boundary in trailing elision.
            Operation::Barrier { .. } => touched.fill(true),
            _ => unreachable!("nonunitary operations excluded above"),
        }
    }
    false
}

/// Count entangling gates, total depth, and entangling-only depth. The final
/// common basis may use CX or CZ; choosing by CX spelling alone is incorrect.
fn output_cost(circuit: &Circuit) -> Option<[usize; 3]> {
    if super::validate_circuit_structure(circuit).is_err()
        || super::validate_static_semantics(circuit).is_err()
    {
        return None;
    }
    let mut two_qubit_count = 0;
    let mut two_depth = vec![0usize; circuit.num_qubits];
    for op in &circuit.operations {
        match op {
            Operation::Gate { name, qubits, .. } => {
                if matches!(name, GateType::Custom(_)) || qubits.len() > 2 {
                    return None;
                }
                if qubits.len() == 2 {
                    two_qubit_count += 1;
                    let next = 1 + two_depth[qubits[0]].max(two_depth[qubits[1]]);
                    two_depth[qubits[0]] = next;
                    two_depth[qubits[1]] = next;
                }
            }
            Operation::Barrier { qubits } => {
                let next = qubits.iter().map(|&q| two_depth[q]).max().unwrap_or(0);
                for &q in qubits {
                    two_depth[q] = next;
                }
            }
            // Terminal readout contributes to Circuit::depth(), but is not
            // an entangling gate and cannot add entangling-only depth.
            Operation::Measure { .. } => {}
            _ => return None,
        }
    }
    Some([
        two_qubit_count,
        circuit.depth(),
        two_depth.into_iter().max().unwrap_or(0),
    ])
}

fn prefer_cost(candidate: [usize; 3], current: [usize; 3], incumbent: [usize; 3]) -> bool {
    candidate.iter().zip(incumbent).all(|(a, b)| *a <= b) && candidate < current
}

/// Use the original canonical output as a fixed nonregression reference.
/// Every chosen output must be no worse on any of the three metrics. Among
/// admissible alternatives prefer fewer entanglers, then lower total depth,
/// then lower entangling depth. Exact ties retain the current circuit.
///
/// The reference must not be replaced by an intermediate winner: two outputs
/// may each improve the canonical circuit while trading depth against count.
pub(crate) fn prefer_output(candidate: &Circuit, current: &Circuit, incumbent: &Circuit) -> bool {
    let (Some(candidate), Some(current), Some(incumbent)) = (
        output_cost(candidate),
        output_cost(current),
        output_cost(incumbent),
    ) else {
        return false;
    };
    prefer_cost(candidate, current, incumbent)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gate(name: GateType, wires: &[usize]) -> Operation {
        Operation::Gate {
            name,
            qubits: wires.to_vec(),
            params: vec![],
        }
    }

    fn circuit(width: usize, ops: Vec<Operation>) -> Circuit {
        let mut result = Circuit::new(width, 0);
        result.operations = ops;
        result
    }

    #[test]
    fn local_certificate_requires_an_interior_swap_and_checks_every_interaction() {
        let backend = Backend::linear(4);
        let local = circuit(
            3,
            vec![gate(GateType::SWAP, &[0, 1]), gate(GateType::CX, &[1, 2])],
        );
        assert!(has_connected_interior_swaps(&local, &backend));
        let mut directed = Backend::new("one-way line", 4);
        directed.set_coupling_map([(0, 1), (1, 2), (2, 3)]);
        assert!(has_connected_interior_swaps(&local, &directed));
        for ops in [
            vec![gate(GateType::CX, &[0, 1])],
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::SWAP, &[0, 1])],
            vec![gate(GateType::SWAP, &[0, 1]), gate(GateType::SWAP, &[1, 2])],
            vec![gate(GateType::SWAP, &[0, 1]), gate(GateType::H, &[2])],
            vec![gate(GateType::SWAP, &[0, 1]), gate(GateType::CX, &[0, 2])],
            vec![
                gate(GateType::CX, &[0, 2]),
                gate(GateType::SWAP, &[0, 1]),
                gate(GateType::H, &[0]),
            ],
        ] {
            assert!(!has_connected_interior_swaps(&circuit(3, ops), &backend));
        }
        let blocked = circuit(
            3,
            vec![
                gate(GateType::SWAP, &[0, 1]),
                Operation::Barrier { qubits: vec![2] },
            ],
        );
        assert!(has_connected_interior_swaps(&blocked, &backend));
    }

    #[test]
    fn certificate_fails_closed_for_invalid_unsupported_or_disconnected_inputs() {
        let backend = Backend::linear(3);
        let local = vec![gate(GateType::SWAP, &[0, 1]), gate(GateType::H, &[0])];
        assert!(!has_connected_interior_swaps(
            &circuit(4, local.clone()),
            &backend
        ));
        assert!(!has_connected_interior_swaps(
            &circuit(3, local.clone()),
            &Backend::new("isolated", 3)
        ));
        for op in [
            gate(GateType::SWAP, &[0]),
            gate(GateType::CX, &[0, 99]),
            gate(GateType::SWAP, &[1, 1]),
            gate(GateType::CCX, &[0, 1, 2]),
            gate(GateType::Custom("missing".into()), &[0, 1]),
            Operation::Reset { qubit: 0 },
            Operation::Measure { qubit: 0, cbit: 0 },
            Operation::Gate {
                name: GateType::RZ,
                qubits: vec![0],
                params: vec![f64::NAN],
            },
        ] {
            let mut ops = local.clone();
            ops.push(op);
            assert!(!has_connected_interior_swaps(&circuit(3, ops), &backend));
        }
    }

    #[test]
    fn selection_requires_a_strict_pareto_improvement_and_counts_cz() {
        let cx = circuit(
            4,
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::CX, &[2, 3])],
        );
        let cz = circuit(
            4,
            vec![gate(GateType::CZ, &[0, 1]), gate(GateType::CZ, &[2, 3])],
        );
        assert!(!prefer_output(&cx, &cx, &cx));
        assert!(!prefer_output(&cz, &cx, &cx));
        assert!(!prefer_output(&cx, &cz, &cz));
        let fewer = circuit(4, vec![gate(GateType::CZ, &[0, 1])]);
        assert!(prefer_output(&fewer, &cx, &cx));
        let deeper = circuit(
            4,
            vec![
                gate(GateType::CX, &[0, 1]),
                gate(GateType::H, &[0]),
                gate(GateType::H, &[0]),
            ],
        );
        assert!(
            !prefer_output(&deeper, &cx, &cx),
            "fewer entanglers cannot pay for deeper output"
        );
        let more = circuit(
            4,
            vec![
                gate(GateType::CX, &[0, 1]),
                gate(GateType::CX, &[2, 3]),
                gate(GateType::CX, &[0, 1]),
            ],
        );
        assert!(!prefer_output(&more, &cx, &cx));
        let malformed = circuit(4, vec![gate(GateType::CX, &[0, 99])]);
        assert!(!prefer_output(&malformed, &cx, &cx));
        assert!(!prefer_output(&cx, &malformed, &malformed));
    }

    #[test]
    fn two_qubit_depth_is_independent_of_local_gate_depth_and_respects_barriers() {
        let parallel = circuit(
            4,
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::CX, &[2, 3])],
        );
        let mut padded = parallel.clone();
        for _ in 0..4 {
            padded.add_op(gate(GateType::H, &[0]));
        }
        let serial = circuit(
            4,
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::CX, &[1, 2])],
        );
        assert_eq!(output_cost(&padded), Some([2, 5, 1]));
        assert_eq!(output_cost(&serial), Some([2, 2, 2]));
        assert!(
            !prefer_output(&serial, &padded, &padded),
            "less total depth cannot pay for more entangling depth"
        );
        let with_barrier = circuit(
            4,
            vec![
                gate(GateType::CX, &[0, 1]),
                Operation::Barrier {
                    qubits: vec![0, 1, 2, 3],
                },
                gate(GateType::CX, &[2, 3]),
            ],
        );
        assert_eq!(output_cost(&parallel), Some([2, 1, 1]));
        assert_eq!(output_cost(&with_barrier), Some([2, 2, 2]));
        assert!(prefer_output(&parallel, &with_barrier, &with_barrier));
    }

    #[test]
    fn terminal_readout_counts_in_total_depth_but_not_entangling_depth() {
        let mut measured = circuit(
            4,
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::CX, &[2, 3])],
        );
        measured.num_cbits = 1;
        measured.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        assert_eq!(output_cost(&measured), Some([2, 2, 1]));
        let mut fewer_but_more_readout_depth = circuit(4, vec![gate(GateType::CX, &[0, 1])]);
        fewer_but_more_readout_depth.num_cbits = 1;
        for _ in 0..2 {
            fewer_but_more_readout_depth.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        }
        assert_eq!(output_cost(&fewer_but_more_readout_depth), Some([1, 3, 1]));
        assert!(!prefer_output(
            &fewer_but_more_readout_depth,
            &measured,
            &measured
        ));
        let mut dynamic = measured.clone();
        dynamic.add_op(gate(GateType::H, &[2]));
        assert_eq!(output_cost(&dynamic), None);
        measured.add_op(Operation::Reset { qubit: 0 });
        assert_eq!(output_cost(&measured), None);
    }

    #[test]
    fn fixed_reference_keeps_incomparable_improvements_from_blocking_lower_count() {
        let incumbent = [29, 40, 24];
        let depth_candidate = [29, 32, 19];
        let count_candidate = [25, 36, 19];
        assert!(prefer_cost(depth_candidate, incumbent, incumbent));
        assert!(prefer_cost(count_candidate, depth_candidate, incumbent));
        assert!(prefer_cost(count_candidate, incumbent, incumbent));
        assert!(!prefer_cost(depth_candidate, count_candidate, incumbent));
        assert!(!prefer_cost([24, 41, 18], count_candidate, incumbent));
        assert!(!prefer_cost([24, 36, 25], count_candidate, incumbent));
        assert!(!prefer_cost(count_candidate, count_candidate, incumbent));
        for candidates in [
            [depth_candidate, count_candidate],
            [count_candidate, depth_candidate],
        ] {
            let mut winner = incumbent;
            for candidate in candidates {
                if prefer_cost(candidate, winner, incumbent) {
                    winner = candidate;
                }
            }
            assert_eq!(winner, count_candidate);
        }
    }

    #[test]
    fn public_cost_selection_preserves_stable_ties_and_nonregression_reference() {
        let incumbent = circuit(
            3,
            vec![gate(GateType::CX, &[0, 1]), gate(GateType::CX, &[1, 2])],
        );
        let candidate = circuit(3, vec![gate(GateType::CZ, &[0, 1])]);
        assert!(prefer_output(&candidate, &incumbent, &incumbent));
        assert!(!prefer_output(&candidate, &candidate, &incumbent));
        let tradeoff = circuit(
            3,
            vec![
                gate(GateType::H, &[0]),
                gate(GateType::H, &[0]),
                gate(GateType::H, &[0]),
            ],
        );
        assert!(!prefer_output(&tradeoff, &candidate, &incumbent));
    }
}
