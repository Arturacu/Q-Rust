//! Iterative SABRE-based layout discovery.
//!
//! Leverages BeamSABRE-style routing with cheap greedy scoring to discover
//! near-optimal initial qubit placements via repeated forward/backward
//! passes. When the circuit's interaction graph is a subgraph of the hardware
//! topology, this often converges to a zero-SWAP layout.

use crate::backend::Backend;
use crate::ir::{Circuit, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use crate::transpiler::routing::Layout;
use std::collections::HashSet;

#[derive(Clone, Debug)]
struct TwoQGate {
    logical_qubits: [usize; 2],
    successors: Vec<usize>,
}

fn build_2q_graph(circuit: &Circuit) -> Vec<TwoQGate> {
    let mut gates: Vec<TwoQGate> = Vec::new();
    let mut last_on_qubit: Vec<Option<usize>> = vec![None; circuit.num_qubits];
    for op in &circuit.operations {
        if let Operation::Gate { qubits, .. } = op {
            if qubits.len() >= 2 {
                let idx = gates.len();
                gates.push(TwoQGate {
                    logical_qubits: [qubits[0], qubits[1]],
                    successors: Vec::new(),
                });
                for &q in &qubits[..2] {
                    if let Some(prev) = last_on_qubit[q] {
                        gates[prev].successors.push(idx);
                    }
                    last_on_qubit[q] = Some(idx);
                }
            }
        }
    }
    gates
}

fn pred_counts(gates: &[TwoQGate]) -> Vec<u16> {
    let mut counts = vec![0u16; gates.len()];
    for g in gates {
        for &s in &g.successors {
            counts[s] += 1;
        }
    }
    counts
}

fn reverse_graph(gates: &[TwoQGate], num_qubits: usize) -> Vec<TwoQGate> {
    let mut reversed: Vec<TwoQGate> = Vec::new();
    let mut last_on_qubit: Vec<Option<usize>> = vec![None; num_qubits];
    for gate in gates.iter().rev() {
        let idx = reversed.len();
        reversed.push(TwoQGate {
            logical_qubits: gate.logical_qubits,
            successors: Vec::new(),
        });
        for &q in &gate.logical_qubits {
            if let Some(prev) = last_on_qubit[q] {
                reversed[prev].successors.push(idx);
            }
            last_on_qubit[q] = Some(idx);
        }
    }
    reversed
}

const SWAP_COST: f64 = 3.0;

fn greedy_sabre_cost(
    gates: &[TwoQGate],
    preds: &[u16],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
) -> (f64, Layout) {
    if gates.is_empty() {
        return (0.0, initial_layout.clone());
    }

    let mut layout = initial_layout.clone();
    let mut remaining = preds.to_vec();
    let mut front: Vec<usize> = (0..gates.len()).filter(|&i| remaining[i] == 0).collect();
    let mut cost = 0.0;
    let mut executed = 0usize;
    let total = gates.len();

    let max_swaps = total * dist.len() * 3;
    let mut swap_count = 0usize;

    while executed < total {
        loop {
            let routable: Vec<usize> = front
                .iter()
                .filter(|&&g| {
                    let p0 = layout.l2p(gates[g].logical_qubits[0]);
                    let p1 = layout.l2p(gates[g].logical_qubits[1]);
                    dist[p0][p1] == 1
                })
                .cloned()
                .collect();
            if routable.is_empty() {
                break;
            }
            for g in routable {
                front.retain(|&x| x != g);
                for &s in &gates[g].successors {
                    remaining[s] -= 1;
                    if remaining[s] == 0 {
                        front.push(s);
                    }
                }
                executed += 1;
            }
        }

        if executed >= total {
            break;
        }

        let mut candidates: HashSet<(usize, usize)> = HashSet::new();
        for &g in &front {
            let p0 = layout.l2p(gates[g].logical_qubits[0]);
            let p1 = layout.l2p(gates[g].logical_qubits[1]);
            for n in backend.neighbors(p0) {
                let edge = if p0 < n { (p0, n) } else { (n, p0) };
                candidates.insert(edge);
            }
            for n in backend.neighbors(p1) {
                let edge = if p1 < n { (p1, n) } else { (n, p1) };
                candidates.insert(edge);
            }
        }

        let best_swap = candidates
            .iter()
            .map(|&swap| {
                let delta = relative_score(&layout, swap, &front, gates, dist);
                (delta, swap)
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .map(|(_, s)| s)
            .expect("no SWAP candidates on a connected graph");

        layout.swap_physical(best_swap.0, best_swap.1);
        cost += SWAP_COST;
        swap_count += 1;
        if swap_count > max_swaps {
            break;
        }
    }
    (cost, layout)
}

fn relative_score(
    layout: &Layout,
    swap: (usize, usize),
    front: &[usize],
    gates: &[TwoQGate],
    dist: &[Vec<usize>],
) -> f64 {
    let (pa, pb) = swap;
    let mut delta = 0.0;
    let mut extended_layer = HashSet::new();

    for &g in front {
        let p0 = layout.l2p(gates[g].logical_qubits[0]);
        let p1 = layout.l2p(gates[g].logical_qubits[1]);
        for &succ in &gates[g].successors {
            extended_layer.insert(succ);
        }
        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }
        let old = dist[p0][p1] as f64;
        let np0 = if p0 == pa {
            pb
        } else if p0 == pb {
            pa
        } else {
            p0
        };
        let np1 = if p1 == pa {
            pb
        } else if p1 == pb {
            pa
        } else {
            p1
        };
        delta += dist[np0][np1] as f64 - old;
    }

    let lookahead_weight = 0.5;
    for g in extended_layer {
        let p0 = layout.l2p(gates[g].logical_qubits[0]);
        let p1 = layout.l2p(gates[g].logical_qubits[1]);
        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }
        let old = dist[p0][p1] as f64;
        let np0 = if p0 == pa {
            pb
        } else if p0 == pb {
            pa
        } else {
            p0
        };
        let np1 = if p1 == pa {
            pb
        } else if p1 == pb {
            pa
        } else {
            p1
        };
        delta += (dist[np0][np1] as f64 - old) * lookahead_weight;
    }

    delta
}

/// Seeded xorshift64 PRNG.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(if seed == 0 { 0xDEAD_BEEF } else { seed })
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

/// Layout discovery via iterative SABRE.
#[derive(Debug, Clone)]
pub struct SabreLayoutPass {
    pub backend: Backend,
    pub num_trials: usize,
    pub num_iterations: usize,
}

impl Pass for SabreLayoutPass {
    fn name(&self) -> &str {
        "SabreLayoutPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        let has_2q = circuit.operations.iter().any(|op| {
            if let Operation::Gate { qubits, .. } = op {
                qubits.len() >= 2
            } else {
                false
            }
        });
        if !has_2q {
            return circuit.clone();
        }

        assert!(self.num_trials >= 1, "num_trials must be >= 1");

        let num_logical = circuit.num_qubits;
        let num_physical = self.backend.num_qubits;
        let dist = self.backend.shortest_path_matrix();
        let fwd_gates = build_2q_graph(circuit);
        let fwd_preds = pred_counts(&fwd_gates);
        let rev_gates = reverse_graph(&fwd_gates, num_logical);
        let rev_preds = pred_counts(&rev_gates);

        let trivial = Layout::trivial(num_logical, num_physical);
        let (trivial_cost, _) =
            greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &trivial);
        let mut best_cost = trivial_cost;
        let mut best_layout = trivial.clone();

        for trial in 0..self.num_trials {
            let mut layout = if trial == 0 {
                trivial.clone()
            } else {
                random_layout(num_logical, num_physical, trial as u64)
            };

            for _ in 0..self.num_iterations {
                let (_, fwd_final) =
                    greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &layout);
                let (_, bwd_final) =
                    greedy_sabre_cost(&rev_gates, &rev_preds, &self.backend, &dist, &fwd_final);
                layout = bwd_final;
            }

            let (final_cost, _) =
                greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &layout);
            if final_cost < best_cost {
                best_cost = final_cost;
                best_layout = layout;
            }
        }

        property_set.insert(
            "sabre_initial_layout",
            best_layout.logical_to_physical.clone(),
        );
        circuit.clone()
    }
}

/// Loop 5 §AF-1: migrated to use [`Layout::from_l2p`] for defense-in-depth.
/// Fisher-Yates shuffle followed by truncation produces an injective
/// mapping into `[0, num_physical)`; the validating constructor
/// re-checks this invariant. The `expect` is appropriate — a panic
/// here indicates a Fisher-Yates correctness bug, not user error.
fn random_layout(num_logical: usize, num_physical: usize, seed: u64) -> Layout {
    let mut rng = Rng::new(seed.wrapping_mul(0x517C_C1B7_2722_0A95).wrapping_add(1));
    let mut positions: Vec<usize> = (0..num_physical).collect();
    rng.shuffle(&mut positions);
    let l2p: Vec<usize> = positions.into_iter().take(num_logical).collect();
    Layout::from_l2p(l2p, num_physical)
        .expect("random_layout: Fisher-Yates output must be injective")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_layout_valid() {
        let layout = random_layout(5, 8, 42);
        let mut seen = HashSet::new();
        for i in 0..5 {
            let p = layout.l2p(i);
            assert!(p < 8);
            assert!(seen.insert(p));
            assert_eq!(layout.physical_to_logical[p], i);
        }
    }

    /// Loop 5 §AF-1: a regression net for the migrated test helper
    /// below — confirms the validating constructor accepts the
    /// shuffle output across a range of (logical, physical) sizes.
    #[test]
    fn test_random_layout_validates_across_sizes() {
        for (l, p, seed) in [(1, 4, 1), (4, 4, 7), (3, 16, 99), (8, 12, 2024)] {
            let layout = random_layout(l, p, seed);
            assert_eq!(layout.logical_to_physical.len(), l);
            assert_eq!(layout.physical_to_logical.len(), p);
        }
    }

    #[test]
    fn test_layout_pass_improves_chain() {
        let mut c = Circuit::new(5, 0);
        for i in (1..5).rev() {
            c.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![i, i - 1],
                params: vec![],
            });
        }
        let backend = Backend::linear(5);
        let pass = SabreLayoutPass {
            backend: backend.clone(),
            num_trials: 20,
            num_iterations: 5,
        };
        let mut ps = PropertySet::new();
        let _ = pass.run(&c, &mut ps);
        let discovered: &Vec<usize> = ps.get("sabre_initial_layout").unwrap();
        // Loop 5 §AF-1: migrated from hand-rolled p2l construction
        // to the validating `Layout::from_l2p` constructor. The
        // `unwrap` is justified — `discovered` was just produced by
        // SabreLayoutPass and is guaranteed injective by construction.
        let layout = Layout::from_l2p(discovered.clone(), backend.num_qubits)
            .expect("SabreLayoutPass output must be a valid layout");
        let dist = backend.shortest_path_matrix();
        for i in (1..5).rev() {
            let p0 = layout.l2p(i);
            let p1 = layout.l2p(i - 1);
            assert_eq!(dist[p0][p1], 1);
        }
    }

    #[test]
    fn test_layout_no_2q_gates() {
        let mut c = Circuit::new(3, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let mut ps = PropertySet::new();
        let pass = SabreLayoutPass {
            backend: Backend::linear(3),
            num_trials: 10,
            num_iterations: 2,
        };
        let _ = pass.run(&c, &mut ps);
        assert!(ps.get::<Vec<usize>>("sabre_initial_layout").is_none());
    }

    #[test]
    fn test_layout_already_optimal() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let mut ps = PropertySet::new();
        let pass = SabreLayoutPass {
            backend: Backend::linear(2),
            num_trials: 5,
            num_iterations: 2,
        };
        let _ = pass.run(&c, &mut ps);
        let layout = ps.get::<Vec<usize>>("sabre_initial_layout").unwrap();
        let dist = Backend::linear(2).shortest_path_matrix();
        assert_eq!(dist[layout[0]][layout[1]], 1);
    }
}