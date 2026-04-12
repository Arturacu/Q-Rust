//! SABRE-Iterative Layout: discovers near-optimal initial qubit placements.
//!
//! Instead of solving the NP-Complete subgraph-isomorphism problem directly,
//! this pass exploits our existing BeamSABRE router to *organically* discover
//! good initial layouts.
//!
//! **Algorithm:**
//! 1. Generate `num_trials` random initial layouts (permutations).
//! 2. For each trial, run `num_iterations` forward/backward routing passes
//!    using a lightweight greedy (beam_width=1) router.  Each backward pass's
//!    final layout becomes the next forward pass's initial layout, letting the
//!    qubits "settle" into positions that minimise routing cost.
//! 3. Pick the trial whose final forward pass produced the fewest SWAPs.
//! 4. Store the winning layout in `PropertySet["initial_layout"]` so the
//!    downstream `BeamSabrePass` can use it instead of the trivial mapping.
//!
//! When the circuit's interaction graph is a subgraph of the hardware topology
//! (e.g. a linear chain on a Heavy-Hex device), this procedure converges to a
//! **zero-SWAP** layout — matching VF2 without any graph-isomorphism solver.

use crate::backend::Backend;
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use crate::transpiler::routing::Layout;
use std::collections::HashSet;

// ─── Lightweight gate dependency graph (mirrors routing.rs, kept private) ───

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

// ─── Lightweight greedy SABRE (beam_width=1, no circuit reconstruction) ────

const SWAP_COST: f64 = 3.0;

/// Runs greedy SABRE on `gates` with the given `layout`, returning
/// (total_swap_cost, final_layout).  This is intentionally stripped down
/// to be as fast as possible — we only care about the cost and the
/// resulting qubit permutation, not the action log.
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
    let mut front: Vec<usize> = (0..gates.len())
        .filter(|&i| remaining[i] == 0)
        .collect();
    let mut cost = 0.0;
    let mut executed = 0usize;
    let total = gates.len();

    // Upper bound on SWAPs to prevent infinite loops in degenerate cases
    let max_swaps = total * dist.len() * 3;
    let mut swap_count = 0usize;

    while executed < total {
        // Execute all routable gates
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

        // Pick the best SWAP (greedy — single candidate)
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
            .expect("No SWAP candidates on connected graph");

        layout.swap_physical(best_swap.0, best_swap.1);
        cost += SWAP_COST;
        swap_count += 1;
        if swap_count > max_swaps {
            break;
        }
    }

    (cost, layout)
}

/// Relative score (mirrors routing.rs scoring).
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
        let np0 = if p0 == pa { pb } else if p0 == pb { pa } else { p0 };
        let np1 = if p1 == pa { pb } else if p1 == pb { pa } else { p1 };
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
        let np0 = if p0 == pa { pb } else if p0 == pb { pa } else { p0 };
        let np1 = if p1 == pa { pb } else if p1 == pb { pa } else { p1 };
        delta += (dist[np0][np1] as f64 - old) * lookahead_weight;
    }
    
    delta
}

// ─── Deterministic seeded shuffle (no external crate needed) ───────────────

/// Simple xorshift64 PRNG — fast, tiny, deterministic.
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
    /// Fisher-Yates shuffle of a mutable slice.
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

// ─── SabreLayoutPass ───────────────────────────────────────────────────────

/// Discovers a near-optimal initial qubit layout via iterative SABRE routing.
pub struct SabreLayoutPass {
    pub backend: Backend,
    /// Number of independent random starting layouts to try.
    pub num_trials: usize,
    /// Number of forward/backward iteration pairs per trial.
    pub num_iterations: usize,
}

impl Pass for SabreLayoutPass {
    fn name(&self) -> &str {
        "SabreLayoutPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        // Only meaningful when there are multi-qubit gates
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

        assert!(
            self.num_trials >= 1,
            "SabreLayoutPass: num_trials must be >= 1, got {}",
            self.num_trials
        );

        let num_logical = circuit.num_qubits;
        let num_physical = self.backend.num_qubits;

        // Pre-compute distance matrix and gate graphs
        let dist = self.backend.shortest_path_matrix();
        let fwd_gates = build_2q_graph(circuit);
        let fwd_preds = pred_counts(&fwd_gates);
        let rev_gates = reverse_graph(&fwd_gates, num_logical);
        let rev_preds = pred_counts(&rev_gates);

        // Score the trivial layout as the baseline
        let trivial = Layout::trivial(num_logical, num_physical);
        let (trivial_cost, _) = greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &trivial);
        let mut best_cost = trivial_cost;
        let mut best_layout = trivial.clone();

        // Trial 0 is always the trivial layout (already scored above).
        // Remaining trials use random permutations.
        for trial in 0..self.num_trials {
            let mut layout = if trial == 0 {
                trivial.clone()
            } else {
                random_layout(num_logical, num_physical, trial as u64)
            };

            // Forward/backward iterations let the layout converge
            for _iter in 0..self.num_iterations {
                // Forward
                let (_, fwd_final) =
                    greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &layout);

                // Backward (using forward's final layout as the starting point)
                let (_, bwd_final) =
                    greedy_sabre_cost(&rev_gates, &rev_preds, &self.backend, &dist, &fwd_final);

                // The backward pass's final layout is the next iteration's start
                layout = bwd_final;
            }

            // Score final converged layout with one last forward pass
            let (final_cost, _) =
                greedy_sabre_cost(&fwd_gates, &fwd_preds, &self.backend, &dist, &layout);

            if final_cost < best_cost {
                best_cost = final_cost;
                best_layout = layout;
            }
        }

        // Store the winning layout for the downstream BeamSabrePass
        property_set.insert(
            "sabre_initial_layout",
            best_layout.logical_to_physical.clone(),
        );

        // The layout pass does NOT modify the circuit — it only sets metadata.
        circuit.clone()
    }
}

/// Creates a random layout by shuffling logical→physical assignments.
fn random_layout(num_logical: usize, num_physical: usize, seed: u64) -> Layout {
    let mut rng = Rng::new(seed.wrapping_mul(0x517CC1B727220A95).wrapping_add(1));
    let mut positions: Vec<usize> = (0..num_physical).collect();
    rng.shuffle(&mut positions);

    let mut l2p = vec![0usize; num_logical];
    let mut p2l = vec![usize::MAX; num_physical];
    for i in 0..num_logical {
        l2p[i] = positions[i];
        p2l[positions[i]] = i;
    }
    Layout {
        logical_to_physical: l2p,
        physical_to_logical: p2l,
    }
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_layout_valid() {
        let layout = random_layout(5, 8, 42);
        // Every logical qubit maps to a distinct physical qubit
        let mut seen = HashSet::new();
        for i in 0..5 {
            let p = layout.l2p(i);
            assert!(p < 8, "Physical qubit out of range");
            assert!(seen.insert(p), "Duplicate physical assignment");
            assert_eq!(layout.physical_to_logical[p], i);
        }
    }

    #[test]
    fn test_layout_pass_improves_chain() {
        // Reverse chain CX(4,3), CX(3,2), CX(2,1), CX(1,0) on a linear(5) backend.
        // Trivial layout puts these on adjacent qubits in the WRONG order.
        // SabreLayoutPass should discover a layout that needs 0 SWAPs.
        let mut circuit = Circuit::new(5, 0);
        for i in (1..5).rev() {
            circuit.add_op(Operation::Gate {
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
        let _ = pass.run(&circuit, &mut ps);

        let discovered: &Vec<usize> = ps.get("sabre_initial_layout").unwrap();

        // Verify the layout produces 0 swaps: every CX(i, i-1) must map to
        // adjacent physical qubits.
        let layout = Layout {
            logical_to_physical: discovered.clone(),
            physical_to_logical: {
                let mut p2l = vec![usize::MAX; backend.num_qubits];
                for (l, &p) in discovered.iter().enumerate() {
                    p2l[p] = l;
                }
                p2l
            },
        };

        let dist = backend.shortest_path_matrix();
        for i in (1..5).rev() {
            let p0 = layout.l2p(i);
            let p1 = layout.l2p(i - 1);
            assert_eq!(
                dist[p0][p1], 1,
                "Logical {i} -> phys {p0} and logical {} -> phys {p1} are not adjacent",
                i - 1
            );
        }
    }

    #[test]
    fn test_layout_no_2q_gates() {
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
        
        let mut ps = PropertySet::new();
        let pass = SabreLayoutPass {
            backend: Backend::linear(3),
            num_trials: 10,
            num_iterations: 2,
        };
        let _ = pass.run(&circuit, &mut ps);
        assert!(ps.get::<Vec<usize>>("sabre_initial_layout").is_none());
    }

    #[test]
    fn test_layout_already_optimal() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
        
        let mut ps = PropertySet::new();
        let pass = SabreLayoutPass {
            backend: Backend::linear(2),
            num_trials: 5,
            num_iterations: 2,
        };
        let _ = pass.run(&circuit, &mut ps);
        
        let layout = ps.get::<Vec<usize>>("sabre_initial_layout").unwrap();
        let dist = Backend::linear(2).shortest_path_matrix();
        assert_eq!(dist[layout[0]][layout[1]], 1);
    }
}
