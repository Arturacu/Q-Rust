//! Iterative SABRE-based layout discovery.
//!
//! Leverages BeamSABRE-style routing with cheap greedy scoring to discover
//! candidate low-cost initial qubit placements via repeated forward/backward
//! passes. When the circuit's interaction graph is a subgraph of the hardware
//! topology, this often converges to a zero-SWAP layout.

use crate::backend::Backend;
#[allow(unused_imports)]
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use crate::transpiler::routing::{
    build_gate_graph, extended_layer, reverse_gate_graph, Layout, TwoQGate,
};
#[cfg(test)]
use std::collections::HashSet;

fn pred_counts(gates: &[TwoQGate]) -> Vec<usize> {
    let mut counts = vec![0usize; gates.len()];
    for g in gates {
        for &s in &g.successors {
            counts[s] += 1;
        }
    }
    counts
}

const SWAP_COST: f64 = 3.0;

fn greedy_sabre_cost(
    gates: &[TwoQGate],
    preds: &[usize],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
    lookahead_size: usize,
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
    let neighbors: Vec<_> = (0..backend.num_qubits)
        .map(|p| backend.neighbors(p))
        .collect();

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

        let extended = extended_layer(&front, gates, &remaining, lookahead_size, |g| &g.successors);
        let mut candidates = Vec::new();
        for &g in &front {
            let p0 = layout.l2p(gates[g].logical_qubits[0]);
            let p1 = layout.l2p(gates[g].logical_qubits[1]);
            for &n in &neighbors[p0] {
                let edge = if p0 < n { (p0, n) } else { (n, p0) };
                candidates.push(edge);
            }
            for &n in &neighbors[p1] {
                let edge = if p1 < n { (p1, n) } else { (n, p1) };
                candidates.push(edge);
            }
        }

        candidates.sort_unstable();
        candidates.dedup();
        let best_swap = candidates
            .iter()
            .map(|&swap| {
                let delta = relative_score(
                    &layout,
                    swap,
                    &front,
                    gates,
                    dist,
                    &extended,
                    lookahead_size > 0,
                );
                (delta, swap)
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then_with(|| a.1.cmp(&b.1)))
            .map(|(_, s)| s);
        let Some(best_swap) = best_swap else {
            // A disconnected candidate placement is not a valid layout.
            // Let other trials compete; the checked router reports failure
            // if no complete placement is found.
            return (f64::INFINITY, layout);
        };

        layout.swap_physical(best_swap.0, best_swap.1);
        cost += SWAP_COST;
        swap_count += 1;
        if swap_count > max_swaps {
            return (f64::INFINITY, layout);
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
    extended: &[usize],
    normalize: bool,
) -> f64 {
    let (pa, pb) = swap;
    let mut delta = 0.0;

    for &g in front {
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
        delta += dist[np0][np1] as f64 - old;
    }

    if normalize {
        delta /= front.len().max(1) as f64;
    }
    let lookahead_weight = if normalize {
        0.5 / extended.len().max(1) as f64
    } else {
        0.5
    };
    for &g in extended {
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

    fn try_run(
        &self,
        circuit: &Circuit,
        property_set: &mut PropertySet,
    ) -> crate::error::Result<Circuit> {
        crate::transpiler::validate_circuit_structure(circuit)?;
        if self.num_trials == 0 || self.num_iterations == 0 {
            return Err(crate::error::QRustError::InvalidConfig(
                "SABRE layout requires at least one trial and iteration".into(),
            ));
        }
        if circuit.num_qubits > self.backend.num_qubits {
            return Err(crate::error::QRustError::InsufficientQubits {
                circuit: circuit.num_qubits,
                backend: self.backend.num_qubits,
            });
        }
        property_set.remove::<String>("sabre_layout_error");
        Ok(self.discover(circuit, property_set))
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        match self.try_run(circuit, property_set) {
            Ok(output) => output,
            Err(error) => {
                property_set.remove::<Vec<usize>>("sabre_initial_layout");
                property_set.remove::<usize>("layout_route_evaluations");
                property_set.insert("sabre_layout_error", error.to_string());
                circuit.clone()
            }
        }
    }
}

impl SabreLayoutPass {
    fn discover(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
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
        let connectivity = crate::transpiler::routing::undirected_backend(&self.backend);
        let dist = connectivity.shortest_path_matrix();
        let commute_diagonal = property_set
            .get::<bool>("routing_commute_diagonal")
            .copied()
            .unwrap_or(false);
        let (fwd_gates, _, _) = build_gate_graph(circuit, commute_diagonal);
        let fwd_preds = pred_counts(&fwd_gates);
        let rev_gates = reverse_gate_graph(&fwd_gates);
        let rev_preds = pred_counts(&rev_gates);

        let lookahead_size = property_set
            .get::<usize>("routing_lookahead_size")
            .copied()
            .unwrap_or(20);
        let trivial = Layout::trivial(num_logical, num_physical);
        let (trivial_cost, _) = greedy_sabre_cost(
            &fwd_gates,
            &fwd_preds,
            &connectivity,
            &dist,
            &trivial,
            lookahead_size,
        );
        let mut evaluations = 1usize;
        let mut best_cost = trivial_cost;
        let mut best_layout = trivial.clone();

        for trial in 0..self.num_trials {
            // Zero is a proven lower bound on SWAP count. Searching another
            // placement cannot improve it (and used to waste hundreds of routes).
            if best_cost == 0.0 {
                break;
            }
            let mut layout = if trial == 0 {
                trivial.clone()
            } else {
                random_layout(num_logical, num_physical, trial as u64)
            };

            for _ in 0..self.num_iterations {
                evaluations += 2;
                let (_, fwd_final) = greedy_sabre_cost(
                    &fwd_gates,
                    &fwd_preds,
                    &connectivity,
                    &dist,
                    &layout,
                    lookahead_size,
                );
                let (_, bwd_final) = greedy_sabre_cost(
                    &rev_gates,
                    &rev_preds,
                    &connectivity,
                    &dist,
                    &fwd_final,
                    lookahead_size,
                );
                let converged = layout.logical_to_physical == bwd_final.logical_to_physical;
                layout = bwd_final;
                if converged {
                    break;
                }
            }

            evaluations += 1;
            let (final_cost, _) = greedy_sabre_cost(
                &fwd_gates,
                &fwd_preds,
                &connectivity,
                &dist,
                &layout,
                lookahead_size,
            );
            if final_cost < best_cost {
                best_cost = final_cost;
                best_layout = layout;
            }
        }

        property_set.insert("layout_route_evaluations", evaluations);
        property_set.insert(
            "sabre_initial_layout",
            best_layout.logical_to_physical.clone(),
        );
        circuit.clone()
    }
}

/// Fisher-Yates shuffle followed by truncation produces an injective mapping
/// into `[0, num_physical)`. The `expect` is appropriate here — a panic
/// indicates a Fisher-Yates correctness bug, not user error.
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
    fn direct_layout_run_validates_local_operands_and_search_dimensions() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![99],
            params: vec![0.2, 0.3, 0.4],
        });
        for pass in [
            SabreLayoutPass {
                backend: Backend::linear(2),
                num_trials: 1,
                num_iterations: 1,
            },
            SabreLayoutPass {
                backend: Backend::linear(2),
                num_trials: 0,
                num_iterations: 1,
            },
            SabreLayoutPass {
                backend: Backend::linear(2),
                num_trials: 1,
                num_iterations: 0,
            },
            SabreLayoutPass {
                backend: Backend::linear(1),
                num_trials: 1,
                num_iterations: 1,
            },
        ] {
            let mut properties = PropertySet::new();
            let output = pass.run(&circuit, &mut properties);
            assert_eq!(output.operations, circuit.operations);
            assert!(properties.get::<String>("sabre_layout_error").is_some());
            assert!(pass.try_run(&circuit, &mut properties).is_err());
        }
        circuit.operations.pop();
        for (trials, iterations, width) in [(0, 1, 2), (1, 0, 2), (1, 1, 1)] {
            let pass = SabreLayoutPass {
                backend: Backend::linear(width),
                num_trials: trials,
                num_iterations: iterations,
            };
            let mut properties = PropertySet::new();
            assert_eq!(
                pass.run(&circuit, &mut properties).operations,
                circuit.operations
            );
            assert!(properties.get::<String>("sabre_layout_error").is_some());
        }
    }

    #[test]
    fn zero_swap_lower_bound_stops_all_remaining_layout_trials() {
        let mut circuit = Circuit::new(8, 0);
        for q in 0..7 {
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![q, q + 1],
                params: vec![],
            });
        }
        let pass = SabreLayoutPass {
            backend: Backend::linear(8),
            num_trials: 50,
            num_iterations: 5,
        };
        let mut properties = PropertySet::new();
        pass.try_run(&circuit, &mut properties).unwrap();
        assert_eq!(
            properties.get::<usize>("layout_route_evaluations"),
            Some(&1)
        );
        assert_eq!(
            properties.get::<Vec<usize>>("sabre_initial_layout"),
            Some(&(0..8).collect())
        );
    }

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
