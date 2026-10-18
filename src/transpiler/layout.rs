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
    build_gate_graph, extended_layer, reverse_gate_graph, Layout, RefinedLayoutCandidate, TwoQGate,
};
use std::collections::HashMap;
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

/// One deterministic routing problem, scoped to a single layout-discovery call.
/// Forward and reversed graphs have separate caches: a placement alone is not
/// a key for both problems. No candidate requests or tie decisions are skipped.
struct GreedyRouteCache<'a> {
    gates: &'a [TwoQGate],
    preds: &'a [usize],
    backend: &'a Backend,
    dist: &'a [Vec<usize>],
    lookahead_size: usize,
    memoize: bool,
    routes: HashMap<Vec<usize>, (f64, Layout)>,
    evaluations: usize,
    hits: usize,
}

impl<'a> GreedyRouteCache<'a> {
    fn new(
        gates: &'a [TwoQGate],
        preds: &'a [usize],
        backend: &'a Backend,
        dist: &'a [Vec<usize>],
        lookahead_size: usize,
        memoize: bool,
    ) -> Self {
        Self {
            gates,
            preds,
            backend,
            dist,
            lookahead_size,
            memoize,
            routes: HashMap::new(),
            evaluations: 0,
            hits: 0,
        }
    }

    fn route(&mut self, initial_layout: &Layout) -> (f64, Layout) {
        if self.memoize {
            if let Some(result) = self.routes.get(&initial_layout.logical_to_physical) {
                self.hits += 1;
                return result.clone();
            }
        }
        self.evaluations += 1;
        let result = greedy_sabre_cost(
            self.gates,
            self.preds,
            self.backend,
            self.dist,
            initial_layout,
            self.lookahead_size,
        );
        if self.memoize {
            self.routes
                .insert(initial_layout.logical_to_physical.clone(), result.clone());
        }
        result
    }
}

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
                property_set.remove::<usize>("layout_route_cache_hits");
                property_set.remove::<usize>("layout_route_requests");
                property_set
                    .remove::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates");
                property_set.insert("sabre_layout_error", error.to_string());
                circuit.clone()
            }
        }
    }
}

impl SabreLayoutPass {
    fn discover(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        self.discover_with_cache(circuit, property_set, true)
    }

    // The uncached mode is an internal regression oracle for preserving the
    // candidate search exactly, not a public search-budget option.
    fn discover_with_cache(
        &self,
        circuit: &Circuit,
        property_set: &mut PropertySet,
        memoize: bool,
    ) -> Circuit {
        property_set.remove::<usize>("layout_route_evaluations");
        property_set.remove::<usize>("layout_route_cache_hits");
        property_set.remove::<usize>("layout_route_requests");
        property_set.remove::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates");
        let retain_candidates = property_set
            .get::<bool>("retain_routing_candidates")
            .copied()
            .unwrap_or(false);
        let mut refined_candidates = Vec::<RefinedLayoutCandidate>::new();
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
        let mut forward_routes = GreedyRouteCache::new(
            &fwd_gates,
            &fwd_preds,
            &connectivity,
            &dist,
            lookahead_size,
            memoize,
        );
        let mut reverse_routes = GreedyRouteCache::new(
            &rev_gates,
            &rev_preds,
            &connectivity,
            &dist,
            lookahead_size,
            memoize,
        );
        let trivial = Layout::trivial(num_logical, num_physical);
        let (trivial_cost, _) = forward_routes.route(&trivial);
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
                let (_, fwd_final) = forward_routes.route(&layout);
                let (_, bwd_final) = reverse_routes.route(&fwd_final);
                let converged = layout.logical_to_physical == bwd_final.logical_to_physical;
                layout = bwd_final;
                if converged {
                    break;
                }
            }

            let (final_cost, _) = forward_routes.route(&layout);
            if retain_candidates
                && final_cost.is_finite()
                && !refined_candidates
                    .iter()
                    .any(|candidate| candidate.layout == layout.logical_to_physical)
            {
                refined_candidates.push(RefinedLayoutCandidate {
                    layout: layout.logical_to_physical.clone(),
                    swap_cost: final_cost,
                });
            }
            if final_cost < best_cost {
                best_cost = final_cost;
                best_layout = layout;
            }
        }

        let evaluations = forward_routes.evaluations + reverse_routes.evaluations;
        let hits = forward_routes.hits + reverse_routes.hits;
        property_set.insert("layout_route_evaluations", evaluations);
        property_set.insert("layout_route_cache_hits", hits);
        property_set.insert("layout_route_requests", evaluations + hits);
        if retain_candidates {
            property_set.insert("sabre_refined_layout_candidates", refined_candidates);
        }
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

    fn interacting_circuit() -> Circuit {
        let mut circuit = Circuit::new(4, 0);
        for (a, b) in [(0, 3), (0, 1), (2, 3), (1, 2), (0, 2), (0, 3)] {
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![a, b],
                params: vec![],
            });
        }
        circuit
    }

    #[test]
    fn retaining_refined_placements_preserves_exact_search_and_scoring_work() {
        let pass = SabreLayoutPass {
            backend: Backend::linear(4),
            num_trials: 8,
            num_iterations: 4,
        };
        let circuit = interacting_circuit();
        let mut ordinary = PropertySet::new();
        let mut retained = PropertySet::new();
        retained.insert("retain_routing_candidates", true);
        assert_eq!(
            pass.try_run(&circuit, &mut ordinary).unwrap(),
            pass.try_run(&circuit, &mut retained).unwrap()
        );
        assert_eq!(
            ordinary.get::<Vec<usize>>("sabre_initial_layout"),
            retained.get::<Vec<usize>>("sabre_initial_layout")
        );
        for key in [
            "layout_route_evaluations",
            "layout_route_cache_hits",
            "layout_route_requests",
        ] {
            assert_eq!(ordinary.get::<usize>(key), retained.get::<usize>(key));
        }
        assert!(ordinary
            .get::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates")
            .is_none());
        let candidates = retained
            .get::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates")
            .unwrap();
        assert!(!candidates.is_empty());
        assert!(candidates.len() <= pass.num_trials);
        let mut seen = HashSet::new();
        for candidate in candidates {
            assert!(candidate.swap_cost.is_finite());
            assert!(seen.insert(candidate.layout.clone()));
            assert!(Layout::from_l2p(candidate.layout.clone(), 4).is_ok());
        }
        retained.insert("retain_routing_candidates", false);
        pass.try_run(&circuit, &mut retained).unwrap();
        assert!(retained
            .get::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates")
            .is_none());
        retained.insert("retain_routing_candidates", true);
        pass.try_run(&circuit, &mut retained).unwrap();
        pass.try_run(&Circuit::new(4, 0), &mut retained).unwrap();
        assert!(retained
            .get::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates")
            .is_none());
    }

    #[test]
    fn cache_preserves_search_and_routing_for_varied_connectivity() {
        let circuit = interacting_circuit();
        let mut directed = Backend::new("directed line", 4);
        directed.set_coupling_map([(0, 1), (1, 2), (2, 3)]);
        let mut disconnected = Backend::new("disconnected", 4);
        disconnected.set_coupling_map([(0, 1), (1, 0), (2, 3), (3, 2)]);
        for backend in [
            Backend::linear(4),
            directed,
            Backend::linear(6),
            disconnected,
        ] {
            for lookahead in [0_usize, 20] {
                let pass = SabreLayoutPass {
                    backend: backend.clone(),
                    num_trials: 8,
                    num_iterations: 4,
                };
                let mut cached = PropertySet::new();
                let mut uncached = PropertySet::new();
                cached.insert("routing_lookahead_size", lookahead);
                uncached.insert("routing_lookahead_size", lookahead);
                assert_eq!(
                    pass.discover_with_cache(&circuit, &mut cached, true),
                    pass.discover_with_cache(&circuit, &mut uncached, false)
                );
                assert_eq!(
                    cached.get::<Vec<usize>>("sabre_initial_layout"),
                    uncached.get::<Vec<usize>>("sabre_initial_layout")
                );
                let evaluations = *cached.get::<usize>("layout_route_evaluations").unwrap();
                let hits = *cached.get::<usize>("layout_route_cache_hits").unwrap();
                let requests = *cached.get::<usize>("layout_route_requests").unwrap();
                assert!(hits > 0, "the initial forward route must be reused");
                assert_eq!(evaluations + hits, requests);
                assert_eq!(
                    Some(&requests),
                    uncached.get::<usize>("layout_route_evaluations")
                );
                assert_eq!(
                    Some(&requests),
                    uncached.get::<usize>("layout_route_requests")
                );
                assert_eq!(uncached.get::<usize>("layout_route_cache_hits"), Some(&0));

                // Identical candidate search must preserve downstream routing,
                // including the failure on a disconnected interaction graph.
                let router = crate::transpiler::routing::BeamSabrePass {
                    backend: backend.clone(),
                    beam_width: 2,
                    branch_factor: 2,
                    bidir_iterations: 2,
                    ..Default::default()
                };
                let left = router.try_run(&circuit, &mut cached);
                let right = router.try_run(&circuit, &mut uncached);
                match (left, right) {
                    (Ok(left), Ok(right)) => {
                        assert_eq!(left, right);
                        for key in ["initial_layout", "final_layout"] {
                            assert_eq!(
                                cached.get::<Vec<usize>>(key),
                                uncached.get::<Vec<usize>>(key)
                            );
                        }
                        assert_eq!(
                            cached.get::<usize>("swaps_inserted"),
                            uncached.get::<usize>("swaps_inserted")
                        );
                    }
                    (Err(left), Err(right)) => assert_eq!(left.to_string(), right.to_string()),
                    result => panic!("cache changed routing success: {result:?}"),
                }
            }
        }
    }

    #[test]
    fn forward_and_reverse_cache_preserve_exact_cost_and_final_placement() {
        let circuit = interacting_circuit();
        let backend = Backend::linear(4);
        let dist = backend.shortest_path_matrix();
        let (forward, _, _) = build_gate_graph(&circuit, false);
        let reverse = reverse_gate_graph(&forward);
        let forward_preds = pred_counts(&forward);
        let reverse_preds = pred_counts(&reverse);
        let mut forward_cache =
            GreedyRouteCache::new(&forward, &forward_preds, &backend, &dist, 20, true);
        let mut reverse_cache =
            GreedyRouteCache::new(&reverse, &reverse_preds, &backend, &dist, 20, true);
        let mut different_problems_observed = false;
        for seed in 0..12 {
            let initial = random_layout(4, 4, seed);
            let direct_forward =
                greedy_sabre_cost(&forward, &forward_preds, &backend, &dist, &initial, 20);
            let direct_reverse =
                greedy_sabre_cost(&reverse, &reverse_preds, &backend, &dist, &initial, 20);
            different_problems_observed |= direct_forward.0 != direct_reverse.0
                || direct_forward.1.logical_to_physical != direct_reverse.1.logical_to_physical;
            for (cache, expected) in [
                (&mut forward_cache, direct_forward),
                (&mut reverse_cache, direct_reverse),
            ] {
                for _ in 0..2 {
                    let actual = cache.route(&initial);
                    assert_eq!(actual.0.to_bits(), expected.0.to_bits());
                    assert_eq!(actual.1.logical_to_physical, expected.1.logical_to_physical);
                    assert_eq!(actual.1.physical_to_logical, expected.1.physical_to_logical);
                }
            }
        }
        assert!(different_problems_observed);
        for cache in [&forward_cache, &reverse_cache] {
            assert_eq!(cache.evaluations + cache.hits, 24);
            assert!(cache.hits >= 12);
            assert_eq!(cache.evaluations, cache.routes.len());
        }
    }

    #[test]
    fn caching_a_failed_placement_preserves_infinite_cost() {
        let circuit = interacting_circuit();
        let backend = Backend::new("isolated", 4);
        let dist = backend.shortest_path_matrix();
        let (gates, _, _) = build_gate_graph(&circuit, false);
        let preds = pred_counts(&gates);
        let initial = Layout::trivial(4, 4);
        let mut cache = GreedyRouteCache::new(&gates, &preds, &backend, &dist, 20, true);
        let expected = greedy_sabre_cost(&gates, &preds, &backend, &dist, &initial, 20);
        assert_eq!(expected.0, f64::INFINITY);
        for _ in 0..2 {
            let actual = cache.route(&initial);
            assert_eq!(actual.0.to_bits(), expected.0.to_bits());
            assert_eq!(actual.1.logical_to_physical, expected.1.logical_to_physical);
            assert_eq!(actual.1.physical_to_logical, expected.1.physical_to_logical);
        }
        assert_eq!(cache.evaluations, 1);
        assert_eq!(cache.hits, 1);
    }

    #[test]
    fn route_cache_is_not_reused_between_pass_invocations() {
        let pass = SabreLayoutPass {
            backend: Backend::linear(4),
            num_trials: 8,
            num_iterations: 4,
        };
        let circuit = interacting_circuit();
        let mut properties = PropertySet::new();
        pass.try_run(&circuit, &mut properties).unwrap();
        let first_evaluations = *properties.get::<usize>("layout_route_evaluations").unwrap();
        let first_hits = *properties.get::<usize>("layout_route_cache_hits").unwrap();
        pass.try_run(&circuit, &mut properties).unwrap();
        assert_eq!(
            properties.get::<usize>("layout_route_evaluations"),
            Some(&first_evaluations)
        );
        assert_eq!(
            properties.get::<usize>("layout_route_cache_hits"),
            Some(&first_hits)
        );

        let mut changed = circuit.clone();
        changed.operations.reverse();
        let mut fresh = PropertySet::new();
        pass.try_run(&changed, &mut properties).unwrap();
        pass.try_run(&changed, &mut fresh).unwrap();
        assert_eq!(
            properties.get::<Vec<usize>>("sabre_initial_layout"),
            fresh.get::<Vec<usize>>("sabre_initial_layout")
        );
        for key in [
            "layout_route_evaluations",
            "layout_route_cache_hits",
            "layout_route_requests",
        ] {
            assert_eq!(properties.get::<usize>(key), fresh.get::<usize>(key));
        }
        pass.try_run(&Circuit::new(4, 0), &mut properties).unwrap();
        assert!(properties
            .get::<usize>("layout_route_evaluations")
            .is_none());
        assert!(properties.get::<usize>("layout_route_cache_hits").is_none());
        assert!(properties.get::<usize>("layout_route_requests").is_none());
    }

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
