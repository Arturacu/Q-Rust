//! BeamSABRE: hardware-aware routing with beam search.
//!
//! ## Lookahead heuristic
//!
//! Routing uses a SABRE-style lookahead score (Li et al. 2019, ASPLOS).
//! Two strategies are supported via [`LookaheadStrategy`]:
//!
//! - [`LookaheadStrategy::Static`]`(w)` — classical SABRE: the extended-layer
//!   contribution is weighted by a constant `w` (paper default: 0.5).
//! - [`LookaheadStrategy::DynamicV2`] — SABRE-v2 (Li et al. 2023, arXiv:2210.12922):
//!   the score is `max(score_F, score_E · |F| / |E|)`. Note: our
//!   implementation operates on delta costs (post-swap minus pre-swap)
//!   rather than absolute scores per the paper's Eq. 5.
//!
//! Default is `Static(0.5)` for backward compatibility with existing tests.

use crate::backend::Backend;
use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::cmp::Ordering;
use std::collections::HashSet;

/// Default weight for the extended-layer (lookahead) contribution to the
/// SABRE relative-score heuristic. Per Li et al. 2019, ASPLOS, §III.B.
pub const DEFAULT_LOOKAHEAD_WEIGHT: f64 = 0.5;

/// Routing SWAPs use undirected connectivity; native CX direction is repaired
/// separately, after all entangling-gate lowering has completed.
pub(crate) fn undirected_backend(backend: &Backend) -> Backend {
    let mut graph = backend.clone();
    let edges: Vec<_> = (0..backend.num_qubits)
        .flat_map(|a| {
            (0..backend.num_qubits)
                .filter_map(move |b| (a != b && backend.is_adjacent(a, b)).then_some((a, b)))
        })
        .collect();
    graph.set_coupling_map(edges);
    graph
}

/// Strategy for combining front-layer and extended-layer scores during
/// the SABRE swap-selection heuristic.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum LookaheadStrategy {
    /// Classical SABRE: `score = score_F + w · score_E` (delta form).
    Static { weight: f64 },
    /// SABRE-v2 (delta form): `score = max(score_F, score_E · |F| / |E|)`.
    DynamicV2,
}

impl Default for LookaheadStrategy {
    fn default() -> Self {
        Self::Static {
            weight: DEFAULT_LOOKAHEAD_WEIGHT,
        }
    }
}

// ─── Layout ────────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct Layout {
    pub logical_to_physical: Vec<usize>,
    pub physical_to_logical: Vec<usize>,
}

impl Layout {
    pub fn trivial(num_logical: usize, num_physical: usize) -> Self {
        let mut l2p = vec![0usize; num_logical];
        let mut p2l = vec![usize::MAX; num_physical];
        for i in 0..num_logical {
            l2p[i] = i;
            p2l[i] = i;
        }
        Layout {
            logical_to_physical: l2p,
            physical_to_logical: p2l,
        }
    }

    /// Constructs a layout from a logical→physical mapping, validating that
    /// the mapping is injective and that all physical indices are in range.
    ///
    /// Returns [`QRustError::InvalidConfig`] if:
    /// - any `l2p[i] >= num_physical`, or
    /// - the same physical qubit appears twice (non-injective mapping).
    pub fn from_l2p(l2p: Vec<usize>, num_physical: usize) -> Result<Self> {
        let mut p2l = vec![usize::MAX; num_physical];
        for (l, &p) in l2p.iter().enumerate() {
            if p >= num_physical {
                return Err(QRustError::InvalidConfig(format!(
                    "Layout::from_l2p: logical {l} maps to physical {p} \
                     but only {num_physical} physical qubits exist"
                )));
            }
            if p2l[p] != usize::MAX {
                return Err(QRustError::InvalidConfig(format!(
                    "Layout::from_l2p: physical qubit {p} mapped twice \
                     (by logical {} and {l})",
                    p2l[p]
                )));
            }
            p2l[p] = l;
        }
        Ok(Layout {
            logical_to_physical: l2p,
            physical_to_logical: p2l,
        })
    }

    pub fn swap_physical(&mut self, p1: usize, p2: usize) {
        let l1 = self.physical_to_logical[p1];
        let l2 = self.physical_to_logical[p2];
        if l1 != usize::MAX {
            self.logical_to_physical[l1] = p2;
        }
        if l2 != usize::MAX {
            self.logical_to_physical[l2] = p1;
        }
        self.physical_to_logical[p1] = l2;
        self.physical_to_logical[p2] = l1;
    }

    #[inline]
    pub fn l2p(&self, logical: usize) -> usize {
        self.logical_to_physical[logical]
    }
}

// ─── Gate dependency tracking ──────────────────────────────────────────────

#[derive(Clone, Debug)]
pub(crate) struct TwoQGate {
    circuit_index: usize,
    pub(crate) logical_qubits: [usize; 2],
    pub(crate) successors: Vec<usize>,
}

/// A refined placement already visited by layout discovery. The score is a
/// greedy estimate, not a lower bound on the eventual beam route.
#[derive(Clone, Debug)]
pub(crate) struct RefinedLayoutCandidate {
    pub layout: Vec<usize>,
    pub swap_cost: f64,
}

/// A complete, independently lowerable alternative to the unchanged incumbent
/// route. Output-permutation composition and target lowering have not run yet.
/// The caller must compare validated, lowered outputs before accepting one.
#[derive(Clone, Debug)]
pub(crate) struct RoutingOutputCandidate {
    pub circuit: Circuit,
    pub initial_layout: Vec<usize>,
    pub final_layout: Vec<usize>,
    pub swaps_inserted: usize,
    pub origin: &'static str,
}

/// Exact, syntactic diagonal classification. In particular, a nearly-zero U
/// angle is never rounded into this class.
fn is_diagonal(name: &GateType, params: &[f64]) -> bool {
    matches!(
        name,
        GateType::ID
            | GateType::Z
            | GateType::S
            | GateType::Sdg
            | GateType::T
            | GateType::Tdg
            | GateType::RZ
            | GateType::CZ
            | GateType::CU1
            | GateType::CRZ
            | GateType::RZZ
    ) || matches!(name, GateType::U) && params.first() == Some(&0.0)
}

/// Bound each commuting epoch so a long pair of diagonal layers separated by
/// a local fence cannot create a quadratic-size dependency graph. The added
/// ordering is conservative: it reduces search freedom, never changes semantics.
const MAX_DIAGONAL_EPOCH: usize = 64;

/// A diagonal epoch can contain mutually commuting interactions sharing a
/// wire. Every interaction after a non-diagonal local gate depends on *all*
/// interactions before that gate, not merely the last interaction in the
/// original list. This also makes reconstruction's per-wire local-gate cursor
/// safe: it cannot flush a non-diagonal local gate past an unfinished epoch.
pub(crate) fn build_gate_graph(
    circuit: &Circuit,
    commute_diagonal: bool,
) -> (Vec<TwoQGate>, Vec<usize>, Vec<Vec<(usize, Operation)>>) {
    let mut two_q: Vec<TwoQGate> = Vec::new();
    let mut single_q: Vec<Vec<(usize, Operation)>> = vec![vec![]; circuit.num_qubits];
    let mut diagonal_epoch: Vec<Vec<usize>> = vec![vec![]; circuit.num_qubits];
    let mut fence: Vec<Vec<usize>> = vec![vec![]; circuit.num_qubits];

    for (i, op) in circuit.operations.iter().enumerate() {
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } if qubits.len() >= 2 => {
                let idx = two_q.len();
                let diagonal = commute_diagonal && is_diagonal(name, params);
                let mut predecessors = Vec::new();
                for &q in &qubits[..2] {
                    if diagonal {
                        if diagonal_epoch[q].len() == MAX_DIAGONAL_EPOCH {
                            fence[q] = std::mem::take(&mut diagonal_epoch[q]);
                        }
                        predecessors.extend(fence[q].iter().copied());
                        diagonal_epoch[q].push(idx);
                    } else {
                        predecessors.extend(fence[q].iter().copied());
                        predecessors.extend(diagonal_epoch[q].iter().copied());
                        diagonal_epoch[q].clear();
                        fence[q] = vec![idx];
                    }
                }
                predecessors.sort_unstable();
                predecessors.dedup();
                two_q.push(TwoQGate {
                    circuit_index: i,
                    logical_qubits: [qubits[0], qubits[1]],
                    successors: Vec::new(),
                });
                for prev in predecessors {
                    two_q[prev].successors.push(idx);
                }
            }
            Operation::Gate {
                name,
                qubits,
                params,
            } if qubits.len() == 1 => {
                let q = qubits[0];
                single_q[q].push((i, op.clone()));
                if commute_diagonal && !is_diagonal(name, params) && !diagonal_epoch[q].is_empty() {
                    fence[q] = std::mem::take(&mut diagonal_epoch[q]);
                }
            }
            _ => {
                // The router handles barriers by separate segments. The layout
                // cost model sees the whole circuit, so fence each wire here.
                // This does not add cross-wire barrier edges to the layout cost
                // approximation; actual routing segment boundaries enforce the
                // global barrier. The public router rejects unsupported semantics.
                for q in 0..circuit.num_qubits {
                    if !diagonal_epoch[q].is_empty() {
                        fence[q] = std::mem::take(&mut diagonal_epoch[q]);
                    }
                }
            }
        }
    }

    let mut preds = vec![0usize; two_q.len()];
    for g in &two_q {
        for &s in &g.successors {
            preds[s] += 1;
        }
    }
    (two_q, preds, single_q)
}

/// Reverse the actual partial order, including one-qubit fences. Rebuilding a
/// strict per-wire order would silently discard the commutation treatment.
pub(crate) fn reverse_gate_graph(gates: &[TwoQGate]) -> Vec<TwoQGate> {
    let n = gates.len();
    let mut reversed: Vec<_> = gates
        .iter()
        .rev()
        .map(|g| TwoQGate {
            circuit_index: g.circuit_index,
            logical_qubits: g.logical_qubits,
            successors: Vec::new(),
        })
        .collect();
    for (i, gate) in gates.iter().enumerate() {
        for &successor in &gate.successors {
            reversed[n - 1 - successor].successors.push(n - 1 - i);
        }
    }
    for gate in &mut reversed {
        gate.successors.sort_unstable();
    }
    reversed
}

// ─── Routing action & beam ─────────────────────────────────────────────────

#[derive(Clone, Debug)]
enum RoutingAction {
    ExecuteGate(usize),
    InsertSwap(usize, usize),
}

#[derive(Clone)]
struct Beam {
    layout: Layout,
    actions: Vec<RoutingAction>,
    remaining_deps: Vec<usize>,
    front_layer: Vec<usize>,
    cumulative_cost: f64,
    executed_count: usize,
    total_gates: usize,
    complete: bool,
    stuck_swaps: usize,
}

impl Beam {
    fn new(layout: Layout, preds: &[usize], n_gates: usize) -> Self {
        let remaining = preds.to_vec();
        let front: Vec<usize> = (0..n_gates).filter(|&i| remaining[i] == 0).collect();
        Beam {
            layout,
            actions: Vec::new(),
            remaining_deps: remaining,
            front_layer: front,
            cumulative_cost: 0.0,
            executed_count: 0,
            total_gates: n_gates,
            complete: n_gates == 0,
            stuck_swaps: 0,
        }
    }

    fn execute_gate(&mut self, gate_idx: usize, gates: &[TwoQGate]) {
        self.front_layer.retain(|&g| g != gate_idx);
        self.actions.push(RoutingAction::ExecuteGate(gate_idx));
        for &succ in &gates[gate_idx].successors {
            self.remaining_deps[succ] -= 1;
            if self.remaining_deps[succ] == 0 {
                self.front_layer.push(succ);
            }
        }
        self.executed_count += 1;
        if self.executed_count == self.total_gates {
            self.complete = true;
        }
        self.stuck_swaps = 0;
    }

    fn apply_swap(&mut self, p1: usize, p2: usize) {
        self.layout.swap_physical(p1, p2);
        self.actions.push(RoutingAction::InsertSwap(p1, p2));
        self.cumulative_cost += SWAP_COST;
        self.stuck_swaps += 1;
    }
}

const SWAP_COST: f64 = 3.0;

#[inline]
fn safe_dist(dist: &[Vec<usize>], a: usize, b: usize) -> f64 {
    let d = dist[a][b];
    if d == usize::MAX {
        1e9
    } else {
        d as f64
    }
}

/// Build a deterministic extended layer once per routing step, shared by every
/// candidate score. Zero preserves the immediate-successor legacy policy.
pub(crate) fn extended_layer<T>(
    front: &[usize],
    gates: &[T],
    remaining: &[usize],
    limit: usize,
    successors: fn(&T) -> &[usize],
) -> Vec<usize> {
    if limit == 0 {
        let mut layer: Vec<_> = front
            .iter()
            .flat_map(|&g| successors(&gates[g]).iter().copied())
            .collect();
        layer.sort_unstable();
        layer.dedup();
        return layer;
    }
    let mut deps = remaining.to_vec();
    let mut queue = front.to_vec();
    let mut layer = Vec::with_capacity(limit);
    let mut cursor = 0;
    while cursor < queue.len() && layer.len() < limit {
        for &next in successors(&gates[queue[cursor]]) {
            deps[next] -= 1;
            if deps[next] == 0 {
                queue.push(next);
                layer.push(next);
                if layer.len() == limit {
                    break;
                }
            }
        }
        cursor += 1;
    }
    layer
}

#[allow(clippy::too_many_arguments)]
fn relative_score(
    layout: &Layout,
    swap: (usize, usize),
    front: &[usize],
    gates: &[TwoQGate],
    dist: &[Vec<usize>],
    strategy: LookaheadStrategy,
    extended: &[usize],
    normalize: bool,
) -> f64 {
    let (pa, pb) = swap;
    let mut delta_f = 0.0;
    let mut delta_e = 0.0;

    for &g in front {
        let p0 = layout.l2p(gates[g].logical_qubits[0]);
        let p1 = layout.l2p(gates[g].logical_qubits[1]);
        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }
        let old = safe_dist(dist, p0, p1);
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
        delta_f += safe_dist(dist, np0, np1) - old;
    }

    for &g in extended {
        let p0 = layout.l2p(gates[g].logical_qubits[0]);
        let p1 = layout.l2p(gates[g].logical_qubits[1]);
        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }
        let old = safe_dist(dist, p0, p1);
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
        delta_e += safe_dist(dist, np0, np1) - old;
    }

    if normalize && matches!(strategy, LookaheadStrategy::Static { .. }) {
        delta_f /= front.len().max(1) as f64;
        delta_e /= extended.len().max(1) as f64;
    }
    match strategy {
        LookaheadStrategy::Static { weight } => delta_f + weight * delta_e,
        LookaheadStrategy::DynamicV2 => {
            if front.is_empty() || extended.is_empty() {
                return delta_f;
            }
            let scaled_e = delta_e * (front.len() as f64 / extended.len() as f64);
            delta_f.max(scaled_e)
        }
    }
}

#[inline]
fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

#[allow(clippy::too_many_arguments)]
fn beam_sabre_forward(
    gates: &[TwoQGate],
    preds: &[usize],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
    beam_width: usize,
    branch_factor: usize,
    strategy: LookaheadStrategy,
    lookahead_size: usize,
    rank_progress: bool,
    retain_completed_ties: bool,
) -> Result<Vec<Beam>> {
    if gates.is_empty() {
        return Ok(vec![Beam::new(initial_layout.clone(), preds, 0)]);
    }
    let mut beams = vec![Beam::new(initial_layout.clone(), preds, gates.len())];
    let neighbors: Vec<_> = (0..backend.num_qubits)
        .map(|p| backend.neighbors(p))
        .collect();

    loop {
        for beam in &mut beams {
            if beam.complete {
                continue;
            }
            execute_routable(beam, gates, dist);
        }
        if beams.iter().all(|b| b.complete) {
            break;
        }

        let mut new_beams: Vec<Beam> = Vec::new();
        for beam in &beams {
            if beam.complete {
                new_beams.push(beam.clone());
                continue;
            }
            let extended = extended_layer(
                &beam.front_layer,
                gates,
                &beam.remaining_deps,
                lookahead_size,
                |g| &g.successors,
            );
            let mut candidates = Vec::new();
            for &gate_idx in &beam.front_layer {
                let g = &gates[gate_idx];
                let p0 = beam.layout.l2p(g.logical_qubits[0]);
                let p1 = beam.layout.l2p(g.logical_qubits[1]);
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
            let mut scored: Vec<(f64, (usize, usize))> = candidates
                .iter()
                .map(|&swap| {
                    let delta = relative_score(
                        &beam.layout,
                        swap,
                        &beam.front_layer,
                        gates,
                        dist,
                        strategy,
                        &extended,
                        lookahead_size > 0,
                    );
                    (beam.cumulative_cost + SWAP_COST + delta, swap)
                })
                .collect();
            scored.sort_by(|a, b| cmp_f64(a.0, b.0).then_with(|| a.1.cmp(&b.1)));
            if scored.is_empty() {
                // Beam is stuck (no path in disconnected topology). Drop it.
                continue;
            }
            let best_delta = scored[0].0 - beam.cumulative_cost - SWAP_COST;
            let effective_branch = if best_delta >= 0.0 {
                1
            } else {
                branch_factor.min(scored.len())
            };
            for i in 0..effective_branch {
                let mut child = beam.clone();
                let (_, swap) = scored[i];
                child.apply_swap(swap.0, swap.1);
                if child.stuck_swaps < 3 * backend.num_qubits {
                    new_beams.push(child);
                }
            }
        }
        if rank_progress {
            for child in &mut new_beams {
                execute_routable(child, gates, dist);
            }
            let distance = |beam: &Beam| -> usize {
                beam.front_layer
                    .iter()
                    .map(|&g| {
                        let [a, b] = gates[g].logical_qubits;
                        dist[beam.layout.l2p(a)][beam.layout.l2p(b)]
                    })
                    .fold(0usize, usize::saturating_add)
            };
            new_beams.sort_by(|a, b| {
                cmp_f64(a.cumulative_cost, b.cumulative_cost)
                    .then_with(|| b.executed_count.cmp(&a.executed_count))
                    .then_with(|| distance(a).cmp(&distance(b)))
                    .then_with(|| {
                        a.layout
                            .logical_to_physical
                            .cmp(&b.layout.logical_to_physical)
                    })
            });
            let mut seen = HashSet::new();
            new_beams.retain(|beam| {
                seen.insert((
                    beam.layout.logical_to_physical.clone(),
                    beam.remaining_deps.clone(),
                    beam.front_layer.clone(),
                ))
            });
            // Every unfinished child has exhausted its free gates. It needs
            // another SWAP, so cannot improve this completed child's cost.
            if let Some(completed) = new_beams.iter().find(|b| b.complete) {
                if retain_completed_ties {
                    let cost = completed.cumulative_cost;
                    return Ok(new_beams
                        .iter()
                        .filter(|b| b.complete && b.cumulative_cost == cost)
                        .cloned()
                        .collect());
                }
                return Ok(vec![completed.clone()]);
            }
        } else {
            new_beams.sort_by(|a, b| cmp_f64(a.cumulative_cost, b.cumulative_cost));
        }

        // Extract witness BEFORE dropping stuck beams.
        let pre_prune_witness: Option<(usize, usize)> = beams
            .iter()
            .find(|b| !b.front_layer.is_empty())
            .and_then(|b| {
                b.front_layer
                    .first()
                    .map(|&gi| (gates[gi].logical_qubits[0], gates[gi].logical_qubits[1]))
            });

        new_beams.truncate(beam_width);
        beams = new_beams;
        if beams.is_empty() {
            let (from, to) = pre_prune_witness.unwrap_or((0, 0));
            return Err(QRustError::DisconnectedTopology { from, to });
        }
    }
    Ok(beams)
}

fn execute_routable(beam: &mut Beam, gates: &[TwoQGate], dist: &[Vec<usize>]) {
    loop {
        let routable: Vec<usize> = beam
            .front_layer
            .iter()
            .filter(|&&g| {
                let p0 = beam.layout.l2p(gates[g].logical_qubits[0]);
                let p1 = beam.layout.l2p(gates[g].logical_qubits[1]);
                dist[p0][p1] == 1
            })
            .cloned()
            .collect();
        if routable.is_empty() {
            break;
        }
        for g in routable {
            beam.execute_gate(g, gates);
        }
    }
}

// ─── Circuit reconstruction ────────────────────────────────────────────────

fn reconstruct_circuit(
    original: &Circuit,
    beam: &Beam,
    initial_layout: &Layout,
    gates: &[TwoQGate],
    num_physical: usize,
    single_q_ops: &[Vec<(usize, Operation)>],
) -> Circuit {
    let mut output = Circuit::new(num_physical, original.num_cbits);
    output.custom_gates = original.custom_gates.clone();
    let mut layout = initial_layout.clone();
    let mut sq_cursors: Vec<usize> = vec![0; original.num_qubits];

    let non_gate_ops: Vec<(usize, &Operation)> = original
        .operations
        .iter()
        .enumerate()
        .filter(|(_, op)| !matches!(op, Operation::Gate { .. }))
        .collect();
    let mut ng_cursor = 0usize;

    for action in &beam.actions {
        match action {
            RoutingAction::ExecuteGate(gate_idx) => {
                let gate = &gates[*gate_idx];
                while ng_cursor < non_gate_ops.len()
                    && non_gate_ops[ng_cursor].0 < gate.circuit_index
                {
                    let (_, op) = non_gate_ops[ng_cursor];
                    emit_non_gate_op(op, &layout, &mut output);
                    ng_cursor += 1;
                }
                for &lq in &gate.logical_qubits {
                    let pq = layout.l2p(lq);
                    while sq_cursors[lq] < single_q_ops[lq].len() {
                        let (orig_idx, ref op) = single_q_ops[lq][sq_cursors[lq]];
                        if orig_idx < gate.circuit_index {
                            if let Operation::Gate { name, params, .. } = op {
                                output.add_op(Operation::Gate {
                                    name: name.clone(),
                                    qubits: vec![pq],
                                    params: params.clone(),
                                });
                            }
                            sq_cursors[lq] += 1;
                        } else {
                            break;
                        }
                    }
                }
                let orig = &original.operations[gate.circuit_index];
                if let Operation::Gate { name, params, .. } = orig {
                    let p0 = layout.l2p(gate.logical_qubits[0]);
                    let p1 = layout.l2p(gate.logical_qubits[1]);
                    output.add_op(Operation::Gate {
                        name: name.clone(),
                        qubits: vec![p0, p1],
                        params: params.clone(),
                    });
                }
            }
            RoutingAction::InsertSwap(p1, p2) => {
                output.add_op(Operation::Gate {
                    name: GateType::SWAP,
                    qubits: vec![*p1, *p2],
                    params: vec![],
                });
                layout.swap_physical(*p1, *p2);
            }
        }
    }

    for lq in 0..original.num_qubits {
        let pq = layout.l2p(lq);
        while sq_cursors[lq] < single_q_ops[lq].len() {
            let (_, ref op) = single_q_ops[lq][sq_cursors[lq]];
            if let Operation::Gate { name, params, .. } = op {
                output.add_op(Operation::Gate {
                    name: name.clone(),
                    qubits: vec![pq],
                    params: params.clone(),
                });
            }
            sq_cursors[lq] += 1;
        }
    }
    while ng_cursor < non_gate_ops.len() {
        let (_, op) = non_gate_ops[ng_cursor];
        emit_non_gate_op(op, &layout, &mut output);
        ng_cursor += 1;
    }
    output
}

fn emit_non_gate_op(op: &Operation, layout: &Layout, output: &mut Circuit) {
    match op {
        Operation::Measure { qubit, cbit } => output.add_op(Operation::Measure {
            qubit: layout.l2p(*qubit),
            cbit: *cbit,
        }),
        Operation::Reset { qubit } => output.add_op(Operation::Reset {
            qubit: layout.l2p(*qubit),
        }),
        Operation::Barrier { qubits } => {
            let phys: Vec<usize> = qubits.iter().map(|&q| layout.l2p(q)).collect();
            output.add_op(Operation::Barrier { qubits: phys });
        }
        _ => {}
    }
}

// ─── BeamSabrePass ─────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct BeamSabrePass {
    pub backend: Backend,
    pub beam_width: usize,
    pub branch_factor: usize,
    pub bidir_iterations: usize,
    /// Heuristic strategy for combining front- and extended-layer scores.
    /// Default: `LookaheadStrategy::Static { weight: 0.5 }` (classical SABRE).
    pub lookahead_strategy: LookaheadStrategy,
}

impl Default for BeamSabrePass {
    fn default() -> Self {
        Self {
            backend: Backend::all_to_all(1),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 2,
            lookahead_strategy: LookaheadStrategy::default(),
        }
    }
}

impl BeamSabrePass {
    pub fn try_run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Result<Circuit> {
        property_set.remove::<Vec<RoutingOutputCandidate>>("routing_output_candidates");
        crate::transpiler::validate_circuit_structure(circuit)?;
        crate::transpiler::validate_static_semantics(circuit)?;
        if circuit.num_qubits > self.backend.num_qubits {
            return Err(QRustError::InsufficientQubits {
                circuit: circuit.num_qubits,
                backend: self.backend.num_qubits,
            });
        }
        if circuit.operations.iter().any(|op| matches!(op,
            Operation::Gate { name, qubits, .. } if qubits.len() > 2 || matches!(name, GateType::Custom(_))
        )) {
            return Err(QRustError::Unsupported("unroll custom and high-arity gates before routing".into()));
        }
        let requested_layout =
            if let Some(l2p) = property_set.get::<Vec<usize>>("sabre_initial_layout") {
                if l2p.len() != circuit.num_qubits {
                    return Err(QRustError::Routing(
                        "initial layout length does not match circuit width".into(),
                    ));
                }
                Layout::from_l2p(l2p.clone(), self.backend.num_qubits)?
            } else {
                Layout::trivial(circuit.num_qubits, self.backend.num_qubits)
            };

        // A barrier is a scheduling boundary. Route each unitary segment in
        // order while carrying its final placement into the next segment.
        // Later segments must not freely choose a new initial placement.
        if circuit.operations.iter().any(Operation::is_barrier) {
            let mut output = Circuit::new(self.backend.num_qubits, circuit.num_cbits);
            output.custom_gates = circuit.custom_gates.clone();
            let mut initial = requested_layout.logical_to_physical.clone();
            let mut current = initial.clone();
            let mut first = true;
            let mut swaps = 0usize;
            let mut start = 0usize;
            for boundary in 0..=circuit.operations.len() {
                if boundary < circuit.operations.len() && !circuit.operations[boundary].is_barrier()
                {
                    continue;
                }
                let mut segment = Circuit::new(circuit.num_qubits, circuit.num_cbits);
                segment.custom_gates = circuit.custom_gates.clone();
                segment.operations = circuit.operations[start..boundary].to_vec();
                let mut properties = PropertySet::new();
                properties.insert("sabre_initial_layout", current.clone());
                if let Some(value) = property_set.get::<usize>("routing_lookahead_size") {
                    properties.insert("routing_lookahead_size", *value);
                }
                if let Some(value) = property_set.get::<bool>("routing_rank_progress") {
                    properties.insert("routing_rank_progress", *value);
                }
                if let Some(value) = property_set.get::<bool>("routing_commute_diagonal") {
                    properties.insert("routing_commute_diagonal", *value);
                }
                let mut segment_pass = self.clone();
                if !first {
                    segment_pass.bidir_iterations = 1;
                }
                let routed = segment_pass.try_run(&segment, &mut properties)?;
                if first {
                    initial = properties
                        .get::<Vec<usize>>("initial_layout")
                        .unwrap()
                        .clone();
                    first = false;
                }
                current = properties
                    .get::<Vec<usize>>("final_layout")
                    .unwrap()
                    .clone();
                swaps += properties
                    .get::<usize>("swaps_inserted")
                    .copied()
                    .unwrap_or(0);
                output.operations.extend(routed.operations);
                if boundary < circuit.operations.len() {
                    if let Operation::Barrier { qubits } = &circuit.operations[boundary] {
                        output.add_op(Operation::Barrier {
                            qubits: qubits.iter().map(|&q| current[q]).collect(),
                        });
                    }
                }
                start = boundary + 1;
            }
            property_set.insert("initial_layout", initial);
            property_set.insert("final_layout", current);
            property_set.insert("swaps_inserted", swaps);
            return Ok(output);
        }
        let has_2q = circuit.operations.iter().any(|op| {
            if let Operation::Gate { qubits, .. } = op {
                qubits.len() >= 2
            } else {
                false
            }
        });
        if self.beam_width < 1 || self.branch_factor < 1 || self.bidir_iterations < 1 {
            return Err(QRustError::InvalidConfig(format!(
                "BeamSabrePass: beam_width={}, branch_factor={}, bidir_iterations={} (all must be >= 1)",
                self.beam_width, self.branch_factor, self.bidir_iterations
            )));
        }
        if !has_2q || self.backend.is_fully_connected() {
            property_set.insert(
                "initial_layout",
                requested_layout.logical_to_physical.clone(),
            );
            property_set.insert("final_layout", requested_layout.logical_to_physical.clone());
            property_set.insert("swaps_inserted", 0_usize);
            let width = requested_layout
                .logical_to_physical
                .iter()
                .copied()
                .max()
                .map_or(circuit.num_qubits, |q| circuit.num_qubits.max(q + 1));
            let mut output = Circuit::new(width, circuit.num_cbits);
            output.custom_gates = circuit.custom_gates.clone();
            for op in &circuit.operations {
                if let Operation::Gate {
                    name,
                    qubits,
                    params,
                } = op
                {
                    output.add_op(Operation::Gate {
                        name: name.clone(),
                        qubits: qubits.iter().map(|&q| requested_layout.l2p(q)).collect(),
                        params: params.clone(),
                    });
                } else {
                    emit_non_gate_op(op, &requested_layout, &mut output);
                }
            }
            return Ok(output);
        }

        let connectivity = undirected_backend(&self.backend);
        let dist = connectivity.shortest_path_matrix();
        let commute_diagonal = property_set
            .get::<bool>("routing_commute_diagonal")
            .copied()
            .unwrap_or(false);
        let (gates, preds, single_q_ops) = build_gate_graph(circuit, commute_diagonal);

        let lookahead_size = property_set
            .get::<usize>("routing_lookahead_size")
            .copied()
            .unwrap_or(20);
        let rank_progress = property_set
            .get::<bool>("routing_rank_progress")
            .copied()
            .unwrap_or(true);
        let retain_candidates = property_set
            .get::<bool>("retain_routing_candidates")
            .copied()
            .unwrap_or(false);
        let mut best_layout = requested_layout;

        let mut best_beam: Option<Beam> = None;
        let mut best_cost = f64::MAX;
        let mut best_initial_layout = best_layout.clone();
        let mut depth_candidate: Option<RoutingOutputCandidate> = None;
        let mut visited_initial_layouts = Vec::new();

        let reversed_gates = reverse_gate_graph(&gates);
        let reversed_preds: Vec<usize> = {
            let mut counts = vec![0usize; reversed_gates.len()];
            for g in &reversed_gates {
                for &s in &g.successors {
                    counts[s] += 1;
                }
            }
            counts
        };

        for iter in 0..self.bidir_iterations {
            let current_initial = best_layout.clone();
            if retain_candidates {
                visited_initial_layouts.push(current_initial.logical_to_physical.clone());
            }
            let fwd = beam_sabre_forward(
                &gates,
                &preds,
                &connectivity,
                &dist,
                &current_initial,
                self.beam_width,
                self.branch_factor,
                self.lookahead_strategy,
                lookahead_size,
                rank_progress,
                retain_candidates,
            )?;
            let fwd_best = fwd
                .iter()
                .min_by(|a, b| cmp_f64(a.cumulative_cost, b.cumulative_cost))
                .ok_or_else(|| {
                    QRustError::Routing(
                        "beam_sabre_forward returned an empty result vector (internal bug)".into(),
                    )
                })?;
            if fwd_best.cumulative_cost < best_cost {
                best_cost = fwd_best.cumulative_cost;
                best_beam = Some(fwd_best.clone());
                best_initial_layout = current_initial.clone();
                depth_candidate = None;
            }
            if retain_candidates && fwd_best.cumulative_cost == best_cost {
                for candidate in &fwd {
                    if candidate.cumulative_cost != best_cost {
                        continue;
                    }
                    let output = reconstruct_circuit(
                        circuit,
                        candidate,
                        &current_initial,
                        &gates,
                        self.backend.num_qubits,
                        &single_q_ops,
                    );
                    if depth_candidate
                        .as_ref()
                        .map_or(true, |previous| output.depth() < previous.circuit.depth())
                    {
                        depth_candidate = Some(RoutingOutputCandidate {
                            circuit: output,
                            initial_layout: current_initial.logical_to_physical.clone(),
                            final_layout: candidate.layout.logical_to_physical.clone(),
                            swaps_inserted: candidate.cumulative_cost as usize / 3,
                            origin: "equal_swap_depth",
                        });
                    }
                }
            }
            if iter + 1 < self.bidir_iterations {
                let bwd = beam_sabre_forward(
                    &reversed_gates,
                    &reversed_preds,
                    &connectivity,
                    &dist,
                    &fwd_best.layout,
                    self.beam_width,
                    self.branch_factor,
                    self.lookahead_strategy,
                    lookahead_size,
                    rank_progress,
                    false,
                )?;
                let bwd_best = bwd
                    .iter()
                    .min_by(|a, b| cmp_f64(a.cumulative_cost, b.cumulative_cost))
                    .ok_or_else(|| {
                        QRustError::Routing(
                            "beam_sabre_forward (bwd) returned an empty result vector".into(),
                        )
                    })?;
                best_layout = bwd_best.layout.clone();
            }
        }

        let winning = best_beam.ok_or_else(|| {
            QRustError::Routing(
                "no routing solution found across all bidirectional iterations".into(),
            )
        })?;
        property_set.insert(
            "initial_layout",
            best_initial_layout.logical_to_physical.clone(),
        );
        property_set.insert("final_layout", winning.layout.logical_to_physical.clone());
        property_set.insert("swaps_inserted", winning.cumulative_cost as usize / 3);

        let output = reconstruct_circuit(
            circuit,
            &winning,
            &best_initial_layout,
            &gates,
            self.backend.num_qubits,
            &single_q_ops,
        );
        if retain_candidates {
            let mut candidates = Vec::with_capacity(2);
            if let Some(candidate) = depth_candidate {
                if candidate.circuit != output
                    || candidate.initial_layout != best_initial_layout.logical_to_physical
                    || candidate.final_layout != winning.layout.logical_to_physical
                {
                    candidates.push(candidate);
                }
            }

            // Spend at most one extra forward search on a refined placement
            // that the greedy score considered close. One SWAP of slack is a
            // bounded exploration policy, not a dominance or optimality proof.
            // The incumbent trajectory and every candidate already visited by
            // it are preserved, so this cannot redirect bidirectional search.
            let alternate = property_set
                .get::<Vec<RefinedLayoutCandidate>>("sabre_refined_layout_candidates")
                .and_then(|layouts| {
                    layouts
                        .iter()
                        .filter(|candidate| {
                            candidate.layout.len() == circuit.num_qubits
                                && candidate.swap_cost.is_finite()
                                && candidate.swap_cost >= 0.0
                                && candidate.swap_cost <= best_cost + SWAP_COST
                                && !visited_initial_layouts.contains(&candidate.layout)
                        })
                        .min_by(|a, b| cmp_f64(a.swap_cost, b.swap_cost))
                });
            if let Some(alternate) = alternate {
                if let Ok(initial) =
                    Layout::from_l2p(alternate.layout.clone(), self.backend.num_qubits)
                {
                    if let Ok(routes) = beam_sabre_forward(
                        &gates,
                        &preds,
                        &connectivity,
                        &dist,
                        &initial,
                        self.beam_width,
                        self.branch_factor,
                        self.lookahead_strategy,
                        lookahead_size,
                        rank_progress,
                        true,
                    ) {
                        let mut alternate_best: Option<RoutingOutputCandidate> = None;
                        for route in routes {
                            let routed = reconstruct_circuit(
                                circuit,
                                &route,
                                &initial,
                                &gates,
                                self.backend.num_qubits,
                                &single_q_ops,
                            );
                            let swaps = route.cumulative_cost as usize / 3;
                            if alternate_best.as_ref().map_or(true, |previous| {
                                (swaps, routed.depth())
                                    < (previous.swaps_inserted, previous.circuit.depth())
                            }) {
                                alternate_best = Some(RoutingOutputCandidate {
                                    circuit: routed,
                                    initial_layout: initial.logical_to_physical.clone(),
                                    final_layout: route.layout.logical_to_physical.clone(),
                                    swaps_inserted: swaps,
                                    origin: "refined_layout_forward",
                                });
                            }
                        }
                        if let Some(candidate) = alternate_best {
                            let distinct = |other: &RoutingOutputCandidate| {
                                other.circuit != candidate.circuit
                                    || other.initial_layout != candidate.initial_layout
                                    || other.final_layout != candidate.final_layout
                            };
                            if (candidate.circuit != output
                                || candidate.initial_layout
                                    != best_initial_layout.logical_to_physical
                                || candidate.final_layout != winning.layout.logical_to_physical)
                                && candidates.iter().all(distinct)
                            {
                                candidates.push(candidate);
                            }
                        }
                    }
                }
            }
            property_set.insert("routing_output_candidates", candidates);
        }
        Ok(output)
    }
}

impl Pass for BeamSabrePass {
    fn name(&self) -> &str {
        "BeamSabrePass"
    }

    fn try_run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Result<Circuit> {
        BeamSabrePass::try_run(self, circuit, property_set)
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        match self.try_run(circuit, property_set) {
            Ok(c) => c,
            Err(e) => {
                property_set.insert("beam_sabre_error", e.to_string());
                crate::transpiler::warn_diagnostic(format_args!(
                    "BeamSabrePass failed: {e}; returning original circuit"
                ));
                circuit.clone()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Circuit, GateType, Operation};
    use crate::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};

    fn default_pass(backend: Backend, beam: usize, branch: usize, bidir: usize) -> BeamSabrePass {
        BeamSabrePass {
            backend,
            beam_width: beam,
            branch_factor: branch,
            bidir_iterations: bidir,
            lookahead_strategy: LookaheadStrategy::default(),
        }
    }

    #[test]
    fn long_diagonal_epochs_bound_edges_without_losing_fence_dependencies() {
        let count = MAX_DIAGONAL_EPOCH * 3 + 5;
        let mut circuit = Circuit::new(3, 0);
        for pair in [[0, 1], [0, 2]] {
            for _ in 0..count {
                circuit.add_op(Operation::Gate {
                    name: GateType::RZZ,
                    qubits: pair.to_vec(),
                    params: vec![0.013],
                });
            }
            circuit.add_op(Operation::Gate {
                name: GateType::H,
                qubits: vec![0],
                params: vec![],
            });
        }
        let (gates, predecessors, _) = build_gate_graph(&circuit, true);
        assert!(predecessors.iter().all(|&n| n <= 4 * MAX_DIAGONAL_EPOCH));
        assert!(predecessors.iter().sum::<usize>() <= 4 * MAX_DIAGONAL_EPOCH * gates.len());
        // Every interaction before H reaches every interaction after H even
        // when the original layer is split into several bounded epochs.
        for source in 0..count {
            let mut reachable = vec![false; gates.len()];
            reachable[source] = true;
            for (i, gate) in gates.iter().enumerate() {
                if reachable[i] {
                    for &successor in &gate.successors {
                        reachable[successor] = true;
                    }
                }
            }
            assert!(reachable[count..].iter().all(|&r| r));
        }
    }

    #[test]
    fn diagonal_epoch_fence_waits_for_every_interaction_and_reverse_preserves_it() {
        let mut circuit = Circuit::new(4, 0);
        for pair in [[0, 1], [0, 2]] {
            circuit.add_op(Operation::Gate {
                name: GateType::RZZ,
                qubits: pair.to_vec(),
                params: vec![0.2],
            });
        }
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        for pair in [[0, 3], [0, 1]] {
            circuit.add_op(Operation::Gate {
                name: GateType::RZZ,
                qubits: pair.to_vec(),
                params: vec![0.3],
            });
        }
        let (gates, predecessors, _) = build_gate_graph(&circuit, true);
        assert_eq!(predecessors, vec![0, 0, 2, 2]);
        assert_eq!(gates[0].successors, vec![2, 3]);
        assert_eq!(gates[1].successors, vec![2, 3]);
        let reversed = reverse_gate_graph(&gates);
        assert_eq!(reversed[0].successors, vec![2, 3]);
        assert_eq!(reversed[1].successors, vec![2, 3]);
    }

    #[test]
    fn controlled_phases_share_diagonal_epochs_but_not_local_fences() {
        let mut circuit = Circuit::new(4, 0);
        for (name, pair, params) in [
            (GateType::CU1, [0, 1], vec![0.2]),
            (GateType::CZ, [2, 0], vec![]),
            (GateType::RZZ, [0, 3], vec![-0.3]),
        ] {
            circuit.add_op(Operation::Gate {
                name,
                qubits: pair.to_vec(),
                params,
            });
        }
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![0.0, 0.7, -0.2],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CU1,
            qubits: vec![1, 0],
            params: vec![1e-15],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        for pair in [[0, 2], [3, 0]] {
            circuit.add_op(Operation::Gate {
                name: GateType::CU1,
                qubits: pair.to_vec(),
                params: vec![-0.4],
            });
        }
        let (gates, predecessors, _) = build_gate_graph(&circuit, true);
        assert_eq!(predecessors, vec![0, 0, 0, 0, 4, 4]);
        for gate in &gates[..4] {
            assert_eq!(gate.successors, vec![4, 5]);
        }
        let reversed = reverse_gate_graph(&gates);
        assert_eq!(reversed[0].successors, vec![2, 3, 4, 5]);
        assert_eq!(reversed[1].successors, vec![2, 3, 4, 5]);
        let (_, strict_predecessors, _) = build_gate_graph(&circuit, false);
        assert_eq!(strict_predecessors, vec![0, 1, 1, 2, 2, 2]);
    }

    #[test]
    fn only_exact_diagonal_locals_allow_routing_reordering() {
        for name in [GateType::RZZ, GateType::CU1] {
            for theta in [0.0, 1e-15, -1e-15, std::f64::consts::TAU] {
                let mut circuit = Circuit::new(3, 0);
                circuit.add_op(Operation::Gate {
                    name: name.clone(),
                    qubits: vec![0, 1],
                    params: vec![0.2],
                });
                circuit.add_op(Operation::Gate {
                    name: GateType::U,
                    qubits: vec![0],
                    params: vec![theta, 0.3, -0.1],
                });
                circuit.add_op(Operation::Gate {
                    name: name.clone(),
                    qubits: vec![0, 2],
                    params: vec![0.4],
                });
                let (_, predecessors, _) = build_gate_graph(&circuit, true);
                assert_eq!(predecessors, vec![0, usize::from(theta != 0.0)]);
                let (_, strict, _) = build_gate_graph(&circuit, false);
                assert_eq!(strict, vec![0, 1]);
            }
        }
    }

    #[test]
    fn extended_layer_sees_beyond_repeated_pair_and_counts_duplicate_dependencies() {
        // CX(0,1), CX(0,1), CX(0,2): the old layer saw only the
        // repeated first pair, hiding the next interaction partner.
        let mut circuit = Circuit::new(3, 0);
        for pair in [[0, 1], [0, 1], [0, 2]] {
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: pair.to_vec(),
                params: vec![],
            });
        }
        let (gates, preds, _) = build_gate_graph(&circuit, false);
        assert_eq!(
            extended_layer(&[0], &gates, &preds, 0, |g| &g.successors),
            vec![1]
        );
        assert_eq!(
            extended_layer(&[0], &gates, &preds, 20, |g| &g.successors),
            vec![1, 2]
        );
        assert_eq!(
            extended_layer(&[0], &gates, &preds, 1, |g| &g.successors),
            vec![1]
        );
    }

    #[test]
    fn tied_routes_are_reproducible() {
        let mut circuit = Circuit::new(5, 0);
        for pair in [[0, 4], [1, 3], [0, 2], [2, 4], [0, 3]] {
            circuit.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: pair.to_vec(),
                params: vec![],
            });
        }
        let pass = default_pass(Backend::linear(5), 4, 3, 2);
        let mut first = PropertySet::new();
        let expected = pass.try_run(&circuit, &mut first).unwrap();
        for _ in 0..8 {
            let mut properties = PropertySet::new();
            let output = pass.try_run(&circuit, &mut properties).unwrap();
            assert_eq!(output.operations, expected.operations);
            assert_eq!(
                properties.get::<Vec<usize>>("initial_layout"),
                first.get::<Vec<usize>>("initial_layout")
            );
            assert_eq!(
                properties.get::<Vec<usize>>("final_layout"),
                first.get::<Vec<usize>>("final_layout")
            );
        }
    }

    #[test]
    fn retained_routes_preserve_incumbent_and_have_valid_complete_mappings() {
        let mut circuit = Circuit::new(5, 0);
        for (i, pair) in [[0, 4], [1, 3], [0, 2], [2, 4], [0, 3]]
            .into_iter()
            .enumerate()
        {
            circuit.add_op(Operation::Gate {
                name: GateType::CU1,
                qubits: pair.to_vec(),
                params: vec![0.13 * (i + 1) as f64],
            });
        }
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 4],
            params: vec![],
        });
        let pass = default_pass(Backend::linear(5), 4, 3, 2);
        let original = circuit_to_unitary(&circuit);
        let mut saw_alternate = false;
        for commute in [false, true] {
            for rank_progress in [false, true] {
                let mut ordinary = PropertySet::new();
                let mut retained = PropertySet::new();
                for properties in [&mut ordinary, &mut retained] {
                    properties.insert("routing_commute_diagonal", commute);
                    properties.insert("routing_rank_progress", rank_progress);
                }
                retained.insert("retain_routing_candidates", true);
                retained.insert(
                    "sabre_refined_layout_candidates",
                    vec![
                        RefinedLayoutCandidate {
                            layout: vec![4, 3, 2, 1, 0],
                            swap_cost: 0.0,
                        },
                        RefinedLayoutCandidate {
                            layout: vec![1, 3, 0, 4, 2],
                            swap_cost: 0.0,
                        },
                        RefinedLayoutCandidate {
                            layout: vec![2, 4, 1, 0, 3],
                            swap_cost: 0.0,
                        },
                    ],
                );
                let incumbent = pass.try_run(&circuit, &mut ordinary).unwrap();
                assert_eq!(incumbent, pass.try_run(&circuit, &mut retained).unwrap());
                for key in ["initial_layout", "final_layout"] {
                    assert_eq!(
                        ordinary.get::<Vec<usize>>(key),
                        retained.get::<Vec<usize>>(key)
                    );
                }
                assert_eq!(
                    ordinary.get::<usize>("swaps_inserted"),
                    retained.get::<usize>("swaps_inserted")
                );
                let candidates = retained
                    .get::<Vec<RoutingOutputCandidate>>("routing_output_candidates")
                    .unwrap();
                assert!(candidates.len() <= 2);
                assert!(
                    candidates
                        .iter()
                        .filter(|c| c.origin == "refined_layout_forward")
                        .count()
                        <= 1
                );
                for candidate in candidates {
                    saw_alternate |= candidate.origin == "refined_layout_forward";
                    assert_eq!(candidate.initial_layout.len(), 5);
                    assert_eq!(candidate.final_layout.len(), 5);
                    assert_eq!(
                        candidate.swaps_inserted,
                        candidate
                            .circuit
                            .operations
                            .iter()
                            .filter(|op| matches!(
                                op,
                                Operation::Gate {
                                    name: GateType::SWAP,
                                    ..
                                }
                            ))
                            .count()
                    );
                    for op in &candidate.circuit.operations {
                        if let Operation::Gate { qubits, .. } = op {
                            if qubits.len() == 2 {
                                assert_eq!(qubits[0].abs_diff(qubits[1]), 1);
                            }
                        }
                    }
                    let logical = extract_logical_unitary(
                        &circuit_to_unitary(&candidate.circuit),
                        5,
                        &candidate.initial_layout,
                        &candidate.final_layout,
                    );
                    // Routing adds no phase: every input column must agree.
                    assert!((&logical - &original).norm() < 1e-10);
                }
            }
        }
        assert!(
            saw_alternate,
            "the alternate-layout branch must be exercised"
        );
    }

    #[test]
    fn routing_candidates_are_cleared_across_disabled_barrier_and_error_runs() {
        let pass = default_pass(Backend::linear(3), 4, 3, 2);
        let stale = RoutingOutputCandidate {
            circuit: Circuit::new(3, 0),
            initial_layout: vec![0, 1, 2],
            final_layout: vec![0, 1, 2],
            swaps_inserted: 0,
            origin: "stale",
        };
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 2],
            params: vec![],
        });
        for mode in 0..3 {
            let mut properties = PropertySet::new();
            properties.insert("retain_routing_candidates", mode != 0);
            properties.insert("routing_output_candidates", vec![stale.clone()]);
            let mut input = circuit.clone();
            if mode == 1 {
                input.add_op(Operation::Barrier {
                    qubits: vec![0, 1, 2],
                });
            } else if mode == 2 {
                input.add_op(Operation::Gate {
                    name: GateType::CX,
                    qubits: vec![0, 9],
                    params: vec![],
                });
            }
            assert_eq!(pass.try_run(&input, &mut properties).is_ok(), mode != 2);
            assert!(properties
                .get::<Vec<RoutingOutputCandidate>>("routing_output_candidates")
                .is_none());
        }
    }

    #[test]
    fn test_layout_trivial() {
        let l = Layout::trivial(3, 5);
        assert_eq!(l.l2p(0), 0);
        assert_eq!(l.physical_to_logical[3], usize::MAX);
    }

    #[test]
    fn test_layout_swap() {
        let mut l = Layout::trivial(3, 3);
        l.swap_physical(0, 2);
        assert_eq!(l.l2p(0), 2);
        assert_eq!(l.l2p(2), 0);
        assert_eq!(l.physical_to_logical[0], 2);
        assert_eq!(l.physical_to_logical[2], 0);
    }

    #[test]
    fn test_layout_from_l2p_accepts_valid() {
        let l = Layout::from_l2p(vec![2, 0, 1], 3).expect("valid mapping");
        assert_eq!(l.l2p(0), 2);
        assert_eq!(l.l2p(1), 0);
        assert_eq!(l.l2p(2), 1);
        assert_eq!(l.physical_to_logical[0], 1);
        assert_eq!(l.physical_to_logical[1], 2);
        assert_eq!(l.physical_to_logical[2], 0);
    }

    #[test]
    fn test_layout_from_l2p_rejects_out_of_range() {
        let r = Layout::from_l2p(vec![0, 5], 3);
        assert!(matches!(r, Err(QRustError::InvalidConfig(_))));
    }

    #[test]
    fn test_layout_from_l2p_rejects_double_mapping() {
        let r = Layout::from_l2p(vec![1, 1, 2], 3);
        assert!(matches!(r, Err(QRustError::InvalidConfig(_))));
    }

    #[test]
    fn test_layout_from_l2p_partial_mapping() {
        let l = Layout::from_l2p(vec![0, 2], 4).expect("valid partial mapping");
        assert_eq!(l.physical_to_logical[0], 0);
        assert_eq!(l.physical_to_logical[1], usize::MAX);
        assert_eq!(l.physical_to_logical[2], 1);
        assert_eq!(l.physical_to_logical[3], usize::MAX);
    }

    #[test]
    fn test_already_routable_noop() {
        let mut c = Circuit::new(3, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let pass = default_pass(Backend::linear(3), 4, 3, 1);
        let r = pass.run(&c, &mut PropertySet::new());
        let swaps = r
            .operations
            .iter()
            .filter(|op| {
                matches!(
                    op,
                    Operation::Gate {
                        name: GateType::SWAP,
                        ..
                    }
                )
            })
            .count();
        assert_eq!(swaps, 0);
    }

    #[test]
    fn test_route_fidelity_basic_cx() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let backend = Backend::linear(3);
        let pass = default_pass(backend, 1, 1, 1);
        let mut ps = PropertySet::new();
        let routed = pass.run(&c, &mut ps);
        let u_orig = circuit_to_unitary(&c);
        let u_phys = circuit_to_unitary(&routed);
        let initial: Vec<usize> = ps.get::<Vec<usize>>("initial_layout").unwrap().clone();
        let final_l: Vec<usize> = ps.get::<Vec<usize>>("final_layout").unwrap().clone();
        let u_log = extract_logical_unitary(&u_phys, 2, &initial, &final_l);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_invalid_config_returns_err() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let pass = default_pass(Backend::linear(2), 0, 1, 1);
        let mut ps = PropertySet::new();
        let res = pass.try_run(&c, &mut ps);
        assert!(matches!(res, Err(QRustError::InvalidConfig(_))));
    }

    #[test]
    fn test_pass_run_graceful_degradation_on_bad_config() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let pass = default_pass(Backend::linear(2), 0, 0, 0);
        let mut ps = PropertySet::new();
        let out = pass.run(&c, &mut ps);
        assert_eq!(out.operations.len(), c.operations.len());
        assert!(
            ps.get::<String>("beam_sabre_error").is_some(),
            "expected beam_sabre_error in PropertySet on failure path"
        );
    }

    #[test]
    fn test_disconnected_topology_reports_actual_pair() {
        let mut backend = Backend::new("split", 4);
        backend.set_coupling_map([(0, 1), (1, 0), (2, 3), (3, 2)]);
        let mut c = Circuit::new(4, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 3],
            params: vec![],
        });
        let pass = default_pass(backend, 4, 3, 1);
        let mut ps = PropertySet::new();
        let res = pass.try_run(&c, &mut ps);
        match res {
            Err(QRustError::DisconnectedTopology { from, to }) => {
                let witnessed = (from == 0 && to == 3) || (from == 3 && to == 0);
                assert!(
                    witnessed,
                    "expected witness pair to involve qubits 0 and 3, got ({from}, {to})"
                );
            }
            other => panic!("expected DisconnectedTopology, got {other:?}"),
        }
    }

    #[test]
    fn test_fully_connected_fast_path_emits_no_swaps() {
        let mut c = Circuit::new(4, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 3],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 3],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![2, 0],
            params: vec![],
        });
        let pass = default_pass(Backend::all_to_all(4), 4, 3, 2);
        let mut ps = PropertySet::new();
        let out = pass.run(&c, &mut ps);
        let swaps = out
            .operations
            .iter()
            .filter(|op| {
                matches!(
                    op,
                    Operation::Gate {
                        name: GateType::SWAP,
                        ..
                    }
                )
            })
            .count();
        assert_eq!(swaps, 0);
        assert_eq!(*ps.get::<usize>("swaps_inserted").unwrap(), 0);
        let initial: &Vec<usize> = ps.get("initial_layout").unwrap();
        assert_eq!(initial, &vec![0, 1, 2, 3]);
        let final_l: &Vec<usize> = ps.get("final_layout").unwrap();
        assert_eq!(final_l, &vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_fully_connected_fast_path_preserves_num_qubits() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let pass = default_pass(Backend::all_to_all(8), 4, 3, 2);
        let mut ps = PropertySet::new();
        let out = pass.run(&c, &mut ps);
        assert_eq!(out.num_qubits, 2);
        let initial: &Vec<usize> = ps.get("initial_layout").unwrap();
        assert_eq!(initial, &vec![0, 1]);
    }

    #[test]
    fn test_lookahead_strategy_dynamic_v2_routes_correctly() {
        let mut c = Circuit::new(5, 0);
        for i in (1..5).rev() {
            c.add_op(Operation::Gate {
                name: GateType::CX,
                qubits: vec![i, i - 1],
                params: vec![],
            });
        }
        let pass = BeamSabrePass {
            backend: Backend::linear(5),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 2,
            lookahead_strategy: LookaheadStrategy::DynamicV2,
        };
        let mut ps = PropertySet::new();
        let routed = pass.run(&c, &mut ps);
        let u_orig = circuit_to_unitary(&c);
        let mut padded = routed.clone();
        padded.num_qubits = 5;
        let u_log = extract_logical_unitary(
            &circuit_to_unitary(&padded),
            5,
            ps.get::<Vec<usize>>("initial_layout").unwrap(),
            ps.get::<Vec<usize>>("final_layout").unwrap(),
        );
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_lookahead_strategy_static_zero_weight_routes_correctly() {
        let mut c = Circuit::new(3, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 2],
            params: vec![],
        });
        let pass = BeamSabrePass {
            backend: Backend::linear(3),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 1,
            lookahead_strategy: LookaheadStrategy::Static { weight: 0.0 },
        };
        let mut ps = PropertySet::new();
        let routed = pass.run(&c, &mut ps);
        let u_orig = circuit_to_unitary(&c);
        let mut padded = routed.clone();
        padded.num_qubits = 3;
        let u_log = extract_logical_unitary(
            &circuit_to_unitary(&padded),
            3,
            ps.get::<Vec<usize>>("initial_layout").unwrap(),
            ps.get::<Vec<usize>>("final_layout").unwrap(),
        );
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_lookahead_strategy_default_is_classical_sabre() {
        let s = LookaheadStrategy::default();
        assert_eq!(
            s,
            LookaheadStrategy::Static {
                weight: DEFAULT_LOOKAHEAD_WEIGHT
            }
        );
    }

    #[test]
    fn test_beam_sabre_default_impl() {
        let p = BeamSabrePass::default();
        assert_eq!(p.beam_width, 4);
        assert_eq!(p.branch_factor, 3);
        assert_eq!(p.bidir_iterations, 2);
        assert_eq!(p.lookahead_strategy, LookaheadStrategy::default());
    }
}
