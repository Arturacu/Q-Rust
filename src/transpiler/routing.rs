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
struct TwoQGate {
    circuit_index: usize,
    logical_qubits: [usize; 2],
    successors: Vec<usize>,
}

fn build_gate_graph(circuit: &Circuit) -> (Vec<TwoQGate>, Vec<u16>, Vec<Vec<(usize, Operation)>>) {
    let mut two_q: Vec<TwoQGate> = Vec::new();
    let mut single_q: Vec<Vec<(usize, Operation)>> = vec![vec![]; circuit.num_qubits];
    let mut last_2q: Vec<Option<usize>> = vec![None; circuit.num_qubits];

    for (i, op) in circuit.operations.iter().enumerate() {
        match op {
            Operation::Gate { qubits, .. } if qubits.len() >= 2 => {
                let idx = two_q.len();
                two_q.push(TwoQGate {
                    circuit_index: i,
                    logical_qubits: [qubits[0], qubits[1]],
                    successors: Vec::new(),
                });
                for &q in &qubits[..2] {
                    if let Some(prev) = last_2q[q] {
                        two_q[prev].successors.push(idx);
                    }
                    last_2q[q] = Some(idx);
                }
            }
            Operation::Gate { qubits, .. } if qubits.len() == 1 => {
                single_q[qubits[0]].push((i, op.clone()));
            }
            _ => {}
        }
    }

    let mut preds = vec![0u16; two_q.len()];
    for g in &two_q {
        for &s in &g.successors {
            preds[s] += 1;
        }
    }
    (two_q, preds, single_q)
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
    remaining_deps: Vec<u16>,
    front_layer: Vec<usize>,
    cumulative_cost: f64,
    executed_count: usize,
    total_gates: usize,
    complete: bool,
    stuck_swaps: usize,
}

impl Beam {
    fn new(layout: Layout, preds: &[u16], n_gates: usize) -> Self {
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

fn relative_score(
    layout: &Layout,
    swap: (usize, usize),
    front: &[usize],
    gates: &[TwoQGate],
    dist: &[Vec<usize>],
    strategy: LookaheadStrategy,
) -> f64 {
    let (pa, pb) = swap;
    let mut delta_f = 0.0;
    let mut delta_e = 0.0;
    let mut extended_layer: HashSet<usize> = HashSet::new();

    for &g in front {
        let p0 = layout.l2p(gates[g].logical_qubits[0]);
        let p1 = layout.l2p(gates[g].logical_qubits[1]);
        for &succ in &gates[g].successors {
            extended_layer.insert(succ);
        }
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

    for &g in &extended_layer {
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

    match strategy {
        LookaheadStrategy::Static { weight } => delta_f + weight * delta_e,
        LookaheadStrategy::DynamicV2 => {
            if front.is_empty() || extended_layer.is_empty() {
                return delta_f;
            }
            let scaled_e = delta_e * (front.len() as f64 / extended_layer.len() as f64);
            delta_f.max(scaled_e)
        }
    }
}

#[inline]
fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}

fn beam_sabre_forward(
    gates: &[TwoQGate],
    preds: &[u16],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
    beam_width: usize,
    branch_factor: usize,
    strategy: LookaheadStrategy,
) -> Result<Vec<Beam>> {
    if gates.is_empty() {
        return Ok(vec![Beam::new(initial_layout.clone(), preds, 0)]);
    }
    let mut beams = vec![Beam::new(initial_layout.clone(), preds, gates.len())];

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
            let mut candidates: HashSet<(usize, usize)> = HashSet::new();
            for &gate_idx in &beam.front_layer {
                let g = &gates[gate_idx];
                let p0 = beam.layout.l2p(g.logical_qubits[0]);
                let p1 = beam.layout.l2p(g.logical_qubits[1]);
                for n in backend.neighbors(p0) {
                    let edge = if p0 < n { (p0, n) } else { (n, p0) };
                    candidates.insert(edge);
                }
                for n in backend.neighbors(p1) {
                    let edge = if p1 < n { (p1, n) } else { (n, p1) };
                    candidates.insert(edge);
                }
            }
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
                    );
                    (beam.cumulative_cost + SWAP_COST + delta, swap)
                })
                .collect();
            scored.sort_by(|a, b| cmp_f64(a.0, b.0));
            if scored.is_empty() {
                // Beam is stuck (no path in disconnected topology). Drop it.
                continue;
            }
            let best_delta = relative_score(
                &beam.layout,
                scored[0].1,
                &beam.front_layer,
                gates,
                dist,
                strategy,
            );
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
        new_beams.sort_by(|a, b| cmp_f64(a.cumulative_cost, b.cumulative_cost));

        // Extract witness BEFORE dropping stuck beams.
        let pre_prune_witness: Option<(usize, usize)> = beams
            .iter()
            .find(|b| !b.front_layer.is_empty())
            .and_then(|b| {
                b.front_layer.first().map(|&gi| {
                    (
                        gates[gi].logical_qubits[0],
                        gates[gi].logical_qubits[1],
                    )
                })
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
    pub fn try_run(
        &self,
        circuit: &Circuit,
        property_set: &mut PropertySet,
    ) -> Result<Circuit> {
        let has_2q = circuit.operations.iter().any(|op| {
            if let Operation::Gate { qubits, .. } = op {
                qubits.len() >= 2
            } else {
                false
            }
        });
        if !has_2q {
            return Ok(circuit.clone());
        }

        if self.beam_width < 1 || self.branch_factor < 1 || self.bidir_iterations < 1 {
            return Err(QRustError::InvalidConfig(format!(
                "BeamSabrePass: beam_width={}, branch_factor={}, bidir_iterations={} (all must be >= 1)",
                self.beam_width, self.branch_factor, self.bidir_iterations
            )));
        }
        if circuit.num_qubits > self.backend.num_qubits {
            return Err(QRustError::InsufficientQubits {
                circuit: circuit.num_qubits,
                backend: self.backend.num_qubits,
            });
        }

        if self.backend.is_fully_connected() {
            let trivial = Layout::trivial(circuit.num_qubits, self.backend.num_qubits);
            property_set.insert("initial_layout", trivial.logical_to_physical.clone());
            property_set.insert("final_layout", trivial.logical_to_physical.clone());
            property_set.insert("swaps_inserted", 0_usize);
            return Ok(circuit.clone());
        }

        let dist = self.backend.shortest_path_matrix();
        let (gates, preds, single_q_ops) = build_gate_graph(circuit);

        let mut best_layout =
            if let Some(l2p) = property_set.get::<Vec<usize>>("sabre_initial_layout") {
                Layout::from_l2p(l2p.clone(), self.backend.num_qubits).map_err(|e| {
                    QRustError::Routing(format!(
                        "invalid sabre_initial_layout from earlier pass: {e}"
                    ))
                })?
            } else {
                Layout::trivial(circuit.num_qubits, self.backend.num_qubits)
            };

        let mut best_beam: Option<Beam> = None;
        let mut best_cost = f64::MAX;
        let mut best_initial_layout = best_layout.clone();

        let reversed_gates: Vec<TwoQGate> = {
            let mut rg = gates.clone();
            rg.reverse();
            let mut reversed: Vec<TwoQGate> = Vec::new();
            let mut last_on: Vec<Option<usize>> = vec![None; circuit.num_qubits];
            for (i, gate) in rg.iter().enumerate() {
                let new_gate = TwoQGate {
                    circuit_index: gate.circuit_index,
                    logical_qubits: gate.logical_qubits,
                    successors: Vec::new(),
                };
                for &q in &gate.logical_qubits {
                    if let Some(prev) = last_on[q] {
                        reversed[prev].successors.push(i);
                    }
                    last_on[q] = Some(i);
                }
                reversed.push(new_gate);
            }
            reversed
        };
        let reversed_preds: Vec<u16> = {
            let mut counts = vec![0u16; reversed_gates.len()];
            for g in &reversed_gates {
                for &s in &g.successors {
                    counts[s] += 1;
                }
            }
            counts
        };

        for iter in 0..self.bidir_iterations {
            let current_initial = best_layout.clone();
            let fwd = beam_sabre_forward(
                &gates,
                &preds,
                &self.backend,
                &dist,
                &current_initial,
                self.beam_width,
                self.branch_factor,
                self.lookahead_strategy,
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
            }
            if iter + 1 < self.bidir_iterations {
                let bwd = beam_sabre_forward(
                    &reversed_gates,
                    &reversed_preds,
                    &self.backend,
                    &dist,
                    &fwd_best.layout,
                    self.beam_width,
                    self.branch_factor,
                    self.lookahead_strategy,
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
            QRustError::Routing("no routing solution found across all bidirectional iterations".into())
        })?;
        property_set.insert(
            "initial_layout",
            best_initial_layout.logical_to_physical.clone(),
        );
        property_set.insert("final_layout", winning.layout.logical_to_physical.clone());
        property_set.insert("swaps_inserted", winning.cumulative_cost as usize / 3);

        Ok(reconstruct_circuit(
            circuit,
            &winning,
            &best_initial_layout,
            &gates,
            self.backend.num_qubits,
            &single_q_ops,
        ))
    }
}

impl Pass for BeamSabrePass {
    fn name(&self) -> &str {
        "BeamSabrePass"
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