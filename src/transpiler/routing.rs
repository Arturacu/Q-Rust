//! BeamSABRE: Hardware-aware qubit routing with beam search.
//!
//! This module implements a novel hybrid routing algorithm that combines
//! LightSABRE's efficient relative scoring heuristic with beam search
//! exploration. At each routing decision point, instead of greedily
//! committing to the single best SWAP, it maintains K candidate paths
//! (beams) and prunes the worst after each step.
//!
//! When `beam_width=1`, this degenerates exactly to greedy LightSABRE.

use crate::backend::Backend;
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashSet;

// ─── Layout ────────────────────────────────────────────────────────────────

/// Bidirectional mapping between logical and physical qubits.
#[derive(Clone, Debug)]
pub struct Layout {
    /// logical_to_physical[logical_qubit] = physical_qubit
    pub logical_to_physical: Vec<usize>,
    /// physical_to_logical[physical_qubit] = logical_qubit (usize::MAX if unused)
    pub physical_to_logical: Vec<usize>,
}

impl Layout {
    /// Creates a trivial layout: logical qubit i → physical qubit i.
    pub fn trivial(num_logical: usize, num_physical: usize) -> Self {
        let mut l2p = vec![0; num_logical];
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

    /// Swaps the logical qubits assigned to physical qubits p1 and p2.
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

    /// Maps a logical qubit to its current physical qubit.
    pub fn l2p(&self, logical: usize) -> usize {
        self.logical_to_physical[logical]
    }
}

// ─── Gate dependency tracking ──────────────────────────────────────────────

/// Lightweight reference to a two-qubit gate for front layer tracking.
#[derive(Clone, Debug)]
struct TwoQGate {
    /// Index into the original circuit's operation list.
    circuit_index: usize,
    /// The two logical qubits this gate acts on.
    logical_qubits: [usize; 2],
    /// Indices of successor 2q gates (gates that depend on this one).
    successors: Vec<usize>,
}

/// Builds the dependency graph for two-qubit gates in the circuit.
/// Returns (list of TwoQGate, initial predecessor counts, single-qubit ops grouped by qubit).
fn build_gate_graph(
    circuit: &Circuit,
) -> (Vec<TwoQGate>, Vec<u16>, Vec<Vec<(usize, Operation)>>) {
    // Identify all two-qubit gates and their positions
    let mut two_q_gates: Vec<TwoQGate> = Vec::new();
    let mut single_q_ops: Vec<Vec<(usize, Operation)>> = vec![vec![]; circuit.num_qubits];
    // Track which 2q gate was the last to touch each qubit
    let mut last_2q_on_qubit: Vec<Option<usize>> = vec![None; circuit.num_qubits];

    for (i, op) in circuit.operations.iter().enumerate() {
        match op {
            Operation::Gate { qubits, .. } if qubits.len() >= 2 => {
                let gate_idx = two_q_gates.len();
                two_q_gates.push(TwoQGate {
                    circuit_index: i,
                    logical_qubits: [qubits[0], qubits[1]],
                    successors: Vec::new(),
                });

                // Record dependencies: this gate depends on the last 2q gate on each qubit
                for &q in &qubits[..2] {
                    if let Some(prev) = last_2q_on_qubit[q] {
                        // prev is a predecessor of gate_idx. Note: this can produce duplicate entries
                        // if two sequential 2q gates share both qubits (e.g. CX(0,1), CX(0,1)).
                        // This is intentional — mirroring these duplicates in pred_counts ensures
                        // the accounting remains correct through the second execution.
                        two_q_gates[prev].successors.push(gate_idx);
                    }
                    last_2q_on_qubit[q] = Some(gate_idx);
                }
            }
            Operation::Gate { qubits, .. } if qubits.len() == 1 => {
                single_q_ops[qubits[0]].push((i, op.clone()));
            }
            _ => {
                // Measures, resets, barriers — handled separately
            }
        }
    }

    // Compute predecessor counts
    let mut pred_counts = vec![0u16; two_q_gates.len()];
    for gate in &two_q_gates {
        for &succ in &gate.successors {
            pred_counts[succ] += 1;
        }
    }

    (two_q_gates, pred_counts, single_q_ops)
}

// ─── Routing Action ────────────────────────────────────────────────────────

/// A single action in the routing decision log.
#[derive(Clone, Debug)]
enum RoutingAction {
    /// Execute the 2q gate at this index in the two_q_gates list.
    ExecuteGate(usize),
    /// Insert a SWAP between these two physical qubits.
    InsertSwap(usize, usize),
}

// ─── Beam ──────────────────────────────────────────────────────────────────

/// A single candidate routing path maintained during beam search.
#[derive(Clone)]
struct Beam {
    layout: Layout,
    actions: Vec<RoutingAction>,
    /// Per-gate remaining predecessor count (decremented as gates execute).
    remaining_deps: Vec<u16>,
    /// Indices into two_q_gates that are ready to execute.
    front_layer: Vec<usize>,
    /// Cumulative routing cost (total distance penalty from SWAPs).
    cumulative_cost: f64,
    /// How many 2q gates have been executed.
    executed_count: usize,
    /// Total number of 2q gates.
    total_gates: usize,
    /// Whether all 2q gates have been executed.
    complete: bool,
}

impl Beam {
    fn new(layout: Layout, pred_counts: &[u16], num_gates: usize) -> Self {
        let remaining_deps = pred_counts.to_vec();
        // Initial front layer: all gates with zero predecessors
        let front_layer: Vec<usize> = (0..num_gates)
            .filter(|&i| remaining_deps[i] == 0)
            .collect();

        Beam {
            layout,
            actions: Vec::new(),
            remaining_deps,
            front_layer,
            cumulative_cost: 0.0,
            executed_count: 0,
            total_gates: num_gates,
            complete: num_gates == 0,
        }
    }

    /// Execute a gate: remove from front layer, update successor deps.
    fn execute_gate(&mut self, gate_idx: usize, gates: &[TwoQGate]) {
        self.front_layer.retain(|&g| g != gate_idx);
        self.actions.push(RoutingAction::ExecuteGate(gate_idx));

        // Update successors
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
    }

    /// Apply a SWAP on two physical qubits.
    fn apply_swap(&mut self, p1: usize, p2: usize) {
        self.layout.swap_physical(p1, p2);
        self.actions.push(RoutingAction::InsertSwap(p1, p2));
        self.cumulative_cost += SWAP_COST;
    }
}

const SWAP_COST: f64 = 3.0; // A SWAP decomposes to 3 CX gates

// ─── Relative scoring ──────────────────────────────────────────────────────

/// Computes the heuristic delta from applying a SWAP on (pa, pb).
/// Only gates in the front layer whose qubits touch pa or pb contribute.
/// This makes scoring effectively O(1) per candidate in practice.
fn relative_score(
    layout: &Layout,
    swap: (usize, usize),
    front_layer: &[usize],
    gates: &[TwoQGate],
    dist: &[Vec<usize>],
) -> f64 {
    let (pa, pb) = swap;
    let mut delta = 0.0;
    
    // Maintain a set of unique successors for the extended lookahead layer
    let mut extended_layer = std::collections::HashSet::new();

    for &gate_idx in front_layer {
        let gate = &gates[gate_idx];
        let p0 = layout.l2p(gate.logical_qubits[0]);
        let p1 = layout.l2p(gate.logical_qubits[1]);
        
        for &succ in &gate.successors {
            extended_layer.insert(succ);
        }

        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }

        let old_dist = dist[p0][p1] as f64;
        let new_p0 = if p0 == pa { pb } else if p0 == pb { pa } else { p0 };
        let new_p1 = if p1 == pa { pb } else if p1 == pb { pa } else { p1 };
        let new_dist = dist[new_p0][new_p1] as f64;
        delta += new_dist - old_dist;
    }
    
    // Incorporate Extended Lookahead Layer (weight = 0.5)
    let lookahead_weight = 0.5;
    for gate_idx in extended_layer {
        let gate = &gates[gate_idx];
        let p0 = layout.l2p(gate.logical_qubits[0]);
        let p1 = layout.l2p(gate.logical_qubits[1]);
        
        if p0 != pa && p0 != pb && p1 != pa && p1 != pb {
            continue;
        }

        let old_dist = dist[p0][p1] as f64;
        let new_p0 = if p0 == pa { pb } else if p0 == pb { pa } else { p0 };
        let new_p1 = if p1 == pa { pb } else if p1 == pb { pa } else { p1 };
        let new_dist = dist[new_p0][new_p1] as f64;
        delta += (new_dist - old_dist) * lookahead_weight;
    }

    delta
}

/// Computes the absolute heuristic score for the current front layer.
#[allow(dead_code)]
fn absolute_score(
    layout: &Layout,
    front_layer: &[usize],
    gates: &[TwoQGate],
    dist: &[Vec<usize>],
) -> f64 {
    let mut score = 0.0;
    for &gate_idx in front_layer {
        let gate = &gates[gate_idx];
        let p0 = layout.l2p(gate.logical_qubits[0]);
        let p1 = layout.l2p(gate.logical_qubits[1]);
        score += dist[p0][p1] as f64;
    }
    score
}

// ─── BeamSABRE core ────────────────────────────────────────────────────────

/// Runs a single forward pass of BeamSABRE.
fn beam_sabre_forward(
    gates: &[TwoQGate],
    pred_counts: &[u16],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
    beam_width: usize,
    branch_factor: usize,
) -> Vec<Beam> {
    if gates.is_empty() {
        return vec![Beam::new(initial_layout.clone(), pred_counts, 0)];
    }

    let mut beams = vec![Beam::new(initial_layout.clone(), pred_counts, gates.len())];

    loop {
        // Execute all routable gates in each beam
        for beam in &mut beams {
            if beam.complete {
                continue;
            }
            execute_routable(beam, gates, dist);
        }

        // Check if all beams are complete
        if beams.iter().all(|b| b.complete) {
            break;
        }

        let mut new_beams: Vec<Beam> = Vec::new();

        for beam in &beams {
            if beam.complete {
                new_beams.push(beam.clone());
                continue;
            }

            // Generate SWAP candidates: edges adjacent to front-layer qubits
            let mut candidates: HashSet<(usize, usize)> = HashSet::new();
            for &gate_idx in &beam.front_layer {
                let gate = &gates[gate_idx];
                let p0 = beam.layout.l2p(gate.logical_qubits[0]);
                let p1 = beam.layout.l2p(gate.logical_qubits[1]);

                for n in backend.neighbors(p0) {
                    let edge = if p0 < n { (p0, n) } else { (n, p0) };
                    candidates.insert(edge);
                }
                for n in backend.neighbors(p1) {
                    let edge = if p1 < n { (p1, n) } else { (n, p1) };
                    candidates.insert(edge);
                }
            }

            // Score each candidate
            let mut scored: Vec<(f64, (usize, usize))> = candidates
                .iter()
                .map(|&swap| {
                    let delta = relative_score(&beam.layout, swap, &beam.front_layer, gates, dist);
                    (beam.cumulative_cost + SWAP_COST + delta, swap)
                })
                .collect();

            scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            if scored.is_empty() {
                // No candidates — shouldn't happen on a connected graph
                new_beams.push(beam.clone());
                continue;
            }

            // Release valve: if no candidate swap improves the heuristic score,
            // restrict branching to 1 to avoid combinatorial explosion on plateaus.
            let best_delta =
                relative_score(&beam.layout, scored[0].1, &beam.front_layer, gates, dist);

            let effective_branch = if best_delta >= 0.0 {
                // Not making progress — use release valve (take only 1 path)
                1
            } else {
                branch_factor.min(scored.len())
            };

            for i in 0..effective_branch {
                let mut child = beam.clone();
                let (_, swap) = scored[i];
                child.apply_swap(swap.0, swap.1);
                new_beams.push(child);
            }
        }

        // Global prune: keep only the best beam_width beams
        new_beams.sort_by(|a, b| {
            a.cumulative_cost
                .partial_cmp(&b.cumulative_cost)
                .unwrap()
        });
        new_beams.truncate(beam_width);

        beams = new_beams;

        if beams.is_empty() {
            panic!("BeamSABRE: all beams pruned — this should not happen.");
        }
    }

    beams
}

/// Execute all immediately-routable gates in a beam.
fn execute_routable(beam: &mut Beam, gates: &[TwoQGate], dist: &[Vec<usize>]) {
    loop {
        let routable: Vec<usize> = beam
            .front_layer
            .iter()
            .filter(|&&gate_idx| {
                let gate = &gates[gate_idx];
                let p0 = beam.layout.l2p(gate.logical_qubits[0]);
                let p1 = beam.layout.l2p(gate.logical_qubits[1]);
                dist[p0][p1] == 1
            })
            .cloned()
            .collect();

        if routable.is_empty() {
            break;
        }

        for gate_idx in routable {
            beam.execute_gate(gate_idx, gates);
        }
    }
}

// ─── Circuit reconstruction ────────────────────────────────────────────────

/// Reconstructs a physical circuit from the winning beam's action log.
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

    // Non-gate ops are emitted just before the first 2q gate whose circuit_index
    // exceeds their original position. If the router executes two parallel 2q gates
    // out of their original order, a non-gate op falling between them will be emitted
    // before both — this is an accepted limitation for circuits with mid-circuit
    // measurements between parallel gates, which the current IR does not support.
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

                let orig_op = &original.operations[gate.circuit_index];
                if let Operation::Gate { name, params, .. } = orig_op {
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
        Operation::Measure { qubit, cbit } => {
            output.add_op(Operation::Measure {
                qubit: layout.l2p(*qubit),
                cbit: *cbit,
            });
        }
        Operation::Reset { qubit } => {
            output.add_op(Operation::Reset {
                qubit: layout.l2p(*qubit),
            });
        }
        Operation::Barrier { qubits } => {
            let phys: Vec<usize> = qubits.iter().map(|&q| layout.l2p(q)).collect();
            output.add_op(Operation::Barrier { qubits: phys });
        }
        _ => {}
    }
}

// ─── BeamSabrePass ─────────────────────────────────────────────────────────

/// Hardware-aware routing pass using the BeamSABRE algorithm.
///
/// Combines LightSABRE's O(1) relative scoring heuristic with beam search
/// exploration. At beam_width=1, this degenerates to greedy LightSABRE.
pub struct BeamSabrePass {
    pub backend: Backend,
    pub beam_width: usize,
    pub branch_factor: usize,
    pub bidir_iterations: usize,
}

impl Pass for BeamSabrePass {
    fn name(&self) -> &str {
        "BeamSabrePass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit {
        // Skip if circuit has no multi-qubit gates
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

        assert!(self.beam_width >= 1, "BeamSabrePass: beam_width must be >= 1, got {}", self.beam_width);
        assert!(self.branch_factor >= 1, "BeamSabrePass: branch_factor must be >= 1, got {}", self.branch_factor);
        assert!(self.bidir_iterations >= 1, "BeamSabrePass: bidir_iterations must be >= 1, got {}", self.bidir_iterations);

        // Validate qubit count
        assert!(
            circuit.num_qubits <= self.backend.num_qubits,
            "Circuit has {} qubits but backend only has {}",
            circuit.num_qubits,
            self.backend.num_qubits
        );

        // Precompute distance matrix
        let dist = self.backend.shortest_path_matrix();

        // Build gate dependency graph
        let (gates, pred_counts, single_q_ops) = build_gate_graph(circuit);

        // Use layout from SabreLayoutPass if available, otherwise trivial
        let mut best_layout =
            if let Some(l2p) = property_set.get::<Vec<usize>>("sabre_initial_layout") {
                let mut p2l = vec![usize::MAX; self.backend.num_qubits];
                for (logical, &physical) in l2p.iter().enumerate() {
                    p2l[physical] = logical;
                }
                Layout {
                    logical_to_physical: l2p.clone(),
                    physical_to_logical: p2l,
                }
            } else {
                Layout::trivial(circuit.num_qubits, self.backend.num_qubits)
            };

        // Bidirectional refinement
        let mut best_beam: Option<Beam> = None;
        let mut best_cost = f64::MAX;
        let mut best_initial_layout = best_layout.clone();

        // For bidirectional: we'll reverse the gate order.
        // NOTE: circuit_index values in reversed_gates retain their original
        // forward-pass positions. This vector is used only for cost estimation,
        // never for circuit reconstruction — do not pass it to reconstruct_circuit.
        let reversed_gates: Vec<TwoQGate> = {
            let mut rg = gates.clone();
            rg.reverse();
            // Rebuild successor links for reversed order
            let mut reversed: Vec<TwoQGate> = Vec::new();
            let mut last_2q_on_qubit: Vec<Option<usize>> = vec![None; circuit.num_qubits];
            for (i, gate) in rg.iter().enumerate() {
                let new_gate = TwoQGate {
                    circuit_index: gate.circuit_index,
                    logical_qubits: gate.logical_qubits,
                    successors: Vec::new(),
                };
                for &q in &gate.logical_qubits {
                    if let Some(prev) = last_2q_on_qubit[q] {
                        reversed[prev].successors.push(i);
                    }
                    last_2q_on_qubit[q] = Some(i);
                }
                reversed.push(new_gate);
            }
            reversed
        };

        let reversed_pred_counts: Vec<u16> = {
            let mut counts = vec![0u16; reversed_gates.len()];
            for gate in &reversed_gates {
                for &succ in &gate.successors {
                    counts[succ] += 1;
                }
            }
            counts
        };

        for iteration in 0..self.bidir_iterations {
            let current_initial = best_layout.clone();

            // Forward pass
            let forward_beams = beam_sabre_forward(
                &gates,
                &pred_counts,
                &self.backend,
                &dist,
                &current_initial,
                self.beam_width,
                self.branch_factor,
            );

            let forward_best = forward_beams
                .iter()
                .min_by(|a, b| a.cumulative_cost.partial_cmp(&b.cumulative_cost).unwrap())
                .unwrap();

            if forward_best.cumulative_cost < best_cost {
                best_cost = forward_best.cumulative_cost;
                best_beam = Some(forward_best.clone());
                best_initial_layout = current_initial.clone();
            }

            // Backward pass: use forward's final layout as backward's initial
            if iteration + 1 < self.bidir_iterations {
                let backward_beams = beam_sabre_forward(
                    &reversed_gates,
                    &reversed_pred_counts,
                    &self.backend,
                    &dist,
                    &forward_best.layout,
                    self.beam_width,
                    self.branch_factor,
                );

                let backward_best = backward_beams
                    .iter()
                    .min_by(|a, b| {
                        a.cumulative_cost.partial_cmp(&b.cumulative_cost).unwrap()
                    })
                    .unwrap();

                // The backward pass's final layout becomes the next forward pass's initial
                best_layout = backward_best.layout.clone();
            }
        }

        let winning_beam = best_beam.expect("No routing solution found");

        // Store layouts in the property set for downstream passes and verification
        property_set.insert(
            "initial_layout",
            best_initial_layout.logical_to_physical.clone(),
        );
        property_set.insert(
            "final_layout",
            winning_beam.layout.logical_to_physical.clone(),
        );
        property_set.insert("swaps_inserted", winning_beam.cumulative_cost as usize / 3);

        // Reconstruct the circuit using the initial layout that produced this beam
        reconstruct_circuit(
            circuit,
            &winning_beam,
            &best_initial_layout,
            &gates,
            self.backend.num_qubits,
            &single_q_ops,
        )
    }
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Circuit, GateType, Operation};
    use crate::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};

    #[test]
    fn test_layout_trivial() {
        let layout = Layout::trivial(3, 5);
        assert_eq!(layout.l2p(0), 0);
        assert_eq!(layout.l2p(1), 1);
        assert_eq!(layout.l2p(2), 2);
        assert_eq!(layout.physical_to_logical[3], usize::MAX);
    }

    #[test]
    fn test_layout_swap() {
        let mut layout = Layout::trivial(3, 3);
        layout.swap_physical(0, 2);
        // Logical 0 was at physical 0, now at physical 2
        assert_eq!(layout.l2p(0), 2);
        // Logical 2 was at physical 2, now at physical 0
        assert_eq!(layout.l2p(2), 0);
        // Logical 1 unchanged
        assert_eq!(layout.l2p(1), 1);
        // Reverse mapping
        assert_eq!(layout.physical_to_logical[0], 2);
        assert_eq!(layout.physical_to_logical[2], 0);
    }

    #[test]
    fn test_layout_swap_consistency() {
        let mut layout = Layout::trivial(5, 5);
        // Perform many swaps and verify invariant
        let swaps = vec![(0, 1), (2, 3), (1, 4), (0, 3), (2, 1)];
        for (a, b) in swaps {
            layout.swap_physical(a, b);
            // Verify bidirectional consistency
            for l in 0..5 {
                let p = layout.l2p(l);
                assert_eq!(layout.physical_to_logical[p], l);
            }
        }
    }

    #[test]
    fn test_already_routable_noop() {
        // CX on adjacent qubits on linear(3): should insert 0 SWAPs
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let backend = Backend::linear(3);
        let pass = BeamSabrePass {
            backend,
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 1,
        };

        let routed = pass.run(&mut circuit, &mut PropertySet::new());

        // Should have exactly 1 gate (the CX), no SWAPs
        let swap_count = routed
            .operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { name: GateType::SWAP, .. }))
            .count();
        assert_eq!(swap_count, 0);
    }

    #[test]
    fn test_route_non_adjacent_linear() {
        // CX(0, 2) on linear(3): qubits 0 and 2 are not adjacent.
        // The router should either:
        // - Insert SWAPs to make them adjacent, OR
        // - Find an initial layout where they're already adjacent (bidir refinement)
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 2],
            params: vec![],
        });

        let backend = Backend::linear(3);
        let pass = BeamSabrePass {
            backend: backend.clone(),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 2,
        };

        let routed = pass.run(&circuit, &mut PropertySet::new());

        // All 2q gates must be on adjacent physical qubits
        for op in &routed.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                if qubits.len() == 2 {
                    assert!(
                        backend.is_adjacent(qubits[0], qubits[1]),
                        "Gate {:?} on qubits {:?} are not adjacent!",
                        name,
                        qubits
                    );
                }
            }
        }

        // Should contain at least one CX gate (possibly remapped)
        let cx_count = routed
            .operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { name: GateType::CX, .. }))
            .count();
        assert!(cx_count >= 1, "Expected at least 1 CX gate in output");
    }

    #[test]
    fn test_route_forces_swap_no_bidir() {
        // CX(0, 2) on linear(3) with NO bidirectional refinement (bidir=1):
        // Trivial layout makes qubits 0 and 2 non-adjacent, so a SWAP is mandatory.
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 2],
            params: vec![],
        });

        let backend = Backend::linear(3);
        let pass = BeamSabrePass {
            backend: backend.clone(),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 1, // No bidir — must use trivial layout
        };

        let routed = pass.run(&circuit, &mut PropertySet::new());

        // With trivial layout on linear(3), CX(0,2) needs a SWAP
        let swap_count = routed
            .operations
            .iter()
            .filter(|op| matches!(op, Operation::Gate { name: GateType::SWAP, .. }))
            .count();
        assert!(swap_count >= 1, "Expected at least 1 SWAP with trivial layout, got {}", swap_count);

        // All 2q gates must be on adjacent physical qubits
        for op in &routed.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                if qubits.len() == 2 {
                    assert!(
                        backend.is_adjacent(qubits[0], qubits[1]),
                        "Gate {:?} on qubits {:?} are not adjacent!",
                        name,
                        qubits
                    );
                }
            }
        }
    }

    #[test]
    fn test_beam_width_1_is_greedy() {
        // beam_width=1 should not crash and should produce a valid result
        let mut circuit = Circuit::new(4, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 3],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 2],
            params: vec![],
        });

        let backend = Backend::linear(4);
        let pass = BeamSabrePass {
            backend: backend.clone(),
            beam_width: 1,
            branch_factor: 1,
            bidir_iterations: 1,
        };

        let routed = pass.run(&circuit, &mut PropertySet::new());

        // All 2q gates should be on adjacent qubits
        for op in &routed.operations {
            if let Operation::Gate { qubits, .. } = op {
                if qubits.len() == 2 {
                    assert!(backend.is_adjacent(qubits[0], qubits[1]));
                }
            }
        }
    }

    #[test]
    fn test_single_qubit_only_passthrough() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![1],
            params: vec![0.5],
        });

        let backend = Backend::linear(2);
        let pass = BeamSabrePass {
            backend,
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 1,
        };

        let routed = pass.run(&circuit, &mut PropertySet::new());

        // Should pass through unchanged (no 2q gates to route)
        assert_eq!(routed.operations.len(), 2);
    }

    #[test]
    fn test_route_fidelity_basic_cx() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
        
        // linear(2) is trivial, but let's use linear(3) and force it to qubits 0,1
        let backend = Backend::linear(3);
        let pass = BeamSabrePass { backend, beam_width: 1, branch_factor: 1, bidir_iterations: 1 };
        
        let mut ps = PropertySet::new();
        let routed = pass.run(&circuit, &mut ps);
        
        let u_orig = circuit_to_unitary(&circuit);
        let u_phys = circuit_to_unitary(&routed);
        
        let initial_layout: Vec<usize> = ps.get::<Vec<usize>>("initial_layout").unwrap().clone();
        let final_layout: Vec<usize> = ps.get::<Vec<usize>>("final_layout").unwrap().clone();
        
        let u_log = extract_logical_unitary(&u_phys, 2, &initial_layout, &final_layout);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_route_fidelity_with_swaps() {
        // CX(0, 2) on linear(3) - requires at least one SWAP
        let mut circuit = Circuit::new(3, 0);
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 2], params: vec![] });
        
        let backend = Backend::linear(3);
        let pass = BeamSabrePass { backend, beam_width: 4, branch_factor: 3, bidir_iterations: 1 };
        
        let mut ps = PropertySet::new();
        let routed = pass.run(&circuit, &mut ps);
        
        let u_orig = circuit_to_unitary(&circuit);
        let u_phys = circuit_to_unitary(&routed);
        
        let initial: &Vec<usize> = ps.get("initial_layout").unwrap();
        let final_l: &Vec<usize> = ps.get("final_layout").unwrap();
        
        let u_log = extract_logical_unitary(&u_phys, 3, initial, final_l);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_route_preserves_non_gates() {
        let mut circuit = Circuit::new(2, 1);
        circuit.add_op(Operation::Gate { name: GateType::H, qubits: vec![0], params: vec![] });
        circuit.add_op(Operation::Barrier { qubits: vec![0, 1] });
        circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        
        let backend = Backend::linear(2);
        let pass = BeamSabrePass { backend, beam_width: 1, branch_factor: 1, bidir_iterations: 1 };
        let mut ps = PropertySet::new();
        let routed = pass.run(&circuit, &mut ps);
        
        // Check for Barrier and Measure presence
        let barrier = routed.operations.iter().any(|op| matches!(op, Operation::Barrier { .. }));
        let measure = routed.operations.iter().any(|op| matches!(op, Operation::Measure { .. }));
        assert!(barrier);
        assert!(measure);
    }

    #[test]
    fn test_route_large_backend_mapping() {
        // 2 logical qubits on 5 physical qubits
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
        
        let backend = Backend::linear(5);
        let pass = BeamSabrePass { backend, beam_width: 2, branch_factor: 2, bidir_iterations: 1 };
        let mut ps = PropertySet::new();
        let routed = pass.run(&circuit, &mut ps);
        
        assert_eq!(routed.num_qubits, 5);
        
        let u_orig = circuit_to_unitary(&circuit);
        let u_phys = circuit_to_unitary(&routed);
        let initial = ps.get::<Vec<usize>>("initial_layout").unwrap();
        let final_l = ps.get::<Vec<usize>>("final_layout").unwrap();
        let u_log = extract_logical_unitary(&u_phys, 2, initial, final_l);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_route_mid_circuit_measurement_ordering() {
        // Circuit: CX(0,1), Measure(0), CX(0,1) on linear(3).
        let mut circuit = Circuit::new(2, 1);
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
        circuit.add_op(Operation::Measure { qubit: 0, cbit: 0 });
        circuit.add_op(Operation::Gate { name: GateType::CX, qubits: vec![0, 1], params: vec![] });
        
        let backend = Backend::linear(3);
        let pass = BeamSabrePass { backend, beam_width: 1, branch_factor: 1, bidir_iterations: 1 };
        let mut ps = PropertySet::new();
        let routed = pass.run(&circuit, &mut ps);
        
        // Verify order: Gate, Measure, Gate
        let ops = &routed.operations;
        assert_eq!(
            ops.len(), 3,
            "Expected exactly 3 operations (CX, Measure, CX), got: {:?}",
            ops
        );
        assert!(matches!(ops[0], Operation::Gate { .. }));
        assert!(matches!(ops[1], Operation::Measure { .. }));
        assert!(matches!(ops[2], Operation::Gate { .. }));
    }
}
