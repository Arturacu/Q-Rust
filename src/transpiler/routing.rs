//! BeamSABRE: hardware-aware routing with beam search.
//!
//! Combines LightSABRE's relative-scoring heuristic with beam-search
//! exploration. At `beam_width = 1` this degenerates to greedy LightSABRE.

use crate::backend::Backend;
use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;
use std::collections::HashSet;

// ─── Layout ────────────────────────────────────────────────────────────────

/// Bidirectional logical↔physical qubit mapping.
#[derive(Clone, Debug)]
pub struct Layout {
    /// `logical_to_physical[l]` = physical qubit holding logical qubit `l`.
    pub logical_to_physical: Vec<usize>,
    /// `physical_to_logical[p]` = logical qubit at physical `p`, or `usize::MAX`.
    pub physical_to_logical: Vec<usize>,
}

impl Layout {
    /// Trivial layout: logical `i` → physical `i` for `0..num_logical`.
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

    /// Swaps the logical qubits placed at physical `p1` and `p2`.
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

    /// Logical → physical lookup.
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
    }

    fn apply_swap(&mut self, p1: usize, p2: usize) {
        self.layout.swap_physical(p1, p2);
        self.actions.push(RoutingAction::InsertSwap(p1, p2));
        self.cumulative_cost += SWAP_COST;
    }
}

const SWAP_COST: f64 = 3.0;

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

fn beam_sabre_forward(
    gates: &[TwoQGate],
    preds: &[u16],
    backend: &Backend,
    dist: &[Vec<usize>],
    initial_layout: &Layout,
    beam_width: usize,
    branch_factor: usize,
) -> Vec<Beam> {
    if gates.is_empty() {
        return vec![Beam::new(initial_layout.clone(), preds, 0)];
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
                    let delta = relative_score(&beam.layout, swap, &beam.front_layer, gates, dist);
                    (beam.cumulative_cost + SWAP_COST + delta, swap)
                })
                .collect();
            scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            if scored.is_empty() {
                new_beams.push(beam.clone());
                continue;
            }
            let best_delta =
                relative_score(&beam.layout, scored[0].1, &beam.front_layer, gates, dist);
            let effective_branch = if best_delta >= 0.0 {
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
        new_beams.sort_by(|a, b| a.cumulative_cost.partial_cmp(&b.cumulative_cost).unwrap());
        new_beams.truncate(beam_width);
        beams = new_beams;
        if beams.is_empty() {
            panic!("BeamSABRE: all beams pruned — unreachable on a connected graph");
        }
    }
    beams
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

/// Hardware-aware routing pass using the BeamSABRE algorithm.
#[derive(Debug, Clone)]
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

        assert!(self.beam_width >= 1);
        assert!(self.branch_factor >= 1);
        assert!(self.bidir_iterations >= 1);
        assert!(circuit.num_qubits <= self.backend.num_qubits);

        let dist = self.backend.shortest_path_matrix();
        let (gates, preds, single_q_ops) = build_gate_graph(circuit);

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
            );
            let fwd_best = fwd
                .iter()
                .min_by(|a, b| a.cumulative_cost.partial_cmp(&b.cumulative_cost).unwrap())
                .unwrap();
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
                );
                let bwd_best = bwd
                    .iter()
                    .min_by(|a, b| a.cumulative_cost.partial_cmp(&b.cumulative_cost).unwrap())
                    .unwrap();
                best_layout = bwd_best.layout.clone();
            }
        }

        let winning = best_beam.expect("no routing solution");
        property_set.insert(
            "initial_layout",
            best_initial_layout.logical_to_physical.clone(),
        );
        property_set.insert("final_layout", winning.layout.logical_to_physical.clone());
        property_set.insert("swaps_inserted", winning.cumulative_cost as usize / 3);

        reconstruct_circuit(
            circuit,
            &winning,
            &best_initial_layout,
            &gates,
            self.backend.num_qubits,
            &single_q_ops,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{Circuit, GateType, Operation};
    use crate::simulator::{circuit_to_unitary, extract_logical_unitary, unitary_fidelity};

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
    fn test_already_routable_noop() {
        let mut c = Circuit::new(3, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let pass = BeamSabrePass {
            backend: Backend::linear(3),
            beam_width: 4,
            branch_factor: 3,
            bidir_iterations: 1,
        };
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
        let pass = BeamSabrePass {
            backend,
            beam_width: 1,
            branch_factor: 1,
            bidir_iterations: 1,
        };
        let mut ps = PropertySet::new();
        let routed = pass.run(&c, &mut ps);
        let u_orig = circuit_to_unitary(&c);
        let u_phys = circuit_to_unitary(&routed);
        let initial: Vec<usize> = ps.get::<Vec<usize>>("initial_layout").unwrap().clone();
        let final_l: Vec<usize> = ps.get::<Vec<usize>>("final_layout").unwrap().clone();
        let u_log = extract_logical_unitary(&u_phys, 2, &initial, &final_l);
        assert!((unitary_fidelity(&u_orig, &u_log) - 1.0).abs() < 1e-9);
    }
}
// /// Beam SABRE routing implementation
