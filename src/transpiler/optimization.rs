//! Optimization passes.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::dag::{DAGCircuit, DAGNode};
use crate::transpiler::pass::Pass;
use crate::transpiler::synthesis::zyz::{u_to_matrix, zyz_decomposition};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Returns true iff `op` is a barrier (passes should not reorder across these).
#[inline]
fn is_barrier(op: &Operation) -> bool {
    matches!(op, Operation::Barrier { .. })
}

/// Fuses adjacent single-qubit `U` gates into one `U`.
#[derive(Debug, Clone, Copy)]
pub struct GateFusionPass;

impl Pass for GateFusionPass {
    fn name(&self) -> &str {
        "GateFusionPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        // Process in segments separated by barriers, so we don't fuse across them.
        run_in_segments(circuit, |seg| fuse_segment(seg))
    }
}

fn fuse_segment(circuit: &Circuit) -> Circuit {
    let mut dag = DAGCircuit::from(circuit);
    let mut progress = true;

    while progress {
        progress = false;
        let edge_indices: Vec<_> = dag.graph.edge_indices().collect();

        for &edge in &edge_indices {
            if dag.graph.edge_weight(edge).is_none() {
                continue;
            }
            let Some((src, dst)) = dag.graph.edge_endpoints(edge) else {
                continue;
            };

            let src_u = match &dag.graph[src] {
                DAGNode::Op(Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                }) if qubits.len() == 1 => Some((qubits[0], params.clone())),
                _ => None,
            };
            let dst_u = match &dag.graph[dst] {
                DAGNode::Op(Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                }) if qubits.len() == 1 => Some((qubits[0], params.clone())),
                _ => None,
            };

            if let (Some((sq, sp)), Some((dq, dp))) = (src_u, dst_u) {
                if sq == dq {
                    let m = u_to_matrix(sp[0], sp[1], sp[2]);
                    let n = u_to_matrix(dp[0], dp[1], dp[2]);
                    let c00 = n[0][0] * m[0][0] + n[0][1] * m[1][0];
                    let c01 = n[0][0] * m[0][1] + n[0][1] * m[1][1];
                    let c10 = n[1][0] * m[0][0] + n[1][1] * m[1][0];
                    let c11 = n[1][0] * m[0][1] + n[1][1] * m[1][1];
                    let c_matrix = [[c00, c01], [c10, c11]];
                    let (theta, phi, lambda, _) = zyz_decomposition(c_matrix);

                    dag.graph[src] = DAGNode::Op(Operation::Gate {
                        name: GateType::U,
                        qubits: vec![sq],
                        params: vec![theta, phi, lambda],
                    });
                    dag.remove_node(dst);
                    progress = true;
                    break;
                }
            }
        }
    }
    Circuit::from(&dag)
}

/// Splits `circuit` into segments separated by `Barrier` ops, runs `f` on each,
/// and re-assembles. Barriers are preserved.
fn run_in_segments<F: Fn(&Circuit) -> Circuit>(circuit: &Circuit, f: F) -> Circuit {
    let mut result = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    result.custom_gates = circuit.custom_gates.clone();
    let mut current = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    current.custom_gates = circuit.custom_gates.clone();

    for op in &circuit.operations {
        if is_barrier(op) {
            let processed = f(&current);
            for p in processed.operations {
                result.add_op(p);
            }
            result.add_op(op.clone());
            current = Circuit::new(circuit.num_qubits, circuit.num_cbits);
            current.custom_gates = circuit.custom_gates.clone();
        } else {
            current.add_op(op.clone());
        }
    }
    let processed = f(&current);
    for p in processed.operations {
        result.add_op(p);
    }
    result
}

/// Cancels pairs of CX gates.
#[derive(Debug, Clone, Copy)]
pub struct CommutationCancellationPass;

impl Pass for CommutationCancellationPass {
    fn name(&self) -> &str {
        "CommutationCancellationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        run_in_segments(circuit, |seg| {
            let mut dag = DAGCircuit::from(seg);
            let mut progress = true;

            while progress {
                progress = false;
                let node_indices: Vec<_> = dag.graph.node_indices().collect();

                for &idx in &node_indices {
                    if dag.graph.node_weight(idx).is_none() {
                        continue;
                    }
                    if let DAGNode::Op(Operation::Gate {
                        name: GateType::CX,
                        qubits,
                        ..
                    }) = &dag.graph[idx]
                    {
                        let ctrl = qubits[0];
                        let target = qubits[1];
                        if let Some(cancel_idx) = find_commuting_cx(&dag, idx, ctrl, target) {
                            dag.remove_node(idx);
                            dag.remove_node(cancel_idx);
                            progress = true;
                            break;
                        }
                    }
                }
            }
            Circuit::from(&dag)
        })
    }
}

fn find_commuting_cx(
    dag: &DAGCircuit,
    start: NodeIndex,
    ctrl_wire: usize,
    target_wire: usize,
) -> Option<NodeIndex> {
    let mut current = start;
    let candidate_idx = loop {
        let next = walk_forward_on_wire(dag, current, ctrl_wire)?;
        match &dag.graph[next] {
            DAGNode::Op(Operation::Gate {
                name: GateType::CX,
                qubits,
                ..
            }) if qubits[0] == ctrl_wire && qubits[1] == target_wire => break next,
            DAGNode::Op(Operation::Gate {
                name: GateType::CX, ..
            }) => return None,
            DAGNode::Op(op) => {
                if !commutes_with_cx(op, ctrl_wire, target_wire) {
                    return None;
                }
            }
            _ => return None,
        }
        current = next;
    };

    current = start;
    loop {
        let next = walk_forward_on_wire(dag, current, target_wire)?;
        if next == candidate_idx {
            return Some(candidate_idx);
        }
        match &dag.graph[next] {
            DAGNode::Op(op) => {
                if !commutes_with_cx(op, ctrl_wire, target_wire) {
                    return None;
                }
            }
            _ => return None,
        }
        current = next;
    }
}

fn walk_forward_on_wire(dag: &DAGCircuit, current: NodeIndex, wire: usize) -> Option<NodeIndex> {
    for edge in dag
        .graph
        .edges_directed(current, petgraph::Direction::Outgoing)
    {
        if edge.weight().index == wire {
            return Some(edge.target());
        }
    }
    None
}

fn commutes_with_cx(op: &Operation, ctrl: usize, target: usize) -> bool {
    match op {
        Operation::Gate {
            name,
            qubits,
            params,
        } => {
            let involves_ctrl = qubits.contains(&ctrl);
            let involves_target = qubits.contains(&target);

            if !involves_ctrl && !involves_target {
                return true;
            }

            if qubits.len() > 1 {
                return false;
            }

            let q = qubits[0];
            if q == ctrl {
                match name {
                    GateType::Z
                    | GateType::S
                    | GateType::Sdg
                    | GateType::T
                    | GateType::Tdg
                    | GateType::RZ
                    | GateType::ID => true,
                    GateType::U => !params.is_empty() && params[0].abs() < 1e-9,
                    _ => false,
                }
            } else if q == target {
                match name {
                    GateType::X | GateType::RX | GateType::ID => true,
                    GateType::U if params.len() >= 3 => {
                        let phi = params[1];
                        let lambda = params[2];
                        let pi_2 = std::f64::consts::FRAC_PI_2;
                        (phi - -pi_2).abs() < 1e-9 && (lambda - pi_2).abs() < 1e-9
                    }
                    _ => false,
                }
            } else {
                true
            }
        }
        _ => false,
    }
}

/// Removes SWAP(a,b)…SWAP(a,b) pairs with no intervening ops on `{a,b}`.
///
/// O(N) via a single forward scan using a per-qubit "last swap index" map.
#[derive(Debug, Clone, Copy)]
pub struct SwapSimplificationPass;

impl Pass for SwapSimplificationPass {
    fn name(&self) -> &str {
        "SwapSimplificationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        run_in_segments(circuit, |seg| swap_simplify_segment(seg))
    }
}

fn swap_simplify_segment(circuit: &Circuit) -> Circuit {
    let n = circuit.operations.len();
    let mut alive = vec![true; n];
    // For each qubit, the index of the last pending SWAP that touches it.
    let mut pending_swap_by_qubit: HashMap<usize, usize> = HashMap::new();

    for i in 0..n {
        if !alive[i] {
            continue;
        }
        let op = &circuit.operations[i];
        match op {
            Operation::Gate {
                name: GateType::SWAP,
                qubits,
                ..
            } if qubits.len() == 2 => {
                let (a, b) = (qubits[0], qubits[1]);
                // Is there a pending swap on BOTH a and b pointing to the same earlier SWAP?
                let pa = pending_swap_by_qubit.get(&a).copied();
                let pb = pending_swap_by_qubit.get(&b).copied();
                if let (Some(pa), Some(pb)) = (pa, pb) {
                    if pa == pb {
                        if let Operation::Gate {
                            name: GateType::SWAP,
                            qubits: pq,
                            ..
                        } = &circuit.operations[pa]
                        {
                            if pq.len() == 2
                                && ((pq[0] == a && pq[1] == b) || (pq[0] == b && pq[1] == a))
                            {
                                // Cancel both.
                                alive[pa] = false;
                                alive[i] = false;
                                pending_swap_by_qubit.remove(&a);
                                pending_swap_by_qubit.remove(&b);
                                continue;
                            }
                        }
                    }
                }
                // Otherwise this SWAP becomes the new pending candidate on both wires.
                pending_swap_by_qubit.insert(a, i);
                pending_swap_by_qubit.insert(b, i);
            }
            _ => {
                // Any op that touches q invalidates the pending swap on q.
                let touched: Vec<usize> = match op {
                    Operation::Gate { qubits, .. } | Operation::Barrier { qubits } => {
                        qubits.clone()
                    }
                    Operation::Measure { qubit, .. } | Operation::Reset { qubit } => vec![*qubit],
                    Operation::Conditional { op, .. } => op.qubits().to_vec(),
                };
                for q in touched {
                    pending_swap_by_qubit.remove(&q);
                }
            }
        }
    }

    let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    out.custom_gates = circuit.custom_gates.clone();
    for (i, op) in circuit.operations.iter().enumerate() {
        if alive[i] {
            out.add_op(op.clone());
        }
    }
    out
}

fn involves_any(op: &Operation, qubits: &[usize]) -> bool {
    match op {
        Operation::Gate { qubits: oq, .. } => oq.iter().any(|q| qubits.contains(q)),
        Operation::Measure { qubit, .. } | Operation::Reset { qubit } => qubits.contains(qubit),
        Operation::Barrier { qubits: bq } => bq.iter().any(|q| qubits.contains(q)),
        Operation::Conditional { op, .. } => involves_any(op, qubits),
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ParameterSimplificationPass {
    pub epsilon: f64,
}

impl Default for ParameterSimplificationPass {
    fn default() -> Self {
        Self { epsilon: 1e-9 }
    }
}

impl Pass for ParameterSimplificationPass {
    fn name(&self) -> &str {
        "ParameterSimplificationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        let two_pi = 2.0 * std::f64::consts::PI;

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } => {
                    let mut np = params.clone();
                    let mut keep = true;
                    if matches!(
                        name,
                        GateType::RX | GateType::RY | GateType::RZ | GateType::U
                    ) {
                        for p in &mut np {
                            *p %= two_pi;
                            if *p < 0.0 {
                                *p += two_pi;
                            }
                        }
                        match name {
                            GateType::RX | GateType::RY | GateType::RZ => {
                                if let Some(theta) = np.first() {
                                    if theta.abs() < self.epsilon
                                        || (theta - two_pi).abs() < self.epsilon
                                    {
                                        keep = false;
                                    }
                                }
                            }
                            GateType::U if np.len() >= 3 => {
                                let theta = np[0];
                                if theta.abs() < self.epsilon {
                                    let phase_sum = (np[1] + np[2]) % two_pi;
                                    if phase_sum.abs() < self.epsilon
                                        || (phase_sum - two_pi).abs() < self.epsilon
                                    {
                                        keep = false;
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                    if keep {
                        out.add_op(Operation::Gate {
                            name: name.clone(),
                            qubits: qubits.clone(),
                            params: np,
                        });
                    }
                }
                other => out.add_op(other.clone()),
            }
        }
        out
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GateCrystallizationPass {
    pub epsilon: f64,
}

impl Default for GateCrystallizationPass {
    fn default() -> Self {
        Self { epsilon: 1e-9 }
    }
}

impl Pass for GateCrystallizationPass {
    fn name(&self) -> &str {
        "GateCrystallizationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        let pi = std::f64::consts::PI;

        for op in &circuit.operations {
            if let Operation::Gate {
                name,
                qubits,
                params,
            } = op
            {
                let mut new_name = name.clone();
                let mut new_params = params.clone();

                if params.len() == 1 {
                    let t = params[0];
                    let is_pi = (t - pi).abs() < self.epsilon || (t + pi).abs() < self.epsilon;
                    let is_pi2 = (t - pi / 2.0).abs() < self.epsilon;
                    let is_mpi2 = (t + pi / 2.0).abs() < self.epsilon;
                    let is_pi4 = (t - pi / 4.0).abs() < self.epsilon;
                    let is_mpi4 = (t + pi / 4.0).abs() < self.epsilon;

                    match name {
                        GateType::RX if is_pi => {
                            new_name = GateType::X;
                            new_params.clear();
                        }
                        GateType::RY if is_pi => {
                            new_name = GateType::Y;
                            new_params.clear();
                        }
                        GateType::RZ => {
                            if is_pi {
                                new_name = GateType::Z;
                                new_params.clear();
                            } else if is_pi2 {
                                new_name = GateType::S;
                                new_params.clear();
                            } else if is_mpi2 {
                                new_name = GateType::Sdg;
                                new_params.clear();
                            } else if is_pi4 {
                                new_name = GateType::T;
                                new_params.clear();
                            } else if is_mpi4 {
                                new_name = GateType::Tdg;
                                new_params.clear();
                            }
                        }
                        _ => {}
                    }
                } else if params.len() == 3 && name == &GateType::U {
                    let theta = params[0];
                    let phi = params[1];
                    let lambda = params[2];
                    if (theta - pi / 2.0).abs() < self.epsilon
                        && phi.abs() < self.epsilon
                        && (lambda - pi).abs() < self.epsilon
                    {
                        new_name = GateType::H;
                        new_params.clear();
                    } else if theta.abs() < self.epsilon && phi.abs() < self.epsilon {
                        let is_pi = (lambda - pi).abs() < self.epsilon
                            || (lambda + pi).abs() < self.epsilon;
                        let is_pi2 = (lambda - pi / 2.0).abs() < self.epsilon;
                        let is_mpi2 = (lambda + pi / 2.0).abs() < self.epsilon;
                        let is_pi4 = (lambda - pi / 4.0).abs() < self.epsilon;
                        let is_mpi4 = (lambda + pi / 4.0).abs() < self.epsilon;
                        if is_pi {
                            new_name = GateType::Z;
                            new_params.clear();
                        } else if is_pi2 {
                            new_name = GateType::S;
                            new_params.clear();
                        } else if is_mpi2 {
                            new_name = GateType::Sdg;
                            new_params.clear();
                        } else if is_pi4 {
                            new_name = GateType::T;
                            new_params.clear();
                        } else if is_mpi4 {
                            new_name = GateType::Tdg;
                            new_params.clear();
                        } else {
                            new_name = GateType::RZ;
                            new_params = vec![lambda];
                        }
                    }
                }
                out.add_op(Operation::Gate {
                    name: new_name,
                    qubits: qubits.clone(),
                    params: new_params,
                });
            } else {
                out.add_op(op.clone());
            }
        }
        out
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RotationMergePass;

impl Pass for RotationMergePass {
    fn name(&self) -> &str {
        "RotationMergePass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        run_in_segments(circuit, rotation_merge_segment)
    }
}

fn rotation_merge_segment(circuit: &Circuit) -> Circuit {
    let mut dag = DAGCircuit::from(circuit);
    let mut progress = true;
    while progress {
        progress = false;
        let edge_indices: Vec<_> = dag.graph.edge_indices().collect();
        for &edge in &edge_indices {
            if dag.graph.edge_weight(edge).is_none() {
                continue;
            }
            let Some((src, dst)) = dag.graph.edge_endpoints(edge) else {
                continue;
            };
            let src_info = match &dag.graph[src] {
                DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) if matches!(name, GateType::RX | GateType::RY | GateType::RZ)
                    && qubits.len() == 1 =>
                {
                    Some((name.clone(), qubits.clone(), params.clone()))
                }
                DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) if matches!(name, GateType::CRX | GateType::CRY | GateType::CRZ)
                    && qubits.len() == 2 =>
                {
                    Some((name.clone(), qubits.clone(), params.clone()))
                }
                _ => None,
            };
            let dst_info = match &dag.graph[dst] {
                DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) if matches!(name, GateType::RX | GateType::RY | GateType::RZ)
                    && qubits.len() == 1 =>
                {
                    Some((name.clone(), qubits.clone(), params.clone()))
                }
                DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) if matches!(name, GateType::CRX | GateType::CRY | GateType::CRZ)
                    && qubits.len() == 2 =>
                {
                    Some((name.clone(), qubits.clone(), params.clone()))
                }
                _ => None,
            };
            if let (Some((sn, sq, sp)), Some((dn, dq, dp))) = (src_info, dst_info) {
                if sq == dq && sn == dn {
                    let new_theta = sp[0] + dp[0];
                    dag.graph[src] = DAGNode::Op(Operation::Gate {
                        name: sn,
                        qubits: sq,
                        params: vec![new_theta],
                    });
                    dag.remove_node(dst);
                    progress = true;
                    break;
                }
            }
        }
    }
    Circuit::from(&dag)
}

#[derive(Debug, Clone, Copy)]
pub struct CrossConjugationPass;

impl Pass for CrossConjugationPass {
    fn name(&self) -> &str {
        "CrossConjugationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        run_in_segments(circuit, cross_conjugation_segment)
    }
}

fn cross_conjugation_segment(circuit: &Circuit) -> Circuit {
    let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    out.custom_gates = circuit.custom_gates.clone();
    let mut skip = std::collections::HashSet::new();

    for i in 0..circuit.operations.len() {
        if skip.contains(&i) {
            continue;
        }
        let op = &circuit.operations[i];
        let mut matched = false;
        if let Operation::Gate {
            name: GateType::H,
            qubits: h1,
            ..
        } = op
        {
            let q = h1[0];
            let mut j = i + 1;
            while j < circuit.operations.len() {
                let next = &circuit.operations[j];
                if involves_any(next, &[q]) {
                    if let Operation::Gate {
                        name: GateType::RZ,
                        params,
                        ..
                    } = next
                    {
                        let mut k = j + 1;
                        while k < circuit.operations.len() {
                            let fin = &circuit.operations[k];
                            if involves_any(fin, &[q]) {
                                if let Operation::Gate {
                                    name: GateType::H, ..
                                } = fin
                                {
                                    skip.insert(j);
                                    skip.insert(k);
                                    out.add_op(Operation::Gate {
                                        name: GateType::RX,
                                        qubits: vec![q],
                                        params: params.clone(),
                                    });
                                    matched = true;
                                }
                                break;
                            }
                            k += 1;
                        }
                    }
                    break;
                }
                j += 1;
            }
        }
        if !matched {
            out.add_op(op.clone());
        }
    }
    out
}

#[derive(Debug, Clone, Copy)]
pub struct InverseCancellationPass;

impl Pass for InverseCancellationPass {
    fn name(&self) -> &str {
        "InverseCancellationPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        run_in_segments(circuit, inverse_cancel_segment)
    }
}

fn inverse_cancel_segment(circuit: &Circuit) -> Circuit {
    let mut dag = DAGCircuit::from(circuit);
    let mut progress = true;
    while progress {
        progress = false;
        let edges: Vec<_> = dag.graph.edge_indices().collect();
        for &edge in &edges {
            if dag.graph.edge_weight(edge).is_none() {
                continue;
            }
            let Some((src, dst)) = dag.graph.edge_endpoints(edge) else {
                continue;
            };
            let Some(src_op) = (if let DAGNode::Op(op) = &dag.graph[src] {
                Some(op.clone())
            } else {
                None
            }) else {
                continue;
            };
            let Some(dst_op) = (if let DAGNode::Op(op) = &dag.graph[dst] {
                Some(op.clone())
            } else {
                None
            }) else {
                continue;
            };
            if are_inverses(&src_op, &dst_op) {
                let shared = get_shared_qubits(&src_op, &dst_op);
                if is_strictly_adjacent(&dag, src, dst, &shared) {
                    dag.remove_node(src);
                    dag.remove_node(dst);
                    progress = true;
                    break;
                }
            }
        }
    }
    Circuit::from(&dag)
}

fn get_shared_qubits(op1: &Operation, op2: &Operation) -> Vec<usize> {
    match (op1, op2) {
        (Operation::Gate { qubits: a, .. }, Operation::Gate { qubits: b, .. }) => {
            a.iter().filter(|q| b.contains(q)).cloned().collect()
        }
        _ => Vec::new(),
    }
}

fn is_strictly_adjacent(
    dag: &DAGCircuit,
    src: NodeIndex,
    dst: NodeIndex,
    qubits: &[usize],
) -> bool {
    for &q in qubits {
        let mut found = false;
        for edge in dag.graph.edges_directed(src, petgraph::Direction::Outgoing) {
            if edge.target() == dst && edge.weight().index == q {
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

fn are_inverses(op1: &Operation, op2: &Operation) -> bool {
    match (op1, op2) {
        (
            Operation::Gate {
                name: n1,
                qubits: q1,
                params: p1,
            },
            Operation::Gate {
                name: n2,
                qubits: q2,
                params: p2,
            },
        ) => {
            if q1 != q2 {
                return false;
            }
            match (n1, n2) {
                (GateType::H, GateType::H)
                | (GateType::X, GateType::X)
                | (GateType::Y, GateType::Y)
                | (GateType::Z, GateType::Z)
                | (GateType::CX, GateType::CX)
                | (GateType::CZ, GateType::CZ)
                | (GateType::CY, GateType::CY)
                | (GateType::CH, GateType::CH)
                | (GateType::CCX, GateType::CCX)
                | (GateType::SWAP, GateType::SWAP)
                | (GateType::S, GateType::Sdg)
                | (GateType::Sdg, GateType::S)
                | (GateType::T, GateType::Tdg)
                | (GateType::Tdg, GateType::T) => true,
                (GateType::RX, GateType::RX)
                | (GateType::RY, GateType::RY)
                | (GateType::RZ, GateType::RZ)
                | (GateType::CRX, GateType::CRX)
                | (GateType::CRY, GateType::CRY)
                | (GateType::CRZ, GateType::CRZ) => match (p1.first(), p2.first()) {
                    (Some(a), Some(b)) => {
                        (a + b).abs() < 1e-9 || (a + b - 2.0 * std::f64::consts::PI).abs() < 1e-9
                    }
                    _ => false,
                },
                _ => false,
            }
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn new_props() -> crate::transpiler::property_set::PropertySet {
        crate::transpiler::property_set::PropertySet::new()
    }

    #[test]
    fn test_gate_fusion_simple() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        c.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        let r = GateFusionPass.run(&c, &mut new_props());
        assert_eq!(r.operations.len(), 1);
    }

    #[test]
    fn test_cx_cancellation() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let r = CommutationCancellationPass.run(&c, &mut new_props());
        assert_eq!(r.operations.len(), 0);
    }

    #[test]
    fn test_barrier_stops_fusion() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        c.add_op(Operation::Barrier { qubits: vec![0] });
        c.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        let r = GateFusionPass.run(&c, &mut new_props());
        // Still 2 U gates + 1 barrier (no fusion across barrier).
        assert_eq!(r.operations.len(), 3);
    }

    #[test]
    fn test_swap_cancellation_simple() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        let r = SwapSimplificationPass.run(&c, &mut new_props());
        assert_eq!(r.operations.len(), 0);
    }

    #[test]
    fn test_swap_stops_on_touching_op() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        let r = SwapSimplificationPass.run(&c, &mut new_props());
        assert_eq!(r.operations.len(), 3);
    }

    #[test]
    fn test_parameter_simplification_drops_zero_rx() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![0.0],
        });
        let out = ParameterSimplificationPass::default().run(&c, &mut new_props());
        assert_eq!(out.operations.len(), 0);
    }
}
