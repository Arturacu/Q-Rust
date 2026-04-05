use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::synthesis::zyz::{u_to_matrix, zyz_decomposition};

/// A pass that merges consecutive single-qubit gates on the same qubit.
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
        let mut dag = DAGCircuit::from(circuit);
        let mut progress = true;

        while progress {
            progress = false;
            let edge_indices: Vec<_> = dag.graph.edge_indices().collect();

            for &edge in &edge_indices {
                if dag.graph.edge_weight(edge).is_none() {
                    continue; // Edge was removed during iteration
                }

                let (src, dst) = match dag.graph.edge_endpoints(edge) {
                    Some(endpoints) => endpoints,
                    None => continue,
                };

                // Check if both nodes are single-qubit U gates
                let mut src_op_u = None;
                if let DAGNode::Op(Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                }) = &dag.graph[src]
                {
                    if qubits.len() == 1 {
                        src_op_u = Some((qubits[0], params.clone()));
                    }
                }

                let mut dst_op_u = None;
                if let DAGNode::Op(Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                }) = &dag.graph[dst]
                {
                    if qubits.len() == 1 {
                        dst_op_u = Some((qubits[0], params.clone()));
                    }
                }

                if let (Some((sq, sp)), Some((dq, dp))) = (src_op_u, dst_op_u) {
                    if sq == dq {
                        // We found an edge connecting two U gates on the same qubit: src -> dst.
                        // New matrix = U(dst) * U(src)
                        let m_matrix = u_to_matrix(sp[0], sp[1], sp[2]);
                        let n_matrix = u_to_matrix(dp[0], dp[1], dp[2]);

                        let c00 = n_matrix[0][0] * m_matrix[0][0] + n_matrix[0][1] * m_matrix[1][0];
                        let c01 = n_matrix[0][0] * m_matrix[0][1] + n_matrix[0][1] * m_matrix[1][1];
                        let c10 = n_matrix[1][0] * m_matrix[0][0] + n_matrix[1][1] * m_matrix[1][0];
                        let c11 = n_matrix[1][0] * m_matrix[0][1] + n_matrix[1][1] * m_matrix[1][1];

                        let c_matrix = [[c00, c01], [c10, c11]];
                        let (theta, phi, lambda, _) = zyz_decomposition(c_matrix);

                        // Update src node to the fused gate
                        dag.graph[src] = DAGNode::Op(Operation::Gate {
                            name: GateType::U,
                            qubits: vec![sq],
                            params: vec![theta, phi, lambda],
                        });

                        // Delete dst node. `remove_node` automatically rewires incoming edges
                        // directly to outgoing, meaning the DAG topology remains perfectly sound!
                        dag.remove_node(dst);

                        progress = true;
                        break; // Restart loop to handle updated DAG
                    }
                }
            }
        }

        Circuit::from(&dag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_gate_fusion_simple() {
        let mut circuit = Circuit::new(1, 0);
        // U(pi/2, 0, pi) -> H
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        // U(pi/2, 0, pi) -> H
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });

        let pass = GateFusionPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        // H * H = I
        // Should result in U(0, 0, 0) or similar identity-equivalent
        assert_eq!(new_circuit.operations.len(), 1);

        if let Operation::Gate { name, params, .. } = &new_circuit.operations[0] {
            if let GateType::U = name {
                assert!(params[0].abs() < 1e-6); // Identity has theta=0
            } else {
                panic!("Expected U gate");
            }
        }
    }

    #[test]
    fn test_gate_fusion_blocking() {
        let mut circuit = Circuit::new(2, 0);
        // H q0
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        // CX q0, q1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H q0
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });

        let pass = GateFusionPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        // Should NOT merge H gates because CX blocks them
        assert_eq!(new_circuit.operations.len(), 3);
    }
}

use crate::transpiler::dag::{DAGCircuit, DAGNode};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

/// A pass that cancels CX gates separated by commuting single-qubit gates,
/// or adjacent identical CX gates.
/// Walks the DAGCircuit topologically to find canceling pairs.
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
        let mut dag = DAGCircuit::from(circuit);
        let mut progress = true;

        while progress {
            progress = false;
            let node_indices: Vec<_> = dag.graph.node_indices().collect();

            for &idx in &node_indices {
                if dag.graph.node_weight(idx).is_none() {
                    continue; // Node might have been removed
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
                        // Found a matching CX with perfectly commuting paths!
                        dag.remove_node(idx);
                        dag.remove_node(cancel_idx);
                        progress = true;
                        break; // Restart the pass with the simplified DAG
                    }
                }
            }
        }

        Circuit::from(&dag)
    }
}

fn find_commuting_cx(
    dag: &DAGCircuit,
    start: NodeIndex,
    ctrl_wire: usize,
    target_wire: usize,
) -> Option<NodeIndex> {
    // Traverse ctrl_wire forward
    let mut current = start;

    let candidate_idx = loop {
        let mut next_node = None;
        let edges = dag
            .graph
            .edges_directed(current, petgraph::Direction::Outgoing);
        for edge in edges {
            if edge.weight().index == ctrl_wire {
                next_node = Some(edge.target());
                break;
            }
        }

        let next = next_node?;

        match &dag.graph[next] {
            DAGNode::Op(Operation::Gate {
                name: GateType::CX,
                qubits,
                ..
            }) => {
                if qubits[0] == ctrl_wire && qubits[1] == target_wire {
                    break next;
                }

                // Another gate, check commutation
                if !commutes_with_cx(&dag.graph[next].clone_op().unwrap(), ctrl_wire, target_wire) {
                    return None;
                }
            }
            DAGNode::Op(op) => {
                if !commutes_with_cx(op, ctrl_wire, target_wire) {
                    return None;
                }
            }
            _ => return None, // In/Out block
        }
        current = next;
    };

    // Traverse target_wire forward to verify it reaches candidate cleanly
    current = start;
    loop {
        let mut next_node = None;
        let edges = dag
            .graph
            .edges_directed(current, petgraph::Direction::Outgoing);
        for edge in edges {
            if edge.weight().index == target_wire {
                next_node = Some(edge.target());
                break;
            }
        }

        let next = next_node?;

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

fn commutes_with_cx(op: &Operation, ctrl: usize, target: usize) -> bool {
    match op {
        Operation::Gate {
            name,
            qubits,
            params,
        } => {
            // If the gate doesn't involve the CX qubits, it commutes
            let involves_ctrl = qubits.contains(&ctrl);
            let involves_target = qubits.contains(&target);

            if !involves_ctrl && !involves_target {
                return true;
            }

            // If it involves other qubits (multi-qubit gate), assume no commutation for now
            // unless we implement more complex logic.
            if qubits.len() > 1 {
                return false;
            }

            // Single qubit gate logic
            let q = qubits[0];

            if q == ctrl {
                // Gate on control qubit: Commutes if it is Z-basis (Z, RZ, U1, S, T)
                // Anti-commutes if X-basis (X, RX)
                match name {
                    GateType::Z
                    | GateType::S
                    | GateType::Sdg
                    | GateType::T
                    | GateType::Tdg
                    | GateType::RZ => true,
                    GateType::U => {
                        // U(0, 0, lambda) is RZ(lambda) -> Commutes
                        // U(theta, phi, lambda) generally doesn't unless theta=0
                        if params.len() >= 1 && params[0].abs() < 1e-9 {
                            true
                        } else {
                            false
                        }
                    }
                    _ => false, // X, Y, H, etc.
                }
            } else if q == target {
                // Gate on target qubit: Commutes if it is X-basis (X, RX)
                // Anti-commutes if Z-basis (Z, RZ)
                match name {
                    GateType::X | GateType::RX => true,
                    GateType::U => {
                        // U(pi, -pi/2, pi/2) is X -> Commutes
                        // U(theta, -pi/2, pi/2) is RX(theta) -> Commutes
                        if params.len() >= 3 {
                            let _theta = params[0];
                            let phi = params[1];
                            let lambda = params[2];
                            // Check if it's RX-like: phi = -pi/2, lambda = pi/2
                            let pi_2 = std::f64::consts::FRAC_PI_2;
                            if (phi - -pi_2).abs() < 1e-9 && (lambda - pi_2).abs() < 1e-9 {
                                true
                            } else {
                                false
                            }
                        } else {
                            false
                        }
                    }
                    _ => false,
                }
            } else {
                true // Should be covered by first check, but safe fallback
            }
        }
        _ => false, // Measure/Barrier/etc assume no commutation
    }
}

#[cfg(test)]
mod optimization_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_cx_cancellation() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_commutation_cancellation_rz_control() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // RZ(pi) 0 (Commutes with Control)
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![PI],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        // Should remove CXs, leaving only RZ
        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::RZ);
        } else {
            panic!("Expected RZ");
        }
    }

    #[test]
    fn test_commutation_cancellation_rx_target() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // RX(pi) 1 (Commutes with Target X-component)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![1],
            params: vec![PI],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::RX);
        } else {
            panic!("Expected RX");
        }
    }

    #[test]
    fn test_commutation_blocking() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // X 0 (Anti-commutes with Control Z-component)
        circuit.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        // Should NOT cancel
        assert_eq!(new_circuit.operations.len(), 3);
    }
}

/// A pass that removes redundant SWAP gates (SWAP a,b; SWAP a,b -> ID).
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
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let mut i = 0;

        while i < circuit.operations.len() {
            let op = &circuit.operations[i];

            // Check start of potential cancellation block: SWAP a,b
            if let Operation::Gate {
                name: GateType::SWAP,
                qubits: swap_qubits,
                ..
            } = op
            {
                // Look ahead for matching SWAP
                let mut j = i + 1;
                let mut commute = true;
                let mut intermediate_ops = Vec::new();

                while j < circuit.operations.len() {
                    let next_op = &circuit.operations[j];

                    // If we find the matching SWAP, we can stop
                    if let Operation::Gate {
                        name: GateType::SWAP,
                        qubits: next_qubits,
                        ..
                    } = next_op
                    {
                        // Check if qubits match (order doesn't matter for SWAP)
                        if (next_qubits[0] == swap_qubits[0] && next_qubits[1] == swap_qubits[1])
                            || (next_qubits[0] == swap_qubits[1]
                                && next_qubits[1] == swap_qubits[0])
                        {
                            // Found match!
                            break;
                        }
                    }

                    // Check commutation: SWAP(a,b) commutes with any gate that doesn't involve a or b
                    if involves_any(next_op, &swap_qubits) {
                        commute = false;
                        break;
                    }

                    intermediate_ops.push(next_op.clone());
                    j += 1;
                }

                if commute && j < circuit.operations.len() {
                    // We found a matching SWAP and everything in between commutes.
                    // Skip the first SWAP (i), add intermediate ops, and skip the second SWAP (j).
                    for mid_op in intermediate_ops {
                        new_circuit.add_op(mid_op);
                    }
                    i = j + 1; // Continue after the second SWAP
                    continue;
                }
            }

            new_circuit.add_op(op.clone());
            i += 1;
        }

        new_circuit
    }
}

fn involves_any(op: &Operation, qubits: &[usize]) -> bool {
    match op {
        Operation::Gate {
            qubits: op_qubits, ..
        } => {
            for q in op_qubits {
                if qubits.contains(q) {
                    return true;
                }
            }
            false
        }
        Operation::Measure { qubit, .. } => qubits.contains(qubit),
        Operation::Reset { qubit } => qubits.contains(qubit),
        Operation::Barrier {
            qubits: barrier_qubits,
        } => {
            for q in barrier_qubits {
                if qubits.contains(q) {
                    return true;
                }
            }
            false
        }
    }
}

#[cfg(test)]
mod swap_tests {
    use super::*;

    #[test]
    fn test_swap_cancellation_simple() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_swap_cancellation_reversed() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![1, 0],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_swap_cancellation_commuting() {
        let mut circuit = Circuit::new(3, 0);
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H 2 (Commutes)
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![2],
            params: vec![],
        });
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::H);
        } else {
            panic!("Expected H");
        }
    }

    #[test]
    fn test_swap_blocking() {
        let mut circuit = Circuit::new(3, 0);
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H 0 (Blocks)
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circuit.operations.len(), 3);
    }
}

/// A pass that simplifies gate parameters (e.g., modulo 2pi, removing negligible rotations).
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
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let two_pi = 2.0 * std::f64::consts::PI;

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } => {
                    let mut new_params = params.clone();
                    let mut keep_gate = true;

                    // Normalize parameters and check for identity
                    match name {
                        GateType::RX | GateType::RY | GateType::RZ | GateType::U => {
                            for param in &mut new_params {
                                // Normalize to [0, 2pi)
                                *param = *param % two_pi;
                                if *param < 0.0 {
                                    *param += two_pi;
                                }
                            }

                            // Check if gate is effectively identity
                            match name {
                                GateType::RX | GateType::RY | GateType::RZ => {
                                    if let Some(theta) = new_params.first() {
                                        if theta.abs() < self.epsilon
                                            || (theta - two_pi).abs() < self.epsilon
                                        {
                                            keep_gate = false;
                                        }
                                    }
                                }
                                GateType::U => {
                                    // U(0, 0, 0) is Identity
                                    // U(0, 0, lambda) is RZ(lambda) - keep unless lambda ~ 0
                                    // U(theta, phi, lambda)
                                    if new_params.len() >= 3 {
                                        let theta = new_params[0];
                                        let phi = new_params[1];
                                        let lambda = new_params[2];

                                        // If theta ~ 0, it's diagonal.
                                        // If phi + lambda ~ 0 (mod 2pi), it's global phase only?
                                        // U(0, phi, lambda) = diag(1, e^(i(phi+lambda)))
                                        // If phi+lambda = 0, it's I.

                                        if theta.abs() < self.epsilon {
                                            let phase_sum = (phi + lambda) % two_pi;
                                            if phase_sum.abs() < self.epsilon
                                                || (phase_sum - two_pi).abs() < self.epsilon
                                            {
                                                keep_gate = false;
                                            }
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                        _ => {}
                    }

                    if keep_gate {
                        new_circuit.add_op(Operation::Gate {
                            name: name.clone(),
                            qubits: qubits.clone(),
                            params: new_params,
                        });
                    }
                }
                _ => new_circuit.add_op(op.clone()),
            }
        }

        new_circuit
    }
}

/// A pass that coerces fuzzy float rotations exactly mapping to multiples of PI / 2
/// downward into algebraic physical constants (X, Y, Z, S, T, H), directly exposing
/// their algebraic CommutationSignature to the topology framework.
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
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        new_circuit.custom_gates = circuit.custom_gates.clone();
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
                    let theta = params[0];
                    let is_pi =
                        (theta - pi).abs() < self.epsilon || (theta + pi).abs() < self.epsilon;
                    let is_pi_2 = (theta - pi / 2.0).abs() < self.epsilon;
                    let is_minus_pi_2 = (theta + pi / 2.0).abs() < self.epsilon;
                    let is_pi_4 = (theta - pi / 4.0).abs() < self.epsilon;
                    let is_minus_pi_4 = (theta + pi / 4.0).abs() < self.epsilon;

                    match name {
                        GateType::RX => {
                            if is_pi {
                                new_name = GateType::X;
                                new_params.clear();
                            }
                        }
                        GateType::RY => {
                            if is_pi {
                                new_name = GateType::Y;
                                new_params.clear();
                            }
                        }
                        GateType::RZ => {
                            if is_pi {
                                new_name = GateType::Z;
                                new_params.clear();
                            } else if is_pi_2 {
                                new_name = GateType::S;
                                new_params.clear();
                            } else if is_minus_pi_2 {
                                new_name = GateType::Sdg;
                                new_params.clear();
                            } else if is_pi_4 {
                                new_name = GateType::T;
                                new_params.clear();
                            } else if is_minus_pi_4 {
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

                    // U(pi/2, 0, pi) is H
                    if (theta - pi / 2.0).abs() < self.epsilon
                        && phi.abs() < self.epsilon
                        && (lambda - pi).abs() < self.epsilon
                    {
                        new_name = GateType::H;
                        new_params.clear();
                    }
                    // U(0, 0, lambda) is RZ(lambda)
                    else if theta.abs() < self.epsilon && phi.abs() < self.epsilon {
                        let is_pi = (lambda - pi).abs() < self.epsilon
                            || (lambda + pi).abs() < self.epsilon;
                        let is_pi_2 = (lambda - pi / 2.0).abs() < self.epsilon;
                        let is_minus_pi_2 = (lambda + pi / 2.0).abs() < self.epsilon;
                        let is_pi_4 = (lambda - pi / 4.0).abs() < self.epsilon;
                        let is_minus_pi_4 = (lambda + pi / 4.0).abs() < self.epsilon;

                        if is_pi {
                            new_name = GateType::Z;
                            new_params.clear();
                        } else if is_pi_2 {
                            new_name = GateType::S;
                            new_params.clear();
                        } else if is_minus_pi_2 {
                            new_name = GateType::Sdg;
                            new_params.clear();
                        } else if is_pi_4 {
                            new_name = GateType::T;
                            new_params.clear();
                        } else if is_minus_pi_4 {
                            new_name = GateType::Tdg;
                            new_params.clear();
                        } else {
                            new_name = GateType::RZ;
                            new_params = vec![lambda];
                        }
                    }
                }

                new_circuit.add_op(Operation::Gate {
                    name: new_name,
                    qubits: qubits.clone(),
                    params: new_params,
                });
            } else {
                new_circuit.add_op(op.clone());
            }
        }

        new_circuit
    }
}

/// A pass that targets sum fusion of topological identical continuous rotation axes sequentially.
/// Preserves explicit structural Lie-algebras rather than collapsing into numerical matrices.
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
        let mut dag = DAGCircuit::from(circuit);
        let mut progress = true;

        while progress {
            progress = false;
            let edge_indices: Vec<_> = dag.graph.edge_indices().collect();

            for &edge in &edge_indices {
                if dag.graph.edge_weight(edge).is_none() {
                    continue; // Edge was removed
                }

                let (src, dst) = match dag.graph.edge_endpoints(edge) {
                    Some(endpoints) => endpoints,
                    None => continue,
                };

                let mut src_info = None;
                if let DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) = &dag.graph[src]
                {
                    if qubits.len() == 1 {
                        if matches!(name, GateType::RX | GateType::RY | GateType::RZ) {
                            src_info = Some((name.clone(), qubits[0], params.clone()));
                        }
                    }
                }

                let mut dst_info = None;
                if let DAGNode::Op(Operation::Gate {
                    name,
                    qubits,
                    params,
                }) = &dag.graph[dst]
                {
                    if qubits.len() == 1 {
                        if matches!(name, GateType::RX | GateType::RY | GateType::RZ) {
                            dst_info = Some((name.clone(), qubits[0], params.clone()));
                        }
                    }
                }

                if let (Some((sn, sq, sp)), Some((dn, dq, dp))) = (src_info, dst_info) {
                    // Check if they are on the exact same axis and the exact same qubit
                    if sq == dq && sn == dn {
                        // Geometrically sum the underlying continuous rotational parameters natively
                        let new_theta = sp[0] + dp[0];

                        dag.graph[src] = DAGNode::Op(Operation::Gate {
                            name: sn,
                            qubits: vec![sq],
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
}

/// A pass that identifies and simplifies cross-conjugated rotation operators natively mapping
/// specific bases without triggering numeric intermediate matrix factorization.
/// Implements structural rewrite rules such as `H * RZ(theta) * H -> RX(theta)`.
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
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        new_circuit.custom_gates = circuit.custom_gates.clone();
        let mut skip = std::collections::HashSet::new();

        for i in 0..circuit.operations.len() {
            if skip.contains(&i) {
                continue;
            }

            let op = &circuit.operations[i];
            let mut matched = false;

            // Algebraic pattern: H(q) -> RZ(t) -> H(q) analytically maps to RX(t)
            if let Operation::Gate {
                name: GateType::H,
                qubits: h1_qubits,
                ..
            } = op
            {
                let q = h1_qubits[0];

                // Scan forward to detect the next dependency boundary on wire `q`
                let mut j = i + 1;
                while j < circuit.operations.len() {
                    let next_op = &circuit.operations[j];
                    if involves_any(next_op, &[q]) {
                        if let Operation::Gate {
                            name: GateType::RZ,
                            params,
                            ..
                        } = next_op
                        {
                            // It's RZ! Now look forward again to close the conjugation target
                            let mut k = j + 1;
                            while k < circuit.operations.len() {
                                let final_op = &circuit.operations[k];
                                if involves_any(final_op, &[q]) {
                                    if let Operation::Gate {
                                        name: GateType::H, ..
                                    } = final_op
                                    {
                                        // Matched: H * RZ * H == RX!
                                        skip.insert(j);
                                        skip.insert(k);
                                        new_circuit.add_op(Operation::Gate {
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
                new_circuit.add_op(op.clone());
            }
        }

        new_circuit
    }
}

/// A pass that cancels adjacent inverse operators (e.g. H-H, X-X, S-Sdg) topologically
/// by evaluating direct edge connections in the unified dependency graph.
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
        let mut dag = DAGCircuit::from(circuit);
        let mut progress = true;

        while progress {
            progress = false;
            let edge_indices: Vec<_> = dag.graph.edge_indices().collect();

            for &edge in &edge_indices {
                if dag.graph.edge_weight(edge).is_none() {
                    continue; // Edge was destructively removed iteratively
                }

                let (src, dst) = match dag.graph.edge_endpoints(edge) {
                    Some(endpoints) => endpoints,
                    None => continue,
                };

                let src_op = if let DAGNode::Op(op) = &dag.graph[src] {
                    op.clone()
                } else {
                    continue;
                };
                let dst_op = if let DAGNode::Op(op) = &dag.graph[dst] {
                    op.clone()
                } else {
                    continue;
                };

                // Since we rely on a strictly routed directed topological edge, the adjacency guarantee
                // mathematically asserts that no interfering structural gates sit between them.
                if are_inverses(&src_op, &dst_op) {
                    dag.remove_node(src);
                    dag.remove_node(dst);
                    progress = true;
                    break;
                }
            }
        }

        Circuit::from(&dag)
    }
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
                // Self-inverses
                (GateType::H, GateType::H) => true,
                (GateType::X, GateType::X) => true,
                (GateType::Y, GateType::Y) => true,
                (GateType::Z, GateType::Z) => true,
                (GateType::CX, GateType::CX) => true, // Already handled by CXCancellation, but good to have
                (GateType::SWAP, GateType::SWAP) => true, // Already handled

                // Inverse pairs
                (GateType::S, GateType::Sdg) => true,
                (GateType::Sdg, GateType::S) => true,
                (GateType::T, GateType::Tdg) => true,
                (GateType::Tdg, GateType::T) => true,

                // Rotations with opposite angles
                (GateType::RX, GateType::RX)
                | (GateType::RY, GateType::RY)
                | (GateType::RZ, GateType::RZ) => {
                    if let (Some(a1), Some(a2)) = (p1.first(), p2.first()) {
                        (a1 + a2).abs() < 1e-9
                            || (a1 + a2 - 2.0 * std::f64::consts::PI).abs() < 1e-9
                    } else {
                        false
                    }
                }

                _ => false,
            }
        }
        _ => false,
    }
}

#[cfg(test)]
mod param_peephole_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_param_simplification_identity() {
        let mut circuit = Circuit::new(1, 0);
        // RZ(0) -> Identity
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![0.0],
        });
        // RZ(2pi) -> Identity
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![2.0 * PI],
        });
        // RX(1e-10) -> Identity (negligible)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![1e-10],
        });

        let pass = ParameterSimplificationPass::default();
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_param_simplification_normalization() {
        let mut circuit = Circuit::new(1, 0);
        let pi = std::f64::consts::PI;
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![3.0 * pi], // Should become pi
        });

        let pass = ParameterSimplificationPass::default();
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { params, .. } = &new_circuit.operations[0] {
            assert!((params[0] - pi).abs() < 1e-9);
        }
    }
}

#[cfg(test)]
mod crystallization_and_merge_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_gate_crystallization_coerces_exact_gates() {
        let mut circuit = Circuit::new(1, 0);

        // RX(pi) -> X
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![PI],
        });

        // RZ(pi/2) -> S
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![PI / 2.0],
        });

        let pass = GateCrystallizationPass::default();
        let new_circ = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circ.operations.len(), 2);
        if let Operation::Gate { name, .. } = &new_circ.operations[0] {
            assert_eq!(*name, GateType::X);
        }
        if let Operation::Gate { name, .. } = &new_circ.operations[1] {
            assert_eq!(*name, GateType::S);
        }
    }

    #[test]
    fn test_rotation_merge_combines_contiguous_parameters() {
        let mut circuit = Circuit::new(1, 0);

        // RX(0.3)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![0.3],
        });
        // RX(0.4)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![0.4],
        });

        // These should mathematically fuse to RX(0.7)
        let pass = RotationMergePass;
        let new_circ = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circ.operations.len(), 1);
        if let Operation::Gate { name, params, .. } = &new_circ.operations[0] {
            assert_eq!(*name, GateType::RX);
            assert!((params[0] - 0.7).abs() < 1e-9);
        }
    }
}

#[cfg(test)]
mod peephole_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_peephole_inverses() {
        let mut circuit = Circuit::new(1, 0);
        // H H -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        // S Sdg -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::S,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::Sdg,
            qubits: vec![0],
            params: vec![],
        });

        let pass = InverseCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_peephole_rotation_cancellation() {
        let mut circuit = Circuit::new(1, 0);
        // RX(pi) RX(-pi) -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![PI],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![-PI],
        });

        let pass = InverseCancellationPass;
        let new_circuit = pass.run(
            &circuit,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        assert_eq!(new_circuit.operations.len(), 0);
    }
}
