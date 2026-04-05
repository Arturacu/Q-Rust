use crate::ir::{Circuit, Operation};
use petgraph::stable_graph::{NodeIndex, StableDiGraph};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WireType {
    Qubit,
    Cbit,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Wire {
    pub wire_type: WireType,
    pub index: usize,
}

#[derive(Debug, Clone)]
pub enum DAGNode {
    /// An operation node (gate, measure, reset, etc.)
    Op(Operation),
    /// Input wire
    In(Wire),
    /// Output wire
    Out(Wire),
}

impl DAGNode {
    pub fn clone_op(&self) -> Option<Operation> {
        if let DAGNode::Op(op) = self {
            Some(op.clone())
        } else {
            None
        }
    }
}

/// A Directed Acyclic Graph (DAG) representation of a quantum circuit.
///
/// This structure expresses the topological dependencies between quantum operations.
/// It enables efficient $O(1)$ adjacency checks and topological traversal,
/// which are critical for advanced optimization passes like commutation cancellation
/// and gate fusion, replacing inefficient $O(N^2)$ array lookaheads.
pub struct DAGCircuit {
    pub graph: StableDiGraph<DAGNode, Wire>,
    pub num_qubits: usize,
    pub num_cbits: usize,
    pub custom_gates: crate::ir::registry::GateRegistry,

    // Internal trackers for building the DAG
    current_q_leaves: Vec<NodeIndex>,
    current_c_leaves: Vec<NodeIndex>,
}

impl DAGCircuit {
    /// Constructs an empty DAGCircuit with the given number of qubits and classical bits.
    pub fn new(num_qubits: usize, num_cbits: usize) -> Self {
        let mut graph = StableDiGraph::new();
        let mut current_q_leaves = Vec::with_capacity(num_qubits);
        let mut current_c_leaves = Vec::with_capacity(num_cbits);

        // Add Input nodes for each qubit
        for i in 0..num_qubits {
            let idx = graph.add_node(DAGNode::In(Wire {
                wire_type: WireType::Qubit,
                index: i,
            }));
            current_q_leaves.push(idx);
        }

        // Add Input nodes for each cbit
        for i in 0..num_cbits {
            let idx = graph.add_node(DAGNode::In(Wire {
                wire_type: WireType::Cbit,
                index: i,
            }));
            current_c_leaves.push(idx);
        }

        Self {
            graph,
            num_qubits,
            num_cbits,
            custom_gates: crate::ir::registry::GateRegistry::new(),
            current_q_leaves,
            current_c_leaves,
        }
    }

    /// Appends an operation to the DAG, connecting it to the necessary data dependencies.
    pub fn add_op(&mut self, op: Operation) -> NodeIndex {
        let (q_deps, c_deps) = match &op {
            Operation::Gate { qubits, .. } => (qubits.clone(), vec![]),
            Operation::Measure { qubit, cbit } => (vec![*qubit], vec![*cbit]),
            Operation::Reset { qubit } => (vec![*qubit], vec![]),
            // Handle other ops... Assuming empty dependencies for anything else for now
            _ => (vec![], vec![]),
        };

        let node_idx = self.graph.add_node(DAGNode::Op(op));

        for q in q_deps {
            let src = self.current_q_leaves[q];
            self.graph.add_edge(
                src,
                node_idx,
                Wire {
                    wire_type: WireType::Qubit,
                    index: q,
                },
            );
            self.current_q_leaves[q] = node_idx;
        }

        for c in c_deps {
            let src = self.current_c_leaves[c];
            self.graph.add_edge(
                src,
                node_idx,
                Wire {
                    wire_type: WireType::Cbit,
                    index: c,
                },
            );
            self.current_c_leaves[c] = node_idx;
        }

        node_idx
    }

    /// Finalizes the DAG by appending Output nodes for all open wires.
    fn finalize(&mut self) {
        for i in 0..self.num_qubits {
            let out_idx = self.graph.add_node(DAGNode::Out(Wire {
                wire_type: WireType::Qubit,
                index: i,
            }));
            let src = self.current_q_leaves[i];
            self.graph.add_edge(
                src,
                out_idx,
                Wire {
                    wire_type: WireType::Qubit,
                    index: i,
                },
            );
            self.current_q_leaves[i] = out_idx;
        }

        for i in 0..self.num_cbits {
            let out_idx = self.graph.add_node(DAGNode::Out(Wire {
                wire_type: WireType::Cbit,
                index: i,
            }));
            let src = self.current_c_leaves[i];
            self.graph.add_edge(
                src,
                out_idx,
                Wire {
                    wire_type: WireType::Cbit,
                    index: i,
                },
            );
            self.current_c_leaves[i] = out_idx;
        }
    }

    /// Removes a node from the DAG, rewiring incoming edges directly to outgoing edges.
    pub fn remove_node(&mut self, node_idx: NodeIndex) {
        // Collect incoming and outgoing edges for each wire
        let mut incoming = Vec::new();
        let mut outgoing = Vec::new();

        let mut walker = self
            .graph
            .neighbors_directed(node_idx, petgraph::Direction::Incoming)
            .detach();
        while let Some((edge_idx, src_idx)) = walker.next(&self.graph) {
            let wire = self.graph[edge_idx].clone();
            incoming.push((src_idx, wire));
        }

        let mut walker_out = self
            .graph
            .neighbors_directed(node_idx, petgraph::Direction::Outgoing)
            .detach();
        while let Some((edge_idx, dst_idx)) = walker_out.next(&self.graph) {
            let wire = self.graph[edge_idx].clone();
            outgoing.push((dst_idx, wire));
        }

        // Delete the node (this removes all connected edges too)
        self.graph.remove_node(node_idx);

        // Rewire src -> dst for matching wires
        for (src, in_wire) in &incoming {
            for (dst, out_wire) in &outgoing {
                if in_wire == out_wire {
                    self.graph.add_edge(*src, *dst, in_wire.clone());
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Converters (Circuit <-> DAGCircuit)
// ---------------------------------------------------------------------------

impl From<&Circuit> for DAGCircuit {
    fn from(circuit: &Circuit) -> Self {
        let mut dag = DAGCircuit::new(circuit.num_qubits, circuit.num_cbits);
        dag.custom_gates = circuit.custom_gates.clone();

        for op in &circuit.operations {
            dag.add_op(op.clone());
        }

        dag.finalize();
        dag
    }
}

impl From<&DAGCircuit> for Circuit {
    fn from(dag: &DAGCircuit) -> Self {
        use petgraph::algo::toposort;

        let mut circuit = Circuit::new(dag.num_qubits, dag.num_cbits);
        circuit.custom_gates = dag.custom_gates.clone();

        let sorted = toposort(&dag.graph, None).expect("DAG contains a cycle!");

        for node_idx in sorted {
            if let DAGNode::Op(op) = &dag.graph[node_idx] {
                circuit.add_op(op.clone());
            }
        }

        circuit
    }
}
