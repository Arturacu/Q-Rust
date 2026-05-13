//! Directed-acyclic-graph IR of a quantum circuit.
//!
//! Each node is either an operation ([`DAGNode::Op`]) or an I/O terminal.
//! Edges carry the [`Wire`] they route (a qubit or a classical bit), which
//! makes commutation and local-rewrite analyses O(1) on adjacency.

use crate::ir::{Circuit, Operation};
use petgraph::stable_graph::{NodeIndex, StableDiGraph};

/// Whether a wire carries quantum or classical data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WireType {
    /// A qubit wire.
    Qubit,
    /// A classical bit wire.
    Cbit,
}

/// An edge label identifying which wire this edge routes.
///
/// Each DAG edge corresponds to a single wire (qubit or classical bit) on
/// which a value flows from a producer node to a consumer node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Wire {
    /// Whether the wire is quantum or classical.
    pub wire_type: WireType,
    /// Index into the corresponding register.
    pub index: usize,
}

/// A node in the circuit DAG.
#[derive(Debug, Clone)]
pub enum DAGNode {
    /// An operation node (gate, measure, reset, etc.).
    Op(Operation),
    /// Input terminal for a wire.
    In(Wire),
    /// Output terminal for a wire.
    Out(Wire),
}

impl DAGNode {
    /// Returns a clone of the inner operation, if this node is an op.
    pub fn clone_op(&self) -> Option<Operation> {
        if let DAGNode::Op(op) = self {
            Some(op.clone())
        } else {
            None
        }
    }
}

/// DAG-form quantum circuit.
///
/// Built from a [`Circuit`] via `DAGCircuit::from(&circuit)`, mutated by
/// optimization passes, and converted back via `Circuit::from(&dag)`.
pub struct DAGCircuit {
    /// Underlying directed graph. Stable indices survive node removal.
    pub graph: StableDiGraph<DAGNode, Wire>,
    /// Number of qubit wires.
    pub num_qubits: usize,
    /// Number of classical-bit wires.
    pub num_cbits: usize,
    /// Custom gate definitions carried over from the source circuit.
    pub custom_gates: crate::ir::registry::GateRegistry,

    current_q_leaves: Vec<NodeIndex>,
    current_c_leaves: Vec<NodeIndex>,
}

impl std::fmt::Debug for DAGCircuit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DAGCircuit")
            .field("num_qubits", &self.num_qubits)
            .field("num_cbits", &self.num_cbits)
            .field("node_count", &self.graph.node_count())
            .field("edge_count", &self.graph.edge_count())
            .finish()
    }
}

impl DAGCircuit {
    /// Creates an empty DAG with `In` terminals for each wire.
    pub fn new(num_qubits: usize, num_cbits: usize) -> Self {
        let mut graph = StableDiGraph::new();
        let mut current_q_leaves = Vec::with_capacity(num_qubits);
        let mut current_c_leaves = Vec::with_capacity(num_cbits);

        for i in 0..num_qubits {
            let idx = graph.add_node(DAGNode::In(Wire {
                wire_type: WireType::Qubit,
                index: i,
            }));
            current_q_leaves.push(idx);
        }
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

    /// Appends an operation, wiring each of its qubit/cbit dependencies
    /// to the most recent producer on that wire.
    pub fn add_op(&mut self, op: Operation) -> NodeIndex {
        let (q_deps, c_deps): (Vec<usize>, Vec<usize>) = match &op {
            Operation::Gate { qubits, .. } => (qubits.clone(), vec![]),
            Operation::Measure { qubit, cbit } => (vec![*qubit], vec![*cbit]),
            Operation::Reset { qubit } => (vec![*qubit], vec![]),
            Operation::Barrier { qubits } => (qubits.clone(), vec![]),
            Operation::Conditional { op: inner, .. } => {
                let q = inner.qubits().to_vec();
                let c = match &**inner {
                    Operation::Measure { cbit, .. } => vec![*cbit],
                    _ => vec![],
                };
                (q, c)
            }
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

    /// Appends `Out` terminals to every open wire — should be called once
    /// all operations have been inserted.
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

    /// Removes a node, reconnecting each incoming edge directly to the
    /// outgoing edge on the same wire (so the DAG remains well-formed).
    pub fn remove_node(&mut self, node_idx: NodeIndex) {
        let mut incoming = Vec::new();
        let mut outgoing = Vec::new();

        let mut walker = self
            .graph
            .neighbors_directed(node_idx, petgraph::Direction::Incoming)
            .detach();
        while let Some((edge_idx, src_idx)) = walker.next(&self.graph) {
            incoming.push((src_idx, self.graph[edge_idx].clone()));
        }
        let mut walker_out = self
            .graph
            .neighbors_directed(node_idx, petgraph::Direction::Outgoing)
            .detach();
        while let Some((edge_idx, dst_idx)) = walker_out.next(&self.graph) {
            outgoing.push((dst_idx, self.graph[edge_idx].clone()));
        }

        self.graph.remove_node(node_idx);

        for (src, in_wire) in &incoming {
            for (dst, out_wire) in &outgoing {
                if in_wire == out_wire {
                    self.graph.add_edge(*src, *dst, in_wire.clone());
                }
            }
        }
    }
}

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