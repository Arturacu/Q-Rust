use petgraph::Directed;
use petgraph::graph::{Graph, NodeIndex};
use std::collections::HashSet;

/// Backend Specification
#[derive(Debug, Clone)]
pub struct Backend {
    pub name: String,
    pub num_qubits: usize,
    pub basis_gates: HashSet<String>,
    /// Coupling map representing physical qubit connectivity.
    /// Nodes are physical qubits, edges represent allowed 2-qubit gates.
    pub coupling_map: Graph<(), (), Directed>,
}

impl Backend {
    pub fn new(name: String, num_qubits: usize) -> Self {
        let mut graph = Graph::new();
        // Initialize nodes for each qubit
        for _ in 0..num_qubits {
            graph.add_node(());
        }

        Self {
            name,
            num_qubits,
            basis_gates: HashSet::new(),
            coupling_map: graph,
        }
    }

    pub fn add_basis_gate(&mut self, gate: &str) {
        self.basis_gates.insert(gate.to_string());
    }

    /// Sets the coupling map from a list of edges.
    /// Edges are directed: (source, target).
    pub fn set_coupling_map(&mut self, edges: Vec<(usize, usize)>) {
        self.coupling_map.clear_edges();
        for (u, v) in edges {
            // Ensure indices are within bounds
            if u < self.num_qubits && v < self.num_qubits {
                self.coupling_map
                    .add_edge(NodeIndex::new(u), NodeIndex::new(v), ());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backend_creation() {
        let backend = Backend::new("test_backend".to_string(), 5);
        assert_eq!(backend.num_qubits, 5);
        assert_eq!(backend.coupling_map.node_count(), 5);
        assert_eq!(backend.coupling_map.edge_count(), 0);
    }

    #[test]
    fn test_coupling_map() {
        let mut backend = Backend::new("test_backend".to_string(), 3);
        // Linear connectivity: 0 -> 1 -> 2
        let edges = vec![(0, 1), (1, 2)];
        backend.set_coupling_map(edges);

        assert_eq!(backend.coupling_map.edge_count(), 2);
        assert!(
            backend
                .coupling_map
                .contains_edge(NodeIndex::new(0), NodeIndex::new(1))
        );
        assert!(
            backend
                .coupling_map
                .contains_edge(NodeIndex::new(1), NodeIndex::new(2))
        );
        assert!(
            !backend
                .coupling_map
                .contains_edge(NodeIndex::new(0), NodeIndex::new(2))
        );
    }
}
