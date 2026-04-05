use petgraph::graph::{Graph, NodeIndex};
use petgraph::Directed;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Serialization surrogate for reading standard JSON hardware definitions.
/// Designed to perfectly interface with generic cloud API execution bounds natively.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendConfig {
    pub backend_name: String,
    pub n_qubits: usize,
    pub basis_gates: Vec<String>,
    pub coupling_map: Vec<[usize; 2]>,
}

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
            if u < self.num_qubits && v < self.num_qubits {
                self.coupling_map
                    .add_edge(NodeIndex::new(u), NodeIndex::new(v), ());
            }
        }
    }

    /// Evaluates a standalone JSON payload schema directly into physical graph logic
    pub fn from_config(config: BackendConfig) -> Self {
        let mut backend = Backend::new(config.backend_name, config.n_qubits);
        for gate in config.basis_gates {
            backend.add_basis_gate(&gate);
        }
        let edges = config
            .coupling_map
            .into_iter()
            .map(|arr| (arr[0], arr[1]))
            .collect();
        backend.set_coupling_map(edges);
        backend
    }

    /// Analytical generator for an unconstrained (infinite-routing) Complete Graph.
    /// Safely collapses routing complexity bounds for pure-simulation objectives natively.
    pub fn all_to_all(num_qubits: usize) -> Self {
        let mut backend = Backend::new("all_to_all".to_string(), num_qubits);
        let mut edges = Vec::new();
        for i in 0..num_qubits {
            for j in 0..num_qubits {
                if i != j {
                    edges.push((i, j));
                }
            }
        }
        backend.set_coupling_map(edges);
        backend
    }

    /// Standard generator for localized 1D-chain architectures (e.g., trapped-ion arrays).
    pub fn linear(num_qubits: usize) -> Self {
        let mut backend = Backend::new("linear".to_string(), num_qubits);
        let mut edges = Vec::new();
        for i in 0..num_qubits.saturating_sub(1) {
            edges.push((i, i + 1));
            edges.push((i + 1, i)); // Assuming symmetric bidirectional connectivity
        }
        backend.set_coupling_map(edges);
        backend
    }

    /// Mathematical generator mapping nearest-neighbor interactions across 2D-lattice squares.
    pub fn grid(rows: usize, cols: usize) -> Self {
        let num_qubits = rows * cols;
        let mut backend = Backend::new(format!("grid_{}x{}", rows, cols), num_qubits);
        let mut edges = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                // Track spatial boundaries to the right
                if c + 1 < cols {
                    edges.push((i, i + 1));
                    edges.push((i + 1, i));
                }
                // Track spatial boundaries exactly downwards
                if r + 1 < rows {
                    edges.push((i, i + cols));
                    edges.push((i + cols, i));
                }
            }
        }
        backend.set_coupling_map(edges);
        backend
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
        assert!(backend
            .coupling_map
            .contains_edge(NodeIndex::new(0), NodeIndex::new(1)));
        assert!(backend
            .coupling_map
            .contains_edge(NodeIndex::new(1), NodeIndex::new(2)));
        assert!(!backend
            .coupling_map
            .contains_edge(NodeIndex::new(0), NodeIndex::new(2)));
    }

    #[test]
    fn test_grid_generator() {
        let backend = Backend::grid(2, 2);
        assert_eq!(backend.num_qubits, 4);
        // Each node in a 2x2 grid has 2 neighbors. With 4 nodes and bidirectionality: 4*2 = 8 edges
        assert_eq!(backend.coupling_map.edge_count(), 8);

        // 0 -- 1
        // |    |
        // 2 -- 3
        assert!(backend
            .coupling_map
            .contains_edge(NodeIndex::new(0), NodeIndex::new(1)));
        assert!(backend
            .coupling_map
            .contains_edge(NodeIndex::new(0), NodeIndex::new(2)));
        assert!(!backend
            .coupling_map
            .contains_edge(NodeIndex::new(0), NodeIndex::new(3))); // Diagonal rejected
    }

    #[test]
    fn test_backend_config_deserialization() {
        let json_data = r#"{
            "backend_name": "generic_5q",
            "n_qubits": 7,
            "basis_gates": ["id", "rz", "sx", "x", "cx"],
            "coupling_map": [
                [0, 1], [1, 0], 
                [1, 2], [2, 1],
                [1, 3], [3, 1]
            ]
        }"#;

        let config: BackendConfig = serde_json::from_str(json_data).unwrap();
        let backend = Backend::from_config(config);

        assert_eq!(backend.name, "generic_5q");
        assert_eq!(backend.num_qubits, 7);
        assert_eq!(backend.coupling_map.edge_count(), 6);
        assert!(backend.basis_gates.contains("cx"));
    }
}
