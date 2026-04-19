//! Hardware backend description.
//!
//! Models the physical qubit topology (coupling map), basis gate set, and
//! metadata of a quantum processor. Used by routing and decomposition passes
//! to produce hardware-executable circuits.

use petgraph::graph::{Graph, NodeIndex};
use petgraph::Directed;
use serde::{Deserialize, Serialize};
use std::collections::{HashSet, VecDeque};

/// Serialization surrogate for JSON backend definitions (compatible with
/// common cloud-API hardware descriptions such as IBM's coupling maps).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendConfig {
    pub backend_name: String,
    pub n_qubits: usize,
    pub basis_gates: Vec<String>,
    pub coupling_map: Vec<[usize; 2]>,
}

/// A quantum backend: its qubits, basis gates, and connectivity.
#[derive(Debug, Clone)]
pub struct Backend {
    /// Backend identifier (e.g. "ibm_quito").
    pub name: String,
    /// Number of physical qubits.
    pub num_qubits: usize,
    /// Native basis gate set.
    pub basis_gates: HashSet<String>,
    /// Directed coupling graph; nodes are physical qubits.
    pub coupling_map: Graph<(), (), Directed>,
}

impl Backend {
    /// Creates an empty backend with `num_qubits` isolated nodes.
    pub fn new(name: impl Into<String>, num_qubits: usize) -> Self {
        let mut graph = Graph::new();
        for _ in 0..num_qubits {
            graph.add_node(());
        }
        Self {
            name: name.into(),
            num_qubits,
            basis_gates: HashSet::new(),
            coupling_map: graph,
        }
    }

    /// Registers a basis gate by name.
    pub fn add_basis_gate(&mut self, gate: impl Into<String>) {
        self.basis_gates.insert(gate.into());
    }

    /// Replaces the coupling map with the given directed edges.
    /// Edges referencing out-of-range qubits are silently discarded.
    pub fn set_coupling_map(&mut self, edges: impl IntoIterator<Item = (usize, usize)>) {
        self.coupling_map.clear_edges();
        for (u, v) in edges {
            if u < self.num_qubits && v < self.num_qubits && u != v {
                self.coupling_map
                    .add_edge(NodeIndex::new(u), NodeIndex::new(v), ());
            }
        }
    }

    /// Builds a [`Backend`] from a JSON-deserialized [`BackendConfig`].
    pub fn from_config(config: BackendConfig) -> Self {
        let mut backend = Backend::new(config.backend_name, config.n_qubits);
        for gate in config.basis_gates {
            backend.add_basis_gate(gate);
        }
        backend.set_coupling_map(config.coupling_map.into_iter().map(|arr| (arr[0], arr[1])));
        backend
    }

    /// Fully-connected (all-to-all) topology. Useful for pure simulation.
    pub fn all_to_all(num_qubits: usize) -> Self {
        let mut backend = Backend::new("all_to_all", num_qubits);
        let edges = (0..num_qubits)
            .flat_map(|i| (0..num_qubits).filter_map(move |j| (i != j).then_some((i, j))));
        backend.set_coupling_map(edges);
        backend
    }

    /// 1D linear nearest-neighbor chain.
    pub fn linear(num_qubits: usize) -> Self {
        let mut backend = Backend::new("linear", num_qubits);
        let edges = (0..num_qubits.saturating_sub(1))
            .flat_map(|i| [(i, i + 1), (i + 1, i)]);
        backend.set_coupling_map(edges);
        backend
    }

    /// 2D grid topology.
    pub fn grid(rows: usize, cols: usize) -> Self {
        let num_qubits = rows * cols;
        let mut backend = Backend::new(format!("grid_{rows}x{cols}"), num_qubits);
        let mut edges = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                if c + 1 < cols {
                    edges.push((i, i + 1));
                    edges.push((i + 1, i));
                }
                if r + 1 < rows {
                    edges.push((i, i + cols));
                    edges.push((i + cols, i));
                }
            }
        }
        backend.set_coupling_map(edges);
        backend
    }

    /// Ring topology (cyclic chain).
    pub fn ring(num_qubits: usize) -> Self {
        let mut backend = Backend::new(format!("ring_{num_qubits}"), num_qubits);
        if num_qubits > 1 {
            let edges = (0..num_qubits).flat_map(|i| {
                let n = (i + 1) % num_qubits;
                [(i, n), (n, i)]
            });
            backend.set_coupling_map(edges);
        }
        backend
    }

    /// Star topology (qubit 0 is the hub).
    pub fn star(num_qubits: usize) -> Self {
        let mut backend = Backend::new(format!("star_{num_qubits}"), num_qubits);
        let edges = (1..num_qubits).flat_map(|i| [(0, i), (i, 0)]);
        backend.set_coupling_map(edges);
        backend
    }

    /// Binary-tree topology.
    pub fn tree(num_qubits: usize) -> Self {
        let mut backend = Backend::new(format!("tree_{num_qubits}"), num_qubits);
        let edges = (1..num_qubits).flat_map(|i| {
            let parent = (i - 1) / 2;
            [(parent, i), (i, parent)]
        });
        backend.set_coupling_map(edges);
        backend
    }

    /// Returns `true` iff qubits `q1` and `q2` are directly connected (either direction).
    #[inline]
    pub fn is_adjacent(&self, q1: usize, q2: usize) -> bool {
        self.coupling_map
            .contains_edge(NodeIndex::new(q1), NodeIndex::new(q2))
            || self
                .coupling_map
                .contains_edge(NodeIndex::new(q2), NodeIndex::new(q1))
    }

    /// Returns `true` iff there is a directed edge `u -> v` in the coupling map.
    ///
    /// For undirected backends (where every edge is stored in both directions),
    /// this is equivalent to checking either direction. For directed hardware
    /// (e.g. asymmetric CX-native devices), this differentiates allowed CX
    /// orientation from its reverse.
    #[inline]
    pub fn has_directed_edge(&self, u: usize, v: usize) -> bool {
        if u >= self.num_qubits || v >= self.num_qubits {
            return false;
        }
        self.coupling_map
            .contains_edge(NodeIndex::new(u), NodeIndex::new(v))
    }

    /// Returns the set of physical qubits directly adjacent to `q`.
    pub fn neighbors(&self, q: usize) -> Vec<usize> {
        self.coupling_map
            .neighbors(NodeIndex::new(q))
            .map(|n| n.index())
            .collect()
    }

    /// All-pairs shortest-path distances in hops (BFS).
    /// Disconnected pairs are reported as `usize::MAX`.
    pub fn shortest_path_matrix(&self) -> Vec<Vec<usize>> {
        let n = self.num_qubits;
        let mut dist = vec![vec![usize::MAX; n]; n];
        for start in 0..n {
            dist[start][start] = 0;
            let mut queue = VecDeque::new();
            queue.push_back(start);
            while let Some(current) = queue.pop_front() {
                let d = dist[start][current];
                for neighbor in self.neighbors(current) {
                    if dist[start][neighbor] == usize::MAX {
                        dist[start][neighbor] = d + 1;
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        dist
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backend_creation() {
        let backend = Backend::new("test_backend", 5);
        assert_eq!(backend.num_qubits, 5);
        assert_eq!(backend.coupling_map.node_count(), 5);
        assert_eq!(backend.coupling_map.edge_count(), 0);
    }

    #[test]
    fn test_coupling_map() {
        let mut backend = Backend::new("test_backend", 3);
        backend.set_coupling_map([(0, 1), (1, 2)]);
        assert_eq!(backend.coupling_map.edge_count(), 2);
        assert!(backend.is_adjacent(0, 1));
        assert!(!backend.is_adjacent(0, 2));
    }

    #[test]
    fn test_has_directed_edge() {
        let mut b = Backend::new("t", 3);
        b.set_coupling_map([(0, 1), (1, 2)]);
        assert!(b.has_directed_edge(0, 1));
        assert!(!b.has_directed_edge(1, 0));
        assert!(b.has_directed_edge(1, 2));
        assert!(!b.has_directed_edge(2, 1));
        // Linear is bidirectional.
        let lin = Backend::linear(3);
        assert!(lin.has_directed_edge(0, 1));
        assert!(lin.has_directed_edge(1, 0));
    }

    #[test]
    fn test_grid_generator() {
        let backend = Backend::grid(2, 2);
        assert_eq!(backend.num_qubits, 4);
        assert_eq!(backend.coupling_map.edge_count(), 8);
        assert!(backend.is_adjacent(0, 1));
        assert!(backend.is_adjacent(0, 2));
        assert!(!backend.is_adjacent(0, 3));
    }

    #[test]
    fn test_backend_config_deserialization() {
        let json = r#"{
            "backend_name": "generic_5q",
            "n_qubits": 7,
            "basis_gates": ["id", "rz", "sx", "x", "cx"],
            "coupling_map": [[0,1],[1,0],[1,2],[2,1],[1,3],[3,1]]
        }"#;
        let config: BackendConfig = serde_json::from_str(json).unwrap();
        let backend = Backend::from_config(config);
        assert_eq!(backend.name, "generic_5q");
        assert_eq!(backend.num_qubits, 7);
        assert_eq!(backend.coupling_map.edge_count(), 6);
        assert!(backend.basis_gates.contains("cx"));
    }

    #[test]
    fn test_is_adjacent() {
        let backend = Backend::linear(5);
        assert!(backend.is_adjacent(0, 1));
        assert!(backend.is_adjacent(1, 0));
        assert!(!backend.is_adjacent(0, 4));
    }

    #[test]
    fn test_neighbors() {
        let backend = Backend::linear(5);
        assert_eq!(backend.neighbors(0), vec![1]);
        let mut mid = backend.neighbors(2);
        mid.sort();
        assert_eq!(mid, vec![1, 3]);
    }

    #[test]
    fn test_shortest_path_matrix() {
        let backend = Backend::linear(5);
        let dist = backend.shortest_path_matrix();
        assert_eq!(dist[0][4], 4);
        assert_eq!(dist[1][3], 2);

        let grid = Backend::grid(2, 3);
        let gd = grid.shortest_path_matrix();
        assert_eq!(gd[0][5], 3);
        assert_eq!(gd[0][4], 2);
    }
}
// /// Hardware backend representation




