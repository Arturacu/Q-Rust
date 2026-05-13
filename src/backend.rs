//! Hardware backend description.
//!
//! A [`Backend`] models a quantum device's qubit count, native gate set, and
//! coupling map. Construct one from a JSON [`BackendConfig`] or via the
//! topology constructors ([`Backend::linear`], [`Backend::grid`], etc.).

use crate::error::{QRustError, Result};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::Directed;
use serde::{Deserialize, Serialize};
use std::collections::{HashSet, VecDeque};
use std::path::Path;

/// JSON-deserializable backend description.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackendConfig {
    /// Human-readable backend identifier.
    pub backend_name: String,
    /// Number of physical qubits.
    pub n_qubits: usize,
    /// Native gate names (lower-case OpenQASM-style).
    pub basis_gates: Vec<String>,
    /// Directed coupling map as `[from, to]` pairs.
    pub coupling_map: Vec<[usize; 2]>,
}

/// Hardware backend: qubit count, native gate set, and coupling graph.
#[derive(Debug, Clone)]
pub struct Backend {
    /// Human-readable backend identifier.
    pub name: String,
    /// Number of physical qubits.
    pub num_qubits: usize,
    /// Native gate names (lower-case).
    pub basis_gates: HashSet<String>,
    /// Directed coupling map: an edge `(u, v)` means a 2-qubit gate with
    /// control `u` and target `v` is natively supported.
    pub coupling_map: Graph<(), (), Directed>,
}

impl Backend {
    /// Creates a backend with `num_qubits` qubits, an empty basis, and no
    /// coupling edges.
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

    /// Adds a gate name to the backend's native basis set.
    pub fn add_basis_gate(&mut self, gate: impl Into<String>) {
        self.basis_gates.insert(gate.into());
    }

    /// Replaces the coupling map with the given directed edges.
    /// Self-loops and out-of-range edges are silently dropped.
    pub fn set_coupling_map(&mut self, edges: impl IntoIterator<Item = (usize, usize)>) {
        self.coupling_map.clear_edges();
        for (u, v) in edges {
            if u < self.num_qubits && v < self.num_qubits && u != v {
                self.coupling_map
                    .add_edge(NodeIndex::new(u), NodeIndex::new(v), ());
            }
        }
    }

    /// Constructs a backend from a JSON [`BackendConfig`].
    pub fn from_config(config: BackendConfig) -> Self {
        let mut backend = Backend::new(config.backend_name, config.n_qubits);
        for gate in config.basis_gates {
            backend.add_basis_gate(gate);
        }
        backend.set_coupling_map(config.coupling_map.into_iter().map(|arr| (arr[0], arr[1])));
        backend
    }

    /// [E2E-NEW-FEATURE] Loads a backend from a JSON file path. Convenience
    /// wrapper around `from_config + serde_json::from_str + fs::read_to_string`.
    ///
    /// # Errors
    /// Returns [`QRustError::ParseError`] if the file is missing, unreadable,
    /// or contains invalid JSON.
    pub fn from_json_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| QRustError::ParseError(format!("cannot read {}: {e}", path.display())))?;
        let cfg: BackendConfig = serde_json::from_str(&text)
            .map_err(|e| QRustError::ParseError(format!("invalid backend JSON: {e}")))?;
        Ok(Backend::from_config(cfg))
    }

    /// Constructs IBM Quito's 5-qubit T-shaped heavy-hex topology.
    ///
    /// Coupling: 0—1—2, with 1—3 and 3—4 forming the T:
    /// ```text
    ///   0 — 1 — 2
    ///       |
    ///       3 — 4
    /// ```
    /// Native basis gates: `{id, rz, sx, x, cx}`.
    ///
    /// Note: this is a *programmatic* construction (no fixture files),
    /// so the library has no runtime file dependency.
    pub fn ibm_quito() -> Self {
        let mut backend = Backend::new("ibm_quito", 5);
        for g in ["id", "rz", "sx", "x", "cx"] {
            backend.add_basis_gate(g);
        }
        backend.set_coupling_map([
            (0, 1),
            (1, 0),
            (1, 2),
            (2, 1),
            (1, 3),
            (3, 1),
            (3, 4),
            (4, 3),
        ]);
        backend
    }

    /// Constructs IBM Nairobi's 7-qubit heavy-hex topology.
    ///
    /// Coupling:
    /// ```text
    ///   0 — 1 — 2
    ///       |
    ///       3
    ///       |
    ///   4 — 5 — 6
    /// ```
    /// Native basis gates: `{id, rz, sx, x, cx}`.
    pub fn ibm_nairobi() -> Self {
        let mut backend = Backend::new("ibm_nairobi", 7);
        for g in ["id", "rz", "sx", "x", "cx"] {
            backend.add_basis_gate(g);
        }
        backend.set_coupling_map([
            (0, 1),
            (1, 0),
            (1, 2),
            (2, 1),
            (1, 3),
            (3, 1),
            (3, 5),
            (5, 3),
            (4, 5),
            (5, 4),
            (5, 6),
            (6, 5),
        ]);
        backend
    }

    /// Builds a fully-connected (all-to-all) topology.
    pub fn all_to_all(num_qubits: usize) -> Self {
        let mut backend = Backend::new("all_to_all", num_qubits);
        let edges = (0..num_qubits)
            .flat_map(|i| (0..num_qubits).filter_map(move |j| (i != j).then_some((i, j))));
        backend.set_coupling_map(edges);
        backend
    }

    /// Builds a linear (1-D chain) topology with bidirectional edges.
    pub fn linear(num_qubits: usize) -> Self {
        let mut backend = Backend::new("linear", num_qubits);
        let edges = (0..num_qubits.saturating_sub(1)).flat_map(|i| [(i, i + 1), (i + 1, i)]);
        backend.set_coupling_map(edges);
        backend
    }

    /// Builds a `rows × cols` grid topology with bidirectional edges.
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

    /// Builds a ring topology (linear chain with wrap-around).
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

    /// Builds a star topology with qubit 0 as the central hub.
    pub fn star(num_qubits: usize) -> Self {
        let mut backend = Backend::new(format!("star_{num_qubits}"), num_qubits);
        let edges = (1..num_qubits).flat_map(|i| [(0, i), (i, 0)]);
        backend.set_coupling_map(edges);
        backend
    }

    /// Builds a complete binary tree topology: each qubit `i > 0` has
    /// parent `(i - 1) / 2`.
    pub fn tree(num_qubits: usize) -> Self {
        let mut backend = Backend::new(format!("tree_{num_qubits}"), num_qubits);
        let edges = (1..num_qubits).flat_map(|i| {
            let parent = (i - 1) / 2;
            [(parent, i), (i, parent)]
        });
        backend.set_coupling_map(edges);
        backend
    }

    /// Returns `true` iff a 2-qubit gate can run between `q1` and `q2` in
    /// either direction.
    #[inline]
    pub fn is_adjacent(&self, q1: usize, q2: usize) -> bool {
        self.coupling_map
            .contains_edge(NodeIndex::new(q1), NodeIndex::new(q2))
            || self
                .coupling_map
                .contains_edge(NodeIndex::new(q2), NodeIndex::new(q1))
    }

    /// Returns `true` iff the directed edge `u -> v` exists in the
    /// coupling map.
    #[inline]
    pub fn has_directed_edge(&self, u: usize, v: usize) -> bool {
        if u >= self.num_qubits || v >= self.num_qubits {
            return false;
        }
        self.coupling_map
            .contains_edge(NodeIndex::new(u), NodeIndex::new(v))
    }

    /// Returns the out-neighbors of qubit `q` in the coupling map.
    pub fn neighbors(&self, q: usize) -> Vec<usize> {
        self.coupling_map
            .neighbors(NodeIndex::new(q))
            .map(|n| n.index())
            .collect()
    }

    /// Returns `true` iff every pair of distinct qubits is directly
    /// connected (i.e., the coupling map is complete / all-to-all).
    pub fn is_fully_connected(&self) -> bool {
        let n = self.num_qubits;
        if n <= 1 {
            return true;
        }
        for i in 0..n {
            for j in 0..n {
                if i != j && !self.is_adjacent(i, j) {
                    return false;
                }
            }
        }
        true
    }

    /// Returns an `n × n` matrix of shortest-path distances between every
    /// pair of qubits, computed via BFS. `usize::MAX` indicates no path.
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
    fn test_is_fully_connected() {
        assert!(Backend::all_to_all(4).is_fully_connected());
        assert!(!Backend::linear(4).is_fully_connected());
        assert!(Backend::all_to_all(1).is_fully_connected());
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

    /// [E2E-NEW-FEATURE] Built-in IBM backends are constructed
    /// programmatically (no fixture file dependency).
    #[test]
    fn test_ibm_quito_builtin() {
        let b = Backend::ibm_quito();
        assert_eq!(b.num_qubits, 5);
        assert_eq!(b.name, "ibm_quito");
        assert!(b.basis_gates.contains("cx"));
        assert!(b.is_adjacent(0, 1));
        assert!(b.is_adjacent(1, 2));
        assert!(b.is_adjacent(1, 3));
        assert!(b.is_adjacent(3, 4));
        assert!(!b.is_adjacent(0, 2));
        assert!(!b.is_adjacent(0, 4));
    }

    #[test]
    fn test_ibm_nairobi_builtin() {
        let b = Backend::ibm_nairobi();
        assert_eq!(b.num_qubits, 7);
        assert_eq!(b.name, "ibm_nairobi");
        assert!(b.basis_gates.contains("cx"));
        assert!(b.is_adjacent(1, 3));
        assert!(b.is_adjacent(3, 5));
        assert!(b.is_adjacent(5, 6));
        assert!(!b.is_adjacent(0, 6));
    }

    #[test]
    fn test_from_json_file_roundtrip() {
        let b = Backend::from_json_file("tests/fixtures/ibm_quito_5q.json").unwrap();
        assert_eq!(b.num_qubits, 5);
        assert!(b.basis_gates.contains("cx"));
    }

    #[test]
    fn test_from_json_file_missing_returns_err() {
        let r = Backend::from_json_file("/nonexistent/path/foo.json");
        assert!(matches!(r, Err(QRustError::ParseError(_))));
    }
}
