//! Pass infrastructure: the [`Pass`] trait and [`PassManager`].

use super::property_set::PropertySet;
use crate::ir::Circuit;

/// A transpiler pass: consumes a circuit, produces a (possibly) transformed one.
///
/// Passes share metadata via a [`PropertySet`] that lives on the [`PassManager`].
pub trait Pass {
    /// A short identifier used in diagnostics and debug output.
    fn name(&self) -> &str;

    /// Runs the pass, possibly reading or writing properties.
    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit;
}

/// Sequentially applies a list of [`Pass`]es, sharing a [`PropertySet`].
pub struct PassManager {
    passes: Vec<Box<dyn Pass>>,
    /// Shared metadata — exposed so tests and downstream tools can inspect it.
    pub property_set: PropertySet,
}

impl std::fmt::Debug for PassManager {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PassManager")
            .field("num_passes", &self.passes.len())
            .finish()
    }
}

impl Default for PassManager {
    fn default() -> Self {
        Self::new()
    }
}

impl PassManager {
    /// Creates a new empty manager.
    pub fn new() -> Self {
        Self {
            passes: Vec::new(),
            property_set: PropertySet::new(),
        }
    }

    /// Appends a pass.
    pub fn add_pass(&mut self, pass: Box<dyn Pass>) {
        self.passes.push(pass);
    }

    /// Runs all passes in order.
    pub fn run(&mut self, circuit: &Circuit) -> Circuit {
        let mut current = circuit.clone();
        for pass in &self.passes {
            current = pass.run(&current, &mut self.property_set);
        }
        current
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{GateType, Operation};

    struct MockPass;
    impl Pass for MockPass {
        fn name(&self) -> &str {
            "MockPass"
        }
        fn run(&self, circuit: &Circuit, _props: &mut PropertySet) -> Circuit {
            let mut c = circuit.clone();
            c.add_op(Operation::Gate {
                name: GateType::ID,
                qubits: vec![0],
                params: vec![],
            });
            c
        }
    }

    #[test]
    fn test_pass_manager() {
        let c = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        pm.add_pass(Box::new(MockPass));
        let out = pm.run(&c);
        assert_eq!(out.operations.len(), 1);
    }
}