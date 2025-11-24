use crate::ir::Circuit;

/// A trait for transpiler passes.
///
/// A pass takes a circuit and returns a transformed circuit.
pub trait Pass {
    /// Returns the name of the pass.
    fn name(&self) -> &str;

    /// Runs the pass on the given circuit.
    fn run(&self, circuit: &Circuit) -> Circuit;
}

/// Manages a sequence of transpiler passes.
pub struct PassManager {
    passes: Vec<Box<dyn Pass>>,
}

impl PassManager {
    /// Creates a new empty PassManager.
    pub fn new() -> Self {
        Self { passes: Vec::new() }
    }

    /// Adds a pass to the manager.
    pub fn add_pass(&mut self, pass: Box<dyn Pass>) {
        self.passes.push(pass);
    }

    /// Runs all passes in sequence on the given circuit.
    pub fn run(&self, circuit: &Circuit) -> Circuit {
        let mut current_circuit = circuit.clone();
        for pass in &self.passes {
            current_circuit = pass.run(&current_circuit);
        }
        current_circuit
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

        fn run(&self, circuit: &Circuit) -> Circuit {
            let mut new_circuit = circuit.clone();
            // Add a dummy gate to verify the pass ran
            new_circuit.add_op(Operation::Gate {
                name: GateType::ID,
                qubits: vec![0],
                params: vec![],
            });
            new_circuit
        }
    }

    #[test]
    fn test_pass_manager() {
        let circuit = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        pm.add_pass(Box::new(MockPass));

        let new_circuit = pm.run(&circuit);
        assert_eq!(new_circuit.operations.len(), 1);
        match &new_circuit.operations[0] {
            Operation::Gate { name, .. } => assert_eq!(*name, GateType::ID),
            _ => panic!("Expected ID gate"),
        }
    }
}
