use crate::ir::ast::ParsedStatement;
use std::collections::HashMap;

/// Represents the definition of a custom gate.
#[derive(Debug, Clone, PartialEq)]
pub struct CustomGateDef {
    /// Variables for angle parameters (e.g. "theta", "phi")
    pub params: Vec<String>,
    /// Placeholders for qubit arguments
    pub qubits: Vec<String>,
    /// The body of the gate definition (statements to be expanded)
    pub body: Vec<ParsedStatement>,
}

/// A registry storing all user-defined gates from the QASM program.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GateRegistry {
    /// Map from gate name to its definition
    pub defs: HashMap<String, CustomGateDef>,
}

impl GateRegistry {
    /// Creates a new empty GateRegistry
    pub fn new() -> Self {
        Self {
            defs: HashMap::new(),
        }
    }

    /// Registers a custom gate
    pub fn register(
        &mut self,
        name: String,
        params: Vec<String>,
        qubits: Vec<String>,
        body: Vec<ParsedStatement>,
    ) {
        self.defs.insert(
            name,
            CustomGateDef {
                params,
                qubits,
                body,
            },
        );
    }

    /// Fetches a custom gate definition by name
    pub fn get(&self, name: &str) -> Option<&CustomGateDef> {
        self.defs.get(name)
    }
}
