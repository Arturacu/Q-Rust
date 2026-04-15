//! Custom gate definition registry.

use crate::ir::ast::ParsedStatement;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct CustomGateDef {
    pub params: Vec<String>,
    pub qubits: Vec<String>,
    pub body: Vec<ParsedStatement>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct GateRegistry {
    pub defs: HashMap<String, CustomGateDef>,
}

impl GateRegistry {
    pub fn new() -> Self {
        Self::default()
    }

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

    #[inline]
    pub fn get(&self, name: &str) -> Option<&CustomGateDef> {
        self.defs.get(name)
    }
}