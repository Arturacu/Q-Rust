//! Unified error type for Q-Rust.

use thiserror::Error;

pub type Result<T> = std::result::Result<T, QRustError>;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum QRustError {
    #[error("QASM parse error: {0}")]
    ParseError(String),
    #[error("Undefined identifier: {0}")]
    Undefined(String),
    #[error("Index out of bounds: {name}[{index}] (size = {size})")]
    IndexOutOfBounds {
        name: String,
        index: usize,
        size: usize,
    },
    #[error("Arithmetic error: {0}")]
    Arithmetic(String),
    #[error("Unknown gate: {0}")]
    UnknownGate(String),
    #[error("Register size mismatch: {0}")]
    SizeMismatch(String),
    #[error("Unsupported feature: {0}")]
    Unsupported(String),
    #[error("Synthesis error: {0}")]
    Synthesis(String),
    #[error("Routing error: {0}")]
    Routing(String),
    #[error("Simulation error: {0}")]
    Simulation(String),
    #[error("Decomposition error: {0}")]
    Decomposition(String),
    #[error("Internal invariant violation: {0}")]
    Internal(String),
    /// Raised when `decompose_basis = true` but no target gate set was provided.
    #[error(
        "no target basis specified: set `target_basis` on TranspilerConfig or \
         provide a Backend with non-empty `basis_gates`"
    )]
    NoBasisSpecified,
    /// Raised when the provided gate set is not quantum-universal.
    #[error("provided gate set is not quantum-universal: {reason}")]
    NonUniversalBasisGateSet { reason: String },
    /// Raised when a gate cannot be expressed in the given target basis.
    #[error("gate '{gate}' cannot be expressed in target basis {basis:?}")]
    UntranslatableGate { gate: String, basis: Vec<String> },
}

impl From<&str> for QRustError {
    fn from(s: &str) -> Self {
        QRustError::Internal(s.to_string())
    }
}

impl From<String> for QRustError {
    fn from(s: String) -> Self {
        QRustError::Internal(s)
    }
}
