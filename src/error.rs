//! Unified error type for Q-Rust.
//!
//! Every fallible public API returns [`Result<T>`] (an alias for
//! `std::result::Result<T, QRustError>`). Variants are tagged with
//! `#[non_exhaustive]` so adding new variants is not a breaking change.
//!
//! See `KNOWN_GAPS_UPDATED.md` § `G-DIAG-02` for the planned migration to
//! structured error codes (e.g. `QR0001`).
//!
//! # Display contract
//!
//! Every variant renders to a single-line, human-readable message via the
//! [`std::fmt::Display`] impl synthesized by `thiserror`, suitable for CLI
//! output without further processing. The exact wording (including
//! capitalization) is part of the public API and is regression-tested by
//! `test_every_error_variant_renders_human_readable` below.

use thiserror::Error;

/// Convenience alias: `Result<T, QRustError>`.
pub type Result<T> = std::result::Result<T, QRustError>;

/// The unified Q-Rust error type.
///
/// All variants render to a human-readable single-line message via the
/// [`std::fmt::Display`] impl synthesized by `thiserror`. The wording is
/// stable across patch releases.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum QRustError {
    /// Failure parsing OpenQASM 2.0 source.
    #[error("QASM parse error: {0}")]
    ParseError(String),

    /// A referenced register, qubit, or classical bit identifier was not
    /// declared in the circuit.
    #[error("Undefined identifier: {0}")]
    Undefined(String),

    /// A register access used an index outside the register's declared bounds.
    #[error("Index out of bounds: {name}[{index}] (size = {size})")]
    IndexOutOfBounds {
        /// Register name.
        name: String,
        /// Requested index.
        index: usize,
        /// Declared register size.
        size: usize,
    },

    /// A parameter expression evaluated to an arithmetic error
    /// (e.g. division by zero, `i64::MIN.wrapping_neg()`).
    #[error("Arithmetic error: {0}")]
    Arithmetic(String),

    /// A gate name was used that does not appear in the standard library
    /// or the circuit's custom-gate registry.
    #[error("Unknown gate: {0}")]
    UnknownGate(String),

    /// A register-wide gate application had mismatched lengths.
    #[error("Register size mismatch: {0}")]
    SizeMismatch(String),

    /// A QASM feature was encountered that Q-Rust does not yet support
    /// (e.g. OpenQASM 3, `defcal`, includes other than `qelib1.inc`).
    #[error("Unsupported feature: {0}")]
    Unsupported(String),

    /// Unitary synthesis (ZYZ / KAK / QSD) failed.
    #[error("Synthesis error: {0}")]
    Synthesis(String),

    /// Routing failed (configuration error, beam-search dead-end, etc.).
    #[error("Routing error: {0}")]
    Routing(String),

    /// Simulator could not materialize a unitary or evolve a state vector.
    #[error("Simulation error: {0}")]
    Simulation(String),

    /// Custom-gate unrolling or basis decomposition failed.
    #[error("Decomposition error: {0}")]
    Decomposition(String),

    /// A library-internal invariant was violated. These should never reach
    /// user code under normal operation; please file a bug at
    /// <https://github.com/Arturacu/Q-Rust/issues> when one does.
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
    NonUniversalBasisGateSet {
        /// Human-readable explanation of why the set was rejected.
        reason: String,
    },

    /// Raised when a gate cannot be expressed in the given target basis.
    #[error("gate '{gate}' cannot be expressed in target basis {basis:?}")]
    UntranslatableGate {
        /// Gate name that could not be translated.
        gate: String,
        /// The basis set the gate was being translated into.
        basis: Vec<String>,
    },

    /// Raised by routing when a 2-qubit gate cannot be routed because the
    /// witnessed pair lies in disconnected components of the coupling map.
    #[error("disconnected coupling-map topology: cannot route {from} -> {to}")]
    DisconnectedTopology {
        /// One endpoint of the un-routable pair.
        from: usize,
        /// The other endpoint.
        to: usize,
    },

    /// Raised when a pass receives an invalid configuration parameter.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    /// Raised when the circuit's logical qubit count exceeds the backend's
    /// physical qubit count.
    #[error("insufficient qubits: circuit needs {circuit}, backend has {backend}")]
    InsufficientQubits {
        /// Logical qubit count required by the circuit.
        circuit: usize,
        /// Physical qubit count provided by the backend.
        backend: usize,
    },
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

#[cfg(test)]
mod tests {
    use super::*;

    /// [FINAL-POLISH] Pin every error variant's Display rendering as a
    /// single, human-readable, CLI-ready line. This locks the contract that
    /// downstream tooling (CLI, LSP, IDE plugins) can rely on the format.
    ///
    /// We match against lower-cased substrings so the test is robust to
    /// future capitalization changes — but the *current* wording (which
    /// preserves the pre-0.3.0 capitalization for backward compatibility
    /// with `tests/parser_test.rs`) is what `to_string()` produces today.
    #[test]
    fn test_every_error_variant_renders_human_readable() {
        let cases: Vec<(QRustError, &str)> = vec![
            (QRustError::ParseError("foo".into()), "parse error"),
            (QRustError::Undefined("x".into()), "undefined identifier"),
            (
                QRustError::IndexOutOfBounds {
                    name: "q".into(),
                    index: 5,
                    size: 3,
                },
                "index out of bounds",
            ),
            (QRustError::Arithmetic("div/0".into()), "arithmetic"),
            (QRustError::UnknownGate("foo".into()), "unknown gate"),
            (
                QRustError::SizeMismatch("len mismatch".into()),
                "size mismatch",
            ),
            (QRustError::Unsupported("qasm3".into()), "unsupported"),
            (QRustError::Synthesis("kak".into()), "synthesis"),
            (QRustError::Routing("dead end".into()), "routing"),
            (QRustError::Simulation("oom".into()), "simulation"),
            (QRustError::Decomposition("custom".into()), "decomposition"),
            (QRustError::Internal("bug".into()), "internal"),
            (QRustError::NoBasisSpecified, "no target basis specified"),
            (
                QRustError::NonUniversalBasisGateSet {
                    reason: "clifford only".into(),
                },
                "not quantum-universal",
            ),
            (
                QRustError::UntranslatableGate {
                    gate: "ccx".into(),
                    basis: vec!["rz".into(), "rx".into()],
                },
                "cannot be expressed",
            ),
            (
                QRustError::DisconnectedTopology { from: 0, to: 5 },
                "disconnected",
            ),
            (
                QRustError::InvalidConfig("beam=0".into()),
                "invalid configuration",
            ),
            (
                QRustError::InsufficientQubits {
                    circuit: 5,
                    backend: 3,
                },
                "insufficient qubits",
            ),
        ];
        for (err, expected_substring) in cases {
            let rendered = err.to_string();
            assert!(
                rendered.to_lowercase().contains(expected_substring),
                "variant rendered as {rendered:?}, expected substring {expected_substring:?}"
            );
            // Single-line invariant: no embedded newlines, suitable for CLI.
            assert!(
                !rendered.contains('\n'),
                "Display output must be a single line, got: {rendered:?}"
            );
        }
    }

    /// Pin the existing capitalized wording that `tests/parser_test.rs`
    /// asserts against. This is the second half of the Display contract:
    /// not just "human readable" but specifically the pre-0.3.0 wording.
    #[test]
    fn test_legacy_capitalized_display_wording_preserved() {
        let undef = QRustError::Undefined("r".into()).to_string();
        assert!(
            undef.contains("Undefined"),
            "legacy capitalization regression: {undef:?}"
        );
        let oob = QRustError::IndexOutOfBounds {
            name: "q".into(),
            index: 10,
            size: 1,
        }
        .to_string();
        assert!(
            oob.contains("Index out of bounds"),
            "legacy capitalization regression: {oob:?}"
        );
    }
}
