//! # Q-Rust — A modular quantum transpiler in Rust
//!
//! Q-Rust takes a quantum circuit (parsed from OpenQASM 2.0 or constructed
//! programmatically), runs it through a configurable
//! optimization → layout → routing → synthesis → basis-translation pipeline,
//! and emits a circuit targeting a chosen hardware backend or basis-gate set.
//!
//! See the [README](https://github.com/Arturacu/Q-Rust/blob/main/README.md)
//! for the high-level architecture and feature matrix.
//!
//! ## Quick start
//!
//! ```
//! use q_rust::backend::Backend;
//! use q_rust::parser::parse_qasm;
//! use q_rust::transpiler::{transpile, TranspilerConfig};
//!
//! let qasm = r#"
//!     OPENQASM 2.0;
//!     qreg q[3];
//!     h q[0];
//!     cx q[0], q[1];
//!     cx q[1], q[2];
//! "#;
//!
//! let circuit = parse_qasm(qasm)?;
//! let cfg = TranspilerConfig::builder()
//!     .optimization_level(2)
//!     .decompose_basis(true)
//!     .backend(Backend::linear(3))
//!     .build();
//! let transpiled = transpile(&circuit, Some(cfg))?;
//! assert!(!transpiled.operations.is_empty());
//! # Ok::<(), q_rust::QRustError>(())
//! ```
//!
//! ## Module map
//!
//! - [`backend`] — hardware backend description (coupling map, basis gates).
//! - [`parser`] — OpenQASM 2.0 parser (nom-based).
//! - [`ir`] — strongly-typed circuit IR.
//! - [`transpiler`] — the optimization, layout, routing, and synthesis pipeline.
//! - [`simulator`] — unitary materialization (≤14q) and state-vector evolution (≤24q).
//! - [`verify`] — top-level equivalence verification, with automatic exact /
//!   statistical mode selection.
//! - [`error`] — the unified [`QRustError`] type.
//!
//! ## Crate features
//!
//! - `serde-ir` — derives `Serialize`/`Deserialize` on the IR types
//!   ([`ir::Circuit`], [`ir::Operation`], [`ir::GateType`],
//!   [`ir::ClassicalCondition`]). Off by default to keep the dep graph minimal.
//!
//! ## Diagnostics
//!
//! Library code is silent by default. Set the `Q_RUST_LOG` environment
//! variable to any non-empty value other than `"0"` or `"off"` to enable
//! warning-level diagnostics (e.g. KAK fallback notifications, custom-gate
//! unroll failures).

#![deny(rust_2018_idioms)]
#![warn(missing_debug_implementations)]
#![allow(
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::approx_constant,
    clippy::useless_vec,
    clippy::redundant_closure
)]

pub mod backend;
pub mod error;
pub mod ir;
pub mod parser;
pub mod simulator;
pub mod transpiler;
pub mod verify;

pub use error::{QRustError, Result};
pub use verify::{verify_equivalence, Verdict};

/// Bit-shift guard for `n`-qubit systems: the simulator and embedding routines
/// compute `1 << n` for matrix dimensions, so `n` must stay strictly below
/// this constant to avoid undefined-behavior shifts.
///
/// In practice this is dwarfed by [`simulator::MAX_QUBITS`] (= 14, the
/// hard ceiling for full unitary materialization) and
/// [`simulator::MAX_STATE_VECTOR_QUBITS`] (= 24, for state-vector evolution).
pub const SHIFT_GUARD: usize = 32;
