//! # Q-Rust

#![deny(rust_2018_idioms)]
#![warn(missing_debug_implementations)]
#![allow(
    clippy::needless_range_loop,
    clippy::type_complexity,
    clippy::approx_constant,
    clippy::useless_vec,
    clippy::redundant_closure,
    unused_imports
)]

pub mod backend;
pub mod error;
pub mod ir;
pub mod parser;
pub mod simulator;
pub mod transpiler;

pub use error::{QRustError, Result};
// /// Q-Rust: A quantum compiler
