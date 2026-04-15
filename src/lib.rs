//! # Q-Rust

#![deny(rust_2018_idioms)]
#![warn(missing_debug_implementations)]

pub mod backend;
pub mod error;
pub mod ir;
pub mod parser;
pub mod simulator;
pub mod transpiler;

pub use error::{QRustError, Result};
// /// Q-Rust: A quantum compiler



