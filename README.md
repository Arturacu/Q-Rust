# Q-Rust

A **quantum transpiler** written in Rust, designed to parse, analyze, and optimize quantum circuits.

[![Crates.io](https://img.shields.io/crates/v/q-rust.svg)](https://crates.io/crates/q-rust)
[![Build Status](https://github.com/Arturacu/Q-Rust/workflows/Rust/badge.svg)](https://github.com/Arturacu/Q-Rust/actions)

## Features

- **Robust QASM 2.0 Support**: Parse OpenQASM 2.0 files with full support for standard gates, parameterized rotations, and measurements.
- **Modular Intermediate Representation (IR)**: A flexible, type-safe IR built for easy analysis and transformation.
- **Graph-Based Backend Topology**: Define custom hardware backends with arbitrary connectivity using efficient graph structures (`petgraph`).

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
q-rust = "0.1"
```

## Usage

### Parsing QASM

```rust
use q_rust::parser::parse_qasm;

fn main() {
    let qasm_code = r#"
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];
        creg c[2];
        h q[0];
        cx q[0], q[1];
        measure q[1] -> c[1];
    "#;

    match parse_qasm(qasm_code) {
        Ok(circuit) => {
            println!("Successfully parsed circuit with {} operations.", circuit.operations.len());
            println!("Qubits: {}, Classical Bits: {}", circuit.num_qubits, circuit.num_cbits);
        }
        Err(e) => eprintln!("Error parsing QASM: {}", e),
    }
}
```

### Defining a Backend

```rust
use q_rust::backend::Backend;

fn main() {
    let mut backend = Backend::new("MyQuantumMachine".to_string(), 5);
    
    // Define connectivity (0 -> 1 -> 2)
    backend.set_coupling_map(vec![(0, 1), (1, 2)]);
    
    println!("Backend {} has {} qubits.", backend.name, backend.num_qubits);
}
```

## License

MIT or Apache-2.0
