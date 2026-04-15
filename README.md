# Q-Rust

A high-performance **quantum transpiler** written in Rust, designed to parse, analyze, and optimize quantum circuits with native OpenQASM 2.0 support and a metadata-driven architecture.

[![Crates.io](https://img.shields.io/crates/v/q-rust.svg)](https://crates.io/crates/q-rust)
[![Build Status](https://github.com/Arturacu/Q-Rust/workflows/Rust/badge.svg)](https://github.com/Arturacu/Q-Rust/actions)

## Features

### QASM 2.0 Parser
- Custom gate definitions with parameter expressions
- Mathematical expressions (`pi`, `+`, `-`, `*`, `/`)
- Register broadcasting (e.g., `h q;`)
- Parameterized gates and measurements
- Single-line comments

### Modernized Transpilation Pipeline
- **Dual IR**: Flat `Vec<Operation>` and DAG via `petgraph`
- **PassManager** with **PropertySet** metadata sharing
- Typed error handling via `QRustError`

### Optimization Passes
- `GateFusionPass`, `RotationMergePass`, `CrossConjugationPass`
- `InverseCancellationPass`, `CommutationCancellationPass`
- `SwapSimplificationPass`, `ParameterSimplificationPass`
- `GateCrystallizationPass` (float → exact Clifford/T)

### Routing
- **BeamSABRE** router with bidirectional refinement
- **SabreLayoutPass** discovers near-optimal initial layouts

## Installation

```toml
[dependencies]
q-rust = "0.2"
```

## Usage

```rust
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, TranspilerConfig};

let qasm = r#"
    OPENQASM 2.0;
    qreg q[2];
    h q[0];
    cx q[0], q[1];
"#;

let circuit = parse_qasm(qasm).expect("parse error");
let config = TranspilerConfig::builder()
    .optimization_level(3)
    .decompose_basis(true)
    .build();
let transpiled = transpile(&circuit, Some(config));
println!("Gates: {}", transpiled.operations.len());
```

## License

MIT or Apache-2.0
## Benchmarking
Run the decathlon script to benchmark.