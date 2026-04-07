# Q-Rust

A high-performance **quantum transpiler** written in Rust, designed to parse, analyze, and optimize quantum circuits with native OpenQASM 2.0 support and a metadata-driven architecture.

[![Crates.io](https://img.shields.io/crates/v/q-rust.svg)](https://crates.io/crates/q-rust)
[![Build Status](https://github.com/Arturacu/Q-Rust/workflows/Rust/badge.svg)](https://github.com/Arturacu/Q-Rust/actions)

## Features

### QASM 2.0 Parser
- ✅ **Custom Gate Definitions**: Define parameterized gates with expressions (supports multi-qubit definitions)
- ✅ **Mathematical Expressions**: Support for `pi`, arithmetic operations (`+`, `-`, `*`, `/`)
- ✅ **Register Broadcasting**: Apply gates to entire registers (e.g., `h q;`)
- ✅ **Parameterized Gates**: Rotation gates with variable angles
- ✅ **Measurements & Barriers**: Measurement and barrier support
- ✅ **Comments**: Single-line comments with `//`

### Modernized Transpilation Pipeline
Q-Rust uses a modular, metadata-driven pipeline inspired by modern compiler architecture (e.g., Qiskit 1.x, LLVM):
- ✅ **Dual-Layer IR**: 
  - **Flat IR**: High-performance sequential `Vec<Operation>` for execution and unrolling.
  - **DAG IR**: Directed Acyclic Graph based on `petgraph`, enabling $O(1)$ topological dependency analysis.
- ✅ **PassManager**: Sequence-controlled pipeline execution with unified error handling.
- ✅ **PropertySet**: Centralized metadata storage allowing decoupled passes to share structural analysis without redundant circuit scans.

### Advanced Optimization Passes
- **Topological Optimization**:
  - `GateFusionPass`: Collapses consecutive single-qubit sequences into individual $U$-gates via matrix multiplication.
  - `CommutationCancellationPass`: Identifies and cancels non-adjacent `CX` and rotation gates by analyzing commutation through intermediate operations.
  - `InverseCancellationPass`: Removes adjacent self-inverse pairs ($H$-$H$, $X$-$X$, $CZ$-$CZ$, $CCX$-$CCX$).
- **Algebraic Simplification**:
  - `RotationMergePass`: Geometrically combines adjacent rotation axes ($RZ(\theta_1) \cdot RZ(\theta_2) \to RZ(\theta_1+\theta_2)$), supporting controlled rotations (`CRX`, `CRY`, `CRZ`).
  - `ParameterSimplificationPass`: Normalizes rotation angles (modulo $2\pi$) and removes negligible rotations.
  - `CrossConjugationPass`: Implements structural rewrites like $H \cdot RZ(\theta) \cdot H \to RX(\theta)$.

### Correctness & Performance
- ✅ **Unitary Simulator**: Native Rust simulator to compute $2^n \times 2^n$ unitary matrices for mathematical verification.
- ✅ **Equivalence Testing**: Automated fidelity checks to ensure transpilation maintains unitary correctness (up to global phase).
- ✅ **Rust Advantage**: Zero-latency binary execution with no Python/GIL overhead, enabling real-time streaming transpilation.

## Supported Gates

### Standard Gates
- **Pauli**: `x`, `y`, `z`
- **Hadamard**: `h`
- **Phase**: `s`, `sdg`, `t`, `tdg`
- **Identity**: `id`

### Parametric Gates
- **Rotations**: `rx(θ)`, `ry(θ)`, `rz(θ)`
- **Universal**: `u3(θ, φ, λ)`, `u1(λ)`, `u2(φ, λ)`

### Multi-Qubit & Controlled Gates
- **Two-Qubit**: `cx`, `cz`, `cy`, `ch`, `csx`, `swap`
- **Three-Qubit**: `ccx` (Toffoli)
- **Controlled Rotations**: `crx(θ)`, `cry(θ)`, `crz(θ)`
- **Ising Interactions**: `rxx(θ)`, `ryy(θ)`, `rzz(θ)`

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
q-rust = "0.1.1"
```

## Usage Example

```rust
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[2];
        h q[0];
        cx q[0], q[1];
        // Redundant cancellation
        cz q[0], q[1];
        cz q[0], q[1];
    "#;

    let circuit = parse_qasm(qasm).expect("Parse error");
    
    let config = TranspilerConfig {
        optimization_level: 3, // Enable full DAG-based optimizations
        decompose_basis: true, // Unroll to {U, CX}
        backend: None,         // Use default all-to-all coupling
    };

    let transpiled = transpile(&circuit, Some(config));
    println!("Gates: {}", transpiled.operations.len());
}
```

## Testing Infrastructure

Run the full suite of fixture-driven tests:
```sh
cargo test
```

Fixtures are located in `tests/fixtures/`. Each circuit is verified for unitary fidelity against its original form across multiple optimization levels.

## Limitations

The following OpenQASM 2.0 features are **not yet supported**:

- ❌ **Include Statements**: e.g., `include "qelib.inc";` (inline gate definitions instead)
- ❌ **OpenQASM 3.0**: Only version 2.0 is supported
  - ❌ **Conditional Operations**: `if(c==1) x q[0];`

## References

This project adheres to the **OpenQASM 2.0** specification:

- **Specification**: [OpenQASM 2.0 on GitHub](https://github.com/openqasm/openqasm/tree/OpenQASM2.x)
- **Paper**: Andrew W. Cross, Lev S. Bishop, John A. Smolin, Jay M. Gambetta, *"Open Quantum Assembly Language"*, [arXiv:1707.03429](https://arxiv.org/abs/1707.03429)

## License

MIT or Apache-2.0
