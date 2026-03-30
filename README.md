# Q-Rust

A **quantum transpiler** written in Rust, designed to parse, analyze, and optimize quantum circuits with OpenQASM 2.0 support.

[![Crates.io](https://img.shields.io/crates/v/q-rust.svg)](https://crates.io/crates/q-rust)
[![Build Status](https://github.com/Arturacu/Q-Rust/workflows/Rust/badge.svg)](https://github.com/Arturacu/Q-Rust/actions)

## Features

### QASM 2.0 Parser
- ✅ **Custom Gate Definitions**: Define parameterized gates with expressions
- ✅ **Mathematical Expressions**: upport for `pi`, arithmetic operations (`+`, `-`, `*`, `/`)
- ✅ **Register Broadcasting**: Apply gates to entire registers (e.g., `h q;`)
- ✅ **Parameterized Gates**: Rotation gates with variable angles
- ✅ **Measurements & Barriers**: Measurement and barrier support
- ✅ **Comments**: Single-line comments with `//`

### Transpilation & Optimization
- ✅ **Basis Decomposition**: automatic conversion to `{U, CX}` universal basis
- ✅ **Single-Qubit Synthesis**: ZYZ decomposition into Euler angles
- ✅ **Two-Qubit Synthesis**: KAK decomposition using the magic basis
- ✅ **Optimization Passes**: 
  - Gate Fusion (merging consecutive single-qubit gates)
  - CX/SWAP Cancellation (removing redundant pairs)
  - Parameter Simplification (modulo 2π, negligible rotations)
  - Peephole Optimization (self-inverse cancellation)

### Verifying Correctness
- ✅ **Unitary Simulator**: computes the full 2ⁿ × 2ⁿ unitary matrix for a circuit
- ✅ **Equivalence Testing**: verifies transpilation by comparing unitaries up to a global phase
- ✅ **Fixture-Driven Tests**: a suite of 21 tests running various QASM circuits through different transpiler configurations

## Supported Gates

### Standard Gates
- **Pauli**: `x`, `y`, `z`
- **Hadamard**: `h`
- **Phase**: `s`, `sdg`, `t`, `tdg`
- **Identity**: `id`

### Parametric Gates
- **Rotations**: `rx(θ)`, `ry(θ)`, `rz(θ)`
- **Universal**: `U(θ, φ, λ)`, `u1(λ)`, `u2(φ, λ)`, `u3(θ, φ, λ)`

### Multi-Qubit Gates
- **CNOT**: `cx`
- **Swap**: `swap`
- **Toffoli**: `ccx` (3-qubit controlled-NOT)

### Custom Gates
- User-defined gates via `gate` definitions with parameters and expressions (single-qubit supported)

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
q-rust = "0.1.1"
```

## Usage Examples

### Transpiling to a Basis Set

```rust
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[1];
        h q[0];
    "#;

    let circuit = parse_qasm(qasm).unwrap();
    let config = TranspilerConfig {
        decompose_basis: true,
        optimization_level: 1,
    };

    let transpiled = transpile(&circuit, Some(config));
    // H gate is now decomposed into U(pi/2, 0, pi)
}
```

### Verifying Unitary Equivalence

```rust
use q_rust::simulator::{circuit_to_unitary, unitary_fidelity};
// ... after transpilation ...
let u_orig = circuit_to_unitary(&circuit);
let u_trans = circuit_to_unitary(&transpiled);

let fidelity = unitary_fidelity(&u_orig, &u_trans);
assert!(fidelity > 0.99999999); // Equivalent up to global phase
```

## Testing Infrastructure

The project includes a robust test suite that uses QASM files as fixtures.

To run the full suite:
```sh
cargo test
```

To run the transpiler fixture suite specifically:
```sh
cargo test --test transpiler_suite
```

Fixtures are located in `tests/fixtures/`. Each fixture is run through multiple transpiler configurations (Default, No Decomposition, Decompose Only) to ensure universal invariants like qubit count preservation and unitary fidelity.

## Limitations

The following OpenQASM 2.0 features are **not yet supported**:

- ❌ **Include Statements**: e.g., `include "qelib.inc";` (inline gate definitions instead)
- ❌ **Conditional Operations**: `if(c==1) x q[0];`
- ❌ **Multi-qubit Custom Gates**: Only single-qubit definitions are currently supported
- ❌ **OpenQASM 3.0**: Only version 2.0 is supported

## References

This project adheres to the **OpenQASM 2.0** specification:

- **Specification**: [OpenQASM 2.0 on GitHub](https://github.com/openqasm/openqasm/tree/OpenQASM2.x)
- **Paper**: Andrew W. Cross, Lev S. Bishop, John A. Smolin, Jay M. Gambetta, *"Open Quantum Assembly Language"*, [arXiv:1707.03429](https://arxiv.org/abs/1707.03429)

## License

MIT or Apache-2.0

