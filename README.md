# Q-Rust

A **quantum transpiler** written in Rust, designed to parse, analyze, and optimize quantum circuits with full OpenQASM 2.0 support.

[![Crates.io](https://img.shields.io/crates/v/q-rust.svg)](https://crates.io/crates/q-rust)
[![Build Status](https://github.com/Arturacu/Q-Rust/workflows/Rust/badge.svg)](https://github.com/Arturacu/Q-Rust/actions)

## Features

### QASM 2.0 Parser
- ✅ **Custom Gate Definitions**: Define parameterized gates with expressions
- ✅ **Mathematical Expressions**: Full support for `pi`, arithmetic operations (`+`, `-`, `*`, `/`)
- ✅ **Register Broadcasting**: Apply gates to entire registers (e.g., `h q;`)
- ✅ **Parameterized Gates**: Rotation gates with variable angles
- ✅ **Measurements & Barriers**: Full measurement and barrier support
- ✅ **Comments**: Single-line comments with `//`

### Intermediate Representation (IR)
- ✅ **Type-Safe Operations**: Enum-based gate types with compile-time safety
- ✅ **Circuit Validation**: Automatic validation for measurements and qubit bounds
- ✅ **Flexible Operations**: Gate, Measure, Reset, Barrier operations

### Backend System
- ✅ **Custom Topologies**: Define arbitrary qubit connectivity
- ✅ **Graph-Based**: Efficient coupling maps using `petgraph`

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
- User-defined gates via `gate` definitions with parameters and expressions

## Limitations

The following OpenQASM 2.0 features are **not yet supported**:

- ❌ **Include Statements**: `include "qelib1.inc";` will error (inline gate definitions instead)
- ❌ **Conditional Operations**: `if(c==1) x q[0];` not supported
- ❌ **OpenQASM 3.0**: Only version 2.0 is supported

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
q-rust = "0.1"
```

## Usage Examples

### Basic Circuit with Standard Gates

```rust
use q_rust::parser::parse_qasm;

fn main() {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[2];
        creg c[2];
        h q[0];
        cx q[0], q[1];
        measure q -> c;
    "#;

    let circuit = parse_qasm(qasm).expect("Failed to parse");
    println!("Circuit: {} qubits, {} operations", 
             circuit.num_qubits, circuit.operations.len());
}
```

### Custom Gate Definitions with Expressions

```rust
use q_rust::parser::parse_qasm;

fn main() {
    let qasm = r#"
        OPENQASM 2.0;
        
        // Define custom gate with parameters
        gate my_rotation(theta) q {
            rx(theta) q;
            rz(pi/2) q;
        }
        
        qreg q[1];
        my_rotation(1.57) q[0];
    "#;

    let circuit = parse_qasm(qasm).expect("Failed to parse");
    println!("Custom gate expanded successfully!");
}
```

### Defining a Backend

```rust
use q_rust::backend::Backend;

fn main() {
    let mut backend = Backend::new("MyQuantumMachine".to_string(), 5);
    
    // Linear topology: 0 -> 1 -> 2 -> 3 -> 4
    backend.set_coupling_map(vec![(0, 1), (1, 2), (2, 3), (3, 4)]);
    
    println!("Backend '{}' with {} qubits", backend.name, backend.num_qubits);
}
```

## Error Handling

Q-Rust enforces strict **OpenQASM 2.0** compliance:

- **Version Check**: Only `OPENQASM 2.0;` is accepted (3.0 will error)
- **Required Header**: Files must start with `OPENQASM 2.0;`
- **Validation**: `Circuit::validate()` warns about missing measurements
- **Bounds Checking**: Qubit and classical bit indices are validated

## References

This project adheres to the **OpenQASM 2.0** specification:

- **Specification**: [OpenQASM 2.0 on GitHub](https://github.com/openqasm/openqasm/tree/OpenQASM2.x)
- **Paper**: Andrew W. Cross, Lev S. Bishop, John A. Smolin, Jay M. Gambetta, *"Open Quantum Assembly Language"*, [arXiv:1707.03429](https://arxiv.org/abs/1707.03429)

## License

MIT or Apache-2.0
