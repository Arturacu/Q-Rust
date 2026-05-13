# Q-Rust

A modular quantum transpiler and OpenQASM 2.0 compiler written in Rust.

Q-Rust takes a quantum circuit — parsed from OpenQASM 2.0 or built programmatically — runs it through a configurable **optimization → layout → routing → synthesis → basis-translation** pipeline, and emits a circuit targeting a chosen hardware backend or basis-gate set. It is designed for researchers, compiler engineers, and tooling authors who want a strongly-typed, panic-free Rust alternative to Python-centric transpiler stacks, with first-class support for IBM heavy-hex devices, custom topologies, and analytic synthesis (ZYZ, KAK).

Q-Rust is a library first (`use qrust::...`) and a CLI second (`qrust circuit.qasm`). It ships with a unitary + state-vector simulator and a verification harness so you can check transpiled circuits for equivalence against the input, either exactly (≤14 qubits) or statistically via Haar sampling (≤22 qubits).

---

## Table of contents

- [Quick start](#quick-start)
- [Features](#features)
- [Architecture](#architecture)
- [Directory layout](#directory-layout)
- [Usage examples](#usage-examples)
- [Configuration](#configuration)
- [Testing](#testing)
- [Origin](#origin)
- [Contributing](#contributing)
- [License](#license)
- [Issues & contact](#issues--contact)

---

## Quick start

**Requirements:** Rust 1.75 or newer.

```bash
git clone https://github.com/Arturacu/Q-Rust
cd Q-Rust
cargo build --release
cargo test
```

Run the CLI on a sample circuit:

```bash
cargo run --release --bin qrust -- tests/fixtures/ghz_3.qasm --opt 2 --report
```

Or use it as a library — add to your `Cargo.toml`:

```toml
[dependencies]
q-rust = "0.3"

# Optional: derive Serialize/Deserialize on IR types.
# q-rust = { version = "0.3", features = ["serde-ir"] }
```

---

## Features

### Parsing & IR

- Full OpenQASM 2.0 parser (nom-based): `include`, `gate` definitions, `if(c==v) op` conditionals, barriers, resets, register-wide application.
- Strongly-typed IR: `Circuit`, `Operation`, `GateType` (30+ variants including `ECR` and `iSWAP`), `GateDefinition`, `GateRegistry`, `CommutationSignature`.
- All IR enums are `#[non_exhaustive]`; round-trip QASM emission is stable.
- Optional `serde-ir` feature derives `Serialize`/`Deserialize` on every IR type.

### Optimization passes

An 8-pass, barrier-aware optimization pipeline:

| Pass | Purpose |
|---|---|
| `GateFusionPass` | Merge consecutive single-qubit gates |
| `CommutationCancellationPass` | Cancel commuting inverse pairs |
| `InverseCancellationPass` | Eliminate adjacent inverses |
| `SwapSimplificationPass` | Remove redundant swaps |
| `RotationMergePass` | Combine same-axis rotations |
| `CrossConjugationPass` | Push rotations through Cliffords |
| `ParameterSimplificationPass` | Fold parameter expressions |
| `GateCrystallizationPass` | Collapse to canonical forms |

A `CircuitProfilerPass` (analysis-only) populates a `ProfileReport` for inspection.

### Synthesis (ZYZ, KAK)

- **`ZyzSynthesizer`** — exact analytic 1-qubit synthesis.
- **`KakSynthesizer`** — exact analytic 2-qubit synthesis with Weyl-chamber branching (0/1/3-CX cases) and a fidelity-pinned analytic fast-path for known 2q gates.
- **`QsdSynthesizer`** — dispatcher that routes to ZYZ/KAK; QSD for N ≥ 3 is a stub.
- **`NelderMead1qSynthesizer`** — numerical 1q synthesis with KAK fall-through for 2q.

### Routing & layout (SABRE)

- **`SabreLayoutPass`** — bidirectional iterative layout with Fisher–Yates seed permutations (10 trials × 3 iterations at opt-2; 10 × 5 at opt-3).
- **`BeamSabrePass`** — beam-search SABRE router with two lookahead strategies:
  - `LookaheadStrategy::Static { weight }` — classical SABRE (Li et al. 2019).
  - `LookaheadStrategy::DynamicV2` — SABRE-v2 (Li et al. 2023).
- Fast-path for fully-connected backends (zero SWAPs).
- `Layout::from_l2p` is a validating constructor that rejects non-injective mappings.

### Basis translation (multi-vendor)

- **`TargetBasisPass`** with universality validation — rejects Clifford-only sets.
- **`CxDirectionPass`** — flips CX direction with H sandwiches when needed.
- **`BasisDecompositionPass`** — uses analytic decompositions from `GateDefinition`.
- Built-in backends: `linear-N`, `grid-RxC`, `ring-N`, `star-N`, `tree-N`, `all2all-N`, `ibm_quito`, `ibm_nairobi`, plus `Backend::from_json_file(path)` for custom hardware.

### Testing & validation

- **Unitary simulator** (≤14q exact) and **state-vector evolution** (≤24q).
- **Verification harness** (`verify_equivalence`): auto-selects exact unitary fidelity (≤14q) → Haar sampling (14 < n ≤ 22) → `Verdict::Unverifiable` (>22q).
- **Transpilation report** (`transpile_with_report`) — per-stage circuit metrics.
- 100+ unit + integration tests, including a fidelity-verified algorithm suite (Bell, GHZ, QFT, Grover, Deutsch–Jozsa, Bernstein–Vazirani, VQE, QPE).

---

## Architecture

```
        ┌────────────────────────────────────────────────────────┐
        │                       Q-Rust                           │
        │                                                        │
QASM 2.0│   ┌────────┐   ┌──────┐   ┌─────────────┐             │
  ──────┼──▶│ Parser ├──▶│  IR  ├──▶│ Optimization├─────┐       │
text    │   └────────┘   └──────┘   │  (8 passes) │     │       │
        │                           └─────────────┘     ▼       │
        │                                       ┌─────────────┐ │
        │                                       │   Layout    │ │
        │                                       │   (SABRE)   │ │
        │                                       └──────┬──────┘ │
        │                                              ▼        │
        │                                       ┌─────────────┐ │
        │                                       │   Routing   │ │
        │                                       │ (BeamSABRE) │ │
        │                                       └──────┬──────┘ │
        │                                              ▼        │
        │                                       ┌─────────────┐ │
        │                                       │  Synthesis  │ │
        │                                       │ (ZYZ / KAK) │ │
        │                                       └──────┬──────┘ │
        │                                              ▼        │
        │                                       ┌─────────────┐ │
        │                                       │    Basis    │ │
        │                                       │ Translation │ │  ┌───────────┐
        │                                       └──────┬──────┘ ├─▶│ Simulator │
        │                                              ▼        │  │ + Verify  │
        │                                       ┌─────────────┐ │  └───────────┘
        │                                       │ QASM emit / │ │
        │                                       │  Report     │ │
        │                                       └─────────────┘ │
        └────────────────────────────────────────────────────────┘
```

The pipeline is driven by a `PassManager`. Each pass implements the `Pass` trait and reads/writes a shared `PropertySet` (carrying e.g. `initial_layout`, `final_layout`, `swaps_inserted`).

---

## Directory layout

```
.
├── Cargo.toml
├── src/
│   ├── lib.rs              # crate root, module map, doctest
│   ├── error.rs            # unified QRustError (15 variants)
│   ├── backend.rs          # topology + basis-gate descriptions
│   ├── parser/             # OpenQASM 2.0 (nom)
│   ├── ir/                 # Circuit, Operation, GateType, registry
│   ├── transpiler/
│   │   ├── pass.rs         # Pass trait, PassManager, PropertySet
│   │   ├── optimization.rs # 8 optimization passes
│   │   ├── layout.rs       # SabreLayoutPass
│   │   ├── routing.rs      # BeamSabrePass + lookahead strategies
│   │   ├── synthesis/      # ZYZ, KAK, QSD, numerical, qsearch
│   │   ├── decomposition.rs
│   │   ├── target_basis.rs # multi-vendor basis translation
│   │   ├── dag.rs          # DAG IR + scheduling
│   │   ├── profiler.rs     # CircuitProfilerPass
│   │   └── report.rs       # TranspilationReport
│   ├── simulator.rs        # ≤14q unitary, ≤24q state-vector
│   ├── verify.rs           # verify_equivalence, Verdict
│   └── bin/qrust.rs        # CLI entry point
├── tests/
│   ├── fixtures/           # *.qasm sample circuits
│   ├── parser_test.rs
│   ├── integration_test.rs
│   ├── routing_suite.rs
│   ├── transpiler_suite.rs
│   ├── e2e_known_circuits.rs       # fidelity-verified algorithms
│   ├── e2e_pipeline_integration.rs # parse → transpile → reparse
│   └── cli_smoke_test.rs           # #[ignore]-gated CLI smoke
└── examples/
    ├── transpile_e2e.rs
    ├── routing_demo.rs
    ├── compare_qrust.rs
    ├── debug_kak.rs
    ├── export_qrust_for_qiskit.rs
    └── benchmark_*.rs
```

---

## Usage examples

### Library

```rust
use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, TranspilerConfig};

fn main() -> Result<(), q_rust::QRustError> {
    let qasm = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[1];
        cx q[1], q[2];
    "#;

    let circuit = parse_qasm(qasm)?;

    let cfg = TranspilerConfig::builder()
        .optimization_level(2)
        .decompose_basis(true)
        .backend(Backend::linear(3))
        .build();

    let transpiled = transpile(&circuit, Some(cfg))?;
    println!("{} operations", transpiled.operations.len());
    Ok(())
}
```

### CLI

```bash
# Basic transpilation
qrust circuit.qasm

# Optimize for IBM Quito with a per-stage report
qrust circuit.qasm --backend ibm_quito --opt 3 --report

# Custom backend, custom basis, file output
qrust circuit.qasm --backend my_device.json --basis rz,sx,cx --output out.qasm

# Verify equivalence after transpilation (≤22 qubits)
qrust circuit.qasm --opt 3 --verify
```

CLI flags:

| Flag | Description |
|---|---|
| `--opt N` | Optimization level (0–3) |
| `--backend SPEC` | Backend specifier (see below) |
| `--basis g1,g2,...` | Comma-separated target basis gates |
| `--report` | Print per-stage metrics |
| `--verify` | Run `verify_equivalence` against the input |
| `--no-decompose` | Skip basis decomposition |
| `--output PATH` | Write QASM to PATH (default: stdout) |

Backend specifiers: `linear-N`, `grid-RxC`, `ring-N`, `star-N`, `tree-N`, `all2all-N`, `ibm_quito`, `ibm_nairobi`, or a path to a JSON file.

### End-to-end: QASM in, verified QASM out

```rust
use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, TranspilerConfig};
use q_rust::{verify_equivalence, Verdict};

fn main() -> Result<(), q_rust::QRustError> {
    let input = r#"
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0], q[1];
        cx q[0], q[2];
    "#;

    let original = parse_qasm(input)?;

    let cfg = TranspilerConfig::builder()
        .optimization_level(3)
        .decompose_basis(true)
        .backend(Backend::ibm_quito())
        .build();
    let transpiled = transpile(&original, Some(cfg))?;

    match verify_equivalence(&original, &transpiled)? {
        Verdict::ExactlyEquivalent { fidelity } => {
            println!("✓ exact match, fidelity = {fidelity:.6}");
        }
        Verdict::StatisticallyEquivalent { samples, .. } => {
            println!("✓ statistical match over {samples} Haar samples");
        }
        Verdict::NotEquivalent { reason } => panic!("regression: {reason}"),
        Verdict::Unverifiable { reason } => println!("skipped: {reason}"),
    }

    println!("{}", transpiled.to_qasm());
    Ok(())
}
```

---

## Configuration

### Selecting a backend

```rust
use q_rust::backend::Backend;

let b = Backend::linear(5);                 // 1D chain
let b = Backend::grid(3, 4);                // 3×4 lattice
let b = Backend::ring(8);                   // ring topology
let b = Backend::all_to_all(6);             // fully connected
let b = Backend::ibm_quito();               // 5q heavy-hex
let b = Backend::ibm_nairobi();             // 7q heavy-hex
let b = Backend::from_json_file("dev.json")?; // custom JSON
```

A backend JSON file describes basis gates and a coupling map; see `tests/fixtures/` for the schema.

### Configuring the pipeline

```rust
use q_rust::transpiler::TranspilerConfig;

let cfg = TranspilerConfig::builder()
    .optimization_level(2)      // 0 = none, 3 = aggressive
    .decompose_basis(true)      // run TargetBasisPass
    .backend(Backend::linear(4))
    .build();
```

The main knobs:

| Knob | Range | Notes |
|---|---|---|
| `optimization_level` | 0–3 | Controls pass selection and SABRE trial counts |
| `decompose_basis` | bool | Whether to translate to the backend's basis |
| `backend` | `Backend` | Coupling map + basis gates for layout/routing |
| `lookahead_strategy` | `Static` / `DynamicV2` | Routing heuristic |

### Custom pass pipelines

```rust
use q_rust::transpiler::pass::PassManager;
use q_rust::transpiler::optimization::{GateFusionPass, RotationMergePass};

let mut pm = PassManager::new();
pm.add_pass(Box::new(GateFusionPass::default()));
pm.add_pass(Box::new(RotationMergePass::default()));
let out = pm.run(&circuit)?;
```

### Diagnostics

Library code is silent by default. Set `Q_RUST_LOG` to any non-empty value other than `"0"` or `"off"` to enable warning-level diagnostics (KAK fallback notices, custom-gate unroll failures, etc.):

```bash
Q_RUST_LOG=1 cargo run --bin qrust -- circuit.qasm
```

---

## Testing

```bash
# All unit + integration tests
cargo test

# A specific test file
cargo test --test e2e_known_circuits

# CLI smoke tests (gated behind #[ignore])
cargo test -- --ignored

# Documentation tests, including the lib.rs doctest
cargo test --doc

# Lints
cargo clippy --all-targets -- -D warnings

# Build the docs locally
cargo doc --no-deps --open
```

Test layout:

- `tests/parser_test.rs` — OpenQASM 2.0 surface area + error messages.
- `tests/integration_test.rs` — full pipeline smoke tests.
- `tests/routing_suite.rs`, `tests/transpiler_suite.rs` — pass-level coverage.
- `tests/e2e_known_circuits.rs` — fidelity-verified algorithm suite.
- `tests/e2e_pipeline_integration.rs` — parse → transpile → emit → re-parse.
- `tests/cli_smoke_test.rs` — `qrust` binary, gated `#[ignore]`.

The library has zero `unsafe` blocks, and `clippy::unwrap_used` is denied on the lib target — public APIs are panic-free.

---

## Origin

This codebase was hardened through an automated, multi-pass review-and-improvement pipeline (5 review/fix iterations followed by 3 polish passes) on top of an earlier prototype. Each pass ran a Principal-Engineer review and applied compile-validated patches before convergence.

For the curious, two companion documents go into more detail:

- `CODE_UPDATES_SUMMARY.md` — what changed in the source, by area.
- The sister thesis update document — the methodology and findings.

These are background reading; nothing about using Q-Rust depends on them.

---

## Contributing

Contributions are welcome. Before opening a PR:

1. `cargo fmt --all`
2. `cargo clippy --all-targets -- -D warnings`
3. `cargo test` (and `cargo test -- --ignored` if your change touches the CLI)
4. Add a test for any bug fix or new feature.
5. Update the docstrings on any public API you touch — `cargo doc --no-deps` should build cleanly.

Public-API changes should be flagged in the PR description; the crate is in `0.x.y`, so semantic breaks are allowed but should be deliberate.

---

## License

Dual-licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

---

## Issues & contact

- **Bug reports & feature requests:** <https://github.com/Arturacu/Q-Rust/issues>
- **Maintainer:** Arturo Acuaviva — `arturoacuaviva@gmail.com>`
- **Repository:** <https://github.com/Arturacu/Q-Rust>
- **Docs:** <https://docs.rs/q-rust>