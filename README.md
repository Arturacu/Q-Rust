# Q-Rust

A modular quantum transpiler and OpenQASM 2.0 compiler written in Rust.

Q-Rust takes a quantum circuit — parsed from OpenQASM 2.0 or built programmatically — runs it through a configurable **optimization → layout → routing → synthesis → basis-translation** pipeline, and emits a circuit targeting a chosen hardware backend or basis-gate set. It is designed for researchers, compiler engineers, and tooling authors who want a strongly typed, memory-safe Rust implementation with custom topologies and analytic synthesis (ZYZ, KAK).

Q-Rust is a library first (`use q_rust::...`) and a CLI second (`qrust circuit.qasm`). It ships with a unitary + state-vector simulator and a verification harness: direct comparisons use full floating-point operators through 14 qubits and use seeded Haar-state sampling through 22 qubits; routed comparisons are layout-aware and sample through 20 physical qubits. Both full-operator and sampled verdicts are numerical evidence, not formal proofs. Verification supports static circuits and matching terminal readout maps; reset, mid-circuit measurement, and classical control return errors. Use `VerificationConfig` to lower the dense budget for your available memory.

---

## Table of contents

- [Quick start](#quick-start)
- [Features](#features)
- [Architecture](#architecture)
- [Directory layout](#directory-layout)
- [Usage examples](#usage-examples)
- [Configuration](#configuration)
- [Current capability boundaries](#current-capability-boundaries)
- [Testing](#testing)
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

Composable, barrier-aware optimization passes:

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
| `TwoQubitBlockSynthesisPass` | Resynthesize dependency-safe U/CX blocks after routing at level 3 |
| `RecoverPhaseGadgetsPass` | Recover exact CX–RZ–CX parity rotations before level-3 hardware routing |
| `PostBlockFusionPass` | Fuse local U runs after block synthesis with independent numerical guards |

A `CircuitProfilerPass` (analysis-only) populates a `ProfileReport` for inspection.

### Synthesis (ZYZ, KAK)

- **`ZyzSynthesizer`** — exact analytic 1-qubit synthesis.
- **`KakSynthesizer`** — analytic 2-qubit synthesis using a Cartan/KAK factorization and a joint construction with at most three CX gates. Local, CNOT-equivalent and two-CX families use reduced constructions. Every returned candidate passes a shared-phase Frobenius residual check; a failed reduced construction retries the untruncated three-CX network, and numerical failures return `None`.
- **`QsdSynthesizer`** — dispatcher to ZYZ (N=1) or KAK (N=2); N ≥ 3 is not implemented and returns `None`.
- **`NelderMead1qSynthesizer`** — numerical 1q synthesis via Nelder–Mead over a ZYZ ansatz; falls through to KAK for 2q inputs.

At optimization level 3, `TwoQubitBlockSynthesisPass` collects lowered U/CX
blocks after routing and accepts only strict CX reductions. It preserves wire
dependencies, barriers and layouts, independently checks each emitted block
with a phase-aligned Frobenius residual of at most `1e-10`, and limits the sum
of accepted residuals to `1e-9`. This numerical budget covers this pass, not
later basis conversion or serialization. Disable it with
`.two_qubit_block_synthesis(false)` for a matched ablation. Pure CX blocks use
an exact six-element binary-linear-map lookup instead of numerical synthesis.

### Routing & layout (SABRE)

- **`SabreLayoutPass`** — bidirectional iterative layout with seeded Fisher–Yates permutations, deterministic ties and convergence/zero-SWAP stopping (10 trials × 3 iterations at opt-2 and opt-3; configurable through `RoutingSearchConfig`). O3 retains a larger routing beam and post-routing cleanup.
- **`BeamSabrePass`** — beam-search SABRE router with two lookahead strategies:
  - `LookaheadStrategy::Static { weight }` — classical SABRE (Li et al. 2019).
  - `LookaheadStrategy::DynamicV2` — SABRE-v2 (Li et al. 2023).
- Fast-path for fully-connected backends (zero SWAPs).
- `Layout::from_l2p` is a validating constructor that rejects non-injective mappings.

At level 3 with a backend, exact parity-rotation recovery and diagonal-aware
routing preserve interaction structure through layout search. Non-diagonal
operations and barriers retain their dependencies; diagonal epochs are capped
at 64 interactions per wire to bound dependency growth. Disable recovery with
`.recover_phase_gadgets(false)` and diagonal scheduling through
`RoutingSearchConfig::commute_diagonal` for separate ablations. These established
rewrites can improve routing and cancellation boundaries, but heuristic routing
can also produce worse results on individual circuits.

At level 3 with basis decomposition, bounded candidate retention compares the
canonical output with up to three alternatives after lowering: an equal-SWAP
route, one refined placement, and an already-local input sequence when eligible.
Selection cannot increase two-qubit count, total depth or two-qubit depth relative
to the canonical output. Evaluating alternatives costs compilation time; disable
it with `.retain_routing_candidates(false)`. This guard does not constrain total
gate count or establish optimality.

After level-3 block synthesis, local U runs are fused without crossing a touching
CX or a non-U/CX boundary. Each replacement is checked against an independent
analytic matrix with a phase-aligned Frobenius residual of at most `1e-12` and
a circuit-wide accepted-residual budget of `1e-10`. These are floating-point
guards, not formal rounding-error certificates. Disable this step with
`.post_block_fusion(false)`. Reports expose recovery and local-fusion counters.

### Basis translation

- **`TargetBasisPass`** — applies early named-gate rewrites after a conservative universality preflight.
- **`BasisClosurePass`** — performs final exact lowering and rejects any residual gate outside the requested basis.
- **`CxDirectionPass`** — flips CX direction with H sandwiches when needed.
- **`BasisDecompositionPass`** — uses analytic decompositions from `GateDefinition`.
- Exact output closure is implemented for `{U, CX/CZ}`, `{RZ, RX, CX/CZ}`, `{RZ, SX, CX/CZ}`, and `{RZ, H, CX/CZ}` families. `{H, T, CX}` is mathematically universal but requires an approximation algorithm that Q-Rust does not yet provide, so arbitrary-angle circuits targeting it are rejected rather than silently emitted in another basis.
- Built-in backends: `linear-N`, `grid-RxC`, `ring-N`, `star-N`, `tree-N`, `all2all-N`, `ibm_quito`, `ibm_nairobi`, plus `Backend::from_json_file(path)` for custom hardware.

### Testing & validation

- **Unitary simulator** (≤14q dense floating-point operators) and **state-vector evolution** (≤24q).
- **Verification harness**: `verify_equivalence` auto-selects full-operator numerical process fidelity (≤14q) → sampled output-state fidelity (14 < n ≤ 22) → `Verdict::Unverifiable`; `verify_equivalence_with_layout` handles routed permutations and ancilla padding (sampling capped at 20 physical qubits).
- **Transpilation report** (`transpile_with_report`) — per-stage circuit metrics plus initial/final layouts when routing runs.
- 200+ unit + integration tests, including a fidelity- and basis-closure-checked algorithm suite (Bell, GHZ, QFT, Grover, Deutsch–Jozsa, Bernstein–Vazirani, VQE, QPE).

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

The pipeline architecture follows the same decomposition used by Qiskit's transpiler (optimization → layout → routing → synthesis → basis). The core algorithms — SABRE (Li et al. 2019 ASPLOS), KAK decomposition (Shende et al. 2004), and ZYZ synthesis — are standard published techniques adopted here with a Rust-native implementation.

---

## Directory layout

```
.
├── Cargo.toml
├── src/
│   ├── lib.rs              # crate root, module map, doctest
│   ├── error.rs            # unified QRustError
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
│   │   ├── target_basis.rs # validation, rewrites, exact basis closure
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

# Verify equivalence after transpilation (direct ≤22q; routed ≤20 physical q)
qrust circuit.qasm --opt 3 --verify
```

CLI flags:

| Flag | Description |
|---|---|
| `--opt N` | Optimization level (0–3) |
| `--backend SPEC` | Backend specifier (see below) |
| `--basis g1,g2,...` | Comma-separated target basis gates |
| `--report` | Print per-stage metrics |
| `--verify` | Verify against the input, accounting for routing layouts when present |
| `--no-decompose` | Skip basis decomposition |
| `--output PATH` | Write QASM to PATH (default: stdout) |

Backend specifiers: `linear-N`, `grid-RxC`, `ring-N`, `star-N`, `tree-N`, `all2all-N`, `ibm_quito`, `ibm_nairobi`, or a path to a JSON file.

### End-to-end: QASM in, verified QASM out

```rust
use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile_with_report, TranspilerConfig};
use q_rust::{verify_equivalence_with_layout, Verdict};

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
    let (transpiled, report) = transpile_with_report(&original, Some(cfg))?;
    let initial = report.initial_layout.as_deref().expect("backend layout");
    let final_layout = report.final_layout.as_deref().expect("backend layout");

    match verify_equivalence_with_layout(
        &original,
        &transpiled,
        initial,
        final_layout,
    )? {
        Verdict::ExactlyEquivalent { fidelity } => {
            println!("✓ exact match, fidelity = {fidelity:.6}");
        }
        Verdict::StatisticallyEquivalent { samples, .. } => {
            println!("✓ statistical match over {samples} Haar samples");
        }
        Verdict::NotEquivalent { fidelity, method } => {
            panic!("regression ({method}): fidelity = {fidelity}")
        }
        Verdict::Unverifiable { reason } => println!("skipped: {reason}"),
    }

    println!("{}", transpiled.to_qasm(None));
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
    .decompose_basis(true)      // decompose and enforce target-basis closure
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
| `retain_routing_candidates` | bool, default true | Bounded level-3 quality search; extra compilation work |

When both a backend basis and an explicit `target_basis` are supplied, the
explicit basis wins. Gate aliases are normalized (`cnot` → `cx`, `u3` → `u`,
`u1`/`p` → `rz`). Unsupported exact targets return an error before the circuit
is emitted.

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

## Current capability boundaries

- **Local synthesis does not establish whole-circuit optimality.** KAK uses at most three CX gates on successful two-qubit synthesis. Block collection is conservative, numerical failures retain the original block, and routing can still dominate output cost. The established two-qubit CX bound is not a global optimality guarantee.
- **Basis universality does not imply implemented lowering.** Exact closure currently covers the four one-qubit families listed above with CX or CZ. Approximate Clifford+T synthesis is not implemented.
- **Statistical verification is not a certificate.** It reports the minimum observed output-state fidelity over seeded Haar samples. Full-operator numerical mode uses normalized Hilbert–Schmidt/process fidelity and is invariant under global phase.
- **Routed verification needs layout metadata.** Use `transpile_with_report` with `verify_equivalence_with_layout`, or use the CLI's `--verify` flag, which does this automatically.
- **General synthesis stops at two qubits.** The QSD dispatcher returns `None` for `n ≥ 3`; named gates such as CCX use explicit templates.
- **Backend presets are static descriptions.** The Quito and Nairobi constructors provide topology/basis fixtures; they do not query live calibrations or submit hardware jobs.

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
- `tests/two_qubit_block_regressions.rs` — independent analytic all-input checks for synthesis, dependency boundaries, routed isometries and W/QFT reductions.
- `tests/qft_w_refinement_regressions.rs` — independent Fourier/W operators, directed ancilla routing, commuting fences, long diagonal epochs, readout mappings and local-fusion guards.
- `tests/cu1_routing_regressions.rs`, `tests/logical_swap_regressions.rs` — independent all-input checks for controlled-phase dependencies, internal wire permutations, custom/high-arity gates and routed readout mappings.
- `tests/routing_candidate_regressions.rs` — local SWAP sequences, QFT/PEA quality bounds, independent all-input checks, directed/spare-wire mappings and report/plain parity.
- `tests/cli_smoke_test.rs` — `qrust` binary, gated `#[ignore]`.

The library has zero `unsafe` blocks. The `try_*` variants of public APIs return `Result` and are the preferred entry points; infallible wrappers (`decompose_basis`, `unroll_custom_gates`) fall back to returning the original circuit on error and emit a diagnostic via `Q_RUST_LOG`.

### Reproducing the cross-tool benchmarks

The [candidate-retention follow-up](experiments/EXTERNAL_RETENTION_FOLLOWUP.md)
recovers the Trotter-4, QFT-4 and PEA-5 line regressions and further improves
AE/QPE output quality. The same 100-case comparison has no CX-count, total-depth
or CX-depth regressions, but compilation is about 17% slower overall and three
cases gain one total gate. The report includes controlled diagnostics, complete
input checks for the focal outputs, and a separate replayable evidence archive;
all earlier results remain unchanged.

The [external routing follow-up](experiments/EXTERNAL_ROUTING_FOLLOWUP.md)
investigates controlled-phase dependencies, logical SWAP absorption and repeated
layout work. It reruns the same 100 cases, with exact-source ablations and
additional full-input checks. AE/QPE-10 line outputs improve substantially;
three other line cases regress. The report retains both outcomes, the remaining
Qiskit speed advantage on the focal cases, and the unchanged earlier evidence.

The [external compatibility follow-up](experiments/EXTERNAL_COMPATIBILITY_FOLLOWUP.md)
fixes Q-Rust CU1 support and Cirq adapter closure. All six recipes now pass
the numerical checks on 100/100 cases in both sessions. It reports speed,
CX count and depth on both the original shared cases and the expanded cohort,
with separate raw evidence and unchanged original results.

The [original external benchmark study](experiments/EXTERNAL_BENCHMARKS.md) compares
50 circuits from pinned QASMBench and MQT Bench sources across two topologies,
six recipes and two sessions of seven trials. It retains input/adapter failures
and independent output checks. That version of Q-Rust accepts 82 of 100 cases; speed and quality
comparisons are explicitly conditional on those successful cases. Methods,
recipe settings, per-case results and compressed raw evidence are included.

The earlier [QFT/W refinement](experiments/QFT_W_REFINEMENT.md) reports improved
QFT routing/cancellation boundaries, W-state depth cleanup, independent
correctness checks and the measured quality/runtime trade-offs. It retains
all adverse cases and two separately planned before/after sessions.

The earlier [two-qubit block investigation](experiments/BLOCK_SYNTHESIS_INVESTIGATION.md)
records the preceding block-only implementation and its historical results.

Use the [repeated comparison study](experiments/COMPARISON_STUDY.md) for the
corrected fixtures, pinned current dependencies, matched U/CX targets, seven
interleaved trials and independent validation of every saved output:

The earlier [routing investigation](experiments/ROUTING_INVESTIGATION.md) reports
its frozen before/after measurements, alternative Cirq/TKET recipes,
and remaining quality regressions. Its compact results and primary cell table
are version-controlled; complete raw campaigns accompany the investigation.

```bash
python3.12 -m venv .venv-comparison
.venv-comparison/bin/python -m pip install -r experiments/requirements-comparison.txt
cargo build --release -p qrust-bench --bin comparison_runner
.venv-comparison/bin/python experiments/comparison_study.py time experiments/results/NEW_COMPARISON \
  --binary target/release/comparison_runner
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-comparison/bin/python \
  experiments/comparison_study.py validate experiments/results/NEW_COMPARISON
.venv-comparison/bin/python experiments/analyze_comparison_study.py experiments/results/NEW_COMPARISON
```

Generated JSON, CSV, LaTeX, and figure payloads under `experiments/results/`
are intentionally ignored; the runner code, pinned dependencies, fixture
manifest, and OpenQASM 2.0 fixtures are version-controlled. The older
`bench_baselines.py` results are historical and are superseded by the validated
comparisons above. Compilation speed, CX count and depth are separate metrics;
this study does not establish mathematical optimality or a causal programming-language
advantage.

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
