# Correctness repair and validation — 2026-10-15

This repair gates fresh benchmarking on implementation regressions and independent algorithm checks. Prior benchmark outputs remain historical; QFT and W-state numbers from the old fixtures must not be reused as measurements of the corrected algorithms.

## Repaired behavior

- **Ancillas and layouts:** layout-aware numerical checks keep leakage probability instead of projecting and renormalizing it away. Logical identity versus H on an unused physical qubit now fails. Layout shape, injectivity, bounds, register widths, finite parameters, and gate arities are validated; malformed input returns an error.
- **Nonunitary operations:** low-level unitary simulators reject measurement, reset, and classical control. High-level verification explicitly compares matching terminal classical readout maps and then unitary prefixes. Reset, mid-circuit readout, conditional control, repeated classical writes, and different readout maps are outside that verifier contract. A differing map produces an error, not a general proof of channel inequivalence.
- **Compilation and readout:** compilation rejects unsupported dynamic operations before optimization. Valid terminal measurements are detached, then appended in their original order with the final physical mapping. Barriers remain routing boundaries. Compilation preserves repeated classical-write order, although the verifier has a narrower contract.
- **Routing and targets:** canonical compilation uses fallible pass execution; disconnected or impossible routing cannot silently return an illegal circuit. SWAP routing uses undirected connectivity. CX direction repair runs after entangling-gate lowering, and final output is checked against width, topology, direction, and requested basis.
- **Controlled rotations:** inverse cancellation uses a 4π period for controlled rotations; a controlled 2π rotation is a relative phase on the control. Rotation merging requires adjacency on every operand, so an intervening noncommuting operation cannot be crossed through a different DAG edge.
- **QFT:** CP(θ) is implemented as P(θ/2) on the control followed by CRZ(θ), with an explicit little-endian Fourier convention and final bit reversal. Both benchmark and Rust test fixtures were corrected.
- **W state:** the generator transfers a single excitation with equal positive final amplitudes; it no longer leaves vacuum population or an unused final qubit.
- **Benchmark evidence:** the runner calls the canonical public compiler and configurable high-level verifier. Parsing and total transpilation are timed separately. Every raw trial retains its output QASM, layouts, hashes, verdict, and timings. Failed repeats remain visible and prevent a successful cell status. Summary numbers are successful-trial medians, not a fictitious representative output. The orchestrator checkpoints results and refuses to replace an existing checkpoint.

## Regression coverage and reproduction

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
cargo test --test cli_smoke_test -- --ignored
python3 -m unittest discover -s experiments -p 'test_*.py'
cargo build --release -p qrust-bench
python3 experiments/bench_runner.py --runner target/release/runner \
  --output-dir experiments/results/NEW_RUN_DIRECTORY --reps 1
```

The new suites cover verifier false positives and malformed inputs, routing legality and readout order, controlled-rotation identities and adjacency, analytic Fourier matrices through six qubits, W amplitudes through fourteen qubits, fixture/manifest synchronization, and preservation of failed benchmark repetitions. Representative regressions were run against the old implementation and demonstrated the failures before repair.

`verify_benchmark_qiskit.py` is an optional independent checker using Qiskit and NumPy. It compares the **full physical isometry**, including dirty-ancilla output branches, with one common global phase and explicit terminal readout maps. It also checks small source QFT and W fixtures against analytic expectations. It consumes the exact saved trial outputs, rather than rerunning randomized routing. The run report records actual dependency versions and every checked, failed, or skipped item.

```sh
python3 experiments/verify_benchmark_qiskit.py \
  --records experiments/results/NEW_RUN_DIRECTORY/raw_benchmark_data.json \
  --output experiments/results/NEW_RUN_DIRECTORY/qiskit_validation.json
```

## Interpretation limits

Full operator comparison is floating-point numerical evidence, not a formal proof. Haar sampling is finite-sample evidence, without a universal error-detection guarantee. The benchmark uses dense checks through eight physical qubits, eight samples above that, and two samples from eighteen qubits; layout-aware sampling is capped at twenty physical qubits. Wider outputs are retained as **unverifiable**, never counted as successful validation.

The library keeps its historical default dense threshold of fourteen qubits for compatibility. That is a resource ceiling, not a memory guarantee; callers should use `VerificationConfig` to set a practical lower dense budget, as the benchmark does. Gate implementations still share decomposition code in places. Independent analytic rotation tests and Qiskit output checks reduce common-mode risk but do not formally verify every gate or pass.

Source-level `two_q` counters historically count gates acting on at least two qubits (including CCX); output U/CX counts are native two-qubit costs. Source-to-output reduction percentages are not an optimality measure. Cleanup-disabled runs are independent randomized compilations, not paired observations of one routed circuit. Kernel RSS is captured before verification; optional psutil polling observes the whole process and can include simulator memory. Timings include canonical report construction and are not yet a controlled cross-tool performance claim.

Qiskit is used here as an independent correctness oracle. Qiskit/Cirq/TKET performance baselines were not rerun in this repair. A publication comparison still needs frozen versions and targets, repeated trials, timeout/failure denominators, and independent validation of every output that contributes to a quality claim.

## Recorded results

- Full workspace/all-features suite: **305 passed**. The four normally ignored CLI tests were run separately: **4 passed**. Python fixture/harness tests: **11 passed**. Formatting and strict Clippy checks pass.
- Main sweep: **597 configurations**, one trial each; **584 passed numerical validation** (241 dense operator, 343 sampled), **13 explicitly unverifiable** at widths 22/25. No detected inequivalence or compilation error.
- QFT/W repeat sweep: **33 configurations × 3 trials = 99 trials**, all passed numerical validation.
- Independent Qiskit checks: **241 main outputs + 69 repeated outputs passed**, with zero failures. The other 356 main and 30 repeated outputs exceed the independent checker’s eight-physical-qubit budget and are explicitly skipped. Eight small corrected source fixtures also pass analytic checks. Maximum normalized full-isometry residual was approximately **1.38 × 10⁻¹⁰**.

The committed [run summary](validation/correctness_20261016_summary.json) and [provenance](validation/correctness_20261016_provenance.json) record counts, configuration, versions, and source/fixture hashes. Full raw outputs, checkpoints, CSV, logs, and independent check reports are generated local artifacts under `results/correctness_20261016/` and are not version-controlled; the directory label is the UTC run date. The commands above reproduce those artifacts in a fresh directory. These counts describe tested instances, not a percentage proof of correctness or an optimality guarantee.
