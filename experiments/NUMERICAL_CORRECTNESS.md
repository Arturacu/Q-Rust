# Numerical correctness repair and existing benchmark revalidation

17 October 2026. This report supersedes correctness and output-quality claims
about the current compiler in earlier reports. Their timings remain historical
measurements of their frozen binaries; they are not measurements of this repair.
No external benchmark study or new speed comparison was run here.

## Problem and repair

The former compiler discarded each of 1,000 `RX(9e-10)` rotations. Their
composition is `RX(9e-7)`, whose residual from identity is approximately
`4.5e-7`, exceeding the independent benchmark oracle's unchanged `1e-7` limit.
The repaired serialized output has independent analytic residual
`7.7292247596e-21` for this case.

The audit and fixes cover the following distinct boundaries:

| Boundary | Previous risk | Current behavior |
|---|---|---|
| Parameter simplification | Per-gate `1e-9` deletion accumulated without a circuit budget | Default candidate threshold is zero; small nonzero rotations survive; changed parameters must pass direct residual admission |
| Gate crystallization | Near-Pauli, phase and Hadamard gates snapped to exact names | Default matching uses exact parameter values; explicitly supplied epsilon only proposes candidates, never bypasses residual admission |
| CX commutation | Almost-diagonal/almost-X U gates treated as commuting | Structural equalities replace epsilon decisions |
| Inverse cancellation | Nearly cancelling angle pairs deleted | Exact period candidates plus reconstruction check; controlled rotations retain the 4π convention |
| Rotation merge and U fusion | Floating arithmetic/decomposition accepted without checking the replacement | Direct local matrix reconstruction and shared allowance; overflow or failed checks retain original operations |
| Angle canonicalization | Tiny negative parameters lost precision on adding 2π; large modulo reductions could be inaccurate | Preserve small signed values; guard any changed normalization against the unreduced operator |
| Block, local fusion and KAK fallback | Separate invocations could reset numerical allowances | Share one cumulative ledger across passes and both pipeline stages, while retaining stricter local/pass caps |
| QASM serialization | Fixed ten-decimal formatting rounded parameters after checks | Decimal formatting round-trips finite f64 parameters, including tiny values |
| Public ZYZ and numerical search APIs | Invalid targets or unconverged results could be accepted | Finite/unitary/reconstruction admission; requested residual bound is enforced; unsupported numerical synthesis returns `None` |
| Symbolic rational angle metadata | Saturating period multiplication/wrapping sign normalization changed fractions at i64 boundaries | Wide exact intermediates and representability-aware sign normalization; not used by the measured pipeline |

`NelderMead1qSynthesizer.tol` now means a direct phase-aligned Frobenius
residual, rather than trace-overlap infidelity. Existing callers may require
more iterations; exhausted search returns `None` instead of silently accepting
a weaker tolerance. This API is not used by the normal comparison pipeline.

## Numerical policy and limits

`numerical_policy.rs` records accepted reconstruction residuals and a
floating-point guard against a shared `1e-9` allowance. Legacy local rewrites
must additionally satisfy `1e-12`; block/local synthesis keep their existing
stricter admission rules. Exhaustion retains the original operation sequence.
The plain and reporting entry points carry the same ledger across the
optimization/lowering split; reports and the comparison worker expose it.

For unitary replacements, an unnormalized local Frobenius residual bounds the
operator norm, and the triangle inequality bounds composition by their sum in
exact arithmetic. The numerical ledger is **not a formal forward-error
certificate**: libm, analytic lowering, matrix accumulation and the checker
have floating-point limits. Arbitrarily huge arguments remain outside any
claim of uniform floating-point accuracy. Direct users composing separate
synthesis API results must manage their own combined tolerance.

The audit also classified routing/layout scores as search heuristics, and
reviewed KAK thresholds, exact phase recovery, basis translation and verifier
thresholds. The public verifier's documented infidelity thresholds were not
reinterpreted as amplitude residuals. Benchmark validation still uses its
independent, unchanged shared-phase residual oracle. No tolerance was relaxed.

## Regression evidence

- 398 Rust tests pass, plus four separately enabled CLI tests; all-feature
  Clippy and formatting pass.
- 112 Python tests pass with no skips. CI builds the actual
  public comparison worker before running its numerical regressions.
- 51 adversarial cases have independently constructed analytic full operators:
  **28 fail the frozen old compiler; all 51 pass the repaired compiler**.
  They cover positive/negative tiny rotations, U phases, near special angles,
  noncommuting echoes, controlled relative phases, routing and composition.
  Worst repaired analytic residual is `4.60519471196e-11`.
- Additional regressions cover exhausted shared budgets across repeated passes,
  both pipeline entry points, large argument reduction, overflow, serialization
  bit round-trips, unsupported synthesis and integer rational boundaries.

These are numerical and regression results, not a proof of universal compiler
correctness. The independent oracle and tests do not reuse the new production
residual helper to establish expected semantics.

## Existing comparison corpus

Inputs, normalization, target edges, compiler settings and layout contracts
are frozen to the prior primary campaign: 51 fixtures, 123 cases, nine families.
All **3,444 saved rows / 698 distinct artifacts** pass fresh validation without
changed verdicts or metrics (123 full numerical isometries, 575 sampled).
This checks the actual saved Qiskit, Cirq, TKET and old Q-Rust outputs; those
three baselines were not recompiled.

The repaired Q-Rust compiler was run once per case, retaining its actual QASM
and layouts: **123/123 pass**, with 26 full numerical isometry checks and 97
eight-input sampled checks. Worst residual is
`1.22673914489e-14`. Sampling is not an all-input proof.
All layouts and recorded search/optimization settings match the old run.

| Metric | Improved | Equal | Worse |
|---|---:|---:|---:|
| CX count | 0 | 123 | 0 |
| CX-only depth | 0 | 123 | 0 |
| Total depth | 0 | 117 | 6 |
| Total gate count | 1 | 107 | 15 |

QFT-14 on a line remains **274 CX / depth 157**; W-state-14 on a line remains
**26 CX / depth 53**. Every adverse case is retained:

| Case | Previous depth | Repaired depth |
|---|---:|---:|
| `clifford_5__linear` | 74 | 75 |
| `graphstate_10__all_to_all` | 20 | 21 |
| `graphstate_12__all_to_all` | 24 | 25 |
| `graphstate_14__all_to_all` | 28 | 29 |
| `graphstate_5__all_to_all` | 10 | 11 |
| `graphstate_7__all_to_all` | 14 | 15 |

The graph-state overhead consists of roundoff-sized U operations that are now
conservatively retained. For example, graphstate-5 all-to-all keeps five
`U(1.224646799e-16,-1.224646799e-16,0)` operations: 10→15 gates and depth10→11,
with five CX unchanged. This was not a routing regression. No cleanup was
added merely to recover the old figures after inspecting these results.

No shared-budget, block numerical or local-fusion rejection occurred on this
corpus. Maximum charged error was
`2.52822426152e-12`, below `1e-9`.
Output text changed through full precision serialization; semantic checks
compare operators, not QASM string equality.

## Reproduction and provenance

The companion evidence bundle retains both worker binaries, source snapshots,
the previous primary inputs/outputs, final output QASM/layouts, all adversarial
cases, source hashes and logs. Machine-readable results are in
[`numerical_correctness_results_20261017.json`](numerical_correctness_results_20261017.json).
Corrected worker SHA-256: `004ca5c49481ef7eed018bb33686ee0e8f27e6504758ba0b0d804783ee90d047`.
Independent oracle SHA-256: `55de70f9fe47d9fcaeb1497bd04b377e5582afc9289e6948dc6c32ebe5bda6d4`.

Build and run the repository checks:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
cargo test --test cli_smoke_test -- --ignored
cargo build --locked --release -p qrust-bench --bin comparison_runner
QRUST_COMPARISON_BINARY=target/release/comparison_runner python3 -m unittest discover -s experiments -p 'test_*.py'
```

Use the pinned environment in `experiments/requirements-comparison.txt` for
the Python checks. With the prior primary data extracted into `../previous`,
fresh revalidation directories can be generated with:

```sh
python3 experiments/revalidate_numerical_correction.py saved --previous ../previous --output ../saved_recheck
python3 experiments/revalidate_numerical_correction.py corrected --previous ../previous --output ../corrected_recheck --binary target/release/comparison_runner
```

The driver checks fixture, normalization, row-grid, artifact and oracle hashes.
It works from a source archive without Git metadata and retains all compile or
validation failures. One-shot timing telemetry is explicitly diagnostic; it
must not be substituted for repeated controlled timing trials. New speed
claims and external generalization remain unmeasured.
