# Dense Clifford20: measured causes and limits

Q-Rust's long compilation time on this development fixture is principally routing work. Beam routing accounts for a median 83.69% of the profiled square request and 89.89% of the heavy-hex request. The frozen pipeline routes two prepared continuations and pays for a substantial progress-ranking search policy. Controlled toggles approximately halve runtime, but their quality effects differ between the two targets. They are diagnostic interventions, not changes to the frozen production compiler.

This is one existing Benchpress Clifford circuit, `clifford_20_12345.qasm`, with 20 logical qubits and 4,000 input operations, on two ideal bidirectional coupling graphs: square has 25 physical vertices/40 edges and heavy-hex has 57/64. It is not an independent sample of Clifford circuits or an equal-time optimization tournament. The separate held-out study is the evidence for generalization beyond this development input.

## Freeze, execution and correctness

The production binary is unchanged: SHA-256 `aa6513a05aa814b8b89cf373472cacaef4a0afa7c66a296d76a67ef4f4824e1e`, compiler commit `3306937c0ae71fb2eb78ad86141f2ccd1b636465`. Source references use repository head `3327d6cdcd937f3a3786c756a0279cf99e522239`. The environment contains Qiskit 2.5.2, Cirq 1.7.0, TKET 2.18.5, NumPy 2.5.4 and SciPy 1.18.1; complete source, case, environment and script hashes are in each frozen protocol.

Each case/recipe has one warmup and three measured repetitions, with fixed seed 20261018 and alternating recipe order. All compilation requests ran sequentially in an exclusive diagnostic timing slot. SDK imports and input parsing are outside the compile clocks; serialization and independent validation follow compilation. Instrumentation can add overhead, so profile clocks describe these diagnostic runs, not a replacement production speed estimate. Repetitions estimate execution variability at one seed, not routing variability across seeds.

| Block | Scheduled outcomes, including warmups | Wall time | Independent output checks |
|---|---:|---:|---|
| Native/profile parity baseline | 16 | 198.983 s | 16 passed; 2 unique artifacts |
| Q-Rust controls | 48 | 487.678 s | 48 passed; 11 unique artifacts |
| SDK controls v2 | 64 | 1074.062 s | 24 passed, 32 unsupported recognition, 8 timeouts; 11 unique compiled artifacts |

The frozen map-aware Clifford Choi-isometry checker checks all logical inputs for the recognized Clifford outputs, including the ancilla/output-map contract, rather than only the zero state. Basis and coupling constraints are checked too. Local gate recognition is numerical; its tolerances and the complete validation records are retained. Q-Rust and TKET outputs passed. The Qiskit and compiled Cirq outputs contain local U3 rotations that this gate-by-gate Clifford recognizer cannot classify, despite a Clifford source. A decomposition can contain non-Clifford local gates whose combined operation is Clifford. These outputs remain **unverified**, not demonstrated incorrect. No tolerance or oracle was changed to admit them, and they are excluded from correctness-qualified timing/quality comparisons.

All eight baseline native responses per target match exactly after excluding only native compile/parse clocks and pass diagnostics, plus controller wall-time/cold-worker fields. Available SDK native/instrumented pairs have identical QASM, output width and maps for every repetition, including warmup. Heavy-hex Cirq originals timed out in both forms, so that pair is unavailable rather than counted as a parity pass. Artifact parity establishes that instrumentation did not change the available output; it does not independently prove its correctness.

The original, unexecuted `sdk_ablation` plan is preserved. Before SDK execution, input eligibility inspection found 402 CY and 402 CZ operations unsupported by the restricted native-CX diagnostic helper. The existing general Cirq one-pass control replaced that arm in `sdk_ablation_v2`, retaining 64 slots and the 30-second request cap. `sdk_protocol_revision.json` records this source-driven, pre-execution substitution. No unsupported attempt is implied to have run.

## Q-Rust: where time goes

The standalone baseline native medians are 6.859 s on square and 18.167 s on heavy-hex. Profiled medians are 6.791 s and 18.075 s, with identical outputs. For every request, exclusive trace events are summed by group, divided by that request's compile time, and only then summarized by the median. These group medians must not be added and interpreted as an exact partition of a median runtime.

| Exclusive group | Square share | Heavy-hex share |
|---|---:|---:|
| Beam routing | 83.69% | 89.89% |
| Initial layout search | 6.43% | 3.43% |
| Local Clifford simplification | 1.29% | 0.71% |
| Two-qubit block synthesis | 0.97% | 0.71% |

All recorded events cover 99.81–99.83% of square request time and 99.88–99.89% of heavy-hex request time. The two-qubit synthesis and local Clifford passes themselves are small fractions. The expensive consequence of Clifford rewriting is that the frozen quality guard also routes and lowers the retained original continuation. This distinction matters when attributing cost.

The final selection uses CX count, total depth and two-qubit depth, not only a SWAP proxy. Square selects the rewritten continuation at `[4416,2464,2017]` over retained `[4491,2795,2223]`. Heavy-hex selects retained `[7898,3612,3138]` over rewritten `[8050,3600,3134]`: the rewritten candidate has slightly lower depths but 152 more CX. The guard therefore rejects it. These observations explain this implementation's work and selection; they do not prove two complete continuations are the only possible way to enforce a quality constraint.

## Q-Rust controlled interventions

All rows below passed the unchanged checker. Runtime is the median of three measurements in the same ablation block. The separate baseline block above is not substituted for the within-block reference.

| Target | Diagnostic recipe | Time (s) | Time/reference | CX | Depth | Two-qubit depth |
|---|---|---:|---:|---:|---:|---:|
| Square | Frozen O3, profiled | 7.286 | 1.000 | 4416 | 2464 | 2017 |
| Square | Disable Clifford rewrite/retained comparison | 3.611 | 0.496 | 4491 | 2795 | 2223 |
| Square | Disable progress ranking | 3.737 | 0.513 | 4638 | 2323 | 1904 |
| Square | One layout trial/refinement | 6.711 | 0.921 | 4463 | 2458 | 2010 |
| Square | Lookahead 5 instead of 20 | 7.170 | 0.984 | 4619 | 2585 | 2135 |
| Square | Disable block synthesis/dependent fusion | 6.913 | 0.949 | 4451 | 2476 | 2039 |
| Heavy-hex | Frozen O3, profiled | 17.776 | 1.000 | 7898 | 3612 | 3138 |
| Heavy-hex | Disable Clifford rewrite/retained comparison | 8.565 | 0.482 | 7898 | 3612 | 3138 |
| Heavy-hex | Disable progress ranking | 8.898 | 0.501 | 7664 | 3252 | 2878 |
| Heavy-hex | One layout trial/refinement | 15.844 | 0.891 | 8076 | 3673 | 3221 |
| Heavy-hex | Lookahead 5 instead of 20 | 16.780 | 0.944 | 7787 | 3398 | 2957 |
| Heavy-hex | Disable block synthesis/dependent fusion | 17.718 | 0.997 | 7959 | 3633 | 3152 |

Disabling Clifford rewriting removes both its local work and the alternate retained-input comparison. It reduces square time by 50.4%, but costs 75 CX and 331 depth. On heavy-hex it reduces time by 51.8% and preserves the exact QASM and maps, because the default selected the same retained original. This demonstrates substantial avoidable work for this particular heavy-hex input if one knew which branch to skip; it does not supply a general safe branch-prediction rule.

Disabling progress ranking reduces square time by 48.7% and heavy-hex time by 49.9%. Heavy-hex improves all three structural costs: 234 fewer CX, 360 less depth and 260 less two-qubit depth. Square instead has 222 more CX while improving both depths. Progress ranking therefore does not dominate the simpler policy on these cases. The toggle changes free-progress execution, ranking/deduplication and subsequent search decisions together; it does not isolate memory allocation cost or prove the same change is beneficial on other circuits.

Reducing layout budget saves about 7.9%/10.9% but increases CX on both targets. Shorter lookahead improves heavy-hex cost and time but worsens all three square costs. These are heuristic search outcomes, so more search or farther lookahead does not monotonically improve every final circuit. Disabling blocks saves only 0.3% on heavy-hex and 5.1% on square in these samples, while worsening CX and both depths on both targets. The direct block-pass shares are below 1%; the whole-pipeline timing change is not the isolated direct cost of that pass, and small deltas are not a reliable general speed result.

The source explanation is in `SOURCE_NOTES.md`: O3 uses beam width 8, conditional branching up to 5, four forward and three backward searches per ordinary route, candidate-state copies and ranking, and up to an additional refined-layout forward search. Candidate scoring scans front/lookahead interactions against a distance table precomputed once per route. It does not run a shortest-path search for every candidate. Child copying/allocation may contribute, but no allocation-level profile measured its share. The larger physical graph and different routing trajectory change work together; these two cases do not isolate a graph-size scaling law.

## TKET: synthesis budget, not routing, dominates this recipe

The frozen TKET recipe applies `FullPeepholeOptimise(allow_swaps=True, target_2qb_gate=CX)` before mapping, followed by the existing mapping, SWAP decomposition, synthesis and U3/CX rebase. The scoped instrumented version applies the identical passes to the same CompilationUnit and reproduces the native artifact/maps exactly. Its initial peephole pass occupies a median 95.53% of square compile time and 93.93% of heavy-hex time; mapping occupies 1.77% and 3.05%. The dominant TKET work therefore differs from Q-Rust's routing-heavy profile on this input.

| Target | TKET recipe | Time (s) | CX | Depth | Two-qubit depth |
|---|---|---:|---:|---:|---:|
| Square | Frozen native | 22.255 | 4813 | 2271 | 2025 |
| Square | Same recipe, stage clocks | 22.490 | 4813 | 2271 | 2025 |
| Square | Replace initial peephole pass with SynthesiseTket | 3.086 | 7870 | 3566 | 3152 |
| Heavy-hex | Frozen native | 23.180 | 8723 | 3478 | 3248 |
| Heavy-hex | Same recipe, stage clocks | 22.842 | 8723 | 3478 | 3248 |
| Heavy-hex | Replace initial peephole pass with SynthesiseTket | 3.554 | 13256 | 5427 | 5021 |

The replacement cuts the observed runtime to 13.9%/15.3% of the frozen native recipe, but adds 3057/4533 CX and 1295/1949 depth. All outputs passed. This is strong evidence for a real synthesis-budget/output-quality trade-off in the selected TKET recipe, rather than an inherent language speed ranking. It does not reveal which internal subroutine of FullPeepholeOptimise consumes that time. Production settings were not changed.

Q-Rust's verified native baseline uses fewer CX than this verified TKET recipe on both targets, while TKET has lower total depth on both; neither dominates all reported structural metrics. Timings are separate diagnostic blocks, so their cross-tool ratio is not promoted to a new formal benchmark score.

## Qiskit and Cirq: informative execution traces, incomplete equivalence evidence

The following clocks and costs are retained observations only. Their compiled outputs were not admitted by the frozen local-Clifford recognizer, so these rows cannot establish a correctness-qualified winner. In particular, their reported CX/depth counts do not overcome the missing equivalence evidence.

| Target | Diagnostic recipe | Observed time (s) | Observed CX / depth | Status |
|---|---|---:|---:|---|
| Square | Qiskit frozen native | 0.365 | 4817 / 2357 | Unsupported recognition |
| Square | Qiskit same recipe, callbacks | 0.371 | 4817 / 2357 | Unsupported recognition; native parity passed |
| Heavy-hex | Qiskit frozen native | 0.514 | 7917 / 3227 | Unsupported recognition |
| Heavy-hex | Qiskit same recipe, callbacks | 0.534 | 7917 / 3227 | Unsupported recognition; native parity passed |
| Square | Cirq frozen native | 24.745 | 7917 / 6068 | Unsupported recognition |
| Square | Cirq same recipe, stage clocks | 24.514 | 7917 / 6068 | Unsupported recognition; native parity passed |
| Square | Cirq one synthesis pass per stage | 12.532 | 7835 / 5428 | Unsupported recognition |
| Heavy-hex | Cirq frozen native and instrumented originals | >30 s request cap each | Unavailable | All 8 requests timed out |
| Heavy-hex | Cirq one synthesis pass per stage | 26.639 | 21585 / 15615 | Unsupported recognition |

Qiskit's frozen level-3 recipe uses a seeded SABRE-based layout/routing flow and exact approximation setting, as specified by the saved adapter. Callback traces locate most reported pass time in SabreLayout (median pass elapsed 0.320 s square/0.446 s heavy-hex), versus 0.371/0.534 s whole profiled compile clocks. The trace locates the work; it does not by itself isolate why Qiskit's search is faster than Q-Rust's different beam policy. No matched internal-search-budget or representation experiment was run for Qiskit here. Callback count/depth collection is outside the reported pass duration but inside the whole profiled clock.

Cirq's selected adapter performs target-gateset synthesis before and after RouteCQC, then closes the representation to CX and fuses one-qubit gates. On square, postrouting synthesis takes a median 66.54% of the instrumented request, prerouting synthesis 21.07%, and graph/routing 1.33%. Reducing synthesis from two passes to one at each stage produces a 12.532 s observed median and changes the output. On heavy-hex the one-pass recipe completes while both originals exceed the cap; its postrouting synthesis still takes 60.77% and CX closure/fusion 27.43% of the request. The timeout is a request-wall cap, not a measured 30-second native compile duration. Neither the apparent timing improvement nor changed structural counts resolves the unsupported correctness check.

These controls explain why the selected TKET and Cirq recipes can spend time in synthesis while Q-Rust spends it in routing. They do not establish all causes of every inter-compiler gap, or attribute the gaps to implementation language. A separately specified, representation-independent equivalence method is needed before these dense Qiskit/Cirq outputs can join correctness-qualified comparisons. That would require a new validation study, not relabeling the preserved unsupported outcomes.

## Audit and retained evidence

`causal_findings.json` binds all three studies' protocols, receipts, raw/validated rows and summaries, along with the derivation and retention-interpretation helper hashes. `baseline_accounting.json`, `qrust_ablation_accounting.json` and `sdk_ablation_accounting.json` retain per-request exclusive accounting. Every failure, timeout, warmup and cold restart remains in the raw ledger. Cold-after-failure, missing-repeat and unsupported cells cannot enter qualified timing comparisons. No failed request was silently retried.

All 13 evidence regression tests passed. They cover disjoint timing groups, per-request fraction calculation, strict native parity rejection for changed output/maps/counters, timing admission rules, protocol-receipt binding, and independent checking of final retention-cost interpretation. These are evidence-integrity checks; the production compiler tests and production code were not changed in this frozen investigation.

The results identify important causal mechanisms and useful future hypotheses, with clear counterexamples to universal toggle benefits. They do not prove optimal output quality, general scaling, noise-aware hardware performance, a Rust-specific speed advantage, or a complete explanation of all compiler behavior.

[Complete pinned source inspection and baseline accounting](./dense_source_notes.md).
