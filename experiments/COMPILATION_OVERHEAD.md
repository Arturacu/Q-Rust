# Reducing compilation overhead while retaining circuit outputs

The compiler removes repeated routing searches for identical prepared inputs and avoids repeated work in the Clifford pass. All 1,664 formal old/new Q-Rust response pairs have identical native results apart from compiler time, parser time and pass diagnostics: emitted QASM, maps, numerical ledgers, selected statistics and retention results are unchanged. This preserves the circuit improvements in [the preceding quality followup](CLIFFORD_QUALITY_RETENTION.md).

Across the 100 external cases, the geometric mean of per-case compile-time ratios is 0.9394 / 0.9372 in sessions 1 / 2, or 6.1% / 6.3% lower. Across the eight matched cases it is 0.7338 / 0.7373, or 26.6% / 26.3% lower. These are equal-weight relative changes per circuit, not a claim of that much less total time for a batch. The two dense 20-qubit cases dominate absolute matched compilation time and barely improve. Summing one compile median per case gives only 1.5% / 0.8% less matched time. These sums are descriptive calculations, not measured batch wall times.

Speed ratios compare the configured compiler pipelines on this arm64 macOS host. They do not causally measure Rust versus Python or establish universal compiler rankings.

## Paired results and unchanged circuit costs

Ratios below are candidate / previous. Each session uses the same stable old/new verified cohort; every timed result and warmup must pass. Sessions are not pooled.

| Cohort | Session | Compile geometric ratio | Request-wall geometric ratio | Sum of compile medians, ms | Ratio of those sums |
|---|---|---|---|---|---|
| External, 100 cases | 1 | 0.9394 | 0.9411 | 955.722 → 697.672 | 0.7300 |
| External, 100 cases | 2 | 0.9372 | 0.9382 | 956.645 → 697.777 | 0.7294 |
| Matched, 8 cases | 1 | 0.7338 | 0.7400 | 26024.434 → 25631.178 | 0.9849 |
| Matched, 8 cases | 2 | 0.7373 | 0.7422 | 25489.915 → 25288.838 | 0.9921 |

A geometric mean gives each relative case change equal weight. A ratio of sums weights long cases by their absolute durations. Neither includes cold startup or measures the wall time of compiling a batch. The immediate previous repository baseline verifies all eight matched cases; the six-case old/new cohort in the earlier report used a different baseline and must not be combined with these ratios.

Every CX count, total gate count, total depth and two-qubit depth is unchanged in the old/new formal comparisons. The retained anchors are:

| Circuit / target | CX | Total depth | Two-qubit depth | Total gates |
|---|---|---|---|---|
| EC-5 / complete | 9 | 11 | 6 | 26 |
| EC-5 / line | 17 | 19 | 12 | 39 |
| QEC encoder-5 / complete | 8 | 10 | 8 | 14 |
| QEC encoder-5 / line | 20 | 21 | 16 | 32 |

The earlier QEC encoder trade-off remains: its line output has 21 total-depth layers and 16 two-qubit-depth layers. This runtime change does not alter that choice or remove its previously documented extra total-depth layer.

The following per-case medians retain small changes and slowdowns. Complete request wall time includes parsing, serialization, diagnostics and IPC as well as compilation. It is measured using the existing warm workers, not a cold command invocation.

| Case | Compile ms, S1 previous → candidate | Compile ms, S2 previous → candidate | Compile ratio, S1 / S2 | Request-wall ratio, S1 / S2 |
|---|---|---|---|---|
| EC-5 complete | 1.538 → 1.315 | 1.618 → 1.377 | 0.8550 / 0.8506 | 0.8746 / 0.8644 |
| EC-5 line | 3.373 → 3.293 | 3.319 → 3.267 | 0.9765 / 0.9845 | 0.9859 / 0.9855 |
| QEC encoder-5 complete | 0.417 → 0.399 | 0.403 → 0.383 | 0.9549 / 0.9498 | 0.9614 / 0.9661 |
| QEC encoder-5 line | 2.397 → 2.384 | 2.424 → 2.319 | 0.9946 / 0.9568 | 1.0082 / 0.9562 |
| AE-10 line | 51.071 → 30.295 | 51.378 → 30.376 | 0.5932 / 0.5912 | 0.5958 / 0.5939 |
| QPE-10 line | 25.256 → 25.296 | 25.177 → 25.108 | 1.0016 / 0.9973 | 1.0014 / 0.9919 |
| QPE-5 complete | 1.064 → 1.087 | 1.078 → 1.067 | 1.0220 / 0.9893 | 1.0192 / 0.9962 |
| Trotter-4 line | 152.794 → 125.732 | 154.653 → 127.287 | 0.8229 / 0.8230 | 0.8247 / 0.8248 |
| bv100 heavy-hex | 174.115 → 89.635 | 176.361 → 90.503 | 0.5148 / 0.5132 | 0.5158 / 0.5142 |
| bv100 square | 149.785 → 78.484 | 149.673 → 78.877 | 0.5240 / 0.5270 | 0.5257 / 0.5286 |
| cat22 heavy-hex | 4.154 → 4.276 | 4.192 → 4.182 | 1.0295 / 0.9977 | 1.0362 / 0.9942 |
| cat22 square | 2.808 → 1.567 | 2.780 → 1.615 | 0.5582 / 0.5809 | 0.5828 / 0.5994 |
| clifford20 heavy-hex | 18464.055 → 18351.117 | 17971.927 → 18016.760 | 0.9939 / 1.0025 | 0.9939 / 1.0025 |
| clifford20 square | 7220.103 → 7098.920 | 7175.574 → 7089.648 | 0.9832 / 0.9880 | 0.9832 / 0.9880 |
| ghz23 heavy-hex | 4.439 → 4.379 | 4.565 → 4.471 | 0.9867 / 0.9794 | 0.9870 / 0.9774 |
| ghz23 square | 4.975 → 2.799 | 4.843 → 2.782 | 0.5626 / 0.5744 | 0.5695 / 0.5870 |

External compile medians are higher on 54/100 cases in session 1 and 46/100 in session 2; complete request-wall medians are higher on 52 and 43 cases, with two ties in session 2. The largest external observed compile increases are 7.0% and 7.5% on different small cases. Matched compile/request-wall medians are higher on one case in each session: cat22 heavy-hex in session 1 and dense Clifford20 heavy-hex in session 2. These counts retain measurement variability and do not assert statistical significance for small per-case shifts. Every case, including slowdowns, is retained in `analysis/reporting/case_timings.json`.

## Why work was removed, and what still costs time

The previous quality-preserving pipeline compiled retained alternatives independently, sometimes repeating layout and routing on exactly the same prepared circuit. The new reuse key compares complete prepared circuit data and the routing context, including exact parameter representations. A compilation keeps at most four routing checkpoints. A reused checkpoint supplies routed output, maps and routing alternatives, while the new branch retains its own prefix metadata, logical permutation and numerical ledger. Each distinct lowering still runs and validates independently. Immutable checkpoint outputs are cloned into branch-owned state; this removes repeated searches, not all allocations. The new routing store exists only within one compilation. Existing Clifford lookup caches already persist in a warm worker and are present in both versions.

The Clifford pass reuses scratch storage for transported signed-Pauli data across endpoint searches, scores the same cost while projecting extracted wire permutations directly, and folds the edited operation stream without an intermediate deep clone. The transport data are cleared and recomputed; sharing that information across cuts was already present in the previous compiler. Final immutable candidate cost is computed lazily once and reused without changing selection ordering. Search bounds, candidate order, fences and acceptance semantics remain unchanged. Regression comparisons exercise the old and new scoring/folding implementations, including cuts, metadata and unsupported-operation fences.

The bounded diagnostic study uses 14 cases, three repetitions per profile mode and one warmup per mode: 112 requests each for the previous compiler, a pass-only variant and the final compiler. These runs occur separately from formal timing. All final diagnostic outputs verify, and all 112 native responses match the previous responses under the same three diagnostic exclusions. The traces show layout/routing call counts changing from 2 to 1 for AE-10 line, cat/GHZ square and both BV100 targets; from 4 to 2 for Trotter-4 line; and from 4 to 3 for EC complete. Explicit `ReuseRoutingCheckpoint` records replace the removed searches. EC line and both QEC cases retain two genuinely distinct routing calls, so they retain substantial route/lowering cost. The dense-20 cases were not included in this targeted profile study; their formal runtimes remain the dominant matched cost.

Summed diagnostic Clifford-pass medians fall from 9.10 to 3.55 ms on BV100 square, 21.69 to 7.70 ms on BV100 heavy-hex, and 6.85 to 2.58 ms on Trotter-4 line. Profile sums include retained and alternate paths and do not partition total compilation time. The separate diagnostic ablation timings are not controlled end-to-end causal estimates; the formal adjacent old/new requests above measure the combined change.

## Contemporaneous SDK comparisons

The unchanged external primary and light comparison cohorts each contain all 100 cases. Compiler/request speed ratios below are SDK / Q-Rust; output-cost ratios are Q-Rust / SDK. Values are geometric means of per-case medians. Both sessions have identical displayed quality ratios.

| Recipe | Compile ratio, S1 / S2 | Request-wall ratio, S1 / S2 | CX ratio | Depth ratio | Two-qubit-depth ratio |
|---|---|---|---|---|---|
| Qiskit | 2.724 / 2.723 | 2.686 / 2.694 | 0.9150 | 0.8690 | 0.8983 |
| Cirq | 63.335 / 62.837 | 61.296 / 60.994 | 0.7147 | 0.6007 | 0.6932 |
| TKET | 168.635 / 168.988 | 139.068 / 139.447 | 0.8951 | 0.9504 | 0.8986 |
| Cirq, one pass | 38.111 / 37.948 | 39.507 / 39.439 | 0.7053 | 0.5976 | 0.6925 |
| TKET, light | 20.378 / 20.558 | 18.578 / 18.704 | 0.6564 | 0.7212 | 0.6663 |

Two zero-CX/two-qubit-depth cases are excluded from the corresponding external geometric ratios and counted explicitly; no pseudocounts are used. All 100 cases remain in timing and total-depth ratios. Aggregate cost ratios do not imply Q-Rust wins every circuit: for example, Qiskit has fewer CX on seven external cases and TKET on six. The unchanged 2,000-draw external bootstrap analysis and its session-agreement results are retained in full.

The stable four-tool matched cohort still contains only the four cat/GHZ cases. It must not be reported as a verified four-tool comparison of all eight cases.

| SDK, four cat/GHZ cases | Compile ratio, S1 / S2 | Request-wall ratio, S1 / S2 | CX ratio | Depth ratio | Two-qubit-depth ratio |
|---|---|---|---|---|---|
| Qiskit | 4.491 / 4.473 | 4.431 / 4.417 | 1.8407 | 1.7637 | 1.7974 |
| Cirq | 28.493 / 28.537 | 29.031 / 29.126 | 1.8407 | 0.9019 | 1.7974 |
| TKET | 111.358 / 110.193 | 107.257 / 106.352 | 1.8407 | 1.7637 | 1.7974 |

Q-Rust uses more CX on all four matched common-cohort cases. Its external aggregate quality advantage does not extend to this cohort. On EC complete it retains 9 CX / 11 depth / 6 two-qubit depth, versus TKET's 9 / 13 / 7; on EC line it retains 17 / 19 / 12, versus TKET's 20 / 25 / 19. Complete per-SDK EC/QEC tables, including compilation time and all gate counts, are retained in `report_metrics.json`.

## Frozen protocol, correctness and limits

The 100-case QASMBench/MQT Bench suite retains six recipes, two sessions, seven timed trials and one initial warmup per cell. The matched extension retains eight circuits/targets, four tools, two sessions, three timed trials and one initial warmup. A previous Q-Rust request runs next to every candidate request; case/trial parity balances the order and session 2 reverses it. The frozen schedule contains 10,040 timed slots and 1,480 initial warmups. Six timed Cirq slots were blocked after failed warmups; scheduled slots are not all successful compilations. Sessions ran sequentially: external 1, external 2, matched 1, matched 2. No agent build, profiling, simulation, validation, compression or other benchmark session overlapped their timing.

The input, normalization, adapter, oracle, schedule, timeout and SDK recipe contracts are unchanged, including Qiskit 2.5.2, Cirq 1.7.0 and pytket 2.18.5. External requests have a 60-second timeout. Matched requests have a 30-second timeout and approximately 100 ms RSS observations enforcing 2 GiB per worker and 3 GiB across workers, including the previous Q-Rust worker. This is sampled limit enforcement, not exact peak memory, descendant accounting or evidence of a peak-memory improvement.

All 11,200 external records verify. In the matched extension, both Q-Rust binaries and TKET verify all eight cases. Qiskit has 32 unsupported-recognition records; Cirq has 24 unsupported-recognition records, two warmup timeouts and six blocked timed slots. These are retained, and unsupported recognition is not called incorrect. No definite semantic counterexample or structural violation was found. Every eligible cell requires all timed outputs and all initial/replacement warmups to pass.

The unchanged external oracle checks complete clean-ancilla isometries through six physical qubits, and eight independently selected inputs above that within its existing limit. Sampling does not prove arbitrary-input equivalence. The matched oracle checks map/ancilla-aware signed-Clifford Choi isometries after bounded numerical recognition of local Clifford matrices; unsupported recognition fails closed. Neither checker has been relaxed. Exact native-result parity is additional evidence: after checking output QASM hashes and binding row/session/schedule identities, all 1,448 timed and 216 warmup pairs match except `compile_ns`, `parse_ns` and `pass_diagnostics`. Scheduler/transport observations are retained separately and are not native compiler fields.

Before timing, the strict 100+8 preflight verifies both binaries and finds exact output parity and unchanged structural costs. Final independent fuzzing verifies 384 compiled outputs from 64 fixtures on two topologies, rejects 128 relative-phase negative controls, checks 91 retention decisions and confirms 128 additional old/new native result pairs. Final checks pass 494 Rust tests (four existing tests intentionally ignored), strict Clippy, formatting and all 249 Python tests with the frozen binary and no skips. Five strict-parity tests, eight formal-pair identity/parity tests, 14 canonical-replay fault-injection tests and nine packaging boundary tests pass.

## Provenance and evidence

The previous repository baseline is `ef9c3552752b819440fad54882f56a77582ba46a`; its timed compiler commit is `65e6c64a66d68da79799c9ea6a344553cb844cc4`, binary SHA-256 `3ae3dcdc258e8b96b41d552964c3863702e17e4112c6e864f13376d33efaf891`. The candidate compiler commit is `3306937c0ae71fb2eb78ad86141f2ccd1b636465`, binary SHA-256 `aa6513a05aa814b8b89cf373472cacaef4a0afa7c66a296d76a67ef4f4824e1e`. The frozen protocol binds 79 Rust/Cargo source files, all adapter/oracle hashes and both binaries. The preservation manifest binds 1,192 earlier files; earlier reports, raw records, protocols and archives remain byte-identical.

The archive retains the baseline/pass-only/final profiles, source/binary provenance, strict preflight, final fuzzing, formal responses, dispatch records and all diagnostics. The first baseline profile attempt omitted the wrapper's required QASM checksum and stopped in instrumentation after compilation; its records and log are preserved, and the corrected wrapper reran the entire study. That failed diagnostic attempt is not a compiler correctness failure or a formal timing result.

The result package is [external_results/20261018_compile_overhead](external_results/20261018_compile_overhead/), with [the evidence archive](external_results/20261018_compile_overhead/evidence.tar.gz) and [checksums](external_results/20261018_compile_overhead/checksums.json). The final campaign is `campaign_release/` locally and `campaign/` in the archive. Saved analysis and a fresh full-output replay agree after normalizing only newly measured `validation.validation_ns` values, preserving field presence and adjusting only their derived row hashes. Compiler/parser/request timings and all remaining evidence remain exact. Both reporting derivations reproduce byte-identical JSON; their receipts are in `analysis/reporting/`.

With the recorded Python dependencies, replay into new paths after extracting the archive:

```sh
python harness/experiments/paired_analysis.py --repo harness --campaign campaign --output replay_analysis
python scripts/canonical_replay.py --reference analysis --replay replay_analysis --output replay_receipt.json
python analysis/reporting/compute_metrics.py --analysis analysis --campaign campaign --output report_metrics_replay.json
python analysis/reporting/derive_case_timings.py --analysis analysis --output case_timings_replay.json
python verification/verify_formal_output_parity.py --campaign campaign --completion receipts/timing_completion.json --allow-relocated --output native_parity_replay.json
```

This reruns the existing matched Q-Rust extension. The earlier native Benchpress workout results remain unchanged; there is no native Q-Rust gym. The bounded Clifford workload, ideal targets and structural metrics do not establish noisy-device execution quality or general large non-Clifford correctness. The [Benchpress pilot](BENCHPRESS_PILOT.md) documents the separate native runs and their different correctness contract.
