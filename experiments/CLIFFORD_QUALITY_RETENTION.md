# Clifford placement and final-output quality retention

The repeated comparison reduces error correction from 11 to 9 CX on the complete graph and from 27 to 17 CX on the line. Q-Rust now matches TKET's complete-graph CX count and uses fewer CX on the line under these recipes. The QEC encoder improves its two-qubit depth from 19 to 16, while total depth increases from 20 to 21 and total gate count from 26 to 32. Compilation is more expensive: across the 100 external cases, the geometric mean of candidate/previous compile-time ratios is 1.239 in session 1 and 1.236 in session 2.

All metrics below come from saved, independently validated outputs. A triple means **CX count / total depth / two-qubit depth** under the unchanged analyzer. The structural results in this table are identical in both sessions.

| Circuit and target | Previous Q-Rust | New Q-Rust | Qiskit | Cirq | TKET |
|---|---:|---:|---:|---:|---:|
| Error correction, complete graph | 11 / 14 / 8 | **9 / 11 / 6** | 35 / 64 / 34 | 35 / 69 / 34 | 9 / 13 / 7 |
| Error correction, line | 27 / 30 / 24 | **17 / 19 / 12** | 51 / 72 / 43 | 52 / 91 / 45 | 20 / 25 / 19 |
| QEC encoder, complete graph | 8 / 10 / 8 | 8 / 10 / 8 | 10 / 14 / 10 | 10 / 20 / 10 | 8 / 10 / 8 |
| QEC encoder, line | 20 / 20 / 19 | **20 / 21 / 16** | 18 / 24 / 17 | 22 / 44 / 22 | 20 / 21 / 20 |

The line QEC result is a trade-off, not improvement on every metric. Qiskit uses two fewer CX there, and TKET uses six fewer total gates than new Q-Rust. The old error-correction line two-qubit depth is 24 under the independent analyzer; a compiler-internal diagnostic value of 22 measures something different and is not substituted here. Nine CX on one fixture is an observed result, not a proof of global optimality.

## Compilation cost

Each entry below is a session median in milliseconds. Old/new values are contemporaneous, adjacent runs with balanced order across sessions.

| Circuit and target | Previous ms, S1 / S2 | New ms, S1 / S2 | New / previous, S1 / S2 |
|---|---:|---:|---:|
| Error correction, complete graph | 0.823 / 0.804 | 1.615 / 1.618 | 1.962 / 2.012 |
| Error correction, line | 2.016 / 2.034 | 3.328 / 3.389 | 1.650 / 1.666 |
| QEC encoder, complete graph | 0.213 / 0.207 | 0.405 / 0.404 | 1.900 / 1.951 |
| QEC encoder, line | 0.958 / 1.063 | 2.373 / 2.597 | 2.476 / 2.444 |

The additional compilation cost is material even where the output is unchanged, as in complete-graph QEC. Error-correction complete-graph time increases by 96.2% / 101.2%; the line increases by 65.0% / 66.6%. QEC line increases by 147.6% / 144.4%.

| Stable verified old/new cohort | Cases | Compile-time new / old, S1 / S2 | Request-wall new / old, S1 / S2 |
|---|---:|---:|---:|
| External suite | 100 | 1.2393 / 1.2363 | 1.2157 / 1.2073 |
| Matched extension | 6 | 1.6044 / 1.6310 | 1.5893 / 1.6169 |

These are geometric means of per-case median ratios, with sessions kept separate. The external compile-time increase is 23.9% / 23.6%; the six-case matched increase is 60.4% / 63.1%. They compare these configured pipelines on this host. They do not isolate Rust versus Python or establish a causal language-speed effect. Compiler time excludes startup, parsing and serialization; request wall time is retained separately.

## Mechanisms

The remaining complete-graph error-correction gap depended on where a legal Clifford rewrite was inserted. The previous pass transported one interaction toward a later interaction and attempted replacement there. Two interactions can instead move toward a shared interior position when each transport is legal. The added search considers these reachable positions using signed-Pauli conjugation, without permitting either interaction to cross a noncommuting operation. Unsupported operations and barriers remain fences. Each pass invocation is bounded to two sweeps, a 96-operation lookahead and 96 full-circuit candidate evaluations; cumulative counters across invocations can exceed 96.

The line-target result also depends on the pipeline representation. Explicit routing SWAPs fence the Clifford pass. After target-basis translation expands them into CX sequences, a physical-wire Clifford pass can expose further cancellations while named local Clifford gates remain available. The promoted placement is after target-basis translation and before general basis decomposition, with SWAP replacements disabled. Subsequent basis closure and direction repair still validate the hardware output. The controlled logical interior-position search reached 21 CX / depth 22 / two-qubit depth 17 on the line; the expanded-routing pass reached 17 / 19 / 12. All 32 input columns passed the independent equivalence check for that control.

The QEC encoder regression had a different cause: a local logical improvement could produce a worse routed circuit. The compiler now retains an independently routed and lowered continuation with cross-block Clifford simplification disabled. When this continuation succeeds, every rewritten candidate must be no worse than that fixed reference in two-qubit gate count, total depth and two-qubit depth. Among eligible candidates it selects lexicographically by those costs. Comparison occurs after direction repair, basis closure and terminal readout restoration.

This guard applies to level 3 compilation with a backend, basis decomposition, output-permutation elision and an actual Clifford rewrite in any candidate. An optional continuation failure is recorded and does not certify nonregression. The reference is the same compiler with the Clifford pass disabled; it is not a guarantee that every metric improves relative to the previous repository baseline. Each continuation owns its numerical ledger and layout/permutation metadata; rejected branches do not spend the selected branch's allowance. The local search's internal bounded CX-family cost is not the final routed hardware cost.

Numerical two-qubit synthesis can express an overall Clifford operation using non-Clifford local factors. That representation prevents the existing efficient Clifford oracle from recognizing the output even when it is mathematically equivalent. If the entire synthesis-pass input is already gatewise Clifford, the compiler now vetoes a synthesized block that leaves that representation, retaining the original block and its unspent numerical budget. Recognition only vetoes a replacement; it never snaps angles or changes accepted gates. Mixed non-Clifford routines retain their previous synthesis policy. The oracle and its tolerances are unchanged.

## External-suite comparison

All 100 cases passed for all six recipes and for the previous Q-Rust control in both sessions, including every warmup and timed output. The cohort is the same 50 external fixtures on complete and line targets. It is not an enlarged or reweighted workload set.

The speed ratios below are SDK compile time divided by new Q-Rust compile time. Values above one favor Q-Rust. Quality ratios are new Q-Rust cost divided by SDK cost; values below one favor Q-Rust. Each ratio is a geometric mean of per-case median ratios. Quality ratios are identical in both sessions.

| Configured SDK recipe | Compile-time ratio, S1 / S2 | CX ratio | Total-depth ratio | Two-qubit-depth ratio | Gate-count ratio |
|---|---:|---:|---:|---:|---:|
| Qiskit, primary | 2.525 / 2.546 | 0.9150 | 0.8690 | 0.8983 | 0.8756 |
| Cirq, primary | 58.532 / 58.434 | 0.7147 | 0.6007 | 0.6932 | 0.5743 |
| TKET, primary | 157.512 / 158.160 | 0.8951 | 0.9504 | 0.8986 | 1.0046 |
| Cirq, one pass | 35.171 / 35.126 | 0.7053 | 0.5976 | 0.6925 | 0.5619 |
| TKET, light | 19.203 / 18.923 | 0.6564 | 0.7212 | 0.6663 | 0.7699 |

All 100 cases enter timing, total-depth and gate-count ratios. The two cases with zero CX and zero two-qubit depth in both compared outputs are counted explicitly and excluded from those geometric ratios, leaving 98 positive ratios; no pseudocount is added. Aggregate advantages do not imply a per-case advantage. For example, new Q-Rust has more CX than Qiskit on seven cases. TKET's primary and light recipes produce materially different speed/quality trade-offs, so the recipe names are retained.

Only five external cases change structurally relative to the previous Q-Rust control: the two error-correction targets, line QEC, and W-state-3 on each graph. W-state depth improves from 13 to 12 on the complete graph and 15 to 14 on the line, with CX unchanged at 6 and 8. The other 95 cases have unchanged CX, gate count, total depth and two-qubit depth, including QPE-5 and Trotter-4 on both targets. CX does not increase on any external case; the one total-depth increase is line QEC.

The original external analyzer was rerun with 2,000 bootstrap draws. Its uncertainty estimates and recipe-sensitivity tables are preserved separately; the descriptive table above adds no new significance claim. All 4,200 same-seed timed artifact comparisons between sessions are byte-identical, across the six recipes. Sessions remain separate, rather than being pooled or treated as extra workload coverage.

## Matched Benchpress extension

The eight unchanged cases cover cat-22, GHZ-23, dense Clifford-20 and BV-100 on square and heavy-hex graphs. New Q-Rust and TKET pass all eight. The previous Q-Rust control passes six, with unsupported local-factor recognition on the two dense cases. Qiskit passes four and has unsupported recognition on four. Cirq passes four, has unsupported recognition on three, and times out during the initial warmup for dense Clifford heavy-hex in each session; the three dependent timed slots in each session are recorded as blocked by that warmup. Unsupported recognition is not a demonstrated semantic counterexample.

Across both sessions, the matched records contain 32 unsupported Qiskit outputs, 24 unsupported Cirq outputs, 16 unsupported old-Q-Rust outputs, two Cirq warmup timeouts and six Cirq timed slots blocked by those warmups. No definite semantic counterexample or structural violation was detected. These outcomes remain in the denominator and the saved evidence; successful subsets cannot qualify an incomplete cell.

The stable four-tool verified cohort contains only the four cat/GHZ cases. On that cohort:

| Primary SDK | Compile-time ratio, S1 / S2 | Q-Rust / SDK CX | Q-Rust / SDK total depth | Q-Rust / SDK two-qubit depth |
|---|---:|---:|---:|---:|
| Qiskit | 3.491 / 3.351 | 1.8407 | 1.7637 | 1.7974 |
| Cirq | 21.604 / 21.617 | 1.8407 | 0.9019 | 1.7974 |
| TKET | 84.515 / 83.040 | 1.8407 | 1.7637 | 1.7974 |

Q-Rust remains faster under these recipes but uses more CX on every case in this cohort. The external-suite aggregate quality advantage therefore does not extend to these larger cat/GHZ cases. This table cannot be presented as a verified four-tool comparison over all eight cases.

| Case | Previous CX / depth / 2Q depth | New CX / depth / 2Q depth | Old/new eligibility |
|---|---:|---:|---|
| BV100, square | 326 / 279 / 265 | 309 / 263 / 247 | Both verified |
| BV100, heavy-hex | 353 / 350 / 336 | 345 / 343 / 328 | Both verified |
| Dense Clifford20, square | 4566 / 2536 / 2046 | 4416 / 2464 / 2017 | Old unsupported; new verified |
| Dense Clifford20, heavy-hex | 7850 / 3543 / 3091 | 7898 / 3612 / 3138 | Old unsupported; new verified |

These structural values are identical in both sessions. The dense rows are retained diagnostics, excluded from verified old/new aggregates. Heavy-hex gains verification coverage but has 48 more CX, 69 more total-depth layers and 47 more two-qubit-depth layers; total gate count decreases from 10,313 to 10,123. Verification coverage is accompanied by a structural cost on this case. These rows do not establish a verified old/new semantic quality comparison or isolate the contribution of each compiler change.

## Protocol and correctness scope

The external study retains 100 cases, six recipes, two sessions, seven timed trials and one initial warmup per cell. The matched extension retains eight cases, four tools, two sessions, three timed trials and one initial warmup. A previous Q-Rust request runs adjacent to every candidate Q-Rust request, with order balanced and reversed across sessions. In total the schedule contains 10,040 timed slots and 1,480 initial warmups. Six matched timed slots were blocked by failed warmups; this is not a claim that all 10,040 slots compiled successfully. Replacement warmups follow the unchanged policies.

Sessions ran sequentially: external 1, external 2, matched 1, matched 2. Their whole-process wall times were 838.45, 829.63, 550.85 and 551.34 seconds. No build, simulation, semantic validation or other benchmark session overlapped timing. External requests retained a 60-second timeout. Matched requests retained a 30-second timeout, approximately 100 ms RSS sampling, a 2 GiB per-worker limit and a 3 GiB aggregate-worker limit including the old-Q-Rust worker. Sampled RSS is neither exact peak memory nor descendant-process accounting.

The unchanged external oracle checks complete clean-ancilla isometries through six qubits and eight independently selected inputs above that within its existing limit. Sampling does not prove arbitrary-input equivalence. The matched oracle checks map/ancilla-aware signed-Clifford Choi isometries after bounded numerical recognition of local Clifford gates. It fails closed when recognition is unsupported. Definite tableau counterexamples and structural violations have separate classifications. Every eligible timing cell requires all its timed outputs and all initial/replacement warmups to pass. The input, adapter, oracle, SDK-version, schedule and resource contracts are unchanged, including Qiskit 2.5.2, Cirq 1.7.0 and pytket 2.18.5.

The new regressions cover interior meeting positions, inverse local conjugation, signed frames, barriers, noncommuting interactions, wire relabeling, directed edges, spare physical wires, repeated terminal classical writes, target-basis changes, numerical-ledger isolation, profiling consistency, retention guards and representation-veto scope. Independent literal matrices compare complete clean-ancilla isometries with one common global phase, without calling the compiler's own verifier.

The final checks passed 484 Rust tests, with four existing tests intentionally ignored; strict Clippy and formatting checks passed. All 249 Python tests ran with the frozen compiler binary and passed with no skips. The 14 canonical-replay fault-injection tests and five packaging boundary tests passed. Final independent fuzzing checked 384 compiled outputs from 64 fixtures on two topologies and 91 retention comparisons, and rejected all 128 relative-phase negative controls. All checks passed; maximum input residual was approximately `4.58e-15`, and 46 candidate outputs selected retained inputs. Provenance binds the exact timed binary and tested sources. These are untimed checks, not speed measurements.

## Provenance and preserved diagnostic history

The previous repository baseline is `1fe70295980fe5e3ac604e903203ddf0bc360b66`. Its timed compiler source was `4a6068c16684d440d3bfcb075025f69d163808b7`, binary SHA-256 `83a820711e00c6abd580ed744bf475b150b490ef7cb30fc91313e4fca08e6c07`. The release compiler commit is `65e6c64a66d68da79799c9ea6a344553cb844cc4`, binary SHA-256 `3ae3dcdc258e8b96b41d552964c3863702e17e4112c6e864f13376d33efaf891`. The paired controller now derives the old compiler identity from both supplied baseline protocols instead of hardcoding an earlier release.

The release `preflight_release/` is bound to this final binary and all 94 frozen Rust/Cargo source files. It verified all 100 external cases and all eight new matched outputs before formal timing. Its structural results match the formal results above. The final campaign is `campaign_release/` locally and `campaign/` inside the evidence archive.

Two earlier frozen configurations were not timed. The cumulative expanded-pass candidate in `campaign/` locally and `preflight_expanded/` introduced unsupported local-factor recognition on BV100 square. The block-local representation veto in `campaign_final/` and `preflight_final/` restored recognition but regressed QPE-5 and Trotter-4. Their configurations, binaries and outputs are preserved as abandoned configurations; neither is substituted for the final campaign. The whole-input scope of the final veto fixes those regressions without changing the checker.

The earlier `postrouting_v1/` test ran while explicit routing SWAPs fenced the search and showed no benefit. It does not establish that optimization after routing is ineffective. The positive `expanded_postrouting_v1/` control runs after their expansion. Both controls retain source, binaries and ablations. `interior_v1/`, `retention_v1/`, `clifford_diagnostics/`, `diagnostics/` and the other verification directories retain development evidence. Their durations are diagnostic only.

A preservation manifest binds 392 earlier files, including the earlier 189-file set, all previous tracked reports/result artifacts and the complete preceding campaign outside redundant packaging scratch and caches. Every file remained byte-identical after this campaign. Earlier results are not replaced.

## Evidence and replay

The result package is [external_results/20261018_clifford_quality](external_results/20261018_clifford_quality/), including [the evidence archive](external_results/20261018_clifford_quality/evidence.tar.gz) and [deliverable checksums](external_results/20261018_clifford_quality/checksums.json). The archive retains raw timed/control outputs, dispatch receipts, frozen protocols, diagnostics, source manifests, binaries, semantic receipts and reporting derivations. The original native Benchpress workout evidence remains unchanged.

The saved paired analysis and a fresh full-output replay agree after normalizing only the newly measured `validation.validation_ns` values. Presence of that field remains mandatory; the comparison updates only the corresponding derived row hashes. Compiler, parser and request timings, row order, output metrics, semantic verdicts and all other fields match exactly. The reporting derivation and local compile-cost derivation also reproduce byte-identical JSON. The reporting metrics SHA-256 is `5eb896490556f09b0079d36bed0224d8bdb1627cbfc0cb550e1076a7983f34d3`.

After extracting the archive, replay against a fresh output directory with the recorded Python dependencies:

```sh
python harness/experiments/paired_analysis.py --repo harness --campaign campaign --output replay_analysis
python scripts/canonical_replay.py --reference analysis --replay replay_analysis --output replay_receipt.json
python analysis/reporting/compute_metrics.py --analysis analysis --campaign campaign --output report_metrics_replay.json
python analysis/reporting/derive_local_costs.py --input analysis/reporting/report_metrics.json --output local_costs_replay.json
```

This campaign reruns the existing matched Q-Rust extension, not a native Q-Rust Benchpress gym. Its bounded Clifford workloads and ideal graph targets do not establish noisy-device performance, general large non-Clifford correctness, or universal compiler rankings. The retained [Benchpress pilot](BENCHPRESS_PILOT.md) documents the separate native execution and its different correctness contract.
