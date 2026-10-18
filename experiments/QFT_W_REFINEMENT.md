# QFT routing and W-state boundary refinement

Frozen compiler `8777f0d4442b8d4844ceb900a99034e518727629`; 16 October 2026. This extends the historical [block-synthesis investigation](BLOCK_SYNTHESIS_INVESTIGATION.md). Exact estimates, ranges, adverse cases and profile traces are in [the compact results](qft_w_results_20261016.json); [the primary cell table](qft_w_cells_20261016.csv) retains every compiler/case cell.

QFT-14 on a line now uses **274 CX at depth 157**, down from **351 CX at depth 200**. W14 on both line and all-to-all targets retains **26 CX** while local gates fall **78 → 66 U** and depth **65 → 53**. The QFT quality improvement costs compilation time: the new/old median-time ratios are **1.60× and 1.27×** in two separately planned sessions. Aggregate warm compile time stays near the previous version; this does not establish negligible overhead on every circuit.

## What changed and why

At level 3 with a backend, exact `CX(c,t) RZ_t(theta) CX(c,t) = RZZ(c,t,theta)` recovery runs before named-angle crystallization. Only disjoint operations can be skipped, and non-gate operations are global fences. This preserves the two-wire interaction as an atomic diagonal operation through layout and routing.

Layout probes and final routing share a dependency graph that permits diagonal interactions to commute. Non-diagonal operations retain dependencies on all earlier unfinished interactions on their wires; reverse layout search reverses the actual edges. The actual router segments at barriers. A fixed 64-interaction cap per wire adds conservative ordering, bounding immediate predecessors at 256 per node. It trades scheduling freedom for bounded dependency construction; it does not make the whole routing search linear.

After two-qubit block synthesis, a linear scan fuses same-wire U runs across disjoint operations. A touching CX closes a run; every non-U/CX or malformed operation is a global boundary. Materialization at the last U preserves the CX sequence and layouts. Independent analytic 2×2 matrices check each candidate: phase-aligned Frobenius residual ≤1e-12 per run, summed accepted residual ≤1e-10, and ≤1e-14 for scalar-identity removal. Failed or over-budget replacements retain the original run. These are floating-point acceptance guards, not certified rounding-error bounds.

In exact arithmetic, a local Frobenius residual upper-bounds its operator norm; embedding preserves that operator norm, and telescoping products of unitary factors bounds the full operator difference by the sum of local bounds, up to one accumulated global phase. Floating-point matrix construction/checking and later QASM serialization require the separate numerical evidence below.

Controls are `.recover_phase_gadgets(bool)`, `RoutingSearchConfig::commute_diagonal`, and `.post_block_fusion(bool)`. They default on for their level-3 stages. The new public struct fields require explicit struct-literal callers to supply the fields or use `..Default::default()`. Reports and the comparison worker retain effective settings and acceptance counters.

### Boundary attribution

The separately instrumented, independently audited QFT14 trace is: old route **83 SWAP → 431 lowered CX → 351 final CX**; new route **78 SWAP → 416 lowered CX → 274 final CX**. Five fewer SWAPs account for 15 CX. The existing inverse-cancellation pass removes 80 CX before and 142 after, accounting for the remaining 62 saved CX. The new route accepts **0 of 98 KAK block candidates**, versus 0 of 134 before. The improvement therefore exposes better routing/cancellation boundaries; it is not a stronger KAK kernel. Recovery alone produces 368 CX. The old fixed-partition local-minimality analysis never ruled out a different route or partition.

Broader search is costly: in the diagnostic profile, layout time rises 7.78 → 22.96 ms and beam-routing time 37.92 → 58.90 ms; layout evaluations rise 57 → 71. Recovery itself takes about 0.095 ms and KAK drops 3.44 → 3.04 ms. These instrumented timings support attribution only. The two uninstrumented sessions establish the observed QFT latency trade-off.

W14 instead accepts 12 local fusions, removes 12 U gates and preserves every CX and layout. Its maximum local residual is 3.91e-16 and summed residual 3.03e-15.

Phase gadgets and optimization-aware routing are described in [Cowtan et al., EPTCS 318 (2020)](https://arxiv.org/abs/1906.01734), [Liu et al., NASSC](https://arxiv.org/abs/2205.10596), and [Qiskit Commuting2qGateRouter](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.Commuting2qGateRouter), which documents positioning RZZ near SWAP to exploit CX cancellations.

## Unchanged experiment and controls

The same 51 fixtures / 123 circuit-target cases, common normalization, U/CX basis, topology/ancilla contract, warm API boundary, one warmup and seven measured trials were used. Versions remain Qiskit 2.5.2, Cirq 1.7.0, pytket 2.18.5 and the pinned Python environment. Host: Apple M3, 8 cores, 8 GiB RAM; thread caps one. Parsing, serialization, IPC and verification are outside compile timing. Builds, tests, simulation and rendering were paused during timing. CPU affinity, frequency and OS background load were uncontrolled.

Two full seven-arm revision sessions were planned before official timing. Each has 6,027 measured rows: `before` (archived binary), `disabled`, `recovery_only`, `commutation_only`, `routing_only` (recovery + commutation), `fusion_only`, and `full`. Block synthesis stays enabled. This contains the recovery × commutation factorial and the selected routing-package × fusion factorial, plus the old/new-disabled control. Outputs and layouts match across sessions; `before` and `disabled` are byte-identical on every case. Timings are kept separate.

The primary four-tool campaign contains 3,444 measured rows. The separately labeled one-pass Cirq / two-qubit TKET sensitivity contains 2,583 measured rows plus 369 warmups. Main measured total: **18,081**. A separate six-case / four-arm profile adds 168 measured rows. Recipe selection and unequal optimization work remain limitations; the secondary harness has coarse stage clocks with unquantified overhead.

## Results

### Aggregate primary comparisons

Equal-cell geometric means of paired medians, with fixture-cluster 95% bootstrap intervals (2,000 draws). Time ratio is baseline / Q-Rust; quality ratios are Q-Rust / baseline. Repeated trials are not independent workloads. Cluster intervals do not capture host/session variation. Zero numerator/denominator cases are counted separately rather than entering log means.

| Baseline | Compile speed ratio [95% CI] | CX ratio [95% CI] | Depth ratio [95% CI] |
|---|---:|---:|---:|
| qiskit | 8.597 [6.651, 11.018] | 0.952 [0.928, 0.975] | 0.863 [0.797, 0.922] |
| cirq | 125.775 [107.517, 148.331] | 0.731 [0.693, 0.770] | 0.519 [0.472, 0.567] |
| tket | 297.388 [242.451, 363.910] | 0.845 [0.791, 0.897] | 0.844 [0.779, 0.907] |

Q-Rust has at least as good a median CX count as the best baseline in **101/123 cases (82.1%)**, is within 5% in **103/123 (83.7%)**, and within 10% in **112/123 (91.1%)**. This is empirical coverage against the named recipes, not a percentage of global optimum. Depth is at or better than the best baseline in 89/123 cases.

Against the alternative recipes, compile-speed geometric means are **72.426×** for one-pass Cirq [62.578, 83.897] and **63.412×** for two-qubit TKET [51.929, 76.914]. Corresponding CX ratios are 0.738 and 0.851. These come from the separate sensitivity campaign.

### Selected primary cells

Medians from the current primary campaign only; compile ranges are retained in the compact data and CSV.

| Case | Compiler | Compile ms | CX | U/CX depth |
|---|---|---:|---:|---:|
| clifford_12__linear | qrust | 167.841 | 730 | 344 |
| clifford_12__linear | qiskit | 96.880 | 637 | 274 |
| clifford_12__linear | cirq | 4036.797 | 930 | 795 |
| clifford_12__linear | tket | 2638.769 | 814 | 392 |
| qft_14__linear | qrust | 113.747 | 274 | 157 |
| qft_14__linear | qiskit | 53.164 | 291 | 186 |
| qft_14__linear | cirq | 4290.473 | 956 | 1301 |
| qft_14__linear | tket | 2634.414 | 501 | 226 |
| wstate_14__all_to_all | qrust | 1.091 | 26 | 53 |
| wstate_14__all_to_all | qiskit | 3.642 | 26 | 53 |
| wstate_14__all_to_all | cirq | 76.274 | 26 | 53 |
| wstate_14__all_to_all | tket | 481.979 | 26 | 53 |
| wstate_14__linear | qrust | 0.959 | 26 | 53 |
| wstate_14__linear | qiskit | 115.269 | 26 | 53 |
| wstate_14__linear | cirq | 272.215 | 47 | 88 |
| wstate_14__linear | tket | 541.174 | 26 | 53 |

QFT14 now has 5.8% fewer CX and 15.6% lower depth than Qiskit’s median, but takes **2.14× longer** in this primary session. W14 matches Qiskit and TKET at 26 CX / depth 53. Clifford12 remains worse than Qiskit in both quality and this session’s latency. The wide timing ranges, particularly Qiskit W14 line and the second revision session, rule out stable per-case latency claims from these runs.

### Controlled revision effects and adverse cases

Old/new compile-speed ratios are **0.98456 [0.94774, 1.02288]** and **0.99077 [0.93748, 1.04586]** in sessions 1 and 2. Both intervals include parity. This supports an observed near-unchanged aggregate cost, not proof of timing equivalence. QFT14 old → new medians are **55.196 → 88.448 ms** and **138.836 → 176.179 ms**. W14 line is **0.818 → 0.843 ms** and **0.805 → 0.827 ms**; W14 all-to-all is **0.777 → 0.811 ms** and **1.186 → 1.999 ms**, with a wide second-session range.

Across all 123 cases, full versus before: CX **19 improved / 102 tied / 2 worsened**; total depth **43 / 78 / 2**; total gates **49 / 71 / 3**. Isolating local fusion on the new route leaves CX identical in all 123 cases, reduces total gates in 32 and depth in 26, with no gate/depth regressions.

Every case with a worse CX, gate count, total depth or CX-only depth is listed below. Entries are before → full.

| Case | CX | Gates | Depth | CX-only depth |
|---|---:|---:|---:|---:|
| graphstate_5__ibm_nairobi | 9 → 11 | 19 → 18 | 13 → 13 | 7 → 9 |
| graphstate_7__ibm_nairobi | 16 → 16 | 23 → 29 | 19 → 13 | 13 → 9 |
| qft_5__ibm_nairobi | 30 → 30 | 61 → 65 | 46 → 42 | 27 → 25 |
| qft_5__linear | 30 → 34 | 61 → 68 | 38 → 43 | 23 → 27 |
| qft_7__linear | 74 → 71 | 138 → 124 | 74 → 74 | 47 → 50 |
| qft_8__linear | 101 → 91 | 186 → 161 | 92 → 97 | 59 → 64 |

The heuristic does not directly minimize final post-cancellation CX or scheduled depth, so extra commutation freedom need not improve every case. Retaining these losses is part of the result. No per-fixture routing policy or fallback was selected.

## Correctness and reproducibility

All **18,081 main measured outputs**, **168 profile outputs**, and **369 sensitivity warmups** passed the unchanged independent analytic NumPy oracle: **18,618 validated rows**, zero failures and zero excluded cells. Repeated semantic artifacts reuse cached checks. Each revision session has 207 unique artifacts; primary 698; sensitivity 369; profiles 16. Primary unique evidence comprises 123 full physical isometries and 575 sampled isometries. Separate primary/revision/profile warmups remain raw records outside that numerical count.

The oracle checks full clean-ancilla physical isometries through six physical qubits; larger outputs use three basis inputs, uniform superposition and four seeded Haar inputs with one shared global phase and all physical amplitudes. It also checks basis, connectivity and layouts. Maximum observed input residual is 3.68434e-10 against the unchanged 1e-7 threshold. Finite samples do not prove all-input correctness at larger widths.

Checks: **378 Rust tests**, four separately enabled CLI tests, **96 Python tests**, formatting and all-feature Clippy. New independent regressions cover analytical Fourier columns, arbitrary-input W operators, directed clean-ancilla routing, barriers/readout, reverse dependency fences, infinitesimal non-diagonal gates, long epochs, malformed direct calls and local-fusion endpoints/budgets. A separate 48-case × 8-control stress study validated **384 outputs / 305 unique artifacts, all checked by full physical isometry** on directed targets, with maximum residual 2.98612e-10. That generated suite is development evidence, not an external holdout.

Compiler CI: [passed](https://github.com/Arturacu/Q-Rust/actions/runs/38974784137). Release binary SHA-256: `9613ed4f819597744711cc755dd6d10f961013b51f1d92dae21bcc83d2a40c2e`. Baseline source/report commit: `1fe218db8eb88262b5d26c99d4b9c63ffff80ec4`; archived baseline worker SHA-256: `1fa371c5967246e93d9e52d3336ed3381e4d9f04e2e05add3c1d5b4a4c100c6a`.

The accompanying `Q_Rust_QFT_W_Refinement_20261016.zip` retains exact binaries, source archives, protocols, all trial QASM/layouts, validation receipts, development records, independent audits, figures and a file-hash manifest. Prior bundles remain separate. Local paths in scripts must be adapted on another machine; changed binaries/environments constitute a new study.

Reproduce a fresh primary campaign using the pinned environment and frozen revision:

```sh
cargo build --release -p qrust-bench --bin comparison_runner
python experiments/comparison_study.py time experiments/results/NEW_PRIMARY --binary target/release/comparison_runner
python experiments/comparison_study.py validate experiments/results/NEW_PRIMARY
python experiments/analyze_comparison_study.py experiments/results/NEW_PRIMARY
```

Set OMP/OPENBLAS/MKL/NUMEXPR/VECLIB/RAYON/QISKIT thread caps to one, QISKIT_PARALLEL=FALSE and PYTHONHASHSEED=20261015 as recorded in the protocol. Use `compare_routing_revisions.py` with the bundled seven-arm `variants.json` for each separate revision session; use the unchanged `diagnose_compiler_tradeoffs.py` recipes for sensitivity. Run validation and analysis only after all timing has completed.

## Scope and further validation

On this development corpus, changes to representation and dependency boundaries improve output quality with measurable runtime costs and some adverse cases. These measurements do not establish whole-circuit optimality, hardware fidelity or a causal speed benefit from Rust. The same corpus guided development. Generalization remains untested on untouched external workloads, larger widths and other topologies. Matched-cost quality/latency curves, controlled repeated sessions and another host would help separate search-budget and machine effects.
