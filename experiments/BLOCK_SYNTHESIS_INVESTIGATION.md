# Two-qubit block synthesis: measured benefit, cost and remaining gaps

The new O3 pass closes the W-state CX gap on all-to-all and line targets, but
it does **not** close the QFT-14 linear-target gap. It is a measured trade-off:
less entangling work, more synthesis time, and sometimes more single-qubit gates
or depth. This investigation extends the immutable
[routing study](ROUTING_INVESTIGATION.md) and uses its same corrected fixtures
and comparator recipes. See [exact compact results](block_results_20261016.json)
and [all primary cells](block_cells_20261016.csv).

## What changed

`TwoQubitBlockSynthesisPass` collects dependency-connected U/CX regions on two
physical wires after routing and existing local cleanup. It stops at a third-wire
interaction; barriers and unsupported operations are conservative boundaries.
Unclaimed leading U runs may be absorbed, while replacement stays anchored at
the first CX so it cannot cross an intervening dependency on the other wire.
Regions with fewer than two CX are skipped. Pure-CX regions use exact finite-group
simplification; other regions use joint Cartan/KAK synthesis with zero-to-three
CX branches. Final direction repair and basis closure still follow this pass.

A candidate is accepted only for a **strict CX decrease** and independently
reconstructed, phase-aligned unnormalized Frobenius residual at most `1e-10`.
The sum of accepted residuals is limited to `1e-9` per pass run. Numerical failure
retains the original region. These local numerical budgets do not include
checker roundoff, earlier passes, later direction/basis closure or serialization.
They are not an end-to-end formal certificate. In exact arithmetic the Frobenius
norm bounds operator norm, and a telescoping product bounds composed replacement
error by the sum; implementation and numerical tests remain necessary.

The shared KAK kernel now synthesizes the Cartan interaction jointly instead
of decomposing its axes separately. Stable ZYZ phase extraction fixes an actual
near-diagonal conditioning failure exposed by an independent W-transfer test.
No acceptance tolerance was loosened. Reduced-CX candidates use a `1e-12`
internal reconstruction limit, with a full three-CX fallback at `1e-10`.
These are implementations of established synthesis results, not new theorems:
[Vatan and Williams](https://arxiv.org/abs/quant-ph/0308006) and
[Shende, Bullock and Markov](https://arxiv.org/abs/quant-ph/0308045).

The pass defaults on for O3 with basis lowering and post-routing cleanup enabled.
`.two_qubit_block_synthesis(false)` supplies a same-binary ablation. The worker
records this setting and candidate/acceptance/CX/residual counters. Default
benchmarking disables detailed pass clocks.

## Frozen experiment and correctness evidence

Compiler/harness commit: `734caeefcb78d1b22c9b69a209d541f1c10b1cfc`.
Release binary SHA-256:
`1fa371c5967246e93d9e52d3336ed3381e4d9f04e2e05add3c1d5b4a4c100c6a`.
The old arm is the previous routing compiler, saved from baseline commit
`1029cfc27fc1f1f6005b7c2633da5e16dd97f4b7`.
Later commits change reports/data only. Inputs, recipes and the oracle were
unchanged; all three revision arms retained identical initial/final layouts.

The pinned environment uses Python 3.12.14, Qiskit 2.5.2, Cirq 1.7.0 and pytket
2.18.5; Rust 1.91.1 release code runs on an Apple M3, 8 cores, 8 GiB RAM.
The corpus has 51 fixtures, nine families, 3–14 logical qubits and 123
circuit/target cases: 51 all-to-all, 51 line and 21 seven-qubit Nairobi cases.
Each cell has one warmup and seven measured repetitions with randomized order,
thread caps and persistent workers. Compilation uses warm public APIs;
parsing, import/startup, serialization, IPC and validation are outside its timer.
Matched repetitions do not mean matched optimizer search budgets.

Tests: **349 Rust tests, four CLI smoke tests and 96 Python tests passed**;
formatting and strict Clippy checks passed. New coverage includes 22 independent
block regressions, arbitrary-input QFT/W behavior, all 0/1/2/3-CX classes,
reversed/nonadjacent operands, third-wire boundaries, directed routing, ancillas,
readout, near-singular matrices, tiny coherent errors and accumulated residuals.
A separate generated stress suite passed **640/640 full-operator checks** on
320 cases with synthesis off/on. Its maximum unnormalized Frobenius residual
was `9.9714e-10`, below its additional `1e-8` stress threshold; no CX increases
or tested class-bound misses occurred. This is generated stress, not an
external application holdout.

| Campaign | Measured outputs checked | Result |
|---|---:|---|
| Old / new pass-off / new pass-on | 2,583 | All pass |
| Primary Q-Rust / Qiskit / Cirq / TKET | 3,444 | All pass |
| Q-Rust / one-pass Cirq / two-qubit TKET | 2,583 | All pass |
| Separately retained revision session repeat | 2,583 | All pass |
| Separate five-case pass profiles | 70 | All pass |

The three planned campaigns account for 8,610 outputs; with the revision repeat
there are 11,193 measured outputs. The 70 diagnostic profile outputs and 369
sensitivity warmups also pass, for **11,632 validated retained rows**. There
are no excluded cells. Validation uses saved output QASM, an independent
analytic NumPy oracle, explicit layouts and full clean-ancilla amplitudes.
Through six physical qubits it checks full isometries; wider circuits use
three basis inputs, a uniform state and four Haar inputs with one shared
phase. The maximum observed per-input residual is `3.6844e-10`, below `1e-7`.
Sampled checks are not all-input proofs for wide circuits; repeated rows share
cached checks when their full semantic artifacts are identical.

## Benefit and cost of enabling block consolidation

The same-binary on/off comparison isolates the additional block pass. Across
123 cases, CX improved/tied/worsened on **33/90/0**. The paired geometric mean
CX ratio was **0.9464** (95% fixture-cluster interval 0.9126–0.9771).
The pass considered 2,370 blocks, attempted numerical synthesis 1,944 times,
used the pure-CX path 426 times, accepted 190 replacements and removed 237 CX
(counted once per case). There were zero numerical rejections; 2,180 blocks
offered no strict CX improvement. The maximum accepted local residual was
`9.3853e-15`, and the largest circuit residual sum was `3.3620e-14`.

Compilation time increased by a geometric factor of **1.3105**
(1.2262–1.4098). A full-corpus repeat, retained separately after unusually wide
first-session timing spread, measured **1.3116** (1.2348–1.3981). The full new
compiler versus the saved old compiler cost 1.2918 and 1.2781 times as much in
those respective sessions. New code with the pass disabled left CX, total gate
count and depth unchanged in every case; its aggregate timing interval includes
parity. Do not attribute all timing changes to shared KAK code or Rust.

Total gate count increased by a geometric factor of **1.0418** (1.0062–1.0832),
and depth by **1.0143** (1.0001–1.0307). Total gates increased on 18 cases;
depth increased on 12. Strict CX improvement alone is therefore not a guarantee
of better circuit duration, total size or hardware fidelity.

| Case, same-binary off → on | CX | One-qubit gates | Depth |
|---|---:|---:|---:|
| W-14, all-to-all or line | 39 → 26 | 27 → 78 | 53 → 65 |
| QFT-14, line | 351 → 351 | 274 → 274 | 200 → 200 |
| Clifford-12, line | 737 → 730 | 240 → 241 | 346 → 345 |

W-14's two-qubit depth improves from 39 to 26, but its extra local rotations
increase overall depth. The W transfer is optimized as a full two-qubit unitary,
not replaced by a circuit that only produces the same state from zero.

## Fresh comparisons with the unchanged primary recipes

| Baseline | Baseline / Q-Rust compile time, geometric mean [95% interval] | Q-Rust / baseline CX [95% interval] |
|---|---:|---:|
| Qiskit O3 | 7.17× [5.33, 9.44] | 0.987 [0.969, 1.003] |
| Cirq original recipe | 136.84× [120.35, 157.76] | 0.758 [0.718, 0.799] |
| TKET original recipe | 328.35× [278.25, 391.64] | 0.876 [0.822, 0.929] |

Q-Rust was faster in 116/123 cells against Qiskit, and 123/123 against each
other primary recipe. Its CX ratio interval against Qiskit spans parity;
this does not establish a general quality advantage. It was at or below the
best baseline CX count on **96/123 cases**, within 5% on **99/123 (80.5%)**,
and within 10% on **110/123 (89.4%)**. These are fractions of this development
corpus within an empirical margin, not percentages of an unknown optimum.

| QFT-14, line: primary campaign | Median compile time | CX | Depth |
|---|---:|---:|---:|
| Q-Rust | 50.18 ms | 351 | 200 |
| Qiskit | 22.65 ms | 291 | 186 |
| Cirq | 3,220.58 ms | 956 | 1,301 |
| TKET | 2,211.60 ms | 501 | 226 |

Qiskit still wins this QFT case in both time and CX. For W-14 all-to-all all
four compilers now use 26 CX, but Q-Rust depth is 65 versus 53 for all three
baselines. On W-14 line Q-Rust/Qiskit/TKET have 26 CX; Cirq has 47. The primary
line timings were 0.818, 3.547, 175.929 and 456.233 ms respectively. W-14
all-to-all showed much larger timing spread and must not be silently replaced
with the line case: its Q-Rust median was 2.817 ms and Qiskit's 36.264 ms.
All medians and ranges remain in the cell table.

The separate alternative-recipe campaign still favors Q-Rust in compilation
time: **74.54×** versus one-pass Cirq (65.72–85.64), and **65.55×** versus
two-qubit TKET (55.00–78.60). CX ratios are **0.7652** (0.7231–0.8062) and
**0.8826** (0.8303–0.9313). These selected recipes are a sensitivity check,
not a per-framework search for the best achievable configuration. Their coarse
stage clocks introduce unquantified instrumentation overhead.

## Why QFT remains unchanged

Independent analytic reconstruction and magic-basis invariants classify the
134 collected QFT-14 line regions as **51 requiring two CX and 83 requiring
three CX** at the stated `1e-9` numerical classification tolerance. Their costs
sum to all 351 CX. Every collected region already meets its individual class
bound. Three-CX cases have a minimum absolute imaginary invariant trace of
3.696, far from the numerical classification threshold; the two-CX cases'
imaginary trace residual is at most `4.45e-16`. Canonical identity, CX, iSWAP,
SWAP and generic Cartan self-checks validate the diagnostic conventions.
The archived `collected_block_classes.json` retains region indices, wire pairs,
margins and source hashes. This is numerical classification using established
results, not a formal implementation proof or a global circuit lower bound.

The block switch preserves the chosen route and dependencies. It cannot move a
third-wire interaction, create a larger region, or choose Qiskit's different
route. Earlier Qiskit tracing showed that its particular routed representation
exposed broader CX savings. Further QFT improvement therefore needs a different
route, ordering or partition, or a larger optimization scope. Simply making
this same two-qubit resynthesis more aggressive cannot reduce these fixed
regions while preserving their full unitary behavior.

Separate diagnostic pass profiles (excluded from headline timing) measured
3.061 ms / 6.66% of total for the QFT-14 block pass, with no accepted replacement;
0.357 ms / 47.82% for W-14 all-to-all, with 13 accepted replacements; and
3.241 ms / 2.88% for Clifford-12 line, with five accepted replacements.
Profiled QASM, layouts and counters exactly match the uninstrumented outputs.

## Timing limits and publication implications

All local builds, simulations, rendering and other agent compute were paused
while timing. The desktop's OS/background load, frequency and scheduling were
not controlled. QFT-14 pass-on medians were 267.76 ms in the first revision run,
50.18 ms in the primary campaign, and 52.52 ms in the retained revision repeat.
The repeat strengthens the aggregate overhead observation, but does not repair
primary session noise. Keep campaigns separate, retain all ranges and avoid
precise deployment latency claims. The 2,000-draw fixture-cluster intervals
resample related topology cases together; they do not model host/session noise.

A defensible result is: on this fixed small-circuit corpus and specified
recipes, bounded block synthesis reduces CX while preserving independent
numerical equivalence checks and Q-Rust's aggregate compilation speed advantage.
It closes a diagnosed W-state gap, quantifies the price in compile time and
local gates, and establishes why this particular local intervention cannot
close the QFT gap. It does not establish whole-circuit optimality, Rust as the
cause of speed, broad hardware quality or a novel synthesis algorithm.

Before publication, prioritize an untouched external circuit holdout, broader
widths/topologies, matched-cost Pareto comparisons, controlled repeated sessions
and another machine. For implementation, investigate local-gate/depth-aware
selection or cleanup for W, and routing/partition objectives that account for
post-synthesis CX for QFT. The current corpus motivated the intervention and
must be described as development evidence. These results strengthen a careful
compiler engineering/evaluation paper; they are not sufficient on their own
for a general algorithmic superiority claim.

## Reproduction and artifacts

Use the environment and primary commands in [COMPARISON_STUDY.md](COMPARISON_STUDY.md).
For the paired experiment, `compare_routing_revisions.py time` accepts a JSON
map of named binaries and request overrides. Use the frozen old binary for
`before`, and the same frozen new binary for `blocks_off` and `blocks_on`, with
`two_qubit_block_synthesis` set to false/true respectively. Then run its
`validate` command and `analyze_routing_revisions.py --reference before`.
The on/off summary uses the same analyzer with `blocks_off` as reference;
preserve both summaries. Recipe sensitivity uses
`diagnose_compiler_tradeoffs.py` and `analyze_recipe_sensitivity.py` as in the
routing study. Do not rerun into existing raw-output directories.

`Q_Rust_Block_Synthesis_Investigation_20261016.zip` retains all five campaigns,
raw QASM and trials, validation, summary audits, exact binaries, frozen source,
design declarations, generated stress suite, classification and profile
diagnostics, figures and the updated publication plan. Its manifest hashes
every payload. The post-results session repeat remains explicitly labeled;
none of the original campaigns or earlier bundles was overwritten.
