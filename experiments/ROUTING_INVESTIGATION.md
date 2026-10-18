# Routing costs and compiler recipe tradeoffs

This investigation follows the [corrected-fixture study](COMPARISON_STUDY.md).
Its purpose is to explain QFT/linear-routing losses and test generic changes
while preserving the same arbitrary-input, clean-ancilla and layout contract.
The original study and original competitor adapters remain reproducible.

The subsequent [two-qubit block investigation](BLOCK_SYNTHESIS_INVESTIGATION.md)
adds guarded resynthesis, fresh comparisons and a retained session repeat.
The measurements below remain the historical routing baseline.

## Measured results, 16 October 2026

The frozen compiler/harness commit is
`a945d1b832bed8cc0a561735516e63a3cea2aba3`. Three separate, sequential timing
campaigns retain the original 51 fixtures and 123 circuit/target cases, with
one warmup and seven measured repetitions per cell. All **7,749 measured
outputs passed** the unchanged independent checker: 1,722 revision comparisons,
3,444 primary comparisons and 2,583 recipe-sensitivity comparisons. Small
outputs use full isometries; larger outputs use sampled arbitrary inputs with
clean ancillas and one shared global phase. Maximum observed residual was
`3.684e-10`, below the declared `1e-7` tolerance. There were no excluded cells.
Numerical success is evidence for these artifacts, not a proof for all inputs.

The fresh interleaved before/after comparison found a **3.48×** compilation
speed ratio (95% fixture-cluster interval 3.11–3.92×) and a CX ratio of
**0.9364** (0.9153–0.9578). On linear targets the corresponding results were
**4.87×** (3.91–6.04×) and **0.8800** (0.8360–0.9236).

| QFT-14 on a line, fresh revision experiment | Median compile time | Median CX | Depth |
|---|---:|---:|---:|
| Saved previous Q-Rust | 401.55 ms | 524 | 282 |
| Revised Q-Rust | 74.57 ms | 351 | 200 |

The primary rerun retains the original competitor recipes:

| QFT-14 on a line, primary campaign | Median compile time | Median CX | Depth |
|---|---:|---:|---:|
| Q-Rust | 47.27 ms | 351 | 200 |
| Qiskit O3 | 26.62 ms | 291 | 186 |
| Cirq original recipe | 3278.91 ms | 956 | 1301 |
| TKET original recipe | 2157.74 ms | 501 | 226 |

Across all primary cases, paired geometric mean baseline/Q-Rust time ratios
were **6.55×** for Qiskit (5.07–8.42), **173×** for Cirq (150–200), and **475×**
for TKET (404–559). Q-Rust/baseline CX ratios were **1.04**, **0.800** and
**0.926**, respectively. Q-Rust was within 10% of the best of these three
baseline CX counts on **92/123 cases (74.8%)**. This is a fraction of cases
within a named empirical margin, not a percentage of optimality or fidelity.
On the 51 linear cases Q-Rust's speed ratio against Qiskit was **4.21×**
(3.00–5.86), but Qiskit still won the difficult QFT-14 case on both objectives.

In the separate full-corpus sensitivity campaign, Q-Rust was **93.39×** faster
than one-pass Cirq (81.42–106.70), with CX ratio **0.8085** (0.7617–0.8546),
and **86.52×** faster than two-qubit-only TKET (69.21–106.85), with CX ratio
**0.9325** (0.8703–1.0026). On QFT-14, its own interleaved Q-Rust result was
63.41 ms/351 CX, versus 1240.11 ms/483 CX and 481.97 ms/469 CX. These settings
were selected after diagnostic inspection; they are not a search for each
framework's best possible pipeline. Coarse stage timers add unquantified
instrumentation cost in this secondary harness. The targeted recipe ablations
also show genuine quality costs when synthesis is removed entirely.

Timing is visibly session-dependent: the same revised QFT artifact took
47.27, 63.41 and 74.57 ms across the three campaigns. All revised Q-Rust QASM
and layouts were identical across all cases/trials/campaigns. The reported ratios
compare results within each campaign; these fixture-cluster intervals do not
capture host/session uncertainty. Repetitions equalize sample counts, not
compiler search budgets.

The revision is not uniformly better: CX improved/tied/worsened on **37/83/3**
cases. The losses were Clifford-10 linear (505→525 CX; depth 256→364),
Clifford-10b linear (517→526 CX), and Ising-7 Nairobi (51→52 CX). Graphstate
depth increased by a geometric factor of **1.143** (1.105–1.190), despite
unchanged CX. W-14 all-to-all retains 39 CX versus 26 for every primary
baseline. These regressions limit the aggregate improvement.

The checked-in [compact results](routing_results_20261016.json) retain exact
estimates and intervals; the [cell table](routing_cells_20261016.csv) retains
all primary medians, ranges and eligibility. Full raw outputs, protocols,
diagnostics, source snapshots and independent audit receipts are preserved in
the accompanying `Q_Rust_Routing_Investigation_20261016.zip` artifact.

## Mechanisms and implementation changes

The original O3 layout budget was 50 trials with five bidirectional iterations:
up to 551 complete greedy routing evaluations before final beam routing. It
continued searching even after finding a zero-SWAP placement. Candidate scoring
also rebuilt the same successor sets repeatedly, and unordered hash-set ties
made otherwise seeded search nondeterministic. The revised implementation
caches shared scoring inputs, orders tied candidates deterministically, stops
at the zero-SWAP lower bound, and stops a deterministic forward/backward cycle
when its layout is unchanged. The zero bound is for inserted routing SWAPs,
not optimal lowered CX count or depth.

The extended-layer policy now includes up to 20 successor interactions and
normalizes its front/extended averages. This is a combined horizon and scoring
policy change, not an isolated experiment in lookahead length. Beam pruning
uses progress and remaining distance to distinguish equal-SWAP candidates and
deduplicates equivalent search states. These are heuristics, with no guarantee
of optimal routing. Common caching/tie/termination changes are a bundled
implementation intervention; the older binary is needed for an exact baseline.

Logical terminal permutations are extracted before placement/routing, then
composed into the returned physical output layout. This avoids routing gates
that can be represented by free output-wire relabeling under the declared
contract. Barriers and terminal measurement mappings remain protected. This
uses an established compiler technique, related to Qiskit's
[ElidePermutations](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.ElidePermutations).

Each SWAP admits two exact decompositions:
`CX(a,b) CX(b,a) CX(a,b)` and `CX(b,a) CX(a,b) CX(b,a)`.
Choosing the orientation next to a matching CX exposes cancellation using
`CX² = I`. The new pass chooses this algebraic orientation; the existing
dependency-aware cancellation pass performs the cancellation. It does not
introduce approximate two-qubit synthesis. Independent full-isometry tests
cover permutations, spare wires, phase preservation and the policy switches.

O3 uses a smaller default placement budget of ten trials and three iterations.
On the complete 123-case development sweep, that budget produced byte-identical
QASM and identical input/output layouts to the 50-by-five configuration.
`RoutingSearchConfig` retains explicit budget and policy overrides. This is an
empirical default choice on this corpus, not a claim about unseen workloads.
O3 retains its larger beam and post-routing cleanup relative to O2.

## Distinct experiments

1. Repeat the original four-compiler protocol on all 123 cases, with one warmup
   and seven measured trials. Keep Qiskit, Cirq and TKET's original recipes.
2. Interleave the saved old and new Q-Rust binaries on all cases, retaining
   every output and validating it after timing. This controls session drift
   better than using the previous study alone as a before/after comparison.
3. Use separately labeled per-pass profiles and routing-policy ablations to
   distinguish permutation, lookahead, progress, orientation and budget effects.
4. Measure competitor recipe sensitivity. The targeted five-case experiment
   compares three Cirq and three TKET recipes. A subsequent whole-corpus run
   interleaves revised Q-Rust with one-pass Cirq and two-qubit-only TKET.
   These alternatives were selected after inspecting the original results.

Timing and numerical simulation must run separately, with no builds or other
CPU-heavy jobs during timing. The original numerical oracle and tolerance are
unchanged. Full small isometries and sampled larger isometries are numerical
evidence, not proofs of every compiler output.

## Why the compiler names alone do not explain the results

Cirq's original recipe performs generic CZ synthesis before and after routing,
then converts CZ to CX and fuses one-qubit gates. Changing its iteration cap
can change the placement/routing problem as well as runtime: the same number
of pre-routing entanglers need not imply the same SWAP count. A cheaper exact
CX-preserving diagnostic supports only the declared 1q/CX/SWAP input subset;
it is not a general replacement for Cirq's compiler.

TKET's FullPeepholeOptimise includes optimization of up to three-qubit blocks.
PeepholeOptimise2Q omits that larger-block work; SynthesiseTket is lighter still.
The [TKET compilation manual](https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html)
documents these choices. More optimization can help synthesis but alter the
later routing problem unfavorably. W-state blocks also illustrate a real
tradeoff: generic two-qubit synthesis can save CX gates that a cheap local
pipeline retains.

Qiskit's pipeline includes pre-layout permutation elision and post-routing
two-qubit optimization. Its [LightSABRE](https://arxiv.org/abs/2409.08368)
implementation already uses Rust. TKET has a native core as well. These runs
therefore do not isolate a causal programming-language effect.

The Qiskit QFT-14 diagnostic makes the distinction between routing cost and
final output cost concrete. In one validated seeded run, its 87 routing SWAPs
gave 443 CX after lowering, reduced to 285 by TwoQubitPeepholeOptimization.
Q-Rust's high-budget diagnostic inserted fewer SWAPs (83), giving 431 CX
before its exact orientation/cancellation and 351 afterward. These are
individual stage observations, not a controlled replacement of one tool's
cleanup with another. On W-14, Qiskit consolidated 13 two-qubit blocks and
reduced 39 to 26 CX. Q-Rust's individual-gate KAK pass does not collect such
blocks and largely passes through an already lowered U/CX circuit.

Common exact input normalization fixes a low-level representation. A separate
native-input versus common-lowered experiment would test its influence; the
current comparisons do not establish behavior on every high-level input form.

## Reproducing the campaigns

Use the pinned environment and release build from [the primary protocol](COMPARISON_STUDY.md).
The measurements used Python 3.12.14, Qiskit 2.5.2, Cirq 1.7.0, pytket 2.18.5,
Rust 1.91.1, and an Apple M3 with 8 cores and 8 GiB RAM. Each timing destination
must be new. Run timing campaigns sequentially, with validation and other
CPU-heavy jobs stopped during each timing campaign.

```sh
# Primary: retain all original compiler recipes.
.venv-comparison/bin/python experiments/comparison_study.py time PRIMARY_DIR \
  --binary target/release/comparison_runner --machine-note 'Record host and controls'
.venv-comparison/bin/python experiments/comparison_study.py validate PRIMARY_DIR
.venv-comparison/bin/python experiments/analyze_comparison_study.py PRIMARY_DIR

# Fresh revision comparison: set before/after binary paths in variants.json.
.venv-comparison/bin/python experiments/compare_routing_revisions.py time REVISION_DIR \
  --variants variants.json --machine-note 'Record host and controls'
.venv-comparison/bin/python experiments/compare_routing_revisions.py validate REVISION_DIR
.venv-comparison/bin/python experiments/analyze_routing_revisions.py REVISION_DIR --reference before

# Separate exploratory sensitivity with Q-Rust interleaved again.
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
  QISKIT_NUM_PROCS=1 QISKIT_PARALLEL=FALSE PYTHONHASHSEED=20261015 \
  .venv-comparison/bin/python experiments/diagnose_compiler_tradeoffs.py time \
  --protocol PRIMARY_DIR/protocol.json --output SENSITIVITY_DIR --all-cases \
  --recipes qrust cirq_one_pass tket_two_qubit \
  --qrust-binary target/release/comparison_runner --trials 7 --probes 0
.venv-comparison/bin/python experiments/diagnose_compiler_tradeoffs.py validate --output SENSITIVITY_DIR
.venv-comparison/bin/python experiments/analyze_recipe_sensitivity.py SENSITIVITY_DIR
```

The old/new binary hashes, exact requests, source hashes and seeds are retained
in the protocols. New Q-Rust uses fixed internal layout starts and deterministic
ties; its record explicitly says the requested trial seed is not used.
Per-pass profiling is disabled in primary timing. All original datasets remain
unchanged. The implementation passed 316 workspace tests plus four CLI tests,
96 Python harness/analysis tests, formatting, strict Clippy and both CI jobs.

## Scope and further validation

The unchanged corpus records both improvements and losses after changes to
costly compiler stages. Stage observations support explanations of the observed
behavior; controlled policy comparisons support causal attribution to the
corresponding intervention. Coupled changes and internal routing details remain
unresolved where no separate ablation isolates them.

Examples of cross-compiler profiling include
[Arline Benchmarks](https://arxiv.org/abs/2202.14025) and
[Benchpress](https://www.nature.com/articles/s43588-025-00792-y).
The internal corpus does not establish behavior on external held-out circuits
or at larger scales. Equal-budget Pareto curves and repeated sessions would help
distinguish compiler search budgets and session effects.
No result here establishes near-optimality, hardware fidelity, memory savings
or intrinsic language speed.
