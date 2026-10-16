# Routing costs and compiler recipe tradeoffs

This investigation follows the [corrected-fixture study](COMPARISON_STUDY.md).
Its purpose is to explain QFT/linear-routing losses and test generic changes
while preserving the same arbitrary-input, clean-ancilla and layout contract.
The original study and original competitor adapters remain reproducible.

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

## Publication interpretation

The useful contribution is a reproducible account of which work pays off:
identify a costly stage, intervene on it, check semantics, and report both
improvements and losses on the unchanged corpus. Stage observations support
mechanisms; only controlled policy comparisons support causal attribution to
the corresponding intervention. Coupled changes and internal routing details
must remain qualified where no separate ablation isolates them.

Cross-compiler profiling itself has precedents, including
[Arline Benchmarks](https://arxiv.org/abs/2202.14025) and
[Benchpress](https://www.nature.com/articles/s43588-025-00792-y).
A paper must explain the reusable insight or software need beyond another
small internal ranking. External held-out circuits, scaling, equal-budget
Pareto curves and repeated sessions remain important before broad claims.
No result here establishes near-optimality, hardware fidelity, memory savings
or intrinsic language speed. Deployment and research usefulness need their
own evidence for a software-focused publication.
