# Exact Clifford optimization and preserved-baseline comparison

The new level-3 passes recognize encoded SWAPs and simplify Clifford
interactions across two-qubit block boundaries. In controlled checks of the
five-qubit error-correction circuit, Q-Rust changes from 47 CX/depth 76 to
11 CX/depth 14 on the complete graph, and from 65 CX/depth 86 to
27 CX/depth 30 on the line. These improvements overlap; they are not additive.

This follow-up preserves the earlier [TKET investigation](TKET_CLIFFORD_INVESTIGATION.md),
[Benchpress pilot](BENCHPRESS_PILOT.md), and all previously recorded results.
It adds separate outputs and contemporaneous runs of the frozen old compiler.

## Implementation and correctness boundary

`EncodedSwapRecognitionPass` searches dependency-local windows of at most
64 operations on two wires. It propagates the signed images of X and Z using
integer algebra. It replaces a window with SWAP only when all four generator
images match SWAP, establishing equality up to one global phase. Disjoint
operations can remain between the selected operations; an interaction with a
third wire stops the window. Barriers and dynamic operations are fences.

The error-correction input contains twelve recognized windows. The earlier
diagnostic found eleven contiguous spellings. The additional window consists
of zero-based source gate indices 93, 95, 96, 98, 99, 100 and 101 on wires 4 and 2;
the intervening CX(0,1) and H(1) are disjoint. It uses the same identity:

```
CX(a,b) H(a) H(b) CX(a,b) H(a) H(b) CX(a,b)
  = CX(a,b) CX(b,a) CX(a,b)
  = SWAP(a,b).
```

Recognizing that additional window explains the improvement from the earlier
14-CX diagnostic with an exchanged output map to 11 CX with identity maps.
The recognizer does not depend on fixture names or literal adjacency.

`CliffordSimplificationPass` transports the signed axes of
`CX = (I + P + Q - PQ)/2`, initially `P=Z_control` and `Q=X_target`, through
named local Clifford gates. It crosses an intervening third-wire interaction
only when that interaction commutes with both axes. At a later CX on the same
pair, an exact finite search considers a local replacement or `L CX R`, and
optional SWAP factors when output permutation absorption is enabled.

The search uses the 24 one-qubit Cliffords and 24×24 local frame pairs, rather
than initializing the full two-qubit Clifford group. Encountered products are
cached. There are two sweeps with 96-operation lookahead and a cap of 96
whole-circuit candidate evaluations per invocation. Exhaustion retains the
valid partial result and is reported explicitly. Each accepted candidate
strictly reduces the pass's CX cost and does not increase total depth or its
internal CX-family weighted depth: CX/CY/CZ have weight one, SWAP has weight
three, and other gates retain dependencies with zero weight. These costs are
projected through permutation absorption when enabled. The final benchmark
outputs use U/CX, where the reported two-qubit depth counts CX layers.

Both passes recognize only named Clifford gates. They introduce no new
floating-point classification tolerance. Existing crystallization remains
subject to the shared numerical error ledger. Exact rewrite identities and
independent numerical tests are distinct evidence: neither makes the entire
floating-point compiler formally verified.

For backend compilation with permutation elision, the order is recognition,
permutation extraction, Clifford simplification, then extraction and
composition of newly introduced permutations. If the first map is `p1` and
the second is `p2`, the composed logical map is `p2[p1[l]]`; it is subsequently
composed with the physical output layout. Without backend permutation
elision, cross-block simplification introduces no SWAP factors.

The new configuration switches `encoded_swap_recognition` and
`clifford_simplification` default to enabled and take effect only at level 3.
The comparison adapter records their effective values and pass counters.

## Controlled ablations

These profiled diagnostics establish mechanisms, not comparative speed.
Every output passes the independent analytic oracle on all 32 input columns
with one shared phase, including layouts, basis and coupling legality. The
largest column residual across the eight outputs is 1.85×10⁻¹⁵.

| Target | Encoded SWAP | Cross-block Clifford | CX | Depth | Two-qubit depth |
|---|---|---|---:|---:|---:|
| Complete | Off | Off | 47 | 76 | 46 |
| Complete | Off | On | 15 | 24 | 13 |
| Complete | On | Off | 11 | 14 | 8 |
| Complete | On | On | 11 | 14 | 8 |
| Line | Off | Off | 65 | 86 | 56 |
| Line | Off | On | 30 | 36 | 26 |
| Line | On | Off | 27 | 32 | 24 |
| Line | On | On | 27 | 30 | 24 |

With both passes enabled, cross-block simplification considers three pairs,
accepts one one-CX reduction and rejects two candidates for depth. Later
block synthesis overlaps with this reduction on the complete graph. On the
line, the final CX count stays at 27 while depth improves from 32 to 30.

TKET's previously established complete-graph result is 9 CX/depth 13. The
pass deliberately restricts search and commutation and uses a conservative
named-gate depth metric. The present controls do not isolate every contribution
to the remaining gap. The rejected candidates do not
establish that this implementation found a valid nine-CX circuit. Nine CX has
not been established as optimal.

## Validation

The implementation adds 22 independent adversarial integration tests and
three inline algebra/budget tests. Coverage includes signed-axis products,
relative phases, interleaved windows, noncommuting boundaries, near-Clifford
angles, composed three-cycles, clean spare wires, directed routing, terminal
readout ordering, W-state isometries and all feature-switch combinations.
At levels 0–2, switching the new options on or off does not change output or
report values.

The all-feature Rust suite passes 461 tests, and strict Clippy checks pass.
All 244 Python regression tests pass with the frozen candidate executable,
including the six public-pipeline numerical tests that require that binary.
Fourteen additional archive-replay fault tests reject changes to verdicts,
metrics, compiler timings, row order, coverage and checksums, while permitting
only newly measured verifier durations to differ.
Independent NumPy checks cover 64 generated fixtures on two topologies with
both compiler versions: all 256 retained outputs pass full-isometry checks,
and all 128 relative-phase negative controls are rejected. Maximum residual
is 4.58×10⁻¹⁵. Initial fuzz-harness limitations and subsequent replay receipts
are preserved separately; the compiler outputs were not changed.

The untimed preflight passes all 100 external cases with both compiler
versions. Of eight larger matched cases, six pass with each version; the two
dense random-Clifford cases remain unsupported by the existing gatewise
Clifford recognizer. Their status is unverified. Additional evidence labels
separate unsupported recognition from actual semantic counterexamples or
structural contract violations without changing historical outcome labels.

## Repeated comparison protocol

The compiler implementation is committed as `36ddc6375d4690f3f16966a67395cfe4708725be`.
The campaign freezes source at `4a6068c` and the candidate executable at
SHA256 `83a820711e00c6abd580ed744bf475b150b490ef7cb30fc91313e4fca08e6c07`.
The preserved old executable is
`f2bfb51d945e4a959450022a4fac4921a30d3a888d8664fb12d6df59bf689a80`,
from implementation commit `837a53ca9a9d4ca3f25cab3e5866b170e582c4a4`.

The original external matrix remains 100 circuit/target cases × six recipes
× seven trials × two sessions, with one warmup per case/recipe/session.
The original matched Benchpress extension remains eight cases × four recipes
× three trials × two sessions, with the same warmups. Old-Q-Rust controls run
immediately adjacent to each new-Q-Rust request; their order is balanced
across sessions. The original schedules remain exact projections of the
expanded dispatch. This prescribes 10,040 timed requests and 1,480 warmups.

Adapters, inputs, SDK versions, recipes, checks and resource policies retain
their frozen contracts. Timing and validation are separate. All warmups and
timed outputs must pass for a cell to qualify; failures and unsupported
outputs remain in the records. Compilation timers exclude parsing,
serialization, process startup and verification. The results describe warm
compilation, including reuse of the new Clifford cache; they do not isolate
the effect of the Rust language or characterize cold process startup.

The native Benchpress results remain preserved. They are not rerun here:
Q-Rust has no native gym, and the unchanged native SDK recipes use different
bases, timers and semantic contracts from this matched extension.

## Measured results

All four sessions completed. Every external output passed its prescribed
independent check: 8,400 primary timed outputs and 1,200 primary warmups,
plus 1,400 old-Q-Rust timed controls and 200 control warmups. The contract
uses full-isometry checks on small inputs and sampled isometry on larger
inputs; a passing sampled check is not an all-input proof.

The matched extension retains 132 verified primary timed outputs, 54 with
unsupported recognition and six blocked by timed-out warmups. Its primary
warmups comprise 44 verified, 18 unsupported and two timeouts. Old-Q-Rust
controls comprise 36 verified and 12 unsupported timed outputs, plus 12
verified and four unsupported warmups. Both Cirq dense-Clifford heavy-hex
warmups timed out; TKET completed both sessions this time. No checked output
produced a definite semantic counterexample or structural contract failure.
Unsupported outputs remain excluded, rather than being counted as correct.

The archived `analysis/reporting/compute_metrics.py` derives the following
tables from the checked records. Sessions are summarized separately using
per-case medians; repeated trials are not treated as additional circuits.

### Change against the frozen old compiler

| Correctness-qualified cohort | Cases | New/old compilation time, session 1 / 2 | Structural changes |
|---|---:|---:|---|
| External | 100 | 1.0254 / 1.0244 | Nine improve CX and depth; one trades depth metrics; 90 unchanged |
| Matched Cat/GHZ/BV | 6 | 1.0060 / 1.0043 | All unchanged |

These are descriptive geometric means of per-case median ratios. The new
passes add about 2.4–2.5% compilation time across the external cohort, not
zero overhead. Individual overhead varies: the largest session-2 increase
is DJ-5 on the complete graph, from 0.163 to 0.195 ms (19.4%). The targeted
error-correction case becomes faster: complete-graph time changes from
0.955/1.001 ms to 0.785/0.828 ms, and line time from 2.279/2.279 ms to
2.046/2.028 ms (session 1/session 2).

The ten external cases with structural changes are the error-correction,
Grover-2, HS4-4 and `iswap_n2` fixtures on both topologies, plus QEC-encoder-5
on both topologies. CX improves in nine cases and never increases; total
depth improves in ten. QEC-encoder-5 on a line retains 20 CX, reduces total
depth 21→20 and gate count 32→26, but increases two-qubit depth 16→19.
AE, QPE, Trotter, QFT, PEA and W-state fixture costs are unchanged.

### External SDK comparison

Versions are Qiskit 2.5.2, Cirq 1.7.0 and pytket 2.18.5. All 100 cases qualify
for each listed recipe in both sessions. Compilation speed is baseline time
divided by Q-Rust time; output costs are Q-Rust divided by baseline. Values
are geometric means of per-case median ratios. Structural ratios are the
same in both sessions. Two zero-CX cases are retained in outcome counts but
excluded from CX-ratio geometric means.

| Baseline recipe | Compilation speed ratio, session 1 / 2 | CX cost ratio | Total-depth ratio |
|---|---:|---:|---:|
| Qiskit primary | 3.09 / 3.10 | 0.921 | 0.876 |
| Cirq primary | 71.77 / 72.05 | 0.720 | 0.606 |
| TKET primary | 194.01 / 196.05 | 0.901 | 0.958 |
| Cirq one-pass | 43.29 / 43.12 | 0.710 | 0.602 |
| TKET light | 22.94 / 23.05 | 0.661 | 0.727 |

The unchanged external analysis also retains 2,000-draw fixture-cluster
bootstrap intervals. Primary compilation speed intervals in session 1 are
2.08–4.40 against Qiskit, 60.7–83.8 against Cirq and 162–233 against TKET;
session 2 gives 2.10–4.43, 61.1–84.1 and 164–235. These intervals reflect
variation in this selected corpus, not corpus-selection bias or hardware
portability. Recipe choice materially changes the ratios, as the light
controls show; equal repetitions do not equalize optimization effort.

Q-Rust is faster than Qiskit in 75/100 cases in session 1 and 78/100 in
session 2, and faster than each Cirq/TKET recipe in all 100 cases. It still
has more CX than primary Qiskit in seven cases and primary TKET in eight.
Against primary TKET on complete graphs alone, the CX ratio is 1.02 and
depth ratio 1.04; on lines they are 0.793 and 0.879. An overall aggregate
therefore does not describe every topology. Across primary recipes, Q-Rust
matches or improves the best observed baseline CX count in 88/100 cases;
it is within 5% in 93/100. These percentages refer to observed baselines,
not mathematical optima.

### Error-correction comparison

CX and depth are identical across sessions; times are separate session
medians in milliseconds. These are the primary recipes.

| Target | Compiler | CX | Depth | Two-qubit depth | Compile ms, session 1 / 2 |
|---|---|---:|---:|---:|---:|
| Complete | Old Q-Rust | 47 | 76 | 46 | 0.955 / 1.001 |
| Complete | New Q-Rust | 11 | 14 | 8 | 0.785 / 0.828 |
| Complete | Qiskit | 35 | 64 | 34 | 2.527 / 2.634 |
| Complete | Cirq | 35 | 69 | 34 | 95.995 / 97.161 |
| Complete | TKET | 9 | 13 | 7 | 270.608 / 277.428 |
| Line | Old Q-Rust | 65 | 86 | 56 | 2.279 / 2.279 |
| Line | New Q-Rust | 27 | 30 | 24 | 2.046 / 2.028 |
| Line | Qiskit | 51 | 72 | 43 | 3.884 / 3.827 |
| Line | Cirq | 52 | 91 | 45 | 218.000 / 220.622 |
| Line | TKET | 20 | 25 | 19 | 277.811 / 280.356 |

### Larger matched cases

The common verified four-SDK cohort contains only Cat-22 and GHZ-23 on the
square and heavy-hex graphs. On those four cases, Q-Rust's compilation speed
ratios are 4.68/4.65 against Qiskit, 29.60/29.86 against Cirq and
115.49/117.58 against TKET. Its CX cost is 1.841 times each baseline; its
depth cost is 1.764 times Qiskit/TKET and 0.902 times Cirq. This routing
weakness is unchanged by the new passes.

BV-100 qualifies for Q-Rust and TKET only. On the square graph Q-Rust has
326 CX/depth 279 versus TKET's 891/444; on heavy-hex it has 353/350 versus
2671/518. Q-Rust compiles in about 70 ms and 77–78 ms respectively, versus
TKET's 8.41–8.45 s and 6.63–6.95 s. These case-specific comparisons are
separate from the four-SDK cohort. Qiskit and Cirq BV outputs remain
unsupported by this checker and cannot be ranked here.

The two dense 20-qubit Q-Rust outputs remain unverified. Their retained
diagnostic costs show square CX 4408→4566, depth 2770→2536 and two-qubit
depth 2176→2046; heavy-hex CX 7729→7850, depth 3568→3543 and two-qubit
depth 3060→3091. These regressions are not hidden, but neither old nor new
outputs qualify for correctness-conditioned speed or quality comparisons.
The limited three-trial, single-machine matched pilot does not establish
general scalability. Native Benchpress results are a separate preserved
track with different contracts, not additional comparable repetitions.

## Scope and routing trade-offs

Logical pass acceptance does not guarantee an improvement after hardware
placement and routing. Changing interactions or intermediate permutations
changes heuristic search and later cleanup. The preflight already exposes a
QEC-encoder line case with unchanged CX count, lower total depth and greater
two-qubit depth. Both dense 20-qubit diagnostic outputs also gain CX gates;
their unsupported semantic checks preclude correctness-qualified rankings.

These trade-offs must be retained when interpreting the repeated results.
Selected structural costs are not hardware fidelity, and improvements on
these circuits do not establish a general ranking across compilers.

## Evidence and reproduction

The separate [evidence directory](external_results/20261018_clifford/) contains
the paired summary, test and replay receipts, preservation manifest and
`evidence.tar.gz`. Its `checksums.json` binds the deliverables; the archive's
`sha256_manifest.json` binds every payload file. The archive includes the raw
requests and outputs, schedules, dispatch receipts, frozen inputs, adapters,
old and new executables, changed compiler source, licenses, validation logs,
controlled ablations and earlier diagnostic attempts.
All 189 files in the historical preservation manifest are checked unchanged
before and after the campaign and during packaging.

After extracting the archive into a new directory, use Python 3.11 or newer
with `harness/experiments/requirements-comparison.txt` in an isolated environment.
From the extracted directory, reproduce the paired analysis without rerunning
compiler timings:

```sh
python harness/experiments/paired_analysis.py \
  --repo harness --campaign campaign --output fresh_analysis
python scripts/canonical_replay.py \
  --reference analysis --replay fresh_analysis --output fresh_replay_receipt.json
```

Both destination paths must be new. The replay checks recorded schedule and
source hashes before validating outputs. Canonical comparison excludes only
the value of `validation.validation_ns` and updates the four derived row-file
hashes accordingly. The timer field must still exist. Compiler, parsing and
request timings, circuit metrics, residuals, verdicts, row order and coverage
must match exactly. Raw validation outputs from both runs remain available.

The frozen executables are host-specific. A new timing campaign requires
building the recorded compiler versions for the destination host and preparing
a new protocol with `harness/experiments/paired_campaign.py`; it must not reuse
or overwrite these measurement directories. The recorded preparation and
timing commands are retained in the archive's `receipts/` directory. Analysis
replay does not execute the compiler binaries or establish cross-host timing
reproducibility.
