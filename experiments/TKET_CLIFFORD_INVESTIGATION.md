# Why TKET reaches nine CX on the five-qubit error-correction fixture

Q-Rust's current 47-CX result misses useful wire permutations encoded as ordinary gates. TKET recognizes equivalent transformations and applies Clifford interaction rewrites across the boundaries of Q-Rust's two-qubit blocks. Controlled tests account for the large reduction and isolate the final 10-to-9 CX improvement. Three-qubit synthesis does not contribute on TKET's actual optimization trajectory for this input.

**This investigation changes no production optimization. Q-Rust still emits 47 CX on the original fixture.** The 14-CX and 17-CX outputs below are diagnostic constructions using an explicitly rewritten input, not newly achieved production benchmark results. Previous benchmark files remain unchanged.

The complete diagnostic inputs, output QASM, native TKET circuits, layouts, phases, pass configurations, scripts, pinned source files and validation records are in the [evidence archive](external_results/20261018_tket_benchpress/evidence.tar.gz). This study supplements the [routing-retention comparison](EXTERNAL_RETENTION_FOLLOWUP.md).

## Input, versions and comparison rules

The case is `qasmbench__error_correctiond3_n5__all_to_all`. Its normalized input has five qubits and 114 operations: 49 CX, 62 H, two Sdg and one identity. Every operation is Clifford. This is a static unitary compilation test; its error-correction name does not mean that measurements, resets, classical feedback or noisy correction cycles were executed.

The Q-Rust executable is the preserved routing-retention binary, SHA-256 `f2bfb51d945e4a959450022a4fac4921a30d3a888d8664fb12d6df59bf689a80`, built from compiler commit `837a53ca9a9d4ca3f25cab3e5866b170e582c4a4`. Its reproduced output matches all 14 historical timed Q-Rust outputs for this case, across both sessions.

TKET is the installed `pytket==2.18.5`. The corresponding upstream tag resolves to commit `4a554cbc2699de2f95da18eb3c30ff0961bcc536`; seven implementation/binding files were checked against their Git blob identifiers in that source tree. No inference from a newer compiler version is needed.

All reported output costs use the same U3/CX basis. Depth is unit-cost ASAP circuit depth, not calibrated hardware duration. Global phase may differ, but validation uses one shared phase for the entire unitary. Initial and final maps are explicit. Timing values captured by diagnostic requests are not used for new speed claims.

## The exact TKET pass trace

Python callbacks expose the outer `FullPeepholeOptimise` standard pass, but not its internal C++ Transform steps. We reconstructed those steps using the installed public Transform bindings and the [pinned implementation](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/tket/src/Transformations/OptimisationPass.cpp#L44). The reconstructed result is an exactly equal TKET Circuit object to the actual full pass, including phase and implicit permutation.

The raw `Transform.ThreeQubitSquash` binding is used for that internal step. The similarly named public pass is a composite containing additional two-qubit and Clifford transformations, so substituting it would confound the trace. See the [pass definition](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/tket/src/Predicates/PassGenerators.cpp#L862) and [Transform binding](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/pytket/binders/transform.cpp#L375).

| Step | CX | U3/CX depth |
|---|---:|---:|
| Original input | 49 | 77 |
| Remove redundancies | 49 | 77 |
| Initial synthesis | 49 | 77 |
| Two-qubit squash, swaps disallowed | 47 | 76 |
| First Clifford simplification | 10 | 15 |
| Synthesis | 10 | 15 |
| Two-qubit squash, swaps allowed | 10 | 15 |
| Raw three-qubit squash | 10 | 15 |
| Second Clifford simplification | **9** | **13** |
| Final synthesis | **9** | **13** |

Each cost is measured from a copy lowered to U3/CX; creating that diagnostic copy does not affect subsequent steps. The nine-CX result has identity initial and final maps and identity implicit permutation. Its advantage is not an uncounted final output SWAP.

Controlled omissions distinguish causal mechanisms from pass names:

| Intervention on the original input | CX | Depth | Interpretation |
|---|---:|---:|---|
| Full pipeline | 9 | 13 | Reference result |
| Two CliffordSimp applications alone | 9 | 13 | Exactly the same final circuit |
| Omit both two-qubit squash stages | 9 | 13 | Exactly the same final circuit |
| Omit raw three-qubit squash | 9 | 13 | Exactly the same final circuit |
| Omit second Clifford simplification | 10 | 15 | Final Clifford round removes one CX |
| Omit both Clifford simplifications | 11 | 13 | Swap-aware two-qubit synthesis provides an alternative route |
| Two-qubit synthesis alone, swaps disallowed | 47 | 76 | Same CX count and depth as current Q-Rust |
| Two-qubit synthesis alone, swaps allowed | 11 | 14 | Permutation freedom matters before later cleanup |
| One CliffordSimp, swaps allowed | 10 | 15 | Large reduction without general three-qubit synthesis |
| One CliffordSimp, swaps disallowed | 34 | 63 | Disallowing swap-producing rewrites removes opportunities |
| Full pipeline, swaps disallowed | 22 | 36 | Other stages still contribute on this different trajectory |

The passes can substitute for one another. Their reductions cannot be added as independent contributions. Raw three-qubit synthesis can improve other inputs or other trajectories; the evidence here establishes that it contributes nothing on the actual reference trajectory for this fixture.

## What Q-Rust misses before routing

The input contains eleven literal seven-operation identities:

```text
CX(a,b); H(a); H(b); CX(a,b); H(a); H(b); CX(a,b)
```

The middle H-conjugated CX equals `CX(b,a)`. The resulting three-CX sequence is exactly `SWAP(a,b)`. This identity holds on arbitrary inputs, not just the all-zero state.

Q-Rust's frozen pass diagnostics show that optimization, layout and routing retain the 49 input CX, with identity maps and zero routing SWAPs. Later two-qubit block synthesis considers 14 blocks: twelve contain three CX and two contain two CX. It accepts the two 2-to-1 reductions and rejects the other twelve for failing the strict CX-improvement condition. None is rejected by the numerical guard. The result is 47 CX at depth 76.

The reconstructed block boundaries match the frozen compiler's counters. They are a read-only transcription of its collection logic, not an instrumented internal block log. Compiling every block separately with synthesis enabled and disabled confirms the same two reductions; all 28 isolated outputs pass complete four-column numerical checks.

The relevant production decisions are visible in the [pipeline order](../src/transpiler/mod.rs), [explicit-SWAP permutation extraction](../src/transpiler/optimization.rs), and [two-qubit block collection and strict CX acceptance rule](../src/transpiler/block_synthesis.rs). Explicit logical SWAP handling runs before layout/routing. Two-wire resynthesis runs later, after lowering, and does not expose a newly recognized SWAP to that earlier permutation pass.

We changed only the diagnostic input representation and ran the unchanged frozen compiler:

| Diagnostic construction | CX | Depth | Final logical-to-physical map |
|---|---:|---:|---|
| Original input | 47 | 76 | `[0,1,2,3,4]` |
| Rewrite any one of the eleven identities to named SWAP | 44 | Varies | Recorded per intervention |
| Reverse the middle CX in all eleven motifs, without naming SWAP | 47 | 54 | `[0,1,2,3,4]` |
| Name all eleven SWAPs, with permutation elision disabled | 47 | 54 | `[0,1,2,3,4]` |
| Name all eleven SWAPs, with the existing permutation handling enabled | **14** | **17** | `[0,1,4,3,2]` |
| Append three CX to restore that final map physically | **17** | **20** | `[0,1,2,3,4]` |

Each of the eleven individual rewrites independently saves three CX. All source rewrites were checked against the original complete unitary before compilation. The 17-CX construction shows that the gain is not explained solely by free output relabeling. Its restoring SWAP was appended without further optimization, so 17 is an achieved upper bound, not a minimum.

This isolates missing encoded-SWAP recognition and pass ordering as major causes. It does not establish that a literal seven-gate matcher is the best production fix. A general fix must account for equivalent local-gate forms, barriers, connectivity, numerical tolerances and interaction with the existing routing alternatives.

## Why TKET goes beyond local block synthesis

TKET's Clifford reduction is more general than this literal SWAP identity. It tracks the Pauli axes and signs of maximally entangling interactions through local Clifford gates and through other gates with which the relevant axis commutes. Matching interactions are replaced by fewer interactions plus local gates. Some replacements include a SWAP, and those SWAPs can then be absorbed into wire permutations. Disallowing swaps excludes those rewrites even when the final net permutation would be identity. See the pinned [interaction replacement rules](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/tket/src/Transformations/CliffordReductionPass.cpp#L52), [propagation and matching](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/tket/src/Transformations/CliffordReductionPass.cpp#L173), and [swap absorption](https://github.com/Quantinuum/tket/blob/4a554cbc2699de2f95da18eb3c30ff0961bcc536/tket/src/Transformations/CliffordReductionPass.cpp#L750).

The final 10-to-9 reduction can be localized further. After the first Clifford simplification, `CX(4,1)` is command 5 and `CX(1,4)` is command 23. Operations involving other qubits lie between them. An identity barrier on wires 1 and 4, placed at every cut from 6 through 23, prevents the second round from reducing below 10 CX. Barriers outside that interval retain the nine-CX result. All 26 fenced circuits remain equivalent to the original full unitary.

A separate local witness, containing those two endpoint interactions and the intervening local rotations on wires 1 and 4, reduces two CX to one. It is validated against its own full unitary. Together with the pinned implementation, the fence interventions support a Clifford interaction rewrite requiring freedom to move across intervening block boundaries. They do not expose individual C++ rewrite callbacks, and the local witness is not claimed to be the entire error-correction circuit.

## Representation sensitivity and output maps

Equivalent inputs do not guarantee the same heuristic output:

| Input supplied to TKET FullPeepholeOptimise | Final CX | Depth |
|---|---:|---:|
| Original fixture | 9 | 13 |
| Q-Rust's 49-CX output with block synthesis disabled | 9 | 13 |
| Q-Rust's 47-CX production output | 10 | 15 |
| Original fixture with the eleven identities written as SWAP | 9 | 13 |
| Q-Rust's diagnostic 14-CX output | 10 | 17 |

The 47-CX representation does not hide the large Clifford opportunity: TKET still removes 37 more CX. Its remaining one-CX difference shows sensitivity to the representation and ordering created by earlier local synthesis. This is not a correctness failure.

For the 14-CX diagnostic input, TKET introduces an implicit swap that composes with Q-Rust's `[0,1,4,3,2]` output map to give identity relative to the original fixture. The checker explicitly composes these maps. It neither ignores the permutation nor materializes it and applies it a second time.

Counting CX before final basis lowering can mislead: a standalone KAK stage on an explicitly SWAP-containing input reports 11 raw CX while retaining eleven explicit SWAP operations. Its lowered U3/CX snapshot contains 44 CX. Every table in this report uses the lowered cost.

## Correctness evidence and limits

The existing independent checker uses analytic NumPy gate matrices and imports none of Q-Rust, Qiskit, Cirq or TKET. For five qubits it compares all 32 computational-basis columns with one shared global phase. This checks the complete unitary, rather than only one state-preparation outcome.

- TKET: 51 trace/omission/representation outputs, 26 fence outputs and one separate local witness pass complete numerical validation. The 51 main outputs have maximum residual approximately `1.66e-15` against the original fixture. The witness is checked against its own input.
- Q-Rust: 26 diagnostic outputs, including the appended output-map restoration, and 14 input rewrites pass complete 32-column checks. Maximum output residual is approximately `1.85e-15`; the source-rewrite residuals are zero in these calculations.
- Isolated Q-Rust blocks: all 28 outputs pass complete four-column checks.
- Seven regression tests cover exact source reconstruction, controlled omissions, Clifford rounds, fences, permutation handling, rejection of a deleted CX, and source integrity.

Nine CX is an achieved heuristic result, **not a proven optimum**. These are numerical equivalence checks, not formal proofs. The investigation explains this fixture and these compiler versions; it does not claim the same mechanisms dominate every circuit, establish a language-level speed advantage, or measure physical hardware fidelity.

## Reproducing the TKET diagnostics

Use the comparison environment with `pytket==2.18.5` and NumPy. From the root of the extracted evidence archive, the following commands read the archived oracle, source files and input artifacts and write new results into previously nonexistent directories:

```sh
python tket/trace_tket.py \
  --repo harness \
  --fixture tket/evidence/input.qasm \
  --source-dir tket \
  --qrust-dir qrust \
  --out replay_tket_trace

python tket/trace_clifford_boundary.py \
  --repo harness \
  --fixture tket/evidence/input.qasm \
  --out replay_tket_boundaries

QRUST_REPO="$PWD/harness" \
QRUST_ERROR_CORRECTION_FIXTURE="$PWD/tket/evidence/input.qasm" \
python -m unittest discover -s tket -p 'test_trace_tket.py' -v
```

The trace script verifies the seven pinned source blobs before running, asserts exact equality with the real full pass, and validates each emitted artifact. Native Circuit JSON retains phases, barriers and implicit permutations that the static QASM snapshots represent through separate maps and phase-insensitive equivalence checks. Negative controls remove a CX or deliberately use the wrong final map and must fail.

The archive also preserves every Q-Rust diagnostic request and response, the frozen executable hash, block boundaries and isolated-block inputs. Original Q-Rust diagnostic scripts retain the workspace paths used for their run. A separate portable replay entry point accepts an explicit frozen binary and archive paths; the TKET commands above do not execute that binary.

After both timing sessions finished, the portable Q-Rust replay executed all **53 archived compiler requests**: eight configuration controls, fourteen input-representation controls, three permutation controls and twenty-eight isolated-block controls. Every emitted QASM and layout matched its preserved counterpart exactly, and the independent numerical checks passed. The appended 17-CX identity-map construction was also recreated and validated. These replay timings are not new benchmark measurements.

Set `QRUST_FROZEN_BINARY` to the preserved executable on a compatible host, then run from the extracted archive root:

```sh
python tket/replay_qrust_controls.py \
  --repo harness \
  --qrust-dir qrust \
  --binary "$QRUST_FROZEN_BINARY" \
  --out replay_qrust_controls
```

This command verifies the executable and oracle hashes before compilation, refuses to overwrite its destination, and independently validates the replayed outputs against their immediate inputs and the original fixture where applicable. A differently rebuilt executable is rejected even if it might be functionally equivalent. The replay result is retained in `qrust_portable_replay/replay_summary.json` in the evidence archive.
