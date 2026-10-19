# External comparison: frozen recipes and interpretation

These are the six measured recipes in the [external study](EXTERNAL_BENCHMARKS.md). They describe specific configurations, rather than every capability or best configuration of their SDKs. The compiler, adapters and numerical thresholds were frozen before timing.

## Common contract

Every recipe receives the same normalized static unitary circuit and targets arbitrary one-qubit U3 gates plus CX. Each target has exactly as many physical as logical qubits. Both the complete graph and the bidirectional line are ideal connectivity models; they carry no measured gate error rates, durations or device calibration data. The native Q-Rust worker receives explicit edges even for the complete graph. The Python adapters use their unconstrained compilation paths for that topology.

Initial and final logical wire permutations are free and are retained separately with every serialized physical circuit. A compiler therefore need not spend gates to restore the original wire order. Correctness is relative to those recorded mappings. Comparing only emitted QASM without its mappings would be incorrect, particularly for TKET's implicit swaps and Q-Rust's output-permutation elimination.

The reference is the retained original unitary prefix, after the corpus frontend has explicitly accounted for terminal measurements, barriers, register canonicalization and gate-definition expansion. The shared CRZ/RZZ normalization happens before timing. This is compilation of the whole unitary on arbitrary logical inputs, rather than testing only the algorithm's usual all-zero initial state. The comparison does not measure dynamic circuits, measurement handling or hardware execution fidelity.

Equivalence permits one whole-circuit global phase. It does not permit a different phase for each input state or basis column. The output serializers omit standalone global phases; they retain relative and controlled phases. Consequently, these checks establish the stated unitary-equivalence convention, not a guarantee that a phase-discarding serialized circuit can later be embedded under an additional quantum control without recovering its phase.

## Six recipes

| Raw recipe ID | Timed compilation configuration | Prespecified difference |
|---|---|---|
| `qrust` | Public `transpile_with_report`, optimization level 3, backend connectivity, target `u/cx`, default routing search, phase-gadget recovery, two-qubit block synthesis and post-block fusion enabled; pass profiling disabled. | One Q-Rust configuration is shared by both analysis profiles. |
| `qiskit` | `transpile(..., basis_gates=["u", "cx"], optimization_level=3, approximation_degree=1.0, num_processes=1, qubits_initially_zero=False, seed_transpiler=trial_seed)`. A bidirectional coupling map is supplied for the line; `coupling_map=None` for the complete graph. | One Qiskit configuration is shared by both profiles. The requested seed is consumed. |
| `cirq` | `CZTargetGateset(atol=1e-10, allow_partial_czs=False, preserve_moment_structure=False, reorder_operations=True)`; target synthesis with `max_num_passes=2`, `ignore_failures=False`; line routing with `RouteCQC`, followed by another target-synthesis call with the same maximum; full CZ-to-CX conversion and one-qubit fusion. | Primary Cirq recipe. There is no invented Cirq “O3” preset. |
| `tket` | `FullPeepholeOptimise(allow_swaps=True, target_2qb_gate=CX)`; on the line, `DefaultMappingPass`, `DecomposeSwapsToCXs(architecture)` and `SynthesiseTket`; finally `AutoRebase({CX,U3})`, applied as one `SequencePass`. | Primary TKET recipe includes two- and three-qubit peephole resynthesis. |
| `cirq_one_pass` | Same CZ target, tolerance, routing and CX closure as primary Cirq, with `max_num_passes=1` for each synthesis call. | Changes the maximum synthesis-iteration count. It does not remove two-qubit synthesis or bypass routing. |
| `tket_light` | Replace the initial `FullPeepholeOptimise` with `SynthesiseTket`; retain line mapping, SWAP decomposition, post-routing `SynthesiseTket` and final `AutoRebase({CX,U3})`. | Changes the initial optimization recipe, not the hardware target or final output basis. |

All six recipes are interleaved within the same two sessions. The primary and lighter profiles reuse exactly the same Q-Rust and Qiskit observations. They do not produce extra independent Q-Rust/Qiskit trials. Seven trial seeds are repeated across the two sessions; Qiskit consumes those seeds, while the other five recipes report no dependence on the requested trial seed. Q-Rust still performs deterministic internal multistart layout search. Repeated deterministic recipes measure run-to-run timing variation, not new randomized searches.

## Q-Rust: bounded search and local replacements

The frozen O3 pipeline first recovers exact parity-rotation patterns, crystallizes recognized gates, simplifies parameters, merges rotations, applies cross-conjugation and inverse/commutation cancellations, and simplifies again. Before routing it extracts output permutations and decomposes gates with more than two operands. It then uses SABRE layout with up to 10 starts and 3 forward/backward refinement iterations per start. The separate beam router has beam width 8, branch factor 5 and 4 bidirectional iterations. The lookahead size is 20; the default static lookahead weight is 0.5. Candidate progress ranking, exact diagonal commutation, SWAP orientation and output-permutation elimination are enabled. The trial seed supplied by the outer harness does not change these deterministic choices.

These are configured search bounds, not a guarantee that all searches execute. For example, the layout implementation can stop on a zero-SWAP lower bound. Circuit structure and connectivity can therefore change how much work the same nominal configuration performs.

After routing, Q-Rust lowers to its target basis, performs KAK synthesis where applicable, fuses gates, simplifies SWAP-related sequences, cancels inverses, and runs two-qubit block consolidation. A consolidated block follows two wires and ends at an interaction with a third wire. Disjoint U/CX operations may be crossed because they commute with the entire block. The pass operates on the already routed sequence and does not redesign the logical placement. It accepts a replacement only when it strictly reduces the CX count and passes reconstructed matrix-residual checks. Pure CX blocks have a small exact-word fast path. Post-block fusion then combines adjacent one-qubit runs while respecting each wire's dependencies. Direction repair and final basis closure follow.

The strict CX-improvement rule is an explicit objective choice: a same-CX replacement that would improve single-qubit count or depth alone is not accepted by this block pass. Likewise, ending a block at a third-wire interaction deliberately bounds synthesis cost and limits which cross-boundary opportunities it can see. Those source facts identify plausible mechanisms for a quality trade-off. They do not, by themselves, prove the cause of a particular benchmark gap.

## Qiskit: the pinned preset, rather than an assumed historical pass list

The adapter delegates to the pinned Qiskit 2.5.2 level-3 preset; it does not construct a custom pass list. The installed implementation has a conditional VF2 layout search and a SABRE fallback on constrained targets. The SABRE layout fallback is configured with four refinement iterations and `_get_trial_count(20)` for both layout and SWAP trials. The default routing plugin delegates to the SABRE routing implementation, using the decay heuristic; when layout has already routed the circuit, the routing stage need not repeat that work. The complete-graph call omits layout and routing stages when there is no coupling map or initial layout.

The pinned default optimization loop contains `TwoQubitPeepholeOptimization`, `RemoveIdentityEquivalent`, `Optimize1qGatesDecomposition`, `CommutativeCancellation`, control-flow idle-wire handling, and conditional basis translation. Its level-3 stopping rule tracks a minimum point in circuit depth and size. This should not be described using an older Qiskit release's `Collect2qBlocks`/`ConsolidateBlocks` loop as if that were the measured code.

The installed two-qubit peephole pass documents its preference order as fewer two-qubit gates, then estimated fidelity, then fewer total gates. That is a broader tie-breaking objective than Q-Rust's strict-CX-only block admission. The benchmark supplies an ideal U/CX target rather than a calibration model, so it does not establish a measured hardware-fidelity benefit from this heuristic.

`approximation_degree=1.0` requests no intentional approximation, and `qubits_initially_zero=False` prevents relying on all-zero logical initialization. Neither setting is a floating-point correctness proof. The unchanged external oracle remains the common acceptance check.

The exact SABRE trial-count helper can be affected by `QISKIT_SABRE_ALL_THREADS` or the user configuration's `sabre_all_threads`. A separate read-only receipt, captured during session 2, found both environment overrides absent and the default `~/.qiskit/settings.conf` absent. This is consistent with the configured default of 20; it is supplemental provenance rather than a snapshot of live worker memory or proof of earlier configuration contents. The harness's frozen worker overrides separately set `QISKIT_NUM_PROCS=1`, `QISKIT_PARALLEL=FALSE` and `RAYON_NUM_THREADS=1`.

## Cirq: target conversion and iteration sensitivity

Cirq synthesizes into arbitrary one-qubit gates and **full** CZ gates. `allow_partial_czs=False` excludes fractional CZs from the chosen output gate set. A full CZ is changed to CX by conjugating its target with Hadamards: `CZ = (I ⊗ H) CX (I ⊗ H)`. Thus this conversion retains one entangling gate per full CZ. The resulting one-qubit gates are fused with `merge_single_qubit_gates_to_phxz(atol=1e-10)` inside the compile timer.

For a line, the pinned `RouteCQC` call uses its default lookahead radius of 8 and `LineInitialMapper`; the adapter then completes mappings for inactive logical wires. The second synthesis phase may jointly simplify routing SWAPs with neighboring same-pair gates. The final map composes the initial placement with the router's returned permutation. All-to-all skips the routing and second synthesis phase.

The conversion is a valid way to reach the shared target, but it is a chosen pipeline. It does not demonstrate that every Cirq workflow should synthesize through CZ before producing CX. The intermediate representation, one-qubit fusion and iteration count can change runtime and block boundaries. The lighter recipe tests the maximum-iteration choice within this pipeline. Because the setting is a maximum, reducing it from two to one need not change every circuit or halve actual synthesis work.

The serializer converts each remaining one-qubit unitary to one U3 representation, and writes one CX per compiled CX. It neither fuses additional gates nor changes the two-qubit sequence. One-qubit matrix-to-Euler conversion is outside `compile_ns`, while CZ-to-CX lowering and subsequent fusion are inside. This representational boundary is disclosed because it differs from directly compiling to an internal U3 representation.

## TKET: richer initial synthesis and implicit permutation handling

The installed `FullPeepholeOptimise` API describes two- and three-qubit sequence resynthesis. `allow_swaps=True` permits implicit wire swaps. The line recipe subsequently uses `DefaultMappingPass`, whose documented placement is `GraphPlacement` and whose architecture handling ignores edge direction. SWAPs are lowered to CX, and `SynthesiseTket` provides post-routing optimization before the final U3/CX rebase. The complete-graph recipe skips the mapping and post-routing section.

`tket_light` replaces the initial full peephole pass with `SynthesiseTket`, which optimizes into CX, TK1 and Phase before the later common rebase. The difference therefore tests a substantial initial synthesis-policy choice. It is not an equal-time-budget comparison and does not establish that either policy is universally preferable. The chosen API calls do not expose a stochastic seed in this adapter, and no explicit numerical tolerance is overridden; internal SDK numerical policies remain those of the pinned version.

The output serializer walks the actual measured commands and converts TKET angle units to radians. It avoids a QASM exporter that might materialize implicit wire permutations. Specifically, final logical locations are computed by applying the inverse implicit output permutation to `CompilationUnit.final_map`. This prevents counting the same permutation twice. Standalone `Phase` commands are omitted under the common global-phase convention.

## Numerical policies are not interchangeable error metrics

| Component | Explicit policy in the measured implementation | Interpretation |
|---|---|---|
| Q-Rust generic guarded local rewrite | Phase-aligned unnormalized Frobenius residual at most `1e-12`, also charged to a shared allowance. | Local acceptance policy, not measured hardware fidelity. |
| Q-Rust two-qubit block rewrite | Residual at most `1e-10`, block-pass residual sum at most `1e-9`, plus shared allowance admission. | A stricter CX objective does not imply a formal whole-circuit proof. |
| Q-Rust post-block one-qubit fusion | Run residual at most `1e-12`; identity removal at most `1e-14`; pass residual sum at most `1e-10`; shared allowance still applies. | Tiny gates can be retained when deletion/fusion is not admissible. |
| Q-Rust shared numerical allowance | `1e-9` across pipeline stages; charges the measured residual plus `32 × machine_epsilon × evaluated_operations`. | Conservative engineering guard. It does not certify all arithmetic, analytic decompositions or transcendental functions. |
| Cirq primary and lighter recipes | `atol=1e-10` for target synthesis and final one-qubit fusion. | A local synthesis/fusion tolerance, not a whole-circuit guarantee. |
| Qiskit recipe | `approximation_degree=1.0`; other numerical thresholds are the pinned SDK defaults. | Requests exact synthesis intent, subject to floating-point arithmetic. |
| TKET recipes | No adapter-level tolerance override. | Uses the pinned SDK's internal policies. |
| Independent comparison oracle | Residual threshold `1e-7`, with one shared phase across all checked inputs. | Full numerical isometry through six physical qubits; eight sampled inputs above six. Sampling is finite evidence, not proof. |

The numerical values above do not all measure the same quantity, so their magnitudes alone cannot rank compiler accuracy. All compared outputs must satisfy the same independent oracle and layout/connectivity checks. A failed or unchecked result must remain visible in coverage and cannot contribute to a fully successful paired cell.

## Timing boundaries

`compile_ns` measures warmed in-process compilation, including circuit copying where the adapter performs it, target/pass configuration, optimization, synthesis, mapping/routing and compiled-circuit retrieval. Q-Rust's `transpile_with_report` work is included; its optional per-pass clocks are disabled. Parsing is recorded separately as `parse_ns`. The recipes retain their existing public adapters; the small pre-timer native request/search-field validation is not silently treated as pass execution.

Imports, process startup, source parsing, output serialization, JSON transport and independent checking are outside `compile_ns`. `request_wall_ns` records the whole worker request/response. `noncompile_nonparse_request_ns` is that wall time minus parse and compile, and combines serialization, layout extraction, diagnostics, transport and other overhead. It is **not** an isolated serialization timer. Parse-plus-compile is a secondary reported boundary and still does not include the full request overhead.

The lighter Cirq and TKET recipes reuse diagnostic functions with small stage-clock and dictionary bookkeeping costs inside their compile timers. Those costs are retained without a speculative subtraction. TKET's lighter measurement uses `SequencePass` with `stage_probe=False`; it does not time separately applied passes as though that were the same API boundary.

Each recipe has its own persistent worker, one initial warmup per circuit/target cell and a common 60-second request deadline. Requests are sequential and their order is frozen before timing. A worker restarted after a timeout must complete a replacement warmup before a measured request is allowed; a failed replacement warmup is retained as a failed scheduled trial. No validation or other benchmark session is permitted concurrently with the timed sessions. Thread environment caps reduce concurrency, but do not establish CPU affinity, fixed frequency or freedom from all operating-system interference.

## What the differences can and cannot explain

Established implementation differences include bounded versus iterative search, conditional layout shortcuts, two-wire versus three-wire synthesis scope, strict versus multi-objective block acceptance, intermediate CZ representation, target-dependent routing, and free output permutations. The primary/lighter Cirq and TKET pairs directly vary two declared recipe choices within the same sessions. They can measure how sensitive this corpus is to those choices while holding the other measured recipes and target contract fixed.

It remains a hypothesis that a particular mechanism causes a particular runtime, CX or depth difference until matching per-circuit evidence supports it. For example, more search can spend additional compilation time without improving an already easy layout; three-qubit synthesis can find transformations unavailable to a two-wire pass; a strict CX-reduction objective can retain a deeper same-CX construction. These possibilities should not be presented as diagnoses of every observed gap. The study does not hold algorithms and implementation languages constant, so it cannot separate language cost from algorithm choice, search effort, data structures, numerical work or integration overhead.

In particular, a Q-Rust runtime advantage would not measure “Rust versus Python.” The pinned Qiskit SABRE and two-qubit peephole implementations call compiled accelerator components, and TKET's passes also run through native components. A causal language comparison would require substantially matched algorithms, data structures, work budgets and timing boundaries. This study measures complete specified tool recipes on a declared corpus.

The defensible output-quality comparison is against the observed outputs of the other recipes under this contract. It does not establish a globally optimal circuit, and a result within a percentage of the best observed baseline is not that percentage of a proven optimum. CX count, unit-duration circuit depth and two-qubit depth characterize structural costs; translating them into execution time or error probability would require a separate hardware/duration/noise model.

## Source anchors

The repository adapters and compiler source are frozen by the study protocol. Source inspection used these files:

- Q-Rust worker and default pass selection: `experiments/src/bin/comparison_runner.rs`; `src/transpiler/mod.rs`, especially `RoutingSearchConfig::default` and `build_pass_manager_for`.
- Q-Rust block boundaries, objective and numerical policy: `src/transpiler/block_synthesis.rs`; `src/transpiler/local_cleanup.rs`; `src/transpiler/numerical_policy.rs`; layout shortcut behavior in `src/transpiler/layout.rs`.
- Primary Python recipes and serialization: `experiments/comparison_qiskit_tket.py`; `experiments/comparison_cirq.py`.
- Lighter recipes: `experiments/diagnose_compiler_tradeoffs.py`, functions `cirq_trial` and `tket_trial`.
- Session, seed, worker and timing contract: `experiments/external_comparison.py`.
- Independent phase, input coverage and metric definitions: `experiments/comparison_validation.py`.
- Installed Qiskit 2.5.2: `qiskit/transpiler/preset_passmanagers/level3.py`; `builtin_plugins.py`, especially `OptimizationPassManager`, `DefaultLayoutPassManager`, `SabreSwapPassManager` and `_get_trial_count`; `qiskit/transpiler/passes/optimization/two_qubit_peephole.py`; `qiskit/user_config.py`.
- Installed Cirq 1.7.0: `cirq/transformers/routing/route_circuit_cqc.py`, `RouteCQC.route_circuit` defaults.
- Installed pytket 2.18.5: `pytket/_tket/passes.pyi`, documentation for `FullPeepholeOptimise`, `SynthesiseTket` and `DefaultMappingPass`.
- Supplemental, explicitly mid-session configuration receipt: `qiskit_environment_config_receipt.json` in the evidence archive. It includes installed Qiskit source hashes and line references and records only relevant environment/configuration values.
