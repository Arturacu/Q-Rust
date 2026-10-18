# Benchpress execution pilot

This pilot executes a pinned portion of Benchpress and extends the comparison to 20–100 logical qubits, square/heavy-hex connectivity, spare physical qubits, and an archived 133-qubit device target. It records compilation time, output structure and sampled process memory. A separate duration calculation uses the device snapshot's instruction timings. Existing benchmark results remain preserved.

There are two deliberately separate execution tracks. Native Benchpress exercises its existing SDK implementations and defaults. The matched extension feeds identical unitary inputs and graphs to the existing Q-Rust, Qiskit, Cirq and TKET adapters. Their gate bases, initialization assumptions and timed work differ, so results from one track cannot be used as a baseline for the other.

## Source and scope

Benchpress is an open-source framework with a peer-reviewed methodology for circuit construction, manipulation and compilation. Its workloads include QASMBench, Feynman and Hamiltonian collections, with both generated connectivity and device snapshots. It provides an established comparison framework; executing a subset is not industry-standard certification. The published methodology does not execute output circuits on quantum hardware. See the [Benchpress study](https://arxiv.org/html/2409.08844v2) and [journal article](https://doi.org/10.1038/s43588-025-00792-y).

The source revision is [`b695f30e83a32bac05b9b4d8e98d37ba9aae5236`](https://github.com/Qiskit/benchpress/tree/b695f30e83a32bac05b9b4d8e98d37ba9aae5236). All 1,069 source blobs were checked against the previously pinned Git tree before execution. The downloaded source archive SHA-256 is `28d37b4230976a97464c90599549598b093a50c7a776a73c2ae2ee5e8e4a0ec8`. The native workouts were not edited.

The machine is an Apple M3 with 8 GiB memory and macOS 14.4.1. SDK versions are Qiskit 2.5.2, Cirq Core 1.7.0 and pytket 2.18.5. Native execution adds qiskit-ibm-runtime 0.50.0, pytket-qiskit 0.79.0, pytest 9.1.1 and Benchpress's pytest-benchmark fork at `e33f5934cb462e719b13b8d352ced107f575b4f7`. These extra dependencies use a separate environment; the previous comparison environment was preserved.

## Native Benchpress results

Five selected workouts were run for each of Qiskit, TKET and Cirq in two sessions. Each supported test has three measured rounds of one iteration, with a warmup. The unmodified benchmark fixture also makes calibration calls and a final invocation whose returned circuit supplies the reported output properties.

The completed campaign contains **20 passing test executions, 10 expected skips and no failures**, with **60 measured compilations**. Qiskit and TKET implement all five selected workouts. The pinned Cirq gym leaves these abstract/device transpilation workouts unimplemented, so its ten executions skip. This is a limitation of that harness integration, not a demonstration that Cirq cannot route circuits.

Q-Rust has no native gym in this pinned Benchpress version. Its results belong
to the separate matched adapter extension below.

Times below are medians of three trials, with the two sessions kept separate. CZ counts and CZ-only depths describe the final returned invocation; they agree between sessions. They do not describe every stochastic timed output, which native Benchpress does not retain.

| Native workload | SDK | Compile ms, S1 / S2 | CZ | CZ depth |
|---|---|---:|---:|---:|
| Cat-22, square | Qiskit | 2.181 / 2.269 | 21 | 21 |
| Cat-22, square | TKET | 359.877 / 342.254 | 30 | 21 |
| Cat-22, heavy-hex | Qiskit | 2.610 / 2.387 | 21 | 21 |
| Cat-22, heavy-hex | TKET | 451.848 / 442.299 | 116 | 90 |
| GHZ-23, square | Qiskit | 2.312 / 2.238 | 22 | 22 |
| GHZ-23, square | TKET | 404.706 / 384.290 | 31 | 28 |
| GHZ-23, heavy-hex | Qiskit | 2.466 / 2.370 | 22 | 22 |
| GHZ-23, heavy-hex | TKET | 447.826 / 459.180 | 64 | 61 |
| BV-100, FakeTorino | Qiskit | 34.951 / 35.719 | 200 | 200 |
| BV-100, FakeTorino | TKET | 1,967.593 / 1,967.676 | 576 | 478 |

These are the native level-2 recipes and `id/sx/x/rz/cz` target basis. Native timing excludes pass-manager construction; TKET's timed function copies the circuit before applying its pass. The matched extension below includes configuration/pass construction in its compile timer. Different timing boundaries and optimization effort prevent direct speed ratios between these tracks.

Native runs completed without timeout or memory-limit termination. The largest sampled pytest-parent RSS was **281,657,344 bytes (268.6 MiB)**. The monitor sampled every 250 ms with a 3.5 GiB termination threshold. This is an observed process high-water value, not the exact allocation peak or the memory attributable solely to a compiler. Tests had a 90-second limit, and each gym process had a 400-second wall-time limit. Compiler work ran serially with recorded thread caps; ordinary operating-system activity and CPU frequency were not controlled.

## Device duration and native correctness contract

The FakeTorino configuration and property files are archived with hashes and provider version. Qiskit's returned BV-100 circuit has a **20.184 microsecond weighted ASAP duration** in both sessions. This calculation uses the snapshot's duration for each emitted instruction, respects qubit dependencies and barriers, and includes measurement. Every required instruction duration was present. It excludes alignment constraints, crosstalk, queueing, repeated shots and classical processing. It is a model prediction, not a measurement of hardware runtime or fidelity.

Upstream native validation checks target gates and connectivity. The additional capture plugin records the returned circuit outside its compile timer. Qiskit's logical maps come from that artifact. TKET's replay through a `CompilationUnit` did not reproduce an identical circuit dictionary, so the replay maps were not attached to its native artifact.

Native workloads permit initialized-state and measurement semantics. In particular, TKET's BV output measures some physical wires and subsequently uses them again. Removing measurements and treating its final wires as a unitary output would change that contract. The preserved `native_dict` contains the full ordered sequence; the measurement-stripped `out_qasm` alone cannot establish its correctness.

**All 20 retained native outputs pass independent joint-readout distribution checks for the specified all-zero input.** The check processes the complete ordered `native_qasm` or `native_dict`, copying each measurement coherently into a fresh record qubit. It then compares signed Z-parity constraints on the last record written to each classical label, tracing out physical wires and older records. This preserves measurement reuse and classical overwrites without relying on unavailable TKET maps. Reset and classical feedback are unsupported and fail closed. Unused classical declarations are excluded from the observed-output contract.

The retained Cat-22 and GHZ-23 outputs have the expected equal mixture of the all-zero and all-one strings. Each BV-100 output deterministically returns one on its 99 measured bits. These checks compare the complete joint distribution, including correlations, rather than individual-bit marginals. The maximum summed local Clifford-recognition residual is `6.925e-13`, below the `1e-7` budget. Signed-Pauli operations are exact integer calculations after numerical gate recognition; the recognition bound is not a certified floating-point proof.

The checker has 12 passing regression tests, including 100 random small-circuit comparisons against independent dense projective-measurement simulation. The archived receipt is `verification/native_readout_validation.json`. This evidence establishes the zero-input classical-readout contract only: it does not establish arbitrary-input unitary/isometry equivalence, quantum-output-state equality, hardware fidelity, or correctness of uncaptured timed invocations.

## Matched four-tool extension

This track is a **Benchpress-sourced adapter extension**, not a native Benchpress score. Four input circuits are frozen before timing, each on two generated topologies:

| Source circuit | Logical qubits | Square physical qubits / edges | Heavy-hex physical qubits / edges |
|---|---:|---:|---:|
| Cat-22 | 22 | 25 / 40 | 57 / 64 |
| GHZ-23 | 23 | 25 / 40 | 57 / 64 |
| Random Clifford-20 | 20 | 25 / 40 | 57 / 64 |
| BV-100 | 100 | 100 / 180 | 115 / 132 |

Edges in this table are undirected; all adapters receive bidirectional connectivity. Graphs come from the pinned `FlexibleBackend` implementation. Only global terminal measurement suffixes and barriers are removed from the source inputs; primitive definitions are checked before normalization. No source optimization is performed. All adapters receive the same normalized unitary and use U3/CX output, clean spare wires and explicit initial/final logical maps. Qiskit's matched recipe disables zero-input-only simplifications.

The recipes are Q-Rust O3, Qiskit optimization level 3, the corrected Cirq CX
adapter with at most two synthesis passes, and TKET's full peephole pipeline
before mapping. They are frozen recipes with unequal search effort. The Q-Rust
binary is unchanged from the retention study (SHA-256
`f2bfb51d945e4a959450022a4fac4921a30d3a888d8664fb12d6df59bf689a80`,
compiler source commit `837a53ca9a9d4ca3f25cab3e5866b170e582c4a4`). The complete
matched protocol SHA-256 is
`486d638da1ec3088aebcdeb8158951148cc0d64adf60a4b3492dacee3afc3682`.

The frozen protocol schedules **192 timed requests and 64 warmups**: eight circuit/topology cases, four tools, three trials and two sessions. Each request has a 30-second limit. Known worker processes are sampled approximately every 100 ms with limits of 2 GiB per worker and 3 GiB in aggregate. These limits include idle resident workers, may be overshot between observations and do not measure descendant processes. Every scheduled outcome, including failures, blocked trials and timeouts, is retained. Validation and profiling run after both timing sessions.

The independent checker uses signed Pauli arithmetic and a Choi-state representation to compare the logical isometry, including spare-wire initialization, without allocating a full statevector. Individual output gates must first match allowed Clifford operations within a recorded numerical bound. Non-Clifford outputs fail closed. Basis and connectivity checks are separate. This contract covers every logical input for accepted Clifford outputs; it is not a general-purpose proof for arbitrary circuits.

The mathematical basis is the polynomial-space stabilizer representation of
Clifford circuits, described by [Aaronson and Gottesman](https://arxiv.org/abs/quant-ph/0406196).
Each logical input is entangled with a reference qubit, while spare physical
wires begin in zero. Equality of the resulting normalized Choi states, up to
one global phase, is equivalent to equality of the isometries on every logical
input. The checker compares their canonical signed stabilizer groups. Its
local numerical recognition residuals are summed: the triangle inequality for
products of unitaries bounds the accumulated operator error by that sum.
The arithmetic on recognized Pauli groups is exact; floating-point gate
recognition is numerical evidence, not a formal floating-point certificate.

The original, frozen gatewise checker certifies **129 of 192 scheduled timed outputs**, plus **43 warmups**. A further **54 timed outputs and 18 warmups compile but remain unverified** because their decompositions contain individual non-Clifford gates. A product of non-Clifford gates can still implement a Clifford operation, so this is a representation limit of the checker, not an observed semantic mismatch. No numerical mismatch was established for these outputs.

Three warmups hit the 30-second limit: Cirq on dense Clifford-20/heavy-hex in both sessions and TKET on that case in session 1. Those failures block nine scheduled timed requests under the frozen policy. TKET completes the same case in session 2: median 22,112.802 ms, 8,723 CX, total depth 3,478 and CX depth 3,248; all its session-2 outputs pass. This is retained as a session-dependent completion result, not a complete two-session comparison.

The following table reports both-session eligible results under the original checker. Each time is the median of three trials. All eligible rows have identical structural counts across their six timed outputs. Dashes retain the unverified or incomplete cases; their actual output metrics and timings remain in the raw records. Minima and maxima, parsing time, request-wall time and resource samples are in `matched_summary.json` and the raw records.

| Matched case | Tool | Compile ms, S1 / S2 | CX | Total depth | CX depth | Original check |
|---|---|---:|---:|---:|---:|---|
| cat22 square | Q-Rust | 1.406 / 1.411 | 28 | 29 | 28 | Pass |
| cat22 square | Qiskit | 16.600 / 16.670 | 21 | 22 | 21 | Pass |
| cat22 square | Cirq | 81.076 / 80.609 | 21 | 43 | 21 | Pass |
| cat22 square | TKET | 316.797 / 313.210 | 21 | 22 | 21 | Pass |
| cat22 heavy-hex | Q-Rust | 4.148 / 4.128 | 51 | 52 | 51 | Pass |
| cat22 heavy-hex | Qiskit | 10.611 / 10.000 | 21 | 22 | 21 | Pass |
| cat22 heavy-hex | Cirq | 81.531 / 81.365 | 21 | 43 | 21 | Pass |
| cat22 heavy-hex | TKET | 311.099 / 308.334 | 21 | 22 | 21 | Pass |
| ghz23 square | Q-Rust | 2.439 / 2.408 | 33 | 31 | 30 | Pass |
| ghz23 square | Qiskit | 16.678 / 16.486 | 22 | 23 | 22 | Pass |
| ghz23 square | Cirq | 85.367 / 84.347 | 22 | 45 | 22 | Pass |
| ghz23 square | TKET | 356.917 / 355.813 | 22 | 23 | 22 | Pass |
| ghz23 heavy-hex | Q-Rust | 4.415 / 4.369 | 52 | 53 | 52 | Pass |
| ghz23 heavy-hex | Qiskit | 11.213 / 11.128 | 22 | 23 | 22 | Pass |
| ghz23 heavy-hex | Cirq | 86.676 / 84.598 | 22 | 45 | 22 | Pass |
| ghz23 heavy-hex | TKET | 337.228 / 334.896 | 22 | 23 | 22 | Pass |
| clifford20 square | Q-Rust | — | — | — | — | Unverified decomposition |
| clifford20 square | Qiskit | — | — | — | — | Unverified decomposition |
| clifford20 square | Cirq | — | — | — | — | Unverified decomposition |
| clifford20 square | TKET | 21,950.229 / 21,833.859 | 4813 | 2271 | 2025 | Pass |
| clifford20 heavy-hex | Q-Rust | — | — | — | — | Unverified decomposition |
| clifford20 heavy-hex | Qiskit | — | — | — | — | Unverified decomposition |
| clifford20 heavy-hex | Cirq | — | — | — | — | Warmup timeouts |
| clifford20 heavy-hex | TKET | — | — | — | — | S1 timeout; S2 pass |
| bv100 square | Q-Rust | 70.013 / 69.362 | 326 | 279 | 265 | Pass |
| bv100 square | Qiskit | — | — | — | — | Unverified decomposition |
| bv100 square | Cirq | — | — | — | — | Unverified decomposition |
| bv100 square | TKET | 8,653.649 / 8,399.150 | 891 | 444 | 427 | Pass |
| bv100 heavy-hex | Q-Rust | 75.355 / 76.575 | 353 | 350 | 336 | Pass |
| bv100 heavy-hex | Qiskit | — | — | — | — | Unverified decomposition |
| bv100 heavy-hex | Cirq | — | — | — | — | Unverified decomposition |
| bv100 heavy-hex | TKET | 9,553.351 / 7,520.405 | 2671 | 518 | 455 | Pass |

Q-Rust is faster on the jointly validated Cat/GHZ cases, but uses 28–52 CX where the other tools use 21–22. On the validated BV-100 comparison with TKET, Q-Rust emits 326/353 CX on square/heavy-hex versus 891/2,671, with lower depth and compilation time. This does not establish a Q-Rust advantage over the two unverified BV decompositions or an overall SDK ranking. The recipes spend unequal optimization effort, and three trials per session support only a bounded feasibility comparison.

Peak sampled worker RSS is 73.3 MiB for Q-Rust, 76.8 MiB for Qiskit, 399.0 MiB for Cirq and 105.1 MiB for TKET; maximum sampled aggregate worker RSS is 440.2 MiB. No memory limit is reached. These include runtime/process overhead and are not isolated allocator measurements.

The original checker was cross-checked against an independent dense simulator on 532 small isometry cases (266 positive and 266 negative), including arbitrary mappings and spare wires, plus 40 generator-basis changes. Its 13 regression tests cover gate recognition, relative phase, accumulated error, maps, basis, connectivity and clean/dirty ancillas. The analysis replay separately checks all 256 scheduled outcomes, frozen source/input/output hashes, eligibility and every per-cell aggregate; its corruption tests prevent failed or unverified outputs from silently entering comparisons.

A separate post-hoc diagnostic tested whether combining gates into strict one- or two-wire blocks could resolve the representation limit. It covered all **17 unique unverified artifacts**, representing the **54 timed outputs and 18 warmups** above. Of 28,962 blocks, 4,558 failed Clifford recognition; no artifact had every block recognized. Candidate Clifford transformations were reconstructed from signed Pauli-generator images and accepted only after a direct full-unitary Frobenius comparison at the unchanged `1e-10` local threshold. Both preserved prototype versions agree on every block assignment and acceptance decision. This bounded partition therefore does not establish global equivalence for these outputs; its failures are not evidence that the circuits are incorrect. No larger-block search or production change was attempted, and all original results and verification statuses remain unchanged. The evidence archive retains both scans and `block_probe_receipt.json`.

## What this expansion establishes

Larger widths can be evaluated reliably when the input family and validation method are stated. This pilot broadens the earlier complete/linear targets, introduces a native CZ basis and a device snapshot, and records memory and a modeled duration. It does not establish performance across arbitrary large quantum programs or hardware fidelity.

The timing differences combine algorithm choice, search effort, implementation
and API overhead. This experiment does not isolate a causal benefit from the
Rust language itself. Likewise, a low entangling-gate count does not by itself
establish lower physical error, shorter calibrated duration or globally optimal
compilation.

Cat-22 and GHZ-23 are both chain GHZ preparations; two sizes do not provide two substantially different application families. Random Clifford-20 and BV-100 add different interaction patterns, but all four matched inputs are Clifford. They were selected to make scalable semantic checking feasible. No performance inference should be extended to non-Clifford programs, larger dense workloads, other machines or untested configurations.

The [upstream resource guidance](https://github.com/Qiskit/benchpress#testing-resource-requirements) recommends 96+ GB for the full multi-SDK suite and notes that serial execution can take a week or more. A bounded subset fits this host; this campaign is not a complete-suite run. The supplied 100-qubit Clifford file contains 9,897 gates and is outside the predeclared pilot selection; it was not run or demonstrated infeasible. Quantum hardware execution would require a separate backend, calibration and sampling protocol.

## Evidence and reproduction

The [evidence archive](external_results/20261018_tket_benchpress/evidence.tar.gz) contains the native protocols and raw outputs under `benchpress/`, matched records under `matched/`, and frozen experiment code under `harness/experiments/`. It includes environment and source manifests, target snapshots, input QASM, timing samples and captured artifacts. The earlier campaigns retain their own unchanged result directories.

The original matched analysis can be replayed from the extracted archive without
running any compiler or importing a compiler SDK. Use the recorded NumPy version
and run this from the archive root:

```sh
python harness/experiments/audit_benchpress_pilot.py --evidence . --harness harness --output matched_archive_replay.json
```

This checks the original protocol and session receipts, every retained outcome,
the semantic verdicts and all per-case aggregates. It preserves the original
unverified and timeout statuses. The native analysis was also recomputed into a
new file and reproduced the original summary byte-for-byte.

The complete Python regression suite passes **227 tests**, including the new
Clifford-isometry, native-readout and archive-corruption regressions. The
five-qubit causal investigation has seven additional TKET regression tests and
a successful portable replay of 53 Q-Rust compiler requests plus the derived
restored-map output. Production compiler code is unchanged. All 64 previously
tracked experiment reports and external-result files match the pre-study hash
manifest.

Two execution-instrumentation problems are retained separately. The first launch failed because the sandbox prohibited the `ps` command; its only selected test was a Cirq skip. The monitor was replaced with `psutil`. A subsequent partial run encountered the capture plugin's default 32-bit QASM export limit on the 100-qubit TKET source after compilation had succeeded. The capture limit was corrected, and all selected native sessions were rerun into `native_results_v2/`. These are instrumentation failures, not SDK failures. The original logs and partial results remain archived and are excluded from the completed-campaign table.

After extraction, install the recorded native requirements in a separate Python 3.12 environment and acquire the exact upstream source:

```sh
python3.12 -m venv benchpress/env
benchpress/env/bin/python -m pip install -r benchpress/native_requirements.txt
curl -L --fail https://codeload.github.com/Qiskit/benchpress/tar.gz/b695f30e83a32bac05b9b4d8e98d37ba9aae5236 -o benchpress/upstream.tar.gz
shasum -a 256 benchpress/upstream.tar.gz
tar -xzf benchpress/upstream.tar.gz -C benchpress
benchpress/env/bin/python benchpress/run_native.py native_rerun
benchpress/env/bin/python benchpress/analyze_native.py --results benchpress/native_rerun --output benchpress/native_rerun_summary.json
```

Check the archive hash above before extraction. Existing evidence directories are never reused as output locations. The full installed environment manifest supplies transitive package versions; the direct requirements file pins the SDKs and benchmark tools.

The retained native readout evidence can be checked without rerunning compilation. From the extracted archive directory, use the frozen harness and a new output filename:

```sh
benchpress/env/bin/python -m unittest discover -s harness/experiments -p 'test_native_clifford_readout.py' -v
benchpress/env/bin/python verification/validate_native_readout.py --repo harness --artifacts benchpress/native_results_v2 --output verification/native_readout_replay.json
```

The replay requires all 20 retained artifacts and writes exclusively; it never replaces the archived receipt or native outputs.

For the matched extension, use a checkout containing this experiment harness and a release build of `comparison_runner`. Set `QRUST_REPO` to that checkout's absolute path, then run the following from the extracted archive directory. Preparation creates a fresh protocol, including the actual compiler/source hashes and installed package versions; newly measured results are a separate campaign.

```sh
cargo build --release --manifest-path "$QRUST_REPO/experiments/Cargo.toml" --bin comparison_runner
benchpress/env/bin/python "$QRUST_REPO/experiments/benchpress_matched.py" prepare --repo "$QRUST_REPO" --binary "$QRUST_REPO/target/release/comparison_runner" --cases benchpress/frozen/cases.json --output matched_rerun
benchpress/env/bin/python "$QRUST_REPO/experiments/benchpress_matched.py" time --repo "$QRUST_REPO" --binary "$QRUST_REPO/target/release/comparison_runner" --output matched_rerun --session 1
benchpress/env/bin/python "$QRUST_REPO/experiments/benchpress_matched.py" time --repo "$QRUST_REPO" --binary "$QRUST_REPO/target/release/comparison_runner" --output matched_rerun --session 2
benchpress/env/bin/python "$QRUST_REPO/experiments/benchpress_matched.py" validate --repo "$QRUST_REPO" --binary "$QRUST_REPO/target/release/comparison_runner" --output matched_rerun
```

Run campaigns serially, without simultaneous compilation, verification or profiling jobs. Timings on another machine are new observations; reproducing archived analysis and validating recorded outputs are separate operations.
