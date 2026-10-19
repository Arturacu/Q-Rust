# External community benchmark methods

This study compares compilation of a fixed common gate representation using Q-Rust, Qiskit, Cirq and TKET. Inputs are drawn from QASMBench and MQT Bench. It measures compilation time and circuit cost on abstract connectivity graphs. It does not execute circuits on hardware, measure application success, or reproduce the complete Benchpress suite.

## Corpus selection and provenance

The selection rule was fixed before compiling the external candidates or observing their compilation time and output quality. The source snapshots are:

| Source | Frozen revision | Role |
| --- | --- | --- |
| [QASMBench](https://github.com/pnnl/QASMBench/tree/357b942396d5c2b7cbc1c229c585a6ef5ccaebac) | `357b942396d5c2b7cbc1c229c585a6ef5ccaebac` | Original, non-transpiled OpenQASM 2 source |
| [MQT Bench 2.3.0](https://github.com/munich-quantum-toolkit/bench/tree/2cd3f3b08d5d8866595948fbce955e701e86a5df) | `2cd3f3b08d5d8866595948fbce955e701e86a5df` | Algorithm-level circuit generators |
| [Benchpress](https://github.com/Qiskit/benchpress/tree/b695f30e83a32bac05b9b4d8e98d37ba9aae5236) | `b695f30e83a32bac05b9b4d8e98d37ba9aae5236` | Methodology reference only |

For QASMBench, the inventory contains every OpenQASM file in the pinned `small` and `medium` directories. The prespecified scope excludes filenames containing `_transpiled`, directory widths above 14 qubits, and raw files larger than 250,000 bytes. Declared directory width is checked against the parsed circuit width. All remaining files are considered; the selection is not a sample ranked by Q-Rust performance.

For MQT Bench, the prespecified algorithm families are `ae`, `dj`, `grover`, `qnn`, `qpeexact`, `qpeinexact`, `qwalk`, `vqe_real_amp` and `vqe_su2`, each at widths 5 and 10. Generation uses `get_benchmark_alg`, the pinned upstream default options, and upstream random parameter binding with seed 10. The target-independent, native-gate and mapped compilation paths are not used.

Both sources must produce a static unitary circuit after the explicitly defined terminal-measurement and barrier removal below. The expanded source is capped at 5,000 gates. Unsupported or malformed candidates are retained in the selection ledger with their source identity and exact reason; no candidate is replaced after its exclusion.

The denominator has several distinct levels:

- The inventory contains **150 entries**: 132 QASMBench paths and 18 prescribed MQT generator configurations.
- **84 entries are excluded by scope before circuit processing:** 64 upstream-transpiled variants, 18 widths above 14, and two files above the raw-size cap.
- The remaining **66 candidates** include 48 QASMBench inputs and 18 MQT configurations.
- **16 candidates are excluded during eligibility processing:** nine dynamic/nonunitary sources, five malformed OpenQASM sources, and two MQT circuits above the expanded-gate cap.
- The retained corpus contains **50 fixtures:** 34 QASMBench and 16 MQT Bench. Applying both target graphs produces **100 circuit-target cases**.

The inventory's 100 exclusions are not 100 failed compiler jobs. They include explicit scope exclusions and duplicate upstream compiled representations. Compiler errors and timeouts, if any, are a separate outcome category over the 100 selected cases.

The corpus manifest retains source URLs and revisions, upstream Git blob identities, original and transformed artifact hashes, generator settings, width, family labels, transformation records, eligibility decisions and numerical normalization checks. Raw QASMBench bytes remain unchanged. Original MQT circuits are retained as QPY; an OpenQASM export is also retained where supported. Qiskit's inability to export the original AE controlled-unitary representation to OpenQASM is recorded as an optional export limitation, rather than an exclusion of the valid source circuit.

The QASMBench license is the upstream Battelle license, retained verbatim with its NOTICE. It should not be relabeled as unmodified BSD-3-Clause text. MQT Bench's MIT license is retained verbatim. Attribution and required notices accompany redistributed derived inputs.

## Common input representation

The comparison uses a static unitary-prefix track. Only a global suffix of terminal measurements may be removed. A quantum operation after the first measurement, any reset, or classical control makes that source ineligible. Barriers are removed and counted. This does not claim to reproduce the measured or dynamic behavior of a complete application.

For QASMBench, the exact pinned `qelib1.inc` is inlined before parsing. This avoids silently substituting an installed header with different definitions or phase conventions. Supported named primitives are retained only after comparing the complete operation matrix against the corresponding analytic standard-gate matrix. Controlled operations are checked as complete controlled matrices: a phase on the uncontrolled base operation is not assumed to be harmless under control.

Other instructions are expanded through source- or Qiskit-provided definitions. Permutations are expanded exactly into SWAPs. This process invokes no compiler optimization pass manager, routing, or placement. However, accessing Qiskit operation definitions can itself invoke internal decomposition or synthesis. The frontend is therefore an explicit dependency of the experiment. These results concern the delivered common gate representation, not each SDK's best source-preserving high-level API path.

The existing exact CRZ/RZZ expansion into RZ/CX then produces the shared compiler input. Every recipe receives identical normalized QASM. No baseline compiler optimizes another compiler's input. Whole-circuit global phase is disregarded, with one shared phase across all checked input states.

Before admission, the original generated/imported unitary prefix is simulated using Qiskit's Statevector implementation and compared with the independent analytic NumPy interpreter on both the canonical and normalized gate text. The check uses every computational-basis input through six qubits and eight inputs at larger widths: three basis states, one uniform superposition and four seeded Haar states. The residual tolerance is `1e-7`. This uses separate numerical implementations, although source parsing and generator definitions remain shared trust boundaries.

## Targets and semantic contract

Each retained circuit is compiled for all-to-all connectivity and an undirected line with the same number of physical and logical qubits. Edges are ideal and bidirectional. The common output basis is arbitrary one-qubit U3 gates and CX gates. Initial and final logical permutations are allowed and must be reported by the compiler adapter.

The required behavior is preservation of the source unitary on arbitrary logical inputs, up to one whole-circuit global phase and the recorded input/output mappings. This is stronger than checking preparation from the all-zero state. The study has no added physical ancillas. It makes no claim about pulse durations, calibrated errors, directed hardware instructions, measurement fidelity or physical execution fidelity.

## Software and execution protocol

The compiler source is frozen at Q-Rust commit `dbc4645f3183b2d44dc076f431b98932510c9a20`. The recorded study executable has SHA256 `004ca5c49481ef7eed018bb33686ee0e8f27e6504758ba0b0d804783ee90d047`. The main SDK versions are Qiskit 2.5.2, Cirq 1.7.0 and pytket 2.18.5. Python dependencies are pinned in `experiments/requirements-comparison.txt`; the protocol additionally records interpreter, platform, source and binary hashes.

The primary comparison retains the existing Q-Rust, Qiskit, Cirq and TKET recipes; [the exact recipe descriptions](EXTERNAL_COMPILER_RECIPES.md) identify their settings and objectives. Two prespecified sensitivity recipes, `cirq_one_pass` and `tket_light`, are interleaved with them. These six recipe identities are recorded independently. The primary and lighter analysis profiles reuse the same Q-Rust and Qiskit observations; they are not extra measurements of those compilers. Recipe settings have unequal search effort and should not be described as identical optimization budgets or each SDK's universal best configuration.

The protocol prescribes two sequential sessions, seven timed trials per case and recipe in each session, and one initial warmup per case and recipe. Cases, trial order and recipe order are shuffled using the frozen schedule. Trial seeds are `20261018 + trial` and are reused across the two sessions. Qiskit consumes the seed; Q-Rust, Cirq and TKET report deterministic recipes. Repeated trials and sessions do not create additional independent circuits. Any worker restart and replacement warmup is recorded. Each request has a 60-second timeout, and all timed outcomes are retained.

Timing uses persistent workers and sequential requests. The primary `compile_ns` measurement covers warm in-process input cloning, configuration, compiler passes and compiled-circuit retrieval. Parsing is recorded separately as `parse_ns`. Imports and process startup, serialization, interprocess communication and correctness checks are outside compilation time. `request_wall_ns` covers the complete worker request and response. The residual `noncompile_nonparse_request_ns` combines serialization, layout reporting, diagnostics, communication and interpreter overhead; it is not an isolated measurement of serialization.

The recorded host is an Apple M3 with 8 GiB RAM, eight physical/logical cores and macOS 14.4.1 on arm64. Native thread pools are capped at one through the recorded environment variables. There is no CPU affinity, frequency or power control. Study sessions, builds and numerical simulations do not overlap with timing; ordinary operating-system background activity remains uncontrolled. Two sessions on this host do not establish performance portability or isolate the effect of Rust as an implementation language.

## Output validation and analysis

Correctness checking starts only after all prescribed timing sessions finish. It examines the retained output produced by each timed request, using the unchanged independent NumPy analytic gate interpreter. The immediate reference is the canonical source gate text, which was separately checked against the original imported/generated unitary prefix before timing. These two numerical checks have separate error tolerances; the original-source frontend remains a trust boundary. Layouts, output basis and topology are checked alongside the complete physical amplitudes. Every input basis state is checked through six physical qubits; the same eight-input strategy is used above that size. The residual threshold remains `1e-7`.

Full numerical checks have complete input coverage for the specified finite-dimensional circuit but are not a symbolic proof of the compiler. Successful sampled checks do not prove equivalence. Failed, unsupported, timed-out and inconclusive outcomes must remain distinct from accepted numerical evidence.

Circuit cost is reported separately as CX count, CX depth and total depth in the common U3/CX basis. These are algebraic circuit metrics, not hardware duration or application fidelity. Best among the tested recipes is not a certified global optimum.

Sessions are analyzed separately using paired comparisons and complete-cell eligibility criteria. The analysis retains every selected case, reports failures and coverage, and stratifies by source suite, family, width and topology. Circuit families and related variants should not be treated as independent applications simply because they have different labels. Paired ratio summaries and fixture-cluster uncertainty estimates describe this bounded corpus; they do not turn it into a representative random sample of quantum applications. The primary and lighter profiles provide sensitivity to recipe effort, without claiming a matched-time-budget comparison.

## Reproduction

The following commands reproduce the workflow in new output directories. The measured results and preserved evidence are linked in [the results report](EXTERNAL_BENCHMARKS.md). Run them from the repository revision containing the frozen external-study harness. Use a new output directory and an otherwise idle host. Source acquisition requires network access; subsequent compilation and checking use the retained local artifacts.

```sh
python3.12 -m venv .venv-external
. .venv-external/bin/activate
python -m pip install -r experiments/requirements-comparison.txt

task_run_dir="$(mktemp -d "${TMPDIR:-/tmp}/qrust-external.XXXXXX")"
python experiments/external_corpus.py \
  --acquire --sources "$task_run_dir/sources"
python experiments/external_corpus.py \
  --sources "$task_run_dir/sources" \
  --output "$task_run_dir/corpus"
```

Build the frozen compiler source in a separate directory so that the current harness and the fixed compiler revision remain distinct:

```sh
task_compiler_dir="$task_run_dir/compiler"
mkdir "$task_compiler_dir"
git archive dbc4645f3183b2d44dc076f431b98932510c9a20 \
  | tar -x -C "$task_compiler_dir"
cargo build --manifest-path "$task_compiler_dir/Cargo.toml" \
  --locked --release -p qrust-bench --bin comparison_runner
task_binary="$task_compiler_dir/target/release/comparison_runner"

python experiments/external_comparison.py prepare \
  "$task_run_dir/corpus/manifest.json" "$task_run_dir/study" \
  --binary "$task_binary" \
  --compiler-commit dbc4645f3183b2d44dc076f431b98932510c9a20 \
  --trials 7 --sessions 2 --timeout 60 --seed 20261018 \
  --machine-note "Describe the actual reproduction host and any uncontrolled activity"
python experiments/external_comparison.py time \
  "$task_run_dir/study" --binary "$task_binary" --session 1
python experiments/external_comparison.py time \
  "$task_run_dir/study" --binary "$task_binary" --session 2
python experiments/external_comparison.py validate "$task_run_dir/study"
python experiments/analyze_external_comparison.py "$task_run_dir/study"
```

The original executable hash identifies the executable used for this study, rather than a promise of bit-identical binaries across toolchains and hosts. A reproduction must retain its own binary hash and actual environment. Do not change the compiler, corpus, adapters, normalization or analysis sources after preparing the protocol; the runner verifies its frozen contract before each timing session.

## Artifact checks

The original study protocol has SHA256 `ddb4fd217613a77638a6fcc20d4b6b72872ed6132f646759683bedfdea821645`; the selected-corpus manifest has SHA256 `2d14203e785b133a9e9d7cacf85b5c8222ce81774988488051beb2cb7288f1ce`. The corpus/harness were committed at `8f72d86` before timing. The compressed [evidence archive](external_results/20261018/evidence.tar.gz) contains the original records and frozen harness; [its manifest](external_results/20261018/evidence_manifest.json) records file checksums. Use a fresh extracted copy if re-running validation/analysis; preserve the original validation and analysis outputs first, because the runner refuses to overwrite them.

## Community references

[QASMBench](https://doi.org/10.1145/3550488) supplies low-level quantum circuits from multiple application areas. [MQT Bench](https://doi.org/10.22331/q-2023-07-20-1062) distinguishes algorithmic, target-independent, native-gate and mapped abstraction levels; this study uses its algorithm-level generators. The [versioned MQT catalog](https://mqt.readthedocs.io/projects/bench/en/v2.3.0/benchmark_selection.html) identifies the supported families.

[Benchpress's methods](https://www.nature.com/articles/s43588-025-00792-y) provide a relevant precedent for common workloads, separate loading/compilation costs, two-qubit count/depth, target variation and explicit unsupported/error/timeout outcomes. Its operational pass status should not be interpreted as an all-input equivalence certificate. The present experiment is a bounded study of externally sourced circuits with its own contract, rather than a complete execution or reproduction of Benchpress.
