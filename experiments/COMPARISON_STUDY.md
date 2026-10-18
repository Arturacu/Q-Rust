# Repeated compiler comparison

This harness compares the repaired public Q-Rust compiler with specified Qiskit,
Cirq and TKET recipes. It supersedes historical performance comparisons that did
not preserve and independently validate each measured artifact. It does not
replace the separate correctness regression suites.

The recorded numerical study below uses the compiler before the subsequent
[routing investigation](ROUTING_INVESTIGATION.md). Current O3 has deterministic
tie handling, a smaller configurable placement budget and new exact cleanup;
its completed rerun is reported separately there, including all measured
outputs' validation and sensitivity to alternative competitor recipes.

The latest [two-qubit block investigation](BLOCK_SYNTHESIS_INVESTIGATION.md)
reports the following synthesis change, its isolated overhead, remaining QFT
gap and fresh four-compiler comparison. Historical data below are unchanged.

## Reproduction

Use Python 3.12 and a clean virtual environment. The complete dependency lock is
`requirements-comparison.txt` (Qiskit 2.5.2, Cirq 1.7.0, pytket 2.18.5).

```sh
python3.12 -m venv .venv-comparison
.venv-comparison/bin/python -m pip install -r experiments/requirements-comparison.txt
.venv-comparison/bin/python -m unittest discover -s experiments -p 'test_comparison*.py'
cargo build --release -p qrust-bench --bin comparison_runner
.venv-comparison/bin/python experiments/comparison_study.py time experiments/results/NEW_COMPARISON \
  --binary target/release/comparison_runner --machine-note 'Describe CPU, RAM and power settings'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-comparison/bin/python experiments/comparison_study.py validate \
  experiments/results/NEW_COMPARISON
.venv-comparison/bin/python experiments/analyze_comparison_study.py experiments/results/NEW_COMPARISON
```

The timing command refuses an existing destination. A `--pilot --trials 1` run
uses a separate 17-case integration subset; pilot observations are not included
in the main results. Timing and simulation must be run separately, with no other
benchmarking, builds or rendering competing for CPU during the timing phase.

## Frozen comparison contract

The main corpus is every existing corrected fixture with at most 14 logical
qubits: 51 fixtures in nine families. All-to-all and linear connectivity use the
logical width; Nairobi connectivity uses seven physical wires where the fixture
fits. This makes 123 fixture/target cases. All edges are bidirectional. The target
is arbitrary one-qubit U3 plus CX, not a calibrated device-native gate set.

Every tool receives the same source, after an exact CRZ/RZZ-to-RZ/CX expansion
needed by the Cirq importer. CCX lowering remains inside compilation. This is a
matched intermediate-representation comparison rather than a test of each
tool's richest native input. The independent checker compares the saved output
against the original unnormalized fixture, so normalization is also checked.

Each tool runs in a persistent process, with one warm-up per case and seven
measured trials. Case order and within-trial tool order are shuffled with seed
20261015. Compilation requests are sequential. Qiskit receives a recorded
transpiler seed. The original Q-Rust study had nondeterministic routing ties;
the revised compiler has deterministic ties and fixed internal layout seeds,
recorded in its worker output. Cirq and TKET use deterministic recipes.
Seven trials do not equalize algorithmic
search effort. Environment variables cap common numerical and native worker
pools at one thread; processor affinity and frequency are not fixed.

The primary timer includes input copying, compiler configuration and passes. It
excludes import/startup, source parsing, serialization, IPC and verification.
Parse time is saved separately; the secondary parse-plus-compile measure adds
the two times **within each trial** before taking medians. It is still not CLI
end-to-end latency. Q-Rust calls `transpile_with_report` with optimization level
3; it does not use a duplicate experimental pass pipeline.

Qiskit uses optimization level 3, exact synthesis, arbitrary-input semantics
(`qubits_initially_zero=False`) and U/CX. TKET uses FullPeepholeOptimise, mapping
where needed, SWAP lowering, post-routing SynthesiseTket and U3/CX rebasing.
Cirq uses CZ synthesis, RouteCQC where needed, post-routing synthesis, an exact
CZ-to-CX conversion and one-qubit fusion. See adapter source and each record's
`pipeline` field for exact settings. These are documented recipes, not a claim
to each framework's best achievable result.

## Evidence and interpretation

`protocol.json` records the case list, original/normalized sources and hashes,
versions, thread environment, compiler revision, source hashes and binary hash.
`warmup.jsonl` and `timed.jsonl` retain outcomes, including errors and timeouts.
Every measured successful output retains QASM, layouts, timing and a hash.
`validated.jsonl` adds independent numerical evidence and common output metrics.
Identical artifacts may reuse validation; source and output hashes are checked
before the artifact cache is consulted.

The NumPy oracle independently implements analytic gate matrices. It compares
full physical amplitudes under input/output layout isometries and a **single**
global phase, without projecting away dirty ancillas. Through six physical
qubits it checks every logical basis input. Above that it checks zero, all-one
and alternating-bit inputs, uniform superposition and four seeded Haar inputs.
Residual tolerance is 1e-7. It also checks basis and coupling legality. Floating
point checks are numerical evidence; sampled checks are neither a proof nor a
process-fidelity estimate. Correctness overlap is not the output-quality score.

Quality means CX count, all-gate dependency depth and CX-only dependency depth,
computed by the same code for all tools. It does not mean hardware success
probability. Runtime summaries use per-case medians and paired baseline/Q-Rust
ratios. Quality ratios use Q-Rust/baseline, with zero denominators reported
separately. Only complete, validated repeated cells enter comparative summaries;
all exclusions and failures remain in the denominator report. Confidence
intervals resample **fixtures**, keeping their topology cases together. They
describe sensitivity within this in-house corpus, not an unseen population.

Descriptive cost margins are measured relative to the best observed baseline
median. No tested baseline establishes the optimum. Runtime measurements
cannot isolate a language effect: Qiskit already contains native Rust code and
TKET has a C++ core. Memory, deployment size, cold start, external holdout
circuits, directed native targets and equal-search-budget comparisons require
separate experiments before supporting broader claims.

## Recorded study on 16 October 2026

On one Apple M3 host with 8 GiB RAM, all 3,444 measured trials compiled and passed
the defined independent checks. This covers 920 distinct semantic artifacts:
728 trial records used full-isometry checks and 2,716 used sampled checks. No
cell was excluded. The 61 Python experiment tests and strict Clippy check for
the new Rust worker passed; the compiler core is unchanged from the prior
305-test correctness repair.

Paired geometric means across 123 cases, using all seven trials' medians:

| Baseline | Baseline / Q-Rust compile time | Q-Rust / baseline CX count | Q-Rust / baseline depth |
| --- | ---: | ---: | ---: |
| Qiskit 2.5.2 | 1.60 | 1.110 | 1.005 |
| Cirq 1.7.0 | 51.7 | 0.851 | 0.604 |
| TKET 2.18.5 | 147 | 0.985 | 0.982 |

Q-Rust was within 10% of the best observed baseline's median CX count in
75/123 cases (61.0%), and within 10% on depth in 95/123 (77.2%). Qiskit was
faster on the linear-target aggregate; QFT14 on a line was 16.9 times slower in
Q-Rust and used 523 versus 291 CX gates. W-state circuits had about 49% more CX
gates than Qiskit on the family geometric mean. These unfavorable cases are
part of the conclusion, not removed from it.

A separate three-repeat timing diagnostic on GHZ8/all-to-all and QFT8/linear
found that over 99% of the selected TKET recipe's time was in pass application,
not configuration construction. This does not establish that the selected
recipe is TKET's fastest useful recipe; lighter recipes remain a follow-up.

The committed [summary](validation/comparison_20261016_summary.json) contains
fixture-cluster intervals, family-balanced estimates, topology/family results
and worst cases. [Provenance](validation/comparison_20261016_provenance.json)
records versions, hashes, numerical coverage and an independent audit. Full
timed/validated QASM records, the frozen source snapshot, CSV and plots are
generated artifacts retained in the local results bundle.
