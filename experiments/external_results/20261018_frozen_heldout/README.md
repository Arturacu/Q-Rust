# Frozen compiler study evidence

Read the [study report](../../FROZEN_HELDOUT_STUDY.md) for results, timing boundaries, coverage and limitations.

- `evidence.tar.gz` contains the frozen inputs, compiler/adaptor sources, native binary, generator mirrors, protocols, raw measurements, diagnostic revisions, validation, analysis, tests and report-generation helpers.
- `checksums.json` binds the archive and its complete `evidence_manifest.json`. The archive contains the same manifest as `sha256_manifest.json`.
- `heldout_analysis.json` and `session_01_cells.csv` / `session_02_cells.csv` expose the aggregate results and individual cells without extraction.
- `failure_ledger.json` retains excluded cells and unqualified diagnostic outcomes. Successful compilation is not treated as independent verification.
- `compiler_freeze.json` records source, binary, runtime, installed-distribution and generator hashes. `original_preservation.json` binds 1,547 historical files, all checked unchanged before packaging.
- `audit/` contains execution logs and receipts. `verification_summary.json` summarizes completed checks. `delivery_checksums.json` binds every delivered file except itself.

## Reading the comparison

The held-out run has 30 prescribed inputs on three targets, two sessions, five timed trials and one warmup per tool/case/session. The four-tool comparison is conditional on 84 of 90 supported, independently checked cases. Six Q-Rust cases reject the retained undeclared `sxdg` input spelling. All compiled held-out outputs passed the unchanged checker; ten- and sixteen-qubit general-circuit checks are sampled, not proofs.

Timing and quality ratios describe the pinned recipes on this host. They do not establish equal optimization effort, a language-level speed effect, best SDK settings or hardware fidelity. The report also retains adverse cases and trial-level TKET variability. Dense outputs outside the frozen local-Clifford recognizer's coverage remain unverified, not incorrect.

## Replaying saved evidence

Use a fresh extraction and the original pinned Python environment. Version pins alone are insufficient for strict replay: the replay helpers additionally check installed-file hashes. Installed distributions are not bundled. The native binary is specific to its original host architecture.

From the extracted archive, the study directory is `working/qrust_frozen_holdout_20261018`. Run:

```sh
/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_frozen_holdout_20261018/replay_controls.py \
  --study /path/to/extracted/working/qrust_frozen_holdout_20261018 \
  --output /path/outside/extracted/new-control-replay

/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_frozen_holdout_20261018/replay_heldout.py \
  --study /path/to/extracted/working/qrust_frozen_holdout_20261018 \
  --output /path/outside/extracted/new-heldout-replay
```

Outputs must be new and outside the extracted archive tree. These commands validate saved circuits and reproduce analysis; they do not execute compilers. Compiler, parsing, request and pass timings, structural costs, maps and verdicts are retained exactly. Only explicitly identified fresh verification durations and independently checked derived validation-file hashes may differ.

The held-out replay began from the first fresh archive extraction while a final prose review corrected one sentence about subgroup summaries. `presentation_revision_binding.json` confirms that all computational inputs and 685 payload files are byte-identical between that extraction and the final archive; only report/completion bookkeeping changed. The first manifest and checksums are retained. Earlier report revisions are inside the final archive. Diagnostic validation and report/figure generation were replayed again from the final archive. All final report exports and plot files reproduce byte-for-byte.
