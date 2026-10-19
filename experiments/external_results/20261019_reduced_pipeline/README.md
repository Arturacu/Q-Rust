# Frozen reduced-pipeline evidence

Read the [technical report](../../REDUCED_PIPELINE_STUDY.md) for the cumulative recipe definitions, timing boundaries, output-quality trade-offs, adverse cases and scope.

- `evidence.tar.gz` contains the complete new study: experiment-only source and native binary, frozen recipes and installed-file inventory, community-generator sources, source preparation revisions, selected input circuits, schedules, raw timed and warmup outputs, independent validation, analysis, tests, diagnostics, and rendering/replay helpers.
- `checksums.json` binds the archive and `evidence_manifest.json`; the same payload manifest is inside the archive.
- `ladder_analysis.json`, `comparisons.csv`, `quality_retention.csv` and the two session cell tables expose results without extraction.
- `stage_interpretation.json` retains every adjacent case effect and adverse case; `variability.json` pairs corresponding trials across sessions.
- `failure_ledger.json` distinguishes the successful formal campaign from preserved development/parser limitations and corrected preparation/audit attempts.
- `verification_summary.json` and `audit/` record completed integrity, replay and export checks. `delivery_checksums.json` binds every delivered file except itself.
- `original_preservation.json` binds 1,956 historical evidence files, checked unchanged. Earlier results remain in their original versioned directories.

## Scope

Seven Q-Rust recipes and three unchanged SDK recipes run on 18 selected MQT Bench instances at 6/10/14 qubits and 54 complete/linear/grid cases. There are two sessions, five timed trials plus one initial warmup per case/recipe/session: all 5,400 timed outputs and 1,080 warmups pass. Six-qubit checks are full numerical isometries; 10/14-qubit checks are sampled. Every comparison cohort contains all 54 cases.

The common executable measures the effect of executing fewer groups, not a smaller binary or Rust-specific causality. R5 preserves optimization groups and placement but reduces routing search. Its favorable mean trade-off includes significant individual regressions; it is not a new tuned default. The selected SDK recipes do not have equal search budgets. This is a custom study using community generators, not execution of the complete native Benchpress workout suite.

Protocol SHA256: `50b54e651176fbe06a68c41d2501a1c71f3db7956424b1834b4c2080c5ab87f9`.

## Replaying saved output validation and analysis

Extract into a fresh directory and use the pinned Python environment. Strict replay checks installed-file hashes; matching package version strings alone is insufficient. Installed distributions are recorded but not bundled. The native executable is specific to its original host architecture.

The study directory after extraction is `working/qrust_reduced_ladder_20261019`. Run:

```sh
/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_reduced_ladder_20261019/replay_heldout.py \
  --study /path/to/extracted/working/qrust_reduced_ladder_20261019 \
  --output /path/outside/extracted/new-replay
```

The output must be a new directory outside the archive tree. This command validates saved circuits and reproduces cell and aggregate analysis; it does not execute compilers. All compilation, parsing, resource, artifact/map, cost and verdict fields remain bound. Only fresh `validation.validation_ns` and the resulting checked validation-file hashes may differ.

## Reproducing report exports

Use fresh output locations. These commands read the archived analysis:

```sh
/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_reduced_ladder_20261019/reporting/export_tables.py \
  --campaign /path/to/extracted/working/qrust_reduced_ladder_20261019/corpus_work/campaign \
  --output /path/outside/extracted/new-tables

/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_reduced_ladder_20261019/reporting/plot_ladder.py \
  --analysis /path/to/extracted/working/qrust_reduced_ladder_20261019/corpus_work/campaign/analysis.json \
  --protocol /path/to/extracted/working/qrust_reduced_ladder_20261019/corpus_work/campaign/protocol.json \
  --output-dir /path/outside/extracted/new-plots --postflight-authorized

/path/to/pinned/venv/bin/python /path/to/extracted/working/qrust_reduced_ladder_20261019/reporting/render_report.py \
  --campaign /path/to/extracted/working/qrust_reduced_ladder_20261019/corpus_work/campaign \
  --tables /path/outside/extracted/new-tables \
  --template /path/to/extracted/working/qrust_reduced_ladder_20261019/reporting/report_template_v2.md \
  --report /path/outside/extracted/new-report.md \
  --receipt /path/outside/extracted/new-report-receipt.json
```

Original absolute paths inside manifests are provenance. The replay helpers remap validated archive-relative inputs and refuse missing or changed bytes. Development width-contract and barrier-check limitations remain recorded; successful held-out verification does not erase them.

