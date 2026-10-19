# Input dialect: retained `sxdg` in the frozen held-out pipeline

The observed `Unknown gate: sxdg` rejection is a frozen pipeline/input-dialect interoperability gap. It is not evidence of a violation of the OpenQASM 2.0 standard or a wrong quantum transformation. No compiler, adapter, normalizer, fixture, timing schedule or validator was changed in response. Exact observed failure counts belong to the completed campaign's failure ledger, not this source-only note.

The source preparation routine `experiments/external_corpus.py:101–176` expands definitions until it reaches a primitive in the independent oracle's `_SIGNATURES`. It checks that primitive's matrix, then retains its name. This is not closure into the intersection of all four front ends or into a strictly portable OpenQASM 2.0 gate vocabulary. The oracle explicitly contains both `sx` and `sxdg` (`comparison_validation.py:39–42`) and supplies their analytic matrices (`:149–152`). A semantic source-eligibility pass therefore does not guarantee that every compiler can parse the resulting source.

The frozen 16-qubit random-circuit input has `sxdg` in its generated raw QASM custom-gate bodies (lines 4, 7 and 10), as well as a direct call at line 204. After source preparation the primitive appears directly, for example at line 32 of both the source and normalized QASM. The normalized text includes `qelib1.inc` but does not declare `sxdg`. This path can expose a dialect assumption that was already present in the generator/exporter; the compiler observations did not cause the input to be rewritten.

The pinned front ends differ as follows:

| Front end | Source-level behavior |
|---|---|
| Q-Rust | `src/ir/gates.rs:62` recognizes `sx`; no `SXdg` variant or `sxdg` case exists. Other spellings become `GateType::Custom` at line 92. `src/parser/mod.rs:408–437` accepts such a name only if a source declaration exists; otherwise it returns `UnknownGate`. `qelib1.inc` is an accepted include marker, not a file from which additional definitions are loaded (`:145–154`). Thus undeclared `sxdg` reaches the reported unknown-gate error before optimization. |
| Qiskit | The adapter explicitly loads with `qasm2.LEGACY_CUSTOM_INSTRUCTIONS` (`comparison_qiskit_tket.py:55`). The installed parser lists both `sx` and `sxdg` with `builtin=True` (`qiskit/qasm2/parse.py:172–173`), accepting them without a source declaration. |
| Cirq | The adapter uses `cirq.contrib.qasm_import.circuit_from_qasm` (`comparison_cirq.py:127`). Its installed QASM gate table maps `sx` to `XPowGate(exponent=0.5)` and `sxdg` to `XPowGate(exponent=-0.5)` (`cirq/contrib/qasm_import/_parser.py:772–776`). |
| TKET | The adapter uses `pytket.qasm.circuit_from_qasm_str` (`comparison_qiskit_tket.py:127`). Its installed no-parameter command table maps `sx` and `sxdg` to `OpType.SX` and `OpType.SXdg` (`pytket/qasm/qasm.py:124–125`). |

The Q-Rust SX decomposition itself is present (`src/ir/gate_def.rs:219`). This inspection did not test a proposed SXdg decomposition or modify support. SDK parser mappings explain why these names are recognized; whether a resulting full circuit passes the study's semantic checks is determined separately by the saved campaign validation.

The failed Q-Rust cells must remain in coverage reporting and cannot contribute Q-Rust speed or output-quality scores. Ratios on eligible paired cases are conditional on this reduced supported corpus. A future portable-basis input policy or additional compiler support would define a new study version and require a new preserved run; neither is substituted into this frozen campaign.

Qiskit's [OpenQASM 2 documentation, Legacy Compatibility](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qasm2) documents the compatibility mechanism. The causal attribution above comes from the exact installed sources and frozen input, not from assuming that the latest documentation matches every package release.

## Inspected source identities

Paths are relative to the original `working/` directory. These hashes identify the source inspected during the frozen run; they do not replace the campaign's own provenance manifest.

| File | SHA-256 |
|---|---|
| `qrust_pr_publish/experiments/external_corpus.py` | `d28d44d07a8ede0ffe5e4e7501bed0d8e46680b0ebb320895b9eb813af1a3104` |
| `qrust_pr_publish/experiments/comparison_validation.py` | `55de70f9fe47d9fcaeb1497bd04b377e5582afc9289e6948dc6c32ebe5bda6d4` |
| `qrust_pr_publish/experiments/comparison_qiskit_tket.py` | `9fa37aaa09208ed2178d17a2c2fdb1e1a52b98a5e75ec5c11a0a9abbbb855d7a` |
| `qrust_pr_publish/experiments/comparison_cirq.py` | `cec7669e6b3f23de3bf305a043a01477da34fe328ba9fd55b77fe6ae82dbfd35` |
| `qrust_pr_publish/src/parser/mod.rs` | `765939f05263d9edbfd1d88a0b5c436d3f16285126f0ad49519b8a856b313421` |
| `qrust_pr_publish/src/ir/gates.rs` | `ad6ba4c7d781b44fdba61175d9aefda7decbc6fbd4ddf1ad8c12e10f2922950e` |
| `qrust_pr_publish/src/ir/gate_def.rs` | `e35fc382e987b8b2dace686d41c56afa4e7c4e90ec7145d8506c4ad77b4ed8e1` |
| `qrust_comparison_env/lib/python3.12/site-packages/qiskit/qasm2/parse.py` | `59d08f1c0bc9841a87a028e27439a1f0350df655ceb60fc6de61759745355616` |
| `qrust_comparison_env/lib/python3.12/site-packages/cirq/contrib/qasm_import/_parser.py` | `b92c65fc0c231efeea7fd6ad8a7862d8f32a1d220fd1f3cb066c3607a77fe7f5` |
| `qrust_comparison_env/lib/python3.12/site-packages/pytket/qasm/qasm.py` | `1dd639f187a257c806913de4ec9af779a953917c7d74bb29cecc89648599269f` |
