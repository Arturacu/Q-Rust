# TKET repeatability in the frozen held-out comparison

The two sessions contain seven TKET matched-request / corresponding-trial pairs with different structural outputs, out of 450 pairs. The remaining 443 pairs have identical artifacts. All changed outputs passed the independent numerical oracle, basis and connectivity checks. Qiskit and Cirq each have 450 identical artifact pairs; Q-Rust has 420 identical pairs and 30 matching input-error pairs.

The frozen protocol expected the non-Qiskit recipes to be deterministic. These observations contradict that expectation for TKET. The protocol and raw outcomes remain unchanged. A recipe ignoring the supplied compiler seed does not establish internal repeatability.

| Fixture / target | Trial (zero-based) | Session 1: CX / depth / gates / two-qubit depth | Session 2: CX / depth / gates / two-qubit depth |
|---|---:|---:|---:|
| bmw_quark_copula_n16_s11__linear | 0 | 1020 / 411 / 1200 / 363 | 874 / 368 / 1054 / 316 |
| bmw_quark_copula_n16_s11__linear | 3 | 669 / 252 / 849 / 213 | 1020 / 411 / 1200 / 363 |
| bmw_quark_copula_n16_s29__linear | 1 | 874 / 368 / 1054 / 316 | 1020 / 411 / 1200 / 363 |
| bmw_quark_copula_n16_s29__linear | 2 | 1020 / 411 / 1200 / 363 | 874 / 368 / 1054 / 316 |
| bmw_quark_copula_n16_s29__linear | 3 | 669 / 252 / 849 / 213 | 1020 / 411 / 1200 / 363 |
| bmw_quark_copula_n16_s29__linear | 4 | 1020 / 411 / 1200 / 363 | 874 / 368 / 1054 / 316 |
| qaoa_n16_s11__linear | 4 | 1175 / 566 / 1347 / 507 | 1066 / 464 / 1238 / 412 |

The changes affect three cells. Their five-trial structural medians are identical between sessions:

| Fixture / target | CX | Depth | Gates | Two-qubit depth |
|---|---:|---:|---:|---:|
| bmw_quark_copula_n16_s11__linear | 1020 | 411 | 1200 | 363 |
| bmw_quark_copula_n16_s29__linear | 1020 | 411 | 1200 | 363 |
| qaoa_n16_s11__linear | 1066 | 464 | 1238 | 412 |

No eligible structural cell median changes between sessions for any recipe: TKET 0/90, Qiskit 0/90, Cirq 0/90 and Q-Rust 0/84. Consequently, the reported cell-median aggregate quality ratios agree between sessions, although TKET trial-level output costs are not repeatable.

The installed API documentation supplies a plausible mechanism, not a causal finding. `GraphPlacement` has a default timeout of 1,000 ms and documents stopping its monomorphism search on elapsed time. `DefaultMappingPass` documents using `GraphPlacement`; the frozen adapter invokes it on constrained targets. Wall-clock-bounded search could therefore yield different candidate sets. We have no trace establishing timeout activation in these seven requests and did not run a placement-budget ablation. Other internal causes have not been excluded.

Exact local source evidence (SHA-256 binds the complete file):

- `qrust_comparison_env/lib/python3.12/site-packages/pytket/_tket/placement.pyi`; SHA-256 `824d27df8d5996dfa1c50a8e95f1f091bf1dbb4602e07b8db0f0a37c6cd14f25`.
  Lines 91–104; see the machine-readable companion for verbatim excerpts.

- `qrust_comparison_env/lib/python3.12/site-packages/pytket/_tket/passes.pyi`; SHA-256 `f9d288a35b77522e5cc3dd93c139663323e3546ae6f833bad2f910e9a57d148d`.
  Lines 481–488; see the machine-readable companion for verbatim excerpts.

- `qrust_pr_publish/experiments/comparison_qiskit_tket.py`; SHA-256 `9fa37aaa09208ed2178d17a2c2fdb1e1a52b98a5e75ec5c11a0a9abbbb855d7a`.
  Lines 137–149; see the machine-readable companion for verbatim excerpts.

`TKET_VARIABILITY.json` stores the exact seven row pairs, output/artifact hashes, oracle evidence, unchanged medians, source excerpts and hashes of the supporting existing files. This note adds no compiler execution and changes no frozen evidence.

[Exact machine-readable companion](./TKET_VARIABILITY.json).
