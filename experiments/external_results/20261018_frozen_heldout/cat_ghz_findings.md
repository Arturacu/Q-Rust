# Cat/GHZ placement controls

The frozen Q-Rust Cat/GHZ CX gap is primarily a placement-search limitation on these four circuits. Every input is a Hadamard followed by a logical CX chain. Each target has a simple path long enough to execute that complete circuit without SWAPs. The other three default recipes discover such paths; Q-Rust's default SABRE search does not. Relabeling each target so a witnessed path occupies its first logical labels lets the unchanged Q-Rust binary attain the same minimum CX counts. This is evidence of heuristic label sensitivity and a missing zero-SWAP placement opportunity; the relabeling intervention is not a fixed-layout API.

## Executed evidence

The final protocol was frozen before compilation. All 124 planned interventions were attempted. They produced 122 compiled outputs and two retained diagnostic errors. Of the 122 compiled intervention outputs, 112 pass the independent all-logical-input Clifford isometry checker and ten are outside its local-Clifford recognition capability. Four additional analytic path witnesses also pass. Thus the 128 artifacts comprise 112 verified compiler interventions, four verified witnesses, ten unsupported compiled controls and two errors. Twelve separate native SDK baseline calls match their instrumented counterparts exactly in output QASM and input/output maps. Five helper regression tests pass. All 348 phase, output-map and missing-edge negative checks associated with the 116 verified artifacts are rejected.

These are diagnostic experiments with a single seed and one execution per intervention. Instrumented timings are retained for trace inspection, not used as comparative speed estimates. The shared compilers and benchmark adapters were not changed.

| Input/target | Input CX | Inserted SWAPs | CX after direct SWAP lowering | Q-Rust final CX | Q-Rust path relabel CX | Qiskit / Cirq / TKET default CX |
|---|---:|---:|---:|---:|---:|---:|
| Cat-22 square | 21 | 3 | 30 | 28 | 21 | 21 / 21 / 21 |
| Cat-22 heavy-hex | 21 | 10 | 51 | 51 | 21 | 21 / 21 / 21 |
| GHZ-23 square | 22 | 5 | 37 | 33 | 22 | 22 / 22 / 22 |
| GHZ-23 heavy-hex | 22 | 10 | 52 | 52 | 22 | 22 / 22 / 22 |

All default and path-relabel rows in this table are verified for arbitrary logical inputs, clean spare wires and the recorded input/output maps. Relabeled outputs are restored to the original physical coordinates before validation. The target transformation preserves every edge and its original position/orientation in the edge list under one bijection. The path controls have identical input/output layouts, no routing SWAPs and exactly the original circuit's one-qubit/CX costs.

The lower-bound argument is separate from that isometry check. On product zero input, H plus the chain produces an n-qubit GHZ state. Under the unitary arbitrary-one-qubit/CX contract, the n logical output wires must lie in a connected interaction component. A connected component containing n initially independent output wires needs at least n−1 two-qubit interactions; clean spare wires cannot reduce that number. The independently validated path realization attains n−1 CX while implementing the complete input unitary. This establishes the CX optimum for these four cases under this contract, not for general GHZ routines, other gate sets or arbitrary circuits.

## Q-Rust: where the extra CX appear

The input and every pre-routing simplification stage retain 21 or 22 CX. Default initial layouts leave 3/7/5/7 logical chain neighbors nonadjacent, with static distance excess 4/10/6/10. Those derived map diagnostics are not the SABRE objective value. BeamSabre introduces 3/10/5/10 SWAPs respectively; the traces show exactly where their CX costs enter. Subsequent Clifford or inverse cancellation recovers two CX on Cat-square and four on GHZ-square, with no recovery on heavy-hex. Two-qubit block synthesis accepts no block on any default case: the considered pure-CX blocks have no strict CX improvement. These observations localize the remaining cost to routing, rather than a gate-basis discrepancy or hidden resynthesis benefit in a competitor.

Disabling Clifford simplification, block synthesis, or both leaves each final default CX/depth result unchanged. On square targets the ordinary inverse-cancellation path recovers the same final CX count when Clifford simplification is disabled. A generic claim that the Clifford pass is responsible for those final savings would therefore be too strong. The default output-retention guard selects among equivalent final costs, which is why a pass can show a reduction inside one trace while its omission leaves the final score unchanged.

Increasing layout trials from 10 to 40 leaves three cases unchanged and reduces GHZ-heavy-hex only from 52 to 50 CX. Increasing layout iterations from three to six changes none. Reducing trials to one worsens Cat-square 28→58 and GHZ-square 33→49; reducing iterations to one worsens them to 40 and 35. Heavy-hex remains at 51/52 for those reductions. More search effort is not an effective substitute for an appropriate path-shaped starting placement in this bounded test.

The path relabel keeps the graph isomorphic and lets the unchanged search encounter a zero-SWAP trivial placement. Source inspection explains why this is effective: Q-Rust first evaluates that trivial placement and zero SWAPs trigger its proven lower-bound early exit. Actual layout-route evaluation counts drop from 56/48/56/53 to one in each case. The only change supplied to Q-Rust is the graph's physical labels; the circuit, compiler, optimization policy and checker remain frozen.

## Why the other recipes get their scores

**Qiskit.** Each default trace executes VF2Layout and no SabreLayout. The returned initial map embeds every logical chain edge directly in the target, and initial/final maps agree. Removing the entire post-routing optimization stage preserves 21/22 CX and depth 22/23. Fixed witnessed layouts also preserve those costs. Thus the observed minimum comes from discovering a valid subgraph embedding, rather than later two-qubit optimization. The pinned preset's source tries VF2 before its SABRE fallback; the official [VF2Layout documentation](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.VF2Layout) describes its subgraph-isomorphism purpose.

**TKET.** The default GraphPlacement-based mapping stage leaves the original 21/22 CX and introduces no SWAP. Omitting FullPeepholeOptimise, SynthesiseTket, or both preserves all final gate/depth costs. Replacing only the placement/mapping construction with explicit identity placement plus RoutingPass gives verified outputs with 63/158/67/165 CX; a witnessed path plus RoutingPass returns to 21/21/22/22. This supports placement as the decisive variable in these cases. [GraphPlacement's documented method](https://docs.quantinuum.com/tket/api-docs/placement.html) matches a weighted circuit-interaction graph to the target graph; the saved pass specifications record the concrete recipe.

**Cirq.** The default RouteCQC trace introduces no SWAP, and its initial/final maps agree. Its [LineInitialMapper](https://quantumai.google/reference/python/cirq/LineInitialMapper) explicitly builds chains in the logical interaction graph and tries to place them along physical paths. Its CZ synthesis stages preserve the 21/22 entangling gates, but leave substantially more one-qubit gates after CX closure: default depths 43/45 and gate counts 65/68, versus the original-chain depths 22/23 and gate counts 22/23. Omitting both synthesis stages while preserving RouteCQC and the existing CX-closure/fusion step attains the original-chain costs, with full isometry verification. Omitting just either synthesis stage does not remove this depth overhead. This is a recipe-specific representation trade-off, not a claim about Cirq's best possible pipeline.

## Preserved limits and unsuccessful controls

The Qiskit forced-SABRE and forced-trivial-layout arms compile, but their arbitrary one-qubit decompositions are not locally Clifford, even when the overall circuit may be Clifford. Eight such outputs are unsupported by the frozen checker. Cirq's identity-mapper square controls produce two similarly unsupported decompositions. Their structural metrics remain in the raw and derived data but cannot support a correctness-qualified output-quality conclusion. They are not classified as incorrect circuits.

Cirq's identity-mapper heavy-hex calls raise `IndexError: list index out of range` in the instrumented diagnostic call. Both outcomes are preserved, and neither is included in verified comparisons. No compiler patch, adapter workaround or checker relaxation was introduced to make these controls pass.

The full saved stage traces, native baseline outputs, layouts, options, failures and verifier details are in `raw.jsonl`, `validated.jsonl` and `native_*.json`. `derived_summary.json` is produced from them by `derive_summary.py`; `protocol.json` binds the frozen binary, fixtures and source hashes. Earlier pre-execution protocol revisions and helper sources remain alongside the final protocol for auditability.
