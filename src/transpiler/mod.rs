//! Transpiler: optimization, decomposition, layout, and routing.
//!
//! See [`transpile`] for the canonical entry point and
//! [`TranspilerConfig::builder`] for configuring the pipeline.

pub mod block_synthesis;
pub mod dag;
pub mod decomposition;
pub mod layout;
pub mod local_cleanup;
pub mod numerical_policy;
pub mod optimization;
pub mod pass;
pub mod pauli_tracker;
pub mod phase_gadgets;
pub mod profiler;
pub mod property_set;
pub mod report;
pub mod routing;
pub mod synthesis;
pub mod target_basis;

pub use report::{PassDiagnostic, StageSnapshot, TranspilationReport};

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateDefinition, GateType, Operation};
use pass::{Pass, PassManager};
use std::collections::HashSet;

/// Independently configurable routing policies for reproducible search ablations.
#[derive(Clone, Debug)]
pub struct RoutingSearchConfig {
    /// Override the optimization-level layout trial budget.
    pub layout_trials: Option<usize>,
    /// Override the optimization-level forward/backward layout iterations.
    pub layout_iterations: Option<usize>,
    /// Number of successor interactions in the extended layer. The static
    /// heuristic normalizes each layer by its size. Zero selects the legacy
    /// immediate-successor layer and unnormalized score for ablations.
    pub lookahead_size: usize,
    /// Rank beam survivors by executed gates and remaining distance, after
    /// duplicate routing states have been collapsed.
    pub rank_progress: bool,
    /// Extract terminal logical permutations before layout and routing.
    pub elide_output_permutations: bool,
    /// Orient equivalent SWAP decompositions to expose adjacent CX cancellation.
    pub orient_swaps: bool,
    /// At level 3, allow diagonal interactions to commute through diagonal
    /// wire segments while retaining non-diagonal and barrier dependencies.
    pub commute_diagonal: bool,
}

impl Default for RoutingSearchConfig {
    fn default() -> Self {
        Self {
            layout_trials: None,
            layout_iterations: None,
            lookahead_size: 20,
            rank_progress: true,
            elide_output_permutations: true,
            orient_swaps: true,
            commute_diagonal: true,
        }
    }
}

/// Top-level transpiler configuration.
///
/// Construct via [`TranspilerConfig::builder`]. All fields are individually
/// addressable on the builder; defaults follow [`Default`].
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct TranspilerConfig {
    /// If `true`, decompose to the internal `{U, CX}` representation and
    /// enforce exact closure over the resolved target basis.
    pub decompose_basis: bool,
    /// Collect per-pass timings and output metrics in `transpile_with_report`.
    /// Disabled for normal compilation and performance measurements.
    pub profile_passes: bool,
    pub routing_search: RoutingSearchConfig,

    /// Run the level-3 post-routing optimization passes. Disable for cleanup
    /// ablations without changing the canonical lowering/validation pipeline.
    pub post_routing_cleanup: bool,
    /// Resynthesize dependency-safe two-qubit U/CX blocks during level-3
    /// post-routing cleanup. Accept only strict CX reductions within a
    /// phase-aligned numerical error budget. Disable for matched ablations.
    pub two_qubit_block_synthesis: bool,
    /// Recover exact lowered CX-RZ-CX parity rotations before level-3
    /// optimization when a backend is present. Disable for routing ablations.
    pub recover_phase_gadgets: bool,
    /// Fuse checked local U runs after level-3 two-qubit block synthesis.
    /// Disable to isolate the local-gate and depth cost of block replacements.
    pub post_block_fusion: bool,
    /// Optimization level (clamped to `[0, 3]` by the builder).
    /// - 0: no optimization passes.
    /// - 1: peephole crystallization + parameter simplification.
    /// - 2: + rotation merging, cross-conjugation, inverse cancellation,
    ///   commutation cancellation, layout (10 trials × 3 iters), routing.
    /// - 3: layout (10 trials × 3 iters), harder routing (beam=8,
    ///   branch=5, 4 bidirectional iterations), post-routing fusion, and
    ///   conservative two-qubit block resynthesis when basis lowering is on.
    ///   Use `routing_search` to increase the layout search budget.
    pub optimization_level: u8,
    /// Optional hardware backend. When set, layout, routing, and CX-direction
    /// passes are added to the pipeline.
    pub backend: Option<crate::backend::Backend>,
    /// Target gate set to translate into after all other passes complete.
    /// When `None`, falls back to the backend's `basis_gates` if any. Exact
    /// closure currently supports `U`, `RZ+RX`, `RZ+SX`, or `RZ+H` with CX
    /// or CZ; unsupported universal sets return an error rather than leaking
    /// gates.
    pub target_basis: Option<HashSet<String>>,
    /// Routing lookahead heuristic. Defaults to classical SABRE
    /// (`Static { weight: 0.5 }`); `DynamicV2` enables the SABRE-v2
    /// max-rule scoring (Li et al. 2023, arXiv:2210.12922).
    pub lookahead_strategy: routing::LookaheadStrategy,
}

impl Default for TranspilerConfig {
    fn default() -> Self {
        Self {
            decompose_basis: true,
            profile_passes: false,
            routing_search: RoutingSearchConfig::default(),
            post_routing_cleanup: true,
            two_qubit_block_synthesis: true,
            recover_phase_gadgets: true,
            post_block_fusion: true,
            optimization_level: 1,
            backend: None,
            target_basis: None,
            lookahead_strategy: routing::LookaheadStrategy::default(),
        }
    }
}

impl TranspilerConfig {
    /// Returns a fresh [`TranspilerConfigBuilder`].
    pub fn builder() -> TranspilerConfigBuilder {
        TranspilerConfigBuilder::default()
    }
}

/// Builder for [`TranspilerConfig`]. Method chaining is the intended idiom.
#[derive(Clone, Debug, Default)]
pub struct TranspilerConfigBuilder {
    decompose_basis: Option<bool>,
    profile_passes: Option<bool>,
    routing_search: Option<RoutingSearchConfig>,
    post_routing_cleanup: Option<bool>,
    two_qubit_block_synthesis: Option<bool>,
    recover_phase_gadgets: Option<bool>,
    post_block_fusion: Option<bool>,
    optimization_level: Option<u8>,
    backend: Option<crate::backend::Backend>,
    target_basis: Option<HashSet<String>>,
    lookahead_strategy: Option<routing::LookaheadStrategy>,
}

impl TranspilerConfigBuilder {
    /// Enable or disable decomposition and target-basis closure. Default:
    /// `true`.
    pub fn decompose_basis(mut self, v: bool) -> Self {
        self.decompose_basis = Some(v);
        self
    }
    /// Configure independent routing search policies and budgets.
    pub fn routing_search(mut self, options: RoutingSearchConfig) -> Self {
        self.routing_search = Some(options);
        self
    }
    /// Enable per-pass diagnostic capture (default: false).
    pub fn profile_passes(mut self, enabled: bool) -> Self {
        self.profile_passes = Some(enabled);
        self
    }
    /// Enable level-3 post-routing cleanup (default: true).
    pub fn post_routing_cleanup(mut self, enabled: bool) -> Self {
        self.post_routing_cleanup = Some(enabled);
        self
    }
    /// Enable two-qubit block resynthesis in level-3 cleanup (default: true).
    pub fn two_qubit_block_synthesis(mut self, enabled: bool) -> Self {
        self.two_qubit_block_synthesis = Some(enabled);
        self
    }
    /// Recover exact phase gadgets before level-3 routing (default: true).
    pub fn recover_phase_gadgets(mut self, enabled: bool) -> Self {
        self.recover_phase_gadgets = Some(enabled);
        self
    }
    /// Fuse adjacent local runs after block synthesis (default: true).
    pub fn post_block_fusion(mut self, enabled: bool) -> Self {
        self.post_block_fusion = Some(enabled);
        self
    }
    /// Set the optimization level. Clamped to `[0, 3]` at `build()` time.
    pub fn optimization_level(mut self, level: u8) -> Self {
        self.optimization_level = Some(level);
        self
    }
    /// Attach a hardware backend (enables layout + routing).
    pub fn backend(mut self, backend: crate::backend::Backend) -> Self {
        self.backend = Some(backend);
        self
    }
    /// Set the target basis explicitly, overriding the backend's `basis_gates`.
    /// See [`target_basis::BasisClosurePass`] for exact supported families.
    pub fn target_basis(mut self, basis: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.target_basis = Some(basis.into_iter().map(|s| s.into()).collect());
        self
    }
    /// Override the routing lookahead strategy.
    pub fn lookahead_strategy(mut self, s: routing::LookaheadStrategy) -> Self {
        self.lookahead_strategy = Some(s);
        self
    }
    /// Finalize the builder.
    pub fn build(self) -> TranspilerConfig {
        let default = TranspilerConfig::default();
        let level = self
            .optimization_level
            .unwrap_or(default.optimization_level)
            .min(3);
        TranspilerConfig {
            decompose_basis: self.decompose_basis.unwrap_or(default.decompose_basis),
            profile_passes: self.profile_passes.unwrap_or(false),
            routing_search: self.routing_search.unwrap_or_default(),
            post_routing_cleanup: self
                .post_routing_cleanup
                .unwrap_or(default.post_routing_cleanup),
            two_qubit_block_synthesis: self
                .two_qubit_block_synthesis
                .unwrap_or(default.two_qubit_block_synthesis),
            recover_phase_gadgets: self
                .recover_phase_gadgets
                .unwrap_or(default.recover_phase_gadgets),
            post_block_fusion: self.post_block_fusion.unwrap_or(default.post_block_fusion),
            optimization_level: level,
            backend: self.backend,
            target_basis: self.target_basis,
            lookahead_strategy: self
                .lookahead_strategy
                .unwrap_or(default.lookahead_strategy),
        }
    }
}

fn diagnostics_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| match std::env::var("Q_RUST_LOG") {
        Ok(v) => !v.is_empty() && v != "0" && !v.eq_ignore_ascii_case("off"),
        Err(_) => false,
    })
}

#[doc(hidden)]
pub(crate) fn warn_diagnostic(msg: std::fmt::Arguments<'_>) {
    if diagnostics_enabled() {
        eprintln!("[Q-Rust][warn] {msg}");
    }
}

#[derive(Debug, Clone, Copy)]
pub struct KakSynthesisPass;

impl Pass for KakSynthesisPass {
    fn name(&self) -> &str {
        "KakSynthesisPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut property_set::PropertySet) -> Circuit {
        use crate::ir::GateDefinition;
        use crate::transpiler::synthesis::kak::KakSynthesizer;
        use crate::transpiler::synthesis::Synthesizer;

        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        let mut kak_invocations: usize = 0;
        let mut analytic_contract_violations: usize = 0;

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } if qubits.len() == 2 && !matches!(name, GateType::CX) => {
                    if let Some(decomp_ops) = name.decompose(qubits, params) {
                        let n = circuit.num_qubits;
                        let invalid = decomp_ops.iter().any(|op| {
                            if let Operation::Gate { qubits: dq, .. } = op {
                                dq.iter().any(|&q| q >= n)
                            } else {
                                false
                            }
                        });
                        if invalid {
                            analytic_contract_violations += 1;
                            warn_diagnostic(format_args!(
                                "KakSynthesisPass: decompose() emitted out-of-range qubit \
                                 indices for gate {name:?} on qubits {qubits:?} (n={n}); \
                                 falling back to KAK"
                            ));
                        } else {
                            for sub in decomp_ops {
                                out.add_op(sub);
                            }
                            continue;
                        }
                    }

                    let u = name.unitary(params);
                    if u.nrows() == 4 && u.ncols() == 4 {
                        if let Some(sub) = KakSynthesizer.synthesize(&u, &[]) {
                            // KAK validates every returned reconstruction at
                            // 1e-10. Reserve that conservative bound across
                            // all invocations instead of resetting per gate.
                            if !numerical_policy::reserve(
                                property_set,
                                1e-10,
                                1e-10,
                                sub.operations.len() + 1,
                            ) {
                                out.add_op(op.clone());
                                continue;
                            }
                            kak_invocations += 1;
                            for sop in sub.operations {
                                if let Operation::Gate {
                                    name: sn,
                                    qubits: sq,
                                    params: sp,
                                } = sop
                                {
                                    let remapped: Vec<usize> =
                                        sq.iter().map(|&q| qubits[q]).collect();
                                    out.add_op(Operation::Gate {
                                        name: sn,
                                        qubits: remapped,
                                        params: sp,
                                    });
                                }
                            }
                            continue;
                        }
                    }
                    out.add_op(op.clone());
                }
                other => out.add_op(other.clone()),
            }
        }

        property_set.insert("kak_fallback_invocations", kak_invocations);
        if analytic_contract_violations > 0 {
            property_set.insert(
                "kak_analytic_contract_violations",
                analytic_contract_violations,
            );
        }

        if kak_invocations > 0 {
            warn_diagnostic(format_args!(
                "KakSynthesisPass: {kak_invocations} expensive KAK fallback(s) invoked"
            ));
        }
        out
    }
}

#[derive(Debug, Clone)]
/// Compatibility wrapper that closes a circuit over a backend's basis.
///
/// New code should prefer [`target_basis::BasisClosurePass::new`], whose
/// constructor is fallible. Because [`Pass::run`] cannot return an error, this
/// wrapper leaves unsupported input unchanged and records the error string as
/// `"native_basis_translation_error"` in the property set.
pub struct NativeBasisTranslationPass {
    pub backend: crate::backend::Backend,
}

impl Pass for NativeBasisTranslationPass {
    fn name(&self) -> &str {
        "NativeBasisTranslationPass"
    }

    fn run(&self, circuit: &Circuit, property_set: &mut property_set::PropertySet) -> Circuit {
        match target_basis::BasisClosurePass::new(self.backend.basis_gates.clone()) {
            Ok(pass) => pass.run(circuit, property_set),
            Err(error) => {
                property_set.insert("native_basis_translation_error", error.to_string());
                circuit.clone()
            }
        }
    }

    fn try_run(
        &self,
        circuit: &Circuit,
        property_set: &mut property_set::PropertySet,
    ) -> Result<Circuit> {
        target_basis::BasisClosurePass::new(self.backend.basis_gates.clone())?
            .try_run(circuit, property_set)
    }
}

/// Stage groups that the pipeline is partitioned into. This split is what
/// lets [`transpile_with_report`] capture an intermediate snapshot without
/// re-running the optimization pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Stage {
    /// Optimization passes (levels 1–2 entries).
    Optimize,
    /// Layout, routing, basis translation, decomposition, and post-routing
    /// (level-3) passes — i.e. everything that depends on a backend or
    /// target-basis being known.
    LayoutAndLower,
}

fn resolved_target_basis(config: &TranspilerConfig) -> Option<HashSet<String>> {
    config.target_basis.clone().or_else(|| {
        config
            .backend
            .as_ref()
            .filter(|backend| !backend.basis_gates.is_empty())
            .map(|backend| backend.basis_gates.clone())
    })
}

/// Validate public, mutable IR before any pass can index operands or use
/// default values for missing parameters. This is also used by simulation.
pub(crate) fn validate_circuit_structure(circuit: &Circuit) -> Result<()> {
    fn check_gate(
        circuit: &Circuit,
        name: &GateType,
        qubits: &[usize],
        params: &[f64],
        active_custom: &mut HashSet<String>,
    ) -> Result<()> {
        let (arity, parameter_count) = match name {
            GateType::Custom(name) => {
                let def = circuit
                    .custom_gates
                    .get(name)
                    .ok_or_else(|| QRustError::UnknownGate(name.clone()))?;
                (def.qubits.len(), def.params.len())
            }
            GateType::Barrier => {
                return Err(QRustError::InvalidConfig(
                    "represent barriers as Operation::Barrier, not as a gate".into(),
                ))
            }
            GateType::U => (1, 3),
            GateType::RX
            | GateType::RY
            | GateType::RZ
            | GateType::CRX
            | GateType::CRY
            | GateType::CRZ
            | GateType::CU1
            | GateType::RXX
            | GateType::RYY
            | GateType::RZZ => (name.num_qubits(), 1),
            _ => (name.num_qubits(), 0),
        };
        if qubits.len() != arity || params.len() != parameter_count {
            return Err(QRustError::InvalidConfig(format!(
                "gate {name} requires {arity} qubits and {parameter_count} parameters; got {} and {}",
                qubits.len(), params.len()
            )));
        }
        if params.iter().any(|v| !v.is_finite()) {
            return Err(QRustError::InvalidConfig(format!(
                "gate {name} has non-finite parameters"
            )));
        }
        if qubits.iter().any(|&q| q >= circuit.num_qubits)
            || qubits.iter().copied().collect::<HashSet<_>>().len() != qubits.len()
        {
            return Err(QRustError::InvalidConfig(format!(
                "gate {name} has invalid or repeated qubits"
            )));
        }
        if let GateType::Custom(custom) = name {
            if !active_custom.insert(custom.clone()) {
                return Err(QRustError::Decomposition(format!(
                    "recursive custom gate {custom}"
                )));
            }
            let def = circuit
                .custom_gates
                .get(custom)
                .expect("definition checked above");
            if def.qubits.iter().collect::<HashSet<_>>().len() != def.qubits.len()
                || def.params.iter().collect::<HashSet<_>>().len() != def.params.len()
            {
                return Err(QRustError::Decomposition(format!(
                    "duplicate formal arguments in custom gate {custom}"
                )));
            }
            let scope = def
                .params
                .iter()
                .cloned()
                .zip(params.iter().copied())
                .collect();
            for statement in &def.body {
                let crate::ir::ast::ParsedStatement::Gate(subname, operands, expressions) =
                    statement
                else {
                    return Err(QRustError::Unsupported(format!(
                        "non-gate statement in custom gate {custom}"
                    )));
                };
                let mut subqubits = Vec::with_capacity(operands.len());
                for (formal, index) in operands {
                    if index.is_some() {
                        return Err(QRustError::Decomposition(format!(
                            "indexed formal qubit in custom gate {custom}"
                        )));
                    }
                    let position = def
                        .qubits
                        .iter()
                        .position(|q| q == formal)
                        .ok_or_else(|| QRustError::Undefined(formal.clone()))?;
                    subqubits.push(qubits[position]);
                }
                let subparams: Vec<f64> = expressions
                    .iter()
                    .map(|e| e.evaluate_with_scope(&scope))
                    .collect::<Result<_>>()?;
                let subgate = subname
                    .parse::<GateType>()
                    .unwrap_or_else(|_| GateType::Custom(subname.clone()));
                check_gate(circuit, &subgate, &subqubits, &subparams, active_custom)?;
            }
            active_custom.remove(custom);
        }
        Ok(())
    }

    fn check_op(circuit: &Circuit, op: &Operation) -> Result<()> {
        if op.qubits().iter().any(|&q| q >= circuit.num_qubits) {
            return Err(QRustError::InvalidConfig(
                "operation qubit index exceeds circuit width".into(),
            ));
        }
        match op {
            Operation::Gate {
                name,
                qubits,
                params,
            } => check_gate(circuit, name, qubits, params, &mut HashSet::new()),
            Operation::Measure { cbit, .. } if *cbit >= circuit.num_cbits => {
                Err(QRustError::InvalidConfig(
                    "measurement classical index exceeds circuit width".into(),
                ))
            }
            Operation::Conditional { op, .. } => check_op(circuit, op),
            _ => Ok(()),
        }
    }
    for op in &circuit.operations {
        check_op(circuit, op)?;
    }
    Ok(())
}

/// Supported compilation domain: gates/barriers followed by terminal readout.
/// Reject unsupported dynamics before optimization can erase their effects.
pub(crate) fn validate_static_semantics(circuit: &Circuit) -> Result<()> {
    let mut measured = false;
    for op in &circuit.operations {
        match op {
            Operation::Measure { .. } => measured = true,
            Operation::Gate { .. } if measured => {
                return Err(QRustError::Unsupported(
                    "mid-circuit measurement: compilation supports terminal measurements only"
                        .into(),
                ))
            }
            Operation::Reset { .. } | Operation::Conditional { .. } => {
                return Err(QRustError::Unsupported(
                    "reset and classical conditional operations are not supported by the compiler"
                        .into(),
                ))
            }
            _ => {}
        }
    }
    Ok(())
}

fn validate_input(circuit: &Circuit, config: &TranspilerConfig) -> Result<()> {
    validate_circuit_structure(circuit)?;
    validate_static_semantics(circuit)?;
    if let Some(backend) = &config.backend {
        if circuit.num_qubits > backend.num_qubits {
            return Err(QRustError::InsufficientQubits {
                circuit: circuit.num_qubits,
                backend: backend.num_qubits,
            });
        }
    }
    Ok(())
}

/// Independent postcondition after every pass, including basis closure.
fn validate_output(circuit: &Circuit, config: &TranspilerConfig) -> Result<()> {
    validate_input(circuit, config)?;
    if let Some(backend) = &config.backend {
        for op in &circuit.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                if qubits.len() > 2 || matches!(name, GateType::Custom(_)) {
                    return Err(QRustError::Unsupported(format!(
                        "hardware output contains unresolved gate {name}"
                    )));
                }
                if qubits.len() == 2 {
                    if !backend.is_adjacent(qubits[0], qubits[1]) {
                        return Err(QRustError::Routing(format!(
                            "output gate {name} uses nonadjacent qubits {} and {}",
                            qubits[0], qubits[1]
                        )));
                    }
                    // Symmetric gates can use either orientation. Directed
                    // controlled gates must match the ordered hardware edge.
                    if !matches!(
                        name,
                        GateType::CZ
                            | GateType::SWAP
                            | GateType::RXX
                            | GateType::RYY
                            | GateType::RZZ
                            | GateType::ISwap
                    ) && !backend.has_directed_edge(qubits[0], qubits[1])
                    {
                        return Err(QRustError::Routing(format!(
                            "output gate {name} uses unsupported direction {} -> {}",
                            qubits[0], qubits[1]
                        )));
                    }
                }
            }
        }
    }
    // A disabled lowering pass is no permission to return an illegal gate set.
    if config.decompose_basis || config.backend.is_some() {
        if let Some(basis) = resolved_target_basis(config) {
            target_basis::validate_circuit_basis(circuit, &basis)?;
        }
    }
    Ok(())
}

fn split_terminal_readout(circuit: &Circuit) -> (Circuit, &[Operation]) {
    let split = circuit
        .operations
        .iter()
        .position(|op| matches!(op, Operation::Measure { .. }))
        .unwrap_or(circuit.operations.len());
    let mut prefix = circuit.clone();
    prefix.operations.truncate(split);
    (prefix, &circuit.operations[split..])
}

fn restore_terminal_readout(
    mut circuit: Circuit,
    terminal: &[Operation],
    layout: Option<&Vec<usize>>,
) -> Circuit {
    let physical = |q: usize| layout.map_or(q, |mapping| mapping[q]);
    for op in terminal {
        match op {
            Operation::Measure { qubit, cbit } => circuit.add_op(Operation::Measure {
                qubit: physical(*qubit),
                cbit: *cbit,
            }),
            Operation::Barrier { qubits } => circuit.add_op(Operation::Barrier {
                qubits: qubits.iter().map(|&q| physical(q)).collect(),
            }),
            _ => unreachable!("terminal readout is validated before compilation"),
        }
    }
    circuit
}

/// Builds and runs a slice of the pipeline. Splitting `build_pass_manager`
/// into stage-aware halves lets `transpile_with_report` produce a
/// mid-pipeline snapshot whose continuation **is bit-for-bit identical to
/// what `transpile()` would have produced** — because we never re-run
/// optimization a second time.
fn build_pass_manager_for(config: &TranspilerConfig, stage: Stage) -> Result<PassManager> {
    let mut pm = PassManager::new();
    pm.property_set.insert(
        "routing_lookahead_size",
        config.routing_search.lookahead_size,
    );
    pm.property_set
        .insert("routing_rank_progress", config.routing_search.rank_progress);
    pm.property_set.insert(
        "routing_commute_diagonal",
        config.optimization_level >= 3 && config.routing_search.commute_diagonal,
    );

    match stage {
        Stage::Optimize => {
            if config.optimization_level >= 3
                && config.backend.is_some()
                && config.recover_phase_gadgets
            {
                // Recover exact parity rotations before crystallization turns
                // RZ angles into named phase gates, retaining their two-wire
                // interaction structure throughout layout and routing.
                pm.add_pass(Box::new(phase_gadgets::RecoverPhaseGadgetsPass));
            }
            if config.optimization_level >= 1 {
                pm.add_pass(Box::new(optimization::GateCrystallizationPass::default()));
                pm.add_pass(Box::new(
                    optimization::ParameterSimplificationPass::default(),
                ));
            }

            if config.optimization_level >= 2 {
                pm.add_pass(Box::new(optimization::RotationMergePass));
                pm.add_pass(Box::new(optimization::CrossConjugationPass));
                pm.add_pass(Box::new(optimization::InverseCancellationPass));
                pm.add_pass(Box::new(optimization::CommutationCancellationPass));
                pm.add_pass(Box::new(
                    optimization::ParameterSimplificationPass::default(),
                ));
            }
        }
        Stage::LayoutAndLower => {
            if let Some(ref backend) = config.backend {
                // Multi-qubit gates (e.g. CCX) must be reduced to ≤2-qubit
                // operations before routing: the SWAP router cannot make a
                // 3-qubit interaction adjacent and would otherwise corrupt the
                // gate. Standard transpiler ordering (decompose, then route).
                if config.routing_search.elide_output_permutations {
                    pm.add_pass(Box::new(optimization::ExtractOutputPermutationPass));
                }
                pm.add_pass(Box::new(decomposition::HighArityDecompositionPass));

                if config.optimization_level >= 2 {
                    // More trials remain available explicitly. On the
                    // development corpus, 10×3 matched 50×5 output exactly;
                    // this is a search-budget choice, not a universal bound.
                    let (num_trials, num_iterations) = (10, 3);
                    pm.add_pass(Box::new(layout::SabreLayoutPass {
                        backend: backend.clone(),
                        num_trials: config.routing_search.layout_trials.unwrap_or(num_trials),
                        num_iterations: config
                            .routing_search
                            .layout_iterations
                            .unwrap_or(num_iterations),
                    }));
                }
                let (beam_width, branch_factor, bidir_iters) = match config.optimization_level {
                    0 | 1 => (1, 1, 1),
                    2 => (4, 3, 2),
                    _ => (8, 5, 4),
                };
                pm.add_pass(Box::new(routing::BeamSabrePass {
                    backend: backend.clone(),
                    beam_width,
                    branch_factor,
                    bidir_iterations: bidir_iters,
                    lookahead_strategy: config.lookahead_strategy,
                }));
                // Fold any trailing SWAPs (e.g. the QFT bit-reversal network)
                // into the output layout rather than emitting them as gates.
                pm.add_pass(Box::new(optimization::TrailingSwapElisionPass));
                if config.routing_search.elide_output_permutations {
                    pm.add_pass(Box::new(optimization::ApplyOutputPermutationPass));
                }
            }

            // Target-basis translation runs BEFORE BasisDecompositionPass so the
            // equivalence rules operate on original named gates (H, CZ, SWAP, etc.)
            // rather than their U-gate expansions.
            if config.decompose_basis {
                let resolved_basis = resolved_target_basis(config);
                if config.optimization_level >= 3
                    && config.post_routing_cleanup
                    && config.routing_search.orient_swaps
                {
                    pm.add_pass(Box::new(optimization::OrientSwapDecompositionPass));
                }
                if let Some(ref basis) = resolved_basis {
                    pm.add_pass(Box::new(target_basis::TargetBasisPass::new(basis.clone())?));
                }

                pm.add_pass(Box::new(decomposition::BasisDecompositionPass));
                pm.add_pass(Box::new(KakSynthesisPass));

                if config.optimization_level >= 3 && config.post_routing_cleanup {
                    pm.add_pass(Box::new(optimization::GateFusionPass));
                    pm.add_pass(Box::new(optimization::SwapSimplificationPass));
                    pm.add_pass(Box::new(optimization::InverseCancellationPass));
                    pm.add_pass(Box::new(
                        optimization::ParameterSimplificationPass::default(),
                    ));
                    if config.two_qubit_block_synthesis {
                        // Operate on the already routed U/CX sequence, keeping
                        // the placement problem fixed. Direction repair and
                        // final target closure still run after synthesis.
                        pm.add_pass(Box::new(block_synthesis::TwoQubitBlockSynthesisPass));
                        if config.post_block_fusion {
                            pm.add_pass(Box::new(local_cleanup::PostBlockFusionPass));
                        }
                    }
                }

                // SWAP/controlled-gate lowering and cleanup can create fresh
                // CX gates. Correct their direction only after these passes.
                if let Some(backend) = &config.backend {
                    pm.add_pass(Box::new(decomposition::CxDirectionPass {
                        backend: backend.clone(),
                    }));
                    // Direction repair introduces H gates; lower only these
                    // remaining named gates to U before final basis closure.
                    pm.add_pass(Box::new(decomposition::BasisDecompositionPass));
                }
                if let Some(basis) = resolved_basis {
                    pm.add_pass(Box::new(target_basis::BasisClosurePass::new(basis)?));
                }
            } else if let Some(backend) = &config.backend {
                pm.add_pass(Box::new(decomposition::CxDirectionPass {
                    backend: backend.clone(),
                }));
            }
        }
    }

    Ok(pm)
}

/// Compile a well-formed circuit and validate the final gate basis and hardware
/// constraints. Supports unitary gates and barriers with optional terminal
/// measurements; reset, conditional operations, and mid-circuit measurement
/// currently return [`QRustError::Unsupported`] before any transformation.
pub fn transpile(circuit: &Circuit, config: Option<TranspilerConfig>) -> Result<Circuit> {
    let config = config.unwrap_or_default();
    validate_input(circuit, &config)?;
    // Optimize only the unitary prefix. Independent measurements must not be
    // pulled earlier by DAG scheduling, and their classical write order must
    // survive routing and trailing-SWAP elision unchanged.
    let (prefix, terminal) = split_terminal_readout(circuit);
    let mut opt_pm = build_pass_manager_for(&config, Stage::Optimize)?;
    let after_opt = opt_pm.try_run(&prefix)?;
    let mut lower_pm = build_pass_manager_for(&config, Stage::LayoutAndLower)?;
    numerical_policy::transfer(&mut opt_pm.property_set, &mut lower_pm.property_set);
    let lowered = lower_pm.try_run(&after_opt)?;
    let output = restore_terminal_readout(
        lowered,
        terminal,
        lower_pm.property_set.get::<Vec<usize>>("final_layout"),
    );
    validate_output(&output, &config)?;
    Ok(output)
}

/// Like [`transpile`], but additionally returns a
/// [`TranspilationReport`] describing per-stage circuit metrics.
///
/// Captures three checkpoints — input, post-optimization, and final —
/// plus the initial/final layouts recorded by the router, when present,
/// and is guaranteed to produce the **same final circuit as
/// [`transpile`]** for the same input/config (the pipeline is split, not
/// re-run). Set `profile_passes` to collect per-pass timing and gate metrics
/// separately from performance measurements.
pub fn transpile_with_report(
    circuit: &Circuit,
    config: Option<TranspilerConfig>,
) -> Result<(Circuit, TranspilationReport)> {
    let config = config.unwrap_or_default();
    validate_input(circuit, &config)?;
    let (prefix, terminal) = split_terminal_readout(circuit);
    let mut report = TranspilationReport::new();
    report.push(StageSnapshot::capture("1. parsed", circuit));

    let mut opt_pm = build_pass_manager_for(&config, Stage::Optimize)?;
    let after_opt = if config.profile_passes {
        let (output, diagnostic) = opt_pm.try_run_profiled(&prefix)?;
        report.pass_diagnostics.extend(diagnostic);
        output
    } else {
        opt_pm.try_run(&prefix)?
    };
    report.push(StageSnapshot::capture(
        "2. optimized",
        &restore_terminal_readout(after_opt.clone(), terminal, None),
    ));

    let mut lower_pm = build_pass_manager_for(&config, Stage::LayoutAndLower)?;
    numerical_policy::transfer(&mut opt_pm.property_set, &mut lower_pm.property_set);
    let lowered = if config.profile_passes {
        let (output, diagnostic) = lower_pm.try_run_profiled(&after_opt)?;
        report.pass_diagnostics.extend(diagnostic);
        output
    } else {
        lower_pm.try_run(&after_opt)?
    };
    let final_circuit = restore_terminal_readout(
        lowered,
        terminal,
        lower_pm.property_set.get::<Vec<usize>>("final_layout"),
    );
    validate_output(&final_circuit, &config)?;
    let final_label = if config.backend.is_some() {
        "3. routed+decomposed"
    } else {
        "3. decomposed"
    };
    report.push(StageSnapshot::capture(final_label, &final_circuit));
    report.initial_layout = lower_pm
        .property_set
        .get::<Vec<usize>>("initial_layout")
        .cloned();
    report.final_layout = lower_pm
        .property_set
        .get::<Vec<usize>>("final_layout")
        .cloned();
    report.two_qubit_block_stats = lower_pm
        .property_set
        .get::<block_synthesis::TwoQubitBlockStats>("two_qubit_block_stats")
        .cloned();
    report.phase_gadgets_recovered = opt_pm
        .property_set
        .get::<usize>("phase_gadgets_recovered")
        .copied();
    report.post_block_fusion_stats = lower_pm
        .property_set
        .get::<local_cleanup::PostBlockFusionStats>("post_block_fusion_stats")
        .cloned();
    report.numerical_budget = numerical_policy::snapshot(&lower_pm.property_set);

    Ok((final_circuit, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::Backend;
    use crate::ir::{GateType, Operation};

    #[test]
    fn test_transpile_default() {
        let mut c = Circuit::new(1, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let t = transpile(&c, None).expect("transpile failed");
        assert_eq!(t.operations.len(), 1);
    }

    #[test]
    fn test_builder_clamps_level() {
        let cfg = TranspilerConfig::builder().optimization_level(255).build();
        assert_eq!(cfg.optimization_level, 3);
    }

    #[test]
    fn test_transpile_no_state_leak_between_calls() {
        let mut c1 = Circuit::new(2, 0);
        c1.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let cfg1 = TranspilerConfig::builder()
            .optimization_level(2)
            .backend(Backend::linear(3))
            .build();
        let _ = transpile(&c1, Some(cfg1));

        let mut c2 = Circuit::new(1, 0);
        c2.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        let cfg2 = TranspilerConfig::default();
        let out2 = transpile(&c2, Some(cfg2)).expect("transpile failed");
        assert!(!out2.operations.is_empty());
    }

    #[test]
    fn test_kak_pass_uses_analytic_path_for_known_gates() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CRZ,
            qubits: vec![0, 1],
            params: vec![1.234],
        });
        c.add_op(Operation::Gate {
            name: GateType::CZ,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        let mut ps = property_set::PropertySet::new();
        let out = KakSynthesisPass.run(&c, &mut ps);
        for op in &out.operations {
            if let Operation::Gate { name, qubits, .. } = op {
                if qubits.len() == 2 {
                    assert!(matches!(name, GateType::CX));
                }
            }
        }
    }

    #[test]
    fn test_builder_threads_lookahead_strategy() {
        let cfg = TranspilerConfig::builder()
            .lookahead_strategy(routing::LookaheadStrategy::DynamicV2)
            .build();
        assert_eq!(
            cfg.lookahead_strategy,
            routing::LookaheadStrategy::DynamicV2
        );
        let cfg2 = TranspilerConfig::builder().build();
        assert_eq!(
            cfg2.lookahead_strategy,
            routing::LookaheadStrategy::default()
        );
    }

    #[test]
    fn test_kak_pass_persists_invocation_count_when_analytic_path_taken() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::CRZ,
            qubits: vec![0, 1],
            params: vec![0.5],
        });
        let mut ps = property_set::PropertySet::new();
        let _ = KakSynthesisPass.run(&c, &mut ps);
        assert_eq!(ps.get::<usize>("kak_fallback_invocations"), Some(&0));
    }

    #[test]
    fn test_kak_pass_uses_analytic_path_for_ecr_and_iswap() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::ECR,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::ISwap,
            qubits: vec![0, 1],
            params: vec![],
        });
        let mut ps = property_set::PropertySet::new();
        let _ = KakSynthesisPass.run(&c, &mut ps);
        assert_eq!(ps.get::<usize>("kak_fallback_invocations"), Some(&0));
    }

    /// Report-producing variant returns a non-empty report with three
    /// stable-named stages.
    #[test]
    fn test_transpile_with_report_produces_report() {
        let mut c = Circuit::new(2, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        let cfg = TranspilerConfig::builder().optimization_level(1).build();
        let (out, report) = transpile_with_report(&c, Some(cfg)).expect("ok");
        assert_eq!(report.stages.len(), 3);
        assert_eq!(report.stages[0].stage, "1. parsed");
        assert_eq!(report.stages[1].stage, "2. optimized");
        assert_eq!(report.stages[2].stage, "3. decomposed");
        assert_eq!(report.stages[0].num_ops, 2);
        assert!(!out.operations.is_empty());
        assert_eq!(report.format_lines().len(), 4);
    }

    /// `transpile_with_report` must produce the
    /// same final circuit as `transpile`. This pins the contract that the
    /// pipeline is split, not run twice (which would yield different
    /// outputs whenever a pass uses non-determinism or a property-set key
    /// from an earlier pass).
    #[test]
    fn test_transpile_with_report_matches_transpile_output() {
        let mut c = Circuit::new(3, 0);
        c.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        c.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![1, 2],
            params: vec![],
        });
        let cfg = TranspilerConfig::builder()
            .optimization_level(2)
            .decompose_basis(true)
            .backend(Backend::linear(3))
            .build();
        let plain = transpile(&c, Some(cfg.clone())).expect("transpile");
        let (with_report, _) = transpile_with_report(&c, Some(cfg)).expect("transpile_with_report");
        assert_eq!(plain.num_qubits, with_report.num_qubits);
        assert_eq!(plain.operations.len(), with_report.operations.len());
        assert_eq!(plain.operations, with_report.operations);
    }
}
