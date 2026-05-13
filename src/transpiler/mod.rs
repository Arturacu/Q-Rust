//! Transpiler: optimization, decomposition, layout, and routing.
//!
//! See [`transpile`] for the canonical entry point and
//! [`TranspilerConfig::builder`] for configuring the pipeline.

pub mod dag;
pub mod decomposition;
pub mod layout;
pub mod optimization;
pub mod pass;
pub mod pauli_tracker;
pub mod profiler;
pub mod property_set;
pub mod report;
pub mod routing;
pub mod synthesis;
pub mod target_basis;

pub use report::{StageSnapshot, TranspilationReport};

use crate::error::Result;
use crate::ir::{Circuit, GateType, Operation};
use pass::{Pass, PassManager};
use std::collections::HashSet;

/// Top-level transpiler configuration.
///
/// Construct via [`TranspilerConfig::builder`]. All fields are individually
/// addressable on the builder; defaults follow [`Default`].
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct TranspilerConfig {
    /// If `true`, decompose every non-basis gate into the target basis as
    /// the final stage of the pipeline.
    pub decompose_basis: bool,
    /// Optimization level (clamped to `[0, 3]` by the builder).
    /// - 0: no optimization passes.
    /// - 1: peephole crystallization + parameter simplification.
    /// - 2: + rotation merging, cross-conjugation, inverse cancellation,
    ///   commutation cancellation, layout (10 trials × 3 iters), routing.
    /// - 3: + harder routing (beam=8) and post-routing fusion.
    pub optimization_level: u8,
    /// Optional hardware backend. When set, layout, routing, and CX-direction
    /// passes are added to the pipeline.
    pub backend: Option<crate::backend::Backend>,
    /// Target gate set to translate into after all other passes complete.
    /// When `None`, falls back to the backend's `basis_gates` if any.
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
    optimization_level: Option<u8>,
    backend: Option<crate::backend::Backend>,
    target_basis: Option<HashSet<String>>,
    lookahead_strategy: Option<routing::LookaheadStrategy>,
}

impl TranspilerConfigBuilder {
    /// Enable or disable basis decomposition. Default: `true`.
    pub fn decompose_basis(mut self, v: bool) -> Self {
        self.decompose_basis = Some(v);
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
    /// Set the target basis gate set explicitly. Overrides the backend's
    /// `basis_gates` when both are present.
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
pub struct NativeBasisTranslationPass {
    pub backend: crate::backend::Backend,
}

impl Pass for NativeBasisTranslationPass {
    fn name(&self) -> &str {
        "NativeBasisTranslationPass"
    }

    fn run(&self, circuit: &Circuit, _property_set: &mut property_set::PropertySet) -> Circuit {
        let native_u = self.backend.basis_gates.contains("u")
            || self.backend.basis_gates.contains("u3")
            || self.backend.basis_gates.contains("U");
        if native_u || self.backend.basis_gates.is_empty() {
            return circuit.clone();
        }

        let has_rz = self.backend.basis_gates.contains("rz");
        let has_sx = self.backend.basis_gates.contains("sx");
        if !has_rz || !has_sx {
            return circuit.clone();
        }

        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();
        let pi = std::f64::consts::PI;

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name: GateType::U,
                    qubits,
                    params,
                } if qubits.len() == 1 && params.len() >= 3 => {
                    let theta = params[0];
                    let phi = params[1];
                    let lambda = params[2];
                    let q = qubits[0];
                    out.add_op(Operation::Gate {
                        name: GateType::RZ,
                        qubits: vec![q],
                        params: vec![lambda - pi],
                    });
                    out.add_op(Operation::Gate {
                        name: GateType::Custom("sx".into()),
                        qubits: vec![q],
                        params: vec![],
                    });
                    out.add_op(Operation::Gate {
                        name: GateType::RZ,
                        qubits: vec![q],
                        params: vec![theta],
                    });
                    out.add_op(Operation::Gate {
                        name: GateType::Custom("sx".into()),
                        qubits: vec![q],
                        params: vec![],
                    });
                    out.add_op(Operation::Gate {
                        name: GateType::RZ,
                        qubits: vec![q],
                        params: vec![phi - pi],
                    });
                }
                other => out.add_op(other.clone()),
            }
        }
        out
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

/// Builds and runs a slice of the pipeline. Splitting `build_pass_manager`
/// into stage-aware halves lets `transpile_with_report` produce a
/// mid-pipeline snapshot whose continuation **is bit-for-bit identical to
/// what `transpile()` would have produced** — because we never re-run
/// optimization a second time.
fn build_pass_manager_for(config: &TranspilerConfig, stage: Stage) -> Result<PassManager> {
    let mut pm = PassManager::new();

    match stage {
        Stage::Optimize => {
            if config.optimization_level >= 1 {
                pm.add_pass(Box::new(optimization::GateCrystallizationPass {
                    epsilon: 1e-9,
                }));
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
                if config.optimization_level >= 2 {
                    let (num_trials, num_iterations) = match config.optimization_level {
                        2 => (10, 3),
                        _ => (50, 5),
                    };
                    pm.add_pass(Box::new(layout::SabreLayoutPass {
                        backend: backend.clone(),
                        num_trials,
                        num_iterations,
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
                pm.add_pass(Box::new(decomposition::CxDirectionPass {
                    backend: backend.clone(),
                }));
            }

            // Target-basis translation runs BEFORE BasisDecompositionPass so the
            // equivalence rules operate on original named gates (H, CZ, SWAP, etc.)
            // rather than their U-gate expansions.
            let resolved_basis: Option<HashSet<String>> =
                config.target_basis.clone().or_else(|| {
                    config
                        .backend
                        .as_ref()
                        .filter(|b| !b.basis_gates.is_empty())
                        .map(|b| b.basis_gates.clone())
                });
            if let Some(ref basis) = resolved_basis {
                let tb_pass = target_basis::TargetBasisPass::new(basis.clone())?;
                pm.add_pass(Box::new(tb_pass));
            }

            if config.decompose_basis {
                pm.add_pass(Box::new(decomposition::BasisDecompositionPass));
                pm.add_pass(Box::new(KakSynthesisPass));
            }

            if config.optimization_level >= 3 {
                pm.add_pass(Box::new(optimization::GateFusionPass));
                pm.add_pass(Box::new(optimization::SwapSimplificationPass));
                pm.add_pass(Box::new(optimization::InverseCancellationPass));
                pm.add_pass(Box::new(
                    optimization::ParameterSimplificationPass::default(),
                ));
            }
        }
    }

    Ok(pm)
}

pub fn transpile(circuit: &Circuit, config: Option<TranspilerConfig>) -> Result<Circuit> {
    let config = config.unwrap_or_default();
    let mut opt_pm = build_pass_manager_for(&config, Stage::Optimize)?;
    let after_opt = opt_pm.run(circuit);
    let mut lower_pm = build_pass_manager_for(&config, Stage::LayoutAndLower)?;
    Ok(lower_pm.run(&after_opt))
}

/// [E2E-NEW-FEATURE] Like [`transpile`], but additionally returns a
/// [`TranspilationReport`] describing per-stage circuit metrics.
///
/// Captures three checkpoints — input, post-optimization, and final —
/// and is guaranteed to produce the **same final circuit as
/// [`transpile`]** for the same input/config (the pipeline is split, not
/// re-run). More granular per-pass tracking would require reworking
/// [`PassManager`] (deferred to G-PM-01).
pub fn transpile_with_report(
    circuit: &Circuit,
    config: Option<TranspilerConfig>,
) -> Result<(Circuit, TranspilationReport)> {
    let config = config.unwrap_or_default();
    let mut report = TranspilationReport::new();
    report.push(StageSnapshot::capture("1. parsed", circuit));

    let mut opt_pm = build_pass_manager_for(&config, Stage::Optimize)?;
    let after_opt = opt_pm.run(circuit);
    report.push(StageSnapshot::capture("2. optimized", &after_opt));

    let mut lower_pm = build_pass_manager_for(&config, Stage::LayoutAndLower)?;
    let final_circuit = lower_pm.run(&after_opt);
    let final_label = if config.backend.is_some() {
        "3. routed+decomposed"
    } else {
        "3. decomposed"
    };
    report.push(StageSnapshot::capture(final_label, &final_circuit));

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

    /// [E2E-NEW-FEATURE] Report-producing variant returns a non-empty report
    /// with three stable-named stages.
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

    /// [E2E-NEW-FEATURE] Critical: `transpile_with_report` must produce the
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
