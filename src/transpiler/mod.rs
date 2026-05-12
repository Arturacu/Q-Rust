//! Transpiler: optimization, decomposition, layout, and routing.

pub mod dag;
pub mod decomposition;
pub mod layout;
pub mod optimization;
pub mod pass;
pub mod pauli_tracker;
pub mod profiler;
pub mod property_set;
pub mod routing;
pub mod synthesis;
pub mod target_basis;

use crate::error::{QRustError, Result};
use crate::ir::{Circuit, GateType, Operation};
use pass::{Pass, PassManager};
use std::collections::HashSet;

#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct TranspilerConfig {
    pub decompose_basis: bool,
    pub optimization_level: u8,
    pub backend: Option<crate::backend::Backend>,
    /// Target gate set to translate into after all other passes complete.
    ///
    /// If `None` and a `Backend` with non-empty `basis_gates` is provided,
    /// those are used automatically. If both are absent, basis translation
    /// is skipped (permissive mode — internal canonical gates are kept).
    pub target_basis: Option<HashSet<String>>,
}

impl Default for TranspilerConfig {
    fn default() -> Self {
        Self {
            decompose_basis: true,
            optimization_level: 1,
            backend: None,
            target_basis: None,
        }
    }
}

impl TranspilerConfig {
    pub fn builder() -> TranspilerConfigBuilder {
        TranspilerConfigBuilder::default()
    }
}

#[derive(Clone, Debug, Default)]
pub struct TranspilerConfigBuilder {
    decompose_basis: Option<bool>,
    optimization_level: Option<u8>,
    backend: Option<crate::backend::Backend>,
    target_basis: Option<HashSet<String>>,
}

impl TranspilerConfigBuilder {
    pub fn decompose_basis(mut self, v: bool) -> Self {
        self.decompose_basis = Some(v);
        self
    }
    pub fn optimization_level(mut self, level: u8) -> Self {
        self.optimization_level = Some(level);
        self
    }
    pub fn backend(mut self, backend: crate::backend::Backend) -> Self {
        self.backend = Some(backend);
        self
    }
    /// Explicitly set the target basis gate set.
    /// Overrides any `basis_gates` coming from the `Backend`.
    pub fn target_basis(mut self, basis: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.target_basis = Some(basis.into_iter().map(|s| s.into()).collect());
        self
    }
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
        }
    }
}

/// KAK-based 2-qubit synthesis pass: finds any 2-qubit gate that is NOT CX
/// and decomposes it (via its unitary) into {U, CX}.
#[derive(Debug, Clone, Copy)]
pub struct KakSynthesisPass;

impl Pass for KakSynthesisPass {
    fn name(&self) -> &str {
        "KakSynthesisPass"
    }

    fn run(&self, circuit: &Circuit, _property_set: &mut property_set::PropertySet) -> Circuit {
        use crate::ir::GateDefinition;
        use crate::transpiler::synthesis::kak::KakSynthesizer;
        use crate::transpiler::synthesis::Synthesizer;

        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } if qubits.len() == 2 && !matches!(name, GateType::CX) => {
                    // Get local 4x4 unitary for this gate and re-synthesize.
                    let u = name.unitary(params);
                    if u.nrows() == 4 && u.ncols() == 4 {
                        if let Some(sub) = KakSynthesizer.synthesize(&u, &[]) {
                            // Remap logical qubits 0,1 of the sub-circuit to physical qubits[0], qubits[1].
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
                    // Fallback: leave untouched.
                    out.add_op(op.clone());
                }
                other => out.add_op(other.clone()),
            }
        }
        out
    }
}

/// Native basis translation: rewrites U(θ,φ,λ) into the backend's native set,
/// e.g. IBM {id, rz, sx, x, cx}: U(θ,φ,λ) = Rz(φ) · Sx · Rz(θ) · Sx · Rz(λ).
///
/// Only fires when the backend basis does NOT already contain {u, cx}. CX
/// gates pass through unchanged (they're in virtually every native set).
#[derive(Debug, Clone)]
pub struct NativeBasisTranslationPass {
    pub backend: crate::backend::Backend,
}

impl Pass for NativeBasisTranslationPass {
    fn name(&self) -> &str {
        "NativeBasisTranslationPass"
    }

    fn run(&self, circuit: &Circuit, _property_set: &mut property_set::PropertySet) -> Circuit {
        // If the backend natively supports the abstract U gate, no translation needed.
        let native_u = self.backend.basis_gates.contains("u")
            || self.backend.basis_gates.contains("u3")
            || self.backend.basis_gates.contains("U");
        if native_u || self.backend.basis_gates.is_empty() {
            return circuit.clone();
        }

        // Check for IBM-style basis: {id, rz, sx, x, cx} (or similar).
        let has_rz = self.backend.basis_gates.contains("rz");
        let has_sx = self.backend.basis_gates.contains("sx");
        if !has_rz || !has_sx {
            // Unknown basis — don't attempt translation.
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
                    // U(θ,φ,λ) ≡ Rz(φ) · Sx · Rz(θ - π) · Sx · Rz(λ - π)   (IBM standard)
                    // Using the cleaner textbook form:
                    //   U(θ,φ,λ) = Rz(φ) · Sx · Rz(θ) · Sx · Rz(λ)   (up to global phase)
                    // The exact identity differs by a global phase; both are valid for
                    // hardware execution, but the Qiskit convention uses the shifted form:
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

pub fn transpile(circuit: &Circuit, config: Option<TranspilerConfig>) -> Result<Circuit> {
    let config = config.unwrap_or_default();
    let mut pm = PassManager::new();

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
        }));
        // After routing, orient CX gates along the coupling map (H-sandwich
        // when only the reverse edge is available).
        pm.add_pass(Box::new(decomposition::CxDirectionPass {
            backend: backend.clone(),
        }));
    }

    if config.decompose_basis {
        // Decompose non-basis gates first (fast path for known identities).
        pm.add_pass(Box::new(decomposition::BasisDecompositionPass));
        // Synthesize any residual 2-qubit gates (that are not CX) via KAK.
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

    // Final target-basis translation.
    // Priority: explicit target_basis on config > backend.basis_gates > skip (permissive).
    let resolved_basis: Option<HashSet<String>> = config.target_basis.clone().or_else(|| {
        config
            .backend
            .as_ref()
            .filter(|b| !b.basis_gates.is_empty())
            .map(|b| b.basis_gates.clone())
    });

    if let Some(basis) = resolved_basis {
        // validate_universality is called inside TargetBasisPass::new.
        let tb_pass = target_basis::TargetBasisPass::new(basis)?;
        pm.add_pass(Box::new(tb_pass));
    }
    // If no basis is resolved: permissive mode — leave canonical gates in place.

    Ok(pm.run(circuit))
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
}
