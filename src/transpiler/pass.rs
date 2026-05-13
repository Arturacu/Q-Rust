//! Pass infrastructure: the [`Pass`] trait and [`PassManager`].
//!
//! Loop 5 review §Finding 3 (Major): the pass manager is the simplest
//! possible "run them all in order" loop. Adding analysis-vs-transform
//! distinction, dependency declaration, or full LLVM-style invalidation
//! is out of scope for this iteration (review §"Out-of-Scope"). What we
//! *do* add here is **conditional pass execution** — a closure-guarded
//! variant of `add_pass` so that downstream callers can use the
//! analysis-pass output of [`crate::transpiler::profiler::CircuitProfilerPass`]
//! to bypass expensive transformation passes when the profile says the
//! transformation has nothing to do.
//!
//! References:
//! - Hoult & Robinson 2018 (LLVM Devmtg), *The New Pass Manager* —
//!   establishes analysis-vs-transform separation and lazy invalidation.
//! - Sivarajah et al. 2020, *t|ket⟩: a retargetable compiler for NISQ
//!   devices*, Quantum Sci. Tech. 6 014003 — `CompilationUnit` /
//!   pass-prerequisite design.

use super::property_set::PropertySet;
use crate::ir::Circuit;

/// A transpiler pass: consumes a circuit, produces a (possibly) transformed one.
///
/// Passes share metadata via a [`PropertySet`] that lives on the [`PassManager`].
pub trait Pass {
    /// A short identifier used in diagnostics and debug output.
    fn name(&self) -> &str;

    /// Runs the pass, possibly reading or writing properties.
    fn run(&self, circuit: &Circuit, property_set: &mut PropertySet) -> Circuit;

    /// Returns `true` if this pass only reads the circuit (does not
    /// transform it). Default: `false` (transformation pass).
    ///
    /// Loop 5 §Finding 3: this is currently advisory only — the pass
    /// manager does not yet use it for invalidation. Pure analysis
    /// passes (e.g. `CircuitProfilerPass`) should override to `true`
    /// so future versions of the manager can skip cloning the circuit.
    ///
    /// **Contract**: a pass overriding this to `true` is asserting that
    /// it will return a circuit whose `operations`, `num_qubits`, and
    /// `num_cbits` are observationally equal to its input. The
    /// `test_circuit_profiler_pass_is_observationally_pure` test
    /// (Loop 5 §NI-3) pins this for `CircuitProfilerPass`. Future
    /// `is_analysis() == true` passes should add a similar pin.
    fn is_analysis(&self) -> bool {
        false
    }
}

/// Internal storage for [`PassManager`] entries: either an unconditional
/// pass or one guarded by a predicate over the [`PropertySet`].
enum PassEntry {
    Always(Box<dyn Pass>),
    Conditional {
        pass: Box<dyn Pass>,
        predicate: Box<dyn Fn(&PropertySet) -> bool>,
    },
}

/// Sequentially applies a list of [`Pass`]es, sharing a [`PropertySet`].
pub struct PassManager {
    entries: Vec<PassEntry>,
    /// Shared metadata — exposed so tests and downstream tools can inspect it.
    pub property_set: PropertySet,
}

impl std::fmt::Debug for PassManager {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PassManager")
            .field("num_passes", &self.entries.len())
            .finish()
    }
}

impl Default for PassManager {
    fn default() -> Self {
        Self::new()
    }
}

impl PassManager {
    /// Creates a new empty manager.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            property_set: PropertySet::new(),
        }
    }

    /// Appends a pass that runs unconditionally.
    pub fn add_pass(&mut self, pass: Box<dyn Pass>) {
        self.entries.push(PassEntry::Always(pass));
    }

    /// Appends a pass that runs only when `predicate(&property_set)` returns `true`.
    ///
    /// The predicate is evaluated at the moment the pass would run,
    /// against the current state of the property set — so passes earlier
    /// in the pipeline can populate analysis data the predicate consults.
    ///
    /// # Example
    /// ```ignore
    /// pm.add_pass(Box::new(profiler::CircuitProfilerPass));
    /// pm.add_conditional(
    ///     Box::new(optimization::InverseCancellationPass),
    ///     |ps| {
    ///         ps.get::<profiler::ProfileReport>("profile_report")
    ///             .map(|r| !r.can_bypass_inverse_cancellation())
    ///             .unwrap_or(true)
    ///     },
    /// );
    /// ```
    pub fn add_conditional<F>(&mut self, pass: Box<dyn Pass>, predicate: F)
    where
        F: Fn(&PropertySet) -> bool + 'static,
    {
        self.entries.push(PassEntry::Conditional {
            pass,
            predicate: Box::new(predicate),
        });
    }

    /// Runs all passes in order, evaluating predicates lazily.
    pub fn run(&mut self, circuit: &Circuit) -> Circuit {
        let mut current = circuit.clone();
        for entry in &self.entries {
            match entry {
                PassEntry::Always(p) => {
                    current = p.run(&current, &mut self.property_set);
                }
                PassEntry::Conditional { pass, predicate } => {
                    if predicate(&self.property_set) {
                        current = pass.run(&current, &mut self.property_set);
                    }
                }
            }
        }
        current
    }

    /// Returns the number of registered entries (conditional + unconditional).
    ///
    /// # Example
    /// ```
    /// use q_rust::transpiler::pass::PassManager;
    /// let pm = PassManager::new();
    /// assert_eq!(pm.num_passes(), 0);
    /// ```
    pub fn num_passes(&self) -> usize {
        self.entries.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{GateType, Operation};

    struct MockPass;
    impl Pass for MockPass {
        fn name(&self) -> &str {
            "MockPass"
        }
        fn run(&self, circuit: &Circuit, _props: &mut PropertySet) -> Circuit {
            let mut c = circuit.clone();
            c.add_op(Operation::Gate {
                name: GateType::ID,
                qubits: vec![0],
                params: vec![],
            });
            c
        }
    }

    /// Pass that records that it ran by setting a flag in the property set.
    struct FlagPass {
        key: &'static str,
    }
    impl Pass for FlagPass {
        fn name(&self) -> &str {
            "FlagPass"
        }
        fn run(&self, circuit: &Circuit, props: &mut PropertySet) -> Circuit {
            props.insert(self.key, true);
            circuit.clone()
        }
    }

    #[test]
    fn test_pass_manager() {
        let c = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        pm.add_pass(Box::new(MockPass));
        let out = pm.run(&c);
        assert_eq!(out.operations.len(), 1);
    }

    /// Loop 5 §Finding 3: conditional pass runs when predicate is true.
    #[test]
    fn test_conditional_pass_runs_when_predicate_true() {
        let c = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        pm.add_conditional(Box::new(FlagPass { key: "ran" }), |_ps| true);
        pm.run(&c);
        assert_eq!(pm.property_set.get::<bool>("ran"), Some(&true));
    }

    /// Loop 5 §Finding 3: conditional pass is skipped when predicate is false.
    #[test]
    fn test_conditional_pass_skipped_when_predicate_false() {
        let c = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        pm.add_conditional(Box::new(FlagPass { key: "ran" }), |_ps| false);
        pm.run(&c);
        assert!(pm.property_set.get::<bool>("ran").is_none());
    }

    /// Loop 5 §Finding 3: predicates can read state set by earlier passes.
    #[test]
    fn test_conditional_pass_can_read_earlier_property_state() {
        let c = Circuit::new(1, 0);
        let mut pm = PassManager::new();
        // First pass sets a flag.
        pm.add_pass(Box::new(FlagPass {
            key: "should_run_next",
        }));
        // Second pass runs only if that flag is set.
        pm.add_conditional(Box::new(FlagPass { key: "second_ran" }), |ps| {
            ps.get::<bool>("should_run_next").copied().unwrap_or(false)
        });
        pm.run(&c);
        assert_eq!(pm.property_set.get::<bool>("second_ran"), Some(&true));
    }

    /// Loop 5 §Finding 3: default `is_analysis()` is `false`.
    #[test]
    fn test_pass_default_is_analysis_is_false() {
        assert!(!MockPass.is_analysis());
        assert!(!FlagPass { key: "x" }.is_analysis());
    }

    /// Loop 5 §NI-2: pin `num_passes()` semantics so it's no longer
    /// dead-by-grep. Cheap counter-style accessor; tracks both
    /// unconditional and conditional entries.
    #[test]
    fn test_num_passes_counts_both_kinds() {
        let mut pm = PassManager::new();
        assert_eq!(pm.num_passes(), 0);
        pm.add_pass(Box::new(MockPass));
        assert_eq!(pm.num_passes(), 1);
        pm.add_conditional(Box::new(FlagPass { key: "x" }), |_| true);
        assert_eq!(pm.num_passes(), 2);
    }

    /// Loop 5 §NI-3: any pass that overrides `is_analysis()` to `true`
    /// asserts a contract that it does not mutate the circuit's
    /// observable state (operations / num_qubits / num_cbits). Pin this
    /// for `CircuitProfilerPass` (currently the only `is_analysis==true`
    /// pass in tree) so future additions cannot quietly violate it.
    /// Without this test, the marker is purely advisory and a future
    /// contributor could set `is_analysis() == true` on a transforming
    /// pass — silently corrupting the optimization in AF-3 once it
    /// becomes load-bearing.
    #[test]
    fn test_circuit_profiler_pass_is_observationally_pure() {
        use crate::transpiler::profiler::CircuitProfilerPass;
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
        // Contract assertion #1: the marker is set.
        assert!(CircuitProfilerPass.is_analysis());
        // Contract assertion #2: the run is observationally pure.
        let mut ps = PropertySet::new();
        let pre = c.clone();
        let post = CircuitProfilerPass.run(&c, &mut ps);
        assert_eq!(pre.operations, post.operations);
        assert_eq!(pre.num_qubits, post.num_qubits);
        assert_eq!(pre.num_cbits, post.num_cbits);
        // And the report was deposited.
        assert!(ps
            .get::<crate::transpiler::profiler::ProfileReport>("profile_report")
            .is_some());
    }
}
