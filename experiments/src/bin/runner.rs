//! Q-Rust benchmark runner — one configuration per invocation.
//!
//! This binary transpiles a single fixture under a single (topology, opt-level,
//! ablation) configuration and prints **one JSON object** describing every
//! metric the evaluation needs.  The Python orchestrator (`bench_runner.py`) sweeps
//! the matrix by invoking this binary once per cell, which also gives us a true
//! per-circuit peak-RSS reading (one process == one transpilation).
//!
//! ## Why a bespoke runner instead of the CLI?
//! The `qrust` CLI emits no timings and only a human-readable fidelity string.
//! Here we drive the transpiler passes *directly*, faithfully reproducing
//! `transpile()`'s two-`PassManager` structure (one `PropertySet` for the
//! optimization stage, one shared across all lowering passes — see
//! `src/transpiler/mod.rs::build_pass_manager_for`).  Driving the passes
//! ourselves buys three things the public API cannot give:
//!   1. a true 5-stage µs breakdown (Parse → LogicOpt → Routing → BasisDecomp → Cleanup),
//!   2. the `initial_layout`/`final_layout` permutation from the property set,
//!      enabling *layout-aware* verification of routed circuits, and
//!   3. a clean Stage-5 (cleanup) bypass toggle for the cleanup-ablation study.
//!
//! ## Faithfulness contract
//! With cleanup enabled, the pass list and ordering here are identical to
//! `build_pass_manager_for`, so the produced circuit matches `transpile()`
//! exactly.  The only intentional deviation is `target_basis = {u, cx}` forced
//! on every run: Q-Rust's rz/sx translator is unwired and it has no
//! Solovay–Kitaev, so {U,CX} is the only basis it can actually emit.

use std::time::Instant;

use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::simulator::{
    equivalence_by_sampling, equivalence_by_sampling_with_layout, extract_logical_unitary,
    try_circuit_to_unitary, unitary_fidelity,
};
use q_rust::transpiler::decomposition::{
    decompose_basis, BasisDecompositionPass, CxDirectionPass, HighArityDecompositionPass,
};
use q_rust::transpiler::layout::SabreLayoutPass;
use q_rust::transpiler::optimization::{
    CommutationCancellationPass, CrossConjugationPass, GateCrystallizationPass, GateFusionPass,
    InverseCancellationPass, ParameterSimplificationPass, RotationMergePass,
    SwapSimplificationPass, TrailingSwapElisionPass,
};
use q_rust::transpiler::pass::Pass;
use q_rust::transpiler::property_set::PropertySet;
use q_rust::transpiler::routing::{BeamSabrePass, LookaheadStrategy};
use q_rust::transpiler::target_basis::TargetBasisPass;
use q_rust::transpiler::KakSynthesisPass;

use serde_json::{json, Value};

/// Full-unitary materialization cap.  At n qubits the unitary is a 2^n × 2^n
/// dense complex matrix and `circuit_to_unitary` is O(g · 8^n), so exact
/// verification is only tractable for small circuits.  Above this we fall back
/// to state-vector sampling (identity layout) or declare the run unverifiable
/// (non-trivial routing layout, which sampling cannot account for).
const EXACT_MAX_QUBITS: usize = 8;

/// Largest width the state-vector sampler can handle (mirrors the library's
/// `MAX_STATE_VECTOR_QUBITS - 2` safety margin used by `verify_equivalence`).
const SAMPLE_MAX_QUBITS: usize = 22;

/// Practical cap for the *layout-aware* state-vector sampler. It evolves the
/// (wider, routed) physical circuit, so it is markedly more expensive than the
/// plain sampler; verifying a routed 22-qubit circuit costs minutes. Capping at
/// 20 keeps re-runs tractable while still certifying the grid (16q) and every
/// circuit up to 18 qubits; wider routed circuits are reported Unverifiable.
const LAYOUT_SAMPLE_MAX_QUBITS: usize = 20;

/// Fixed seed so the statistical verdict is reproducible across runs.
const VERIFY_SEED: u64 = 0x00C0_FFEE_DEAD_BEEF;
const VERIFY_SAMPLES: usize = 8;

// --------------------------------------------------------------------------- //
// Peak RSS via getrusage                                                       //
// --------------------------------------------------------------------------- //
fn peak_rss_bytes() -> u64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if rc != 0 {
        return 0;
    }
    let usage = unsafe { usage.assume_init() };
    let maxrss = usage.ru_maxrss as u64;
    // macOS reports ru_maxrss in bytes; Linux/BSD in kibibytes.
    if cfg!(target_os = "macos") {
        maxrss
    } else {
        maxrss * 1024
    }
}

// --------------------------------------------------------------------------- //
// Circuit metrics                                                              //
// --------------------------------------------------------------------------- //
struct Counts {
    gates: usize,
    cx: usize,
    two_q: usize,
    swaps: usize,
    depth: usize,
}

fn count(circuit: &Circuit) -> Counts {
    let mut cx = 0;
    let mut two_q = 0;
    let mut swaps = 0;
    for op in &circuit.operations {
        if let Operation::Gate { name, qubits, .. } = op {
            if qubits.len() >= 2 {
                two_q += 1;
                match name {
                    GateType::CX => cx += 1,
                    GateType::SWAP => swaps += 1,
                    _ => {}
                }
            }
        }
    }
    Counts {
        gates: circuit.operations.len(),
        cx,
        two_q,
        swaps,
        depth: circuit.depth(),
    }
}

fn pct_reduction(before: usize, after: usize) -> Value {
    if before == 0 {
        Value::Null
    } else {
        json!(((before as f64 - after as f64) / before as f64) * 100.0)
    }
}

// --------------------------------------------------------------------------- //
// Backend resolution                                                          //
// --------------------------------------------------------------------------- //
fn resolve_backend(topology: &str, n_logical: usize) -> Option<Backend> {
    match topology {
        "ibm_quito" => Some(Backend::ibm_quito()),
        "ibm_nairobi" => Some(Backend::ibm_nairobi()),
        "all_to_all" => Some(Backend::all_to_all(n_logical.max(1))),
        "grid_4x4" => Some(Backend::grid(4, 4)),
        "linear" => Some(Backend::linear(n_logical.max(1))),
        "ring" => Some(Backend::ring(n_logical.max(2))),
        "star" => Some(Backend::star(n_logical.max(1))),
        "none" => None,
        other => {
            eprintln!("unknown topology: {other}");
            std::process::exit(2);
        }
    }
}

// --------------------------------------------------------------------------- //
// Stage runners — each mirrors a slice of build_pass_manager_for and threads   //
// the supplied PropertySet (opt stage gets its own; all lowering passes share  //
// one, exactly as transpile() does).                                           //
// --------------------------------------------------------------------------- //
fn run_pass(p: &dyn Pass, c: &Circuit, ps: &mut PropertySet) -> Circuit {
    p.run(c, ps)
}

/// Stage 2: logic optimization (backend-independent).
fn stage_logic_opt(circuit: &Circuit, opt: u8, ps: &mut PropertySet) -> Circuit {
    let mut c = circuit.clone();
    if opt >= 1 {
        c = run_pass(&GateCrystallizationPass { epsilon: 1e-9 }, &c, ps);
        c = run_pass(&ParameterSimplificationPass::default(), &c, ps);
    }
    if opt >= 2 {
        c = run_pass(&RotationMergePass, &c, ps);
        c = run_pass(&CrossConjugationPass, &c, ps);
        c = run_pass(&InverseCancellationPass, &c, ps);
        c = run_pass(&CommutationCancellationPass, &c, ps);
        c = run_pass(&ParameterSimplificationPass::default(), &c, ps);
    }
    c
}

/// Stage 3: layout + routing + CX-direction (only when a backend is present).
fn stage_routing(
    circuit: &Circuit,
    opt: u8,
    backend: &Option<Backend>,
    ps: &mut PropertySet,
) -> Circuit {
    let mut c = circuit.clone();
    if let Some(backend) = backend {
        // Mirror the library: reduce >2-qubit gates before routing.
        c = run_pass(&HighArityDecompositionPass, &c, ps);
        if opt >= 2 {
            let (num_trials, num_iterations) = if opt == 2 { (10, 3) } else { (50, 5) };
            c = run_pass(
                &SabreLayoutPass {
                    backend: backend.clone(),
                    num_trials,
                    num_iterations,
                },
                &c,
                ps,
            );
        }
        let (beam_width, branch_factor, bidir_iterations) = match opt {
            0 | 1 => (1, 1, 1),
            2 => (4, 3, 2),
            _ => (8, 5, 4),
        };
        c = run_pass(
            &BeamSabrePass {
                backend: backend.clone(),
                beam_width,
                branch_factor,
                bidir_iterations,
                lookahead_strategy: LookaheadStrategy::default(),
            },
            &c,
            ps,
        );
        c = run_pass(
            &CxDirectionPass {
                backend: backend.clone(),
            },
            &c,
            ps,
        );
        // Fold trailing SWAPs (e.g. the QFT bit-reversal network) into the
        // output layout (mirrors the library pipeline).
        c = run_pass(&TrailingSwapElisionPass, &c, ps);
    }
    c
}

/// Stage 4: target-basis translation + basis decomposition + KAK synthesis.
/// The benchmark fixes `{u, cx}` as the common cross-tool comparison basis.
fn stage_basis_decomp(circuit: &Circuit, ps: &mut PropertySet) -> Result<Circuit, String> {
    let basis: std::collections::HashSet<String> =
        ["u", "cx"].iter().map(|s| s.to_string()).collect();
    let tb = TargetBasisPass::new(basis).map_err(|e| format!("target basis: {e}"))?;
    let mut c = run_pass(&tb, circuit, ps);
    c = run_pass(&BasisDecompositionPass, &c, ps);
    c = run_pass(&KakSynthesisPass, &c, ps);
    Ok(c)
}

/// Stage 5: post-routing cleanup (level-3 only). Bypassed entirely when
/// `enabled == false` — this is the cleanup-ablation toggle.
fn stage_cleanup(circuit: &Circuit, opt: u8, enabled: bool, ps: &mut PropertySet) -> Circuit {
    let mut c = circuit.clone();
    if opt >= 3 && enabled {
        c = run_pass(&GateFusionPass, &c, ps);
        c = run_pass(&SwapSimplificationPass, &c, ps);
        c = run_pass(&InverseCancellationPass, &c, ps);
        c = run_pass(&ParameterSimplificationPass::default(), &c, ps);
    }
    c
}

// --------------------------------------------------------------------------- //
// Verification (layout-aware)                                                  //
// --------------------------------------------------------------------------- //
struct VerifyResult {
    class: &'static str, // "Exact" | "Statistical" | "Unverifiable" | "NotEquivalent"
    fidelity: Option<f64>,
    verified: bool,
    method: String,
}

fn is_identity(layout: &Option<Vec<usize>>, n: usize) -> bool {
    match layout {
        None => true,
        Some(v) => v.len() == n && v.iter().enumerate().all(|(i, &x)| i == x),
    }
}

fn verify(
    orig: &Circuit,
    transpiled: &Circuit,
    initial_layout: &Option<Vec<usize>>,
    final_layout: &Option<Vec<usize>>,
) -> VerifyResult {
    let n_log = orig.num_qubits;
    let n_phys = transpiled.num_qubits;
    let identity =
        n_phys == n_log && is_identity(initial_layout, n_log) && is_identity(final_layout, n_log);

    if identity {
        // No routing permutation: a direct comparison is valid.
        if n_log <= EXACT_MAX_QUBITS {
            match (
                try_circuit_to_unitary(orig),
                try_circuit_to_unitary(transpiled),
            ) {
                (Ok(u1), Ok(u2)) => {
                    let fid = unitary_fidelity(&u1, &u2);
                    let ok = (fid - 1.0).abs() < 1e-9;
                    VerifyResult {
                        class: if ok { "Exact" } else { "NotEquivalent" },
                        fidelity: Some(fid),
                        verified: ok,
                        method: "exact-unitary".into(),
                    }
                }
                _ => VerifyResult {
                    class: "Unverifiable",
                    fidelity: None,
                    verified: false,
                    method: "exact-unitary-failed".into(),
                },
            }
        } else if n_log <= SAMPLE_MAX_QUBITS {
            // Decompose the source so the sampler never sees a >2-qubit gate it
            // cannot embed at this width (e.g. the 16-qubit adder's Toffolis);
            // use fewer Haar samples for wide circuits (each is 2^{-n/2}-tight).
            let orig2 = decompose_basis(orig);
            let samples = if n_log >= 18 { 2 } else { VERIFY_SAMPLES };
            match equivalence_by_sampling(&orig2, transpiled, samples, VERIFY_SEED) {
                Ok(min_fid) => {
                    let ok = (min_fid - 1.0).abs() < 1e-6;
                    VerifyResult {
                        class: if ok { "Statistical" } else { "NotEquivalent" },
                        fidelity: Some(min_fid),
                        verified: ok,
                        method: format!("statistical-{samples}-haar"),
                    }
                }
                Err(e) => VerifyResult {
                    class: "Unverifiable",
                    fidelity: None,
                    verified: false,
                    method: format!("sampling-error: {e}"),
                },
            }
        } else {
            VerifyResult {
                class: "Unverifiable",
                fidelity: None,
                verified: false,
                method: format!("{n_log} qubits exceeds sampling limit {SAMPLE_MAX_QUBITS}"),
            }
        }
    } else {
        // Routing applied a layout permutation and/or widened the register.
        // The only layout-aware primitive we have needs the full routed
        // unitary, so this is feasible only for small physical widths.
        if n_phys <= EXACT_MAX_QUBITS {
            let init = initial_layout
                .clone()
                .unwrap_or_else(|| (0..n_log).collect());
            let fin = final_layout
                .clone()
                .unwrap_or_else(|| (0..n_phys).collect());
            match (
                try_circuit_to_unitary(orig),
                try_circuit_to_unitary(transpiled),
            ) {
                (Ok(u_orig), Ok(u_routed)) => {
                    let u_log = extract_logical_unitary(&u_routed, n_log, &init, &fin);
                    let fid = unitary_fidelity(&u_orig, &u_log);
                    let ok = (fid - 1.0).abs() < 1e-6;
                    VerifyResult {
                        class: if ok { "Exact" } else { "NotEquivalent" },
                        fidelity: Some(fid),
                        verified: ok,
                        method: "exact-layout-aware".into(),
                    }
                }
                _ => VerifyResult {
                    class: "Unverifiable",
                    fidelity: None,
                    verified: false,
                    method: "layout-aware-unitary-failed".into(),
                },
            }
        } else if n_phys <= LAYOUT_SAMPLE_MAX_QUBITS {
            // Layout-aware state-vector sampling: O(2^n_phys), certifies routed
            // circuits well beyond the exact boundary. Fewer Haar samples for
            // wide circuits (each is already 2^{-n_phys/2} concentrated).
            let init = initial_layout
                .clone()
                .unwrap_or_else(|| (0..n_log).collect());
            let fin = final_layout
                .clone()
                .unwrap_or_else(|| (0..n_phys).collect());
            let samples = if n_phys >= 18 { 2 } else { VERIFY_SAMPLES };
            match equivalence_by_sampling_with_layout(
                orig,
                transpiled,
                &init,
                &fin,
                samples,
                VERIFY_SEED,
            ) {
                Ok(min_fid) => {
                    let ok = (min_fid - 1.0).abs() < 1e-6;
                    VerifyResult {
                        class: if ok { "Statistical" } else { "NotEquivalent" },
                        fidelity: Some(min_fid),
                        verified: ok,
                        method: format!("statistical-layout-aware-{samples}-haar"),
                    }
                }
                Err(e) => VerifyResult {
                    class: "Unverifiable",
                    fidelity: None,
                    verified: false,
                    method: format!("layout-aware-sampling-error: {e}"),
                },
            }
        } else {
            VerifyResult {
                class: "Unverifiable",
                fidelity: None,
                verified: false,
                method: format!(
                    "routed width {n_phys} exceeds layout-aware sampling budget \
                     {LAYOUT_SAMPLE_MAX_QUBITS}"
                ),
            }
        }
    }
}

// --------------------------------------------------------------------------- //
// Argument parsing                                                            //
// --------------------------------------------------------------------------- //
struct Args {
    input: String,
    fixture: String,
    family: String,
    nominal_n: i64,
    topology: String,
    opt: u8,
    cleanup: bool,   // false => Stage-5 bypassed (ablation)
    ablation: bool,  // marks this as an ablation run for tagging
    emit_qasm: bool, // also emit the routed output QASM + layout (for external verification)
    no_verify: bool, // skip equivalence verification (e.g. for CX/depth-only ablation sweeps)
}

fn parse_args() -> Args {
    let mut input = None;
    let mut fixture = String::new();
    let mut family = String::new();
    let mut nominal_n = 0i64;
    let mut topology = String::from("all_to_all");
    let mut opt = 1u8;
    let mut cleanup = true;
    let mut ablation = false;
    let mut emit_qasm = false;
    let mut no_verify = false;

    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--input" => {
                i += 1;
                input = Some(argv[i].clone());
            }
            "--fixture" => {
                i += 1;
                fixture = argv[i].clone();
            }
            "--family" => {
                i += 1;
                family = argv[i].clone();
            }
            "--nominal-n" => {
                i += 1;
                nominal_n = argv[i].parse().unwrap_or(0);
            }
            "--topology" => {
                i += 1;
                topology = argv[i].clone();
            }
            "--opt" => {
                i += 1;
                opt = argv[i].parse().unwrap_or(1);
            }
            "--no-cleanup" => {
                cleanup = false;
            }
            "--ablation" => {
                ablation = true;
            }
            "--emit-qasm" => {
                emit_qasm = true;
            }
            "--no-verify" => {
                no_verify = true;
            }
            other => {
                eprintln!("unknown arg: {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    Args {
        input: input.unwrap_or_else(|| {
            eprintln!("--input required");
            std::process::exit(2);
        }),
        fixture,
        family,
        nominal_n,
        topology,
        opt: opt.min(3),
        cleanup,
        ablation,
        emit_qasm,
        no_verify,
    }
}

// --------------------------------------------------------------------------- //
// Main                                                                         //
// --------------------------------------------------------------------------- //
fn main() {
    let args = parse_args();

    let mut out = json!({
        "fixture": args.fixture,
        "family": args.family,
        "nominal_n": args.nominal_n,
        "topology": args.topology,
        "opt_level": args.opt,
        "basis": ["u", "cx"],
        "cleanup_enabled": args.cleanup,
        "ablation": args.ablation,
    });

    // --- Stage 1: parse (timed) ---
    let src = match std::fs::read_to_string(&args.input) {
        Ok(s) => s,
        Err(e) => {
            emit_error(&mut out, &format!("read {}: {e}", args.input));
            return;
        }
    };
    let t_parse = Instant::now();
    let parsed = match parse_qasm(&src) {
        Ok(c) => c,
        Err(e) => {
            emit_error(&mut out, &format!("parse: {e}"));
            return;
        }
    };
    let us_parse = t_parse.elapsed().as_micros() as u64;

    let n_logical = parsed.num_qubits;
    out["num_qubits_logical"] = json!(n_logical);
    let backend = resolve_backend(&args.topology, n_logical);
    let backend_qubits = backend.as_ref().map(|b| b.num_qubits);
    out["backend_qubits"] = json!(backend_qubits);

    // Width guard: a fixed-size backend cannot host a wider circuit.
    if let Some(bq) = backend_qubits {
        if n_logical > bq {
            out["status"] = json!("skipped_too_wide");
            out["reason"] = json!(format!("circuit {n_logical}q > backend {bq}q"));
            println!("{out}");
            return;
        }
    }

    let pre = count(&parsed);

    // --- Stages 2–5 (timed); two PropertySets, exactly as transpile() ---
    let mut opt_ps = PropertySet::new();
    let t_opt = Instant::now();
    let after_opt = stage_logic_opt(&parsed, args.opt, &mut opt_ps);
    let us_opt = t_opt.elapsed().as_micros() as u64;

    let mut low_ps = PropertySet::new();

    let t_route = Instant::now();
    let after_route = stage_routing(&after_opt, args.opt, &backend, &mut low_ps);
    let us_route = t_route.elapsed().as_micros() as u64;

    let t_decomp = Instant::now();
    let after_decomp = match stage_basis_decomp(&after_route, &mut low_ps) {
        Ok(c) => c,
        Err(e) => {
            emit_error(&mut out, &e);
            return;
        }
    };
    let us_decomp = t_decomp.elapsed().as_micros() as u64;

    // Controlled ablation: metrics of the *same* routed+decomposed circuit
    // before cleanup, so `pre_cleanup` vs `post` isolates the cleanup effect
    // from routing randomness (the two share one routing).
    let pre_cleanup = count(&after_decomp);

    let t_clean = Instant::now();
    let final_circuit = stage_cleanup(&after_decomp, args.opt, args.cleanup, &mut low_ps);
    let us_clean = t_clean.elapsed().as_micros() as u64;

    // Peak RSS captured *before* verification (which materializes unitaries and
    // would otherwise dominate the reading). This is the transpilation peak.
    let peak_rss = peak_rss_bytes();

    let post = count(&final_circuit);
    let us_total = us_parse + us_opt + us_route + us_decomp + us_clean;

    // --- Verification (layout-aware) ---
    let initial_layout = low_ps.get::<Vec<usize>>("initial_layout").cloned();
    let final_layout = low_ps.get::<Vec<usize>>("final_layout").cloned();
    let v = if args.no_verify {
        VerifyResult {
            class: "Skipped",
            fidelity: None,
            verified: false,
            method: "skipped".into(),
        }
    } else {
        verify(&parsed, &final_circuit, &initial_layout, &final_layout)
    };

    out["status"] = json!("ok");
    out["pre"] = json!({
        "gates": pre.gates, "cx": pre.cx, "two_q": pre.two_q, "depth": pre.depth,
    });
    out["post"] = json!({
        "gates": post.gates, "cx": post.cx, "two_q": post.two_q,
        "swaps": post.swaps, "depth": post.depth, "num_qubits": final_circuit.num_qubits,
    });
    out["pre_cleanup"] = json!({
        "gates": pre_cleanup.gates, "cx": pre_cleanup.cx, "depth": pre_cleanup.depth,
    });
    out["cx_reduction_pct"] = pct_reduction(pre.cx, post.cx);
    out["two_q_reduction_pct"] = pct_reduction(pre.two_q, post.two_q);
    out["depth_reduction_pct"] = pct_reduction(pre.depth, post.depth);
    out["timing_us"] = json!({
        "parse": us_parse, "logic_opt": us_opt, "routing": us_route,
        "basis_decomp": us_decomp, "cleanup": us_clean, "total": us_total,
    });
    out["peak_rss_bytes"] = json!(peak_rss);
    out["verification"] = json!({
        "verdict_class": v.class,
        "fidelity": v.fidelity,
        "verified": v.verified,
        "method": v.method,
    });

    // Expose the routing permutation so an external tool can perform
    // layout-aware verification of circuits beyond this harness's simulator
    // reach (e.g. the 16-qubit adders). Identity when no routing occurred.
    out["initial_layout"] = json!(initial_layout);
    out["final_layout"] = json!(final_layout);
    if args.emit_qasm {
        out["out_qasm"] = json!(final_circuit.to_qasm(None));
    }

    println!("{out}");
}

fn emit_error(out: &mut Value, msg: &str) {
    out["status"] = json!("error");
    out["error"] = json!(msg);
    println!("{out}");
}
