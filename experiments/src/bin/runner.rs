//! Q-Rust benchmark runner — one configuration per invocation.
//!
//! This binary transpiles a single fixture under a single (topology, opt-level,
//! ablation) configuration and prints **one JSON object** describing every
//! metric the evaluation needs.  The Python orchestrator (`bench_runner.py`) sweeps
//! the matrix by invoking this binary once per cell, which also gives us a true
//! per-circuit peak-RSS reading (one process == one transpilation).
//!
//! Uses the canonical public compiler pipeline, including its final legality
//! checks. Timing covers parsing and complete transpilation separately. Raw
//! output QASM and layouts can be retained for independent validation.

use std::time::Instant;

use q_rust::backend::Backend;
use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile_with_report, TranspilerConfig};
use q_rust::verify::{verify_equivalence_with_layout_configured, Verdict, VerificationConfig};

use serde_json::{json, Value};

/// Dense checks are numerical operator comparisons, not formal proofs.
const EXACT_MAX_QUBITS: usize = 8;
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

struct VerifyResult {
    class: &'static str,
    fidelity: Option<f64>,
    verified: bool,
    method: String,
}

fn verify(
    orig: &Circuit,
    output: &Circuit,
    initial: &Option<Vec<usize>>,
    final_layout: &Option<Vec<usize>>,
) -> VerifyResult {
    let initial = initial
        .clone()
        .unwrap_or_else(|| (0..orig.num_qubits).collect());
    let final_layout = final_layout
        .clone()
        .unwrap_or_else(|| (0..output.num_qubits).collect());
    let config = VerificationConfig {
        max_exact_qubits: EXACT_MAX_QUBITS,
        samples: if output.num_qubits >= 18 {
            2
        } else {
            VERIFY_SAMPLES
        },
        seed: VERIFY_SEED,
        ..VerificationConfig::default()
    };
    match verify_equivalence_with_layout_configured(orig, output, &initial, &final_layout, &config)
    {
        Ok(Verdict::ExactlyEquivalent { fidelity }) => VerifyResult {
            class: "NumericalOperator",
            fidelity: Some(fidelity),
            verified: true,
            method: "full floating-point operator comparison, including layout and clean ancillas"
                .into(),
        },
        Ok(Verdict::StatisticallyEquivalent {
            min_fidelity,
            samples,
        }) => VerifyResult {
            class: "Sampled",
            fidelity: Some(min_fidelity),
            verified: true,
            method: format!("{samples} seeded Haar inputs; numerical evidence, not proof"),
        },
        Ok(Verdict::NotEquivalent { fidelity, method }) => VerifyResult {
            class: "NotEquivalent",
            fidelity: Some(fidelity),
            verified: false,
            method: method.into(),
        },
        Ok(Verdict::Unverifiable { reason }) => VerifyResult {
            class: "Unverifiable",
            fidelity: None,
            verified: false,
            method: reason,
        },
        Err(error) => VerifyResult {
            class: "Unverifiable",
            fidelity: None,
            verified: false,
            method: error.to_string(),
        },
    }
}

fn verification_status(v: &VerifyResult) -> &'static str {
    if v.verified {
        "ok"
    } else if v.class == "NotEquivalent" {
        "incorrect"
    } else if v.class == "Skipped" {
        "unchecked"
    } else {
        "unverifiable"
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

    let mut builder = TranspilerConfig::builder()
        .optimization_level(args.opt)
        .post_routing_cleanup(args.cleanup)
        .target_basis(["u", "cx"]);
    if let Some(backend) = backend {
        builder = builder.backend(backend);
    }
    let t_transpile = Instant::now();
    let (final_circuit, report) = match transpile_with_report(&parsed, Some(builder.build())) {
        Ok(result) => result,
        Err(error) => {
            emit_error(&mut out, &format!("transpile: {error}"));
            return;
        }
    };
    let us_transpile = t_transpile.elapsed().as_micros() as u64;
    // Capture process peak before verification allocates simulator memory.
    let peak_rss = peak_rss_bytes();
    let post = count(&final_circuit);
    let us_total = us_parse + us_transpile;

    // --- Verification (layout-aware) ---
    let initial_layout = report.initial_layout;
    let final_layout = report.final_layout;
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

    out["schema_version"] = json!(2);
    out["compiled"] = json!(true);
    out["status"] = json!(verification_status(&v));
    out["pre"] = json!({
        "gates": pre.gates, "cx": pre.cx, "two_q": pre.two_q, "depth": pre.depth,
    });
    out["post"] = json!({
        "gates": post.gates, "cx": post.cx, "two_q": post.two_q,
        "swaps": post.swaps, "depth": post.depth, "num_qubits": final_circuit.num_qubits,
    });
    out["cx_reduction_pct"] = pct_reduction(pre.cx, post.cx);
    out["two_q_reduction_pct"] = pct_reduction(pre.two_q, post.two_q);
    out["depth_reduction_pct"] = pct_reduction(pre.depth, post.depth);
    out["timing_us"] = json!({
        "parse": us_parse, "transpile": us_transpile, "total": us_total,
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nonpassing_verification_never_reports_success() {
        for (class, expected) in [
            ("NotEquivalent", "incorrect"),
            ("Skipped", "unchecked"),
            ("Unverifiable", "unverifiable"),
        ] {
            let result = VerifyResult {
                class,
                fidelity: None,
                verified: false,
                method: String::new(),
            };
            assert_eq!(verification_status(&result), expected);
        }
    }
}
