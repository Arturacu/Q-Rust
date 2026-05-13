//! Q-Rust CLI: read a `.qasm` file, transpile it, emit the result.
//!
//! # Usage
//! ```text
//! qrust <input.qasm> [--output <out.qasm>] [--backend <name|path>]
//!       [--opt <0|1|2|3>] [--basis <gate,gate,...>] [--report]
//!       [--verify] [--no-decompose]
//! ```
//!
//! # Examples
//! ```text
//! qrust circuit.qasm
//! qrust circuit.qasm --backend ibm_quito --opt 3 --report
//! qrust circuit.qasm --backend my_device.json --basis rz,sx,cx --output out.qasm
//! qrust circuit.qasm --opt 3 --verify
//! ```
//!
//! Built-in backend names: `linear-N`, `grid-RxC`, `ring-N`, `star-N`,
//! `tree-N`, `all2all-N`, `ibm_quito`, `ibm_nairobi`.
//! Anything else is treated as a path to a JSON config.

use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile, transpile_with_report, TranspilerConfig};
use q_rust::verify::verify_equivalence;
use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "--help" || a == "-h") {
        print_help();
        return ExitCode::SUCCESS;
    }

    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn print_help() {
    println!(
        "Q-Rust quantum transpiler\n\
         \n\
         USAGE:\n\
         \x20   qrust <input.qasm> [OPTIONS]\n\
         \n\
         OPTIONS:\n\
         \x20   --output <path>      Write transpiled circuit (default: stdout)\n\
         \x20   --backend <spec>     Built-in name or JSON path (see below)\n\
         \x20   --opt <0|1|2|3>      Optimization level (default: 1)\n\
         \x20   --basis <g,g,...>    Target basis (e.g. rz,sx,cx)\n\
         \x20   --no-decompose       Disable basis decomposition\n\
         \x20   --report             Print per-stage report to stderr\n\
         \x20   --verify             Verify equivalence (≤22 qubits)\n\
         \x20   --help, -h           Show this help\n\
         \n\
         BACKEND SPECIFIERS:\n\
         \x20   linear-N             Linear chain of N qubits\n\
         \x20   grid-RxC             R×C grid\n\
         \x20   ring-N, star-N       Ring / star topology\n\
         \x20   tree-N, all2all-N    Tree / all-to-all topology\n\
         \x20   ibm_quito            IBM Quito (5q heavy-hex)\n\
         \x20   ibm_nairobi          IBM Nairobi (7q heavy-hex)\n\
         \x20   <path/to.json>       Custom JSON backend config\n"
    );
}

fn run(args: &[String]) -> Result<(), String> {
    let mut input: Option<&str> = None;
    let mut output: Option<&str> = None;
    let mut backend_spec: Option<&str> = None;
    let mut opt_level: u8 = 1;
    let mut basis: Option<Vec<String>> = None;
    let mut decompose = true;
    let mut report = false;
    let mut verify = false;

    let mut i = 0;
    while i < args.len() {
        let arg = args[i].as_str();
        match arg {
            "--output" | "-o" => {
                i += 1;
                output = Some(args.get(i).ok_or("--output requires a value")?);
            }
            "--backend" | "-b" => {
                i += 1;
                backend_spec = Some(args.get(i).ok_or("--backend requires a value")?);
            }
            "--opt" => {
                i += 1;
                let v = args.get(i).ok_or("--opt requires a value")?;
                opt_level = v.parse().map_err(|_| format!("invalid opt level: {v}"))?;
            }
            "--basis" => {
                i += 1;
                let v = args.get(i).ok_or("--basis requires a value")?;
                basis = Some(v.split(',').map(|s| s.trim().to_string()).collect());
            }
            "--no-decompose" => decompose = false,
            "--report" => report = true,
            "--verify" => verify = true,
            other if !other.starts_with("--") && input.is_none() => {
                input = Some(other);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }

    let input = input.ok_or("missing <input.qasm>")?;
    let qasm = std::fs::read_to_string(input)
        .map_err(|e| format!("cannot read {input}: {e}"))?;

    let circuit = parse_qasm(&qasm).map_err(|e| format!("parse error: {e}"))?;
    if report {
        eprintln!("=== Q-Rust Transpilation Report ===");
        eprintln!(
            "Input:  {} qubits, {} cbits, {} ops, depth {}",
            circuit.num_qubits,
            circuit.num_cbits,
            circuit.operations.len(),
            circuit.depth()
        );
    }

    let mut builder = TranspilerConfig::builder()
        .optimization_level(opt_level)
        .decompose_basis(decompose);

    if let Some(spec) = backend_spec {
        let backend = parse_backend_spec(spec)?;
        if report {
            eprintln!(
                "Backend: {} ({} qubits, {} edges)",
                backend.name,
                backend.num_qubits,
                backend.coupling_map.edge_count()
            );
        }
        builder = builder.backend(backend);
    }
    if let Some(b) = basis {
        builder = builder.target_basis(b);
    }
    let cfg = builder.build();

    let (out_circ, stage_report) = if report {
        let (c, r) = transpile_with_report(&circuit, Some(cfg))
            .map_err(|e| format!("transpile failed: {e}"))?;
        (c, Some(r))
    } else {
        let c = transpile(&circuit, Some(cfg)).map_err(|e| format!("transpile failed: {e}"))?;
        (c, None)
    };

    if let Some(r) = stage_report {
        eprintln!("\n--- Stage outputs ---");
        for line in r.format_lines() {
            eprintln!("{line}");
        }
        eprintln!(
            "\nFinal:  {} qubits, {} ops, depth {}",
            out_circ.num_qubits,
            out_circ.operations.len(),
            out_circ.depth()
        );
    }

    if verify {
        match verify_equivalence(&circuit, &out_circ) {
            Ok(verdict) => {
                if report {
                    eprintln!("\nEquivalence: {}", verdict.describe());
                }
                if !verdict.is_equivalent() {
                    return Err(format!("verification failed: {}", verdict.describe()));
                }
            }
            Err(e) => return Err(format!("verification error: {e}")),
        }
    }

    let qasm_out = out_circ.to_qasm(None);
    match output {
        Some(path) => std::fs::write(path, &qasm_out)
            .map_err(|e| format!("cannot write {path}: {e}"))?,
        None => print!("{qasm_out}"),
    }
    Ok(())
}

fn parse_backend_spec(spec: &str) -> Result<Backend, String> {
    if let Some(n) = spec.strip_prefix("linear-") {
        return Ok(Backend::linear(parse_usize(n)?));
    }
    if let Some(n) = spec.strip_prefix("ring-") {
        return Ok(Backend::ring(parse_usize(n)?));
    }
    if let Some(n) = spec.strip_prefix("star-") {
        return Ok(Backend::star(parse_usize(n)?));
    }
    if let Some(n) = spec.strip_prefix("tree-") {
        return Ok(Backend::tree(parse_usize(n)?));
    }
    if let Some(n) = spec.strip_prefix("all2all-") {
        return Ok(Backend::all_to_all(parse_usize(n)?));
    }
    if let Some(rest) = spec.strip_prefix("grid-") {
        let parts: Vec<&str> = rest.split('x').collect();
        if parts.len() != 2 {
            return Err(format!("grid spec must be RxC, got: {spec}"));
        }
        return Ok(Backend::grid(parse_usize(parts[0])?, parse_usize(parts[1])?));
    }
    if spec == "ibm_quito" {
        return Ok(Backend::ibm_quito());
    }
    if spec == "ibm_nairobi" {
        return Ok(Backend::ibm_nairobi());
    }
    if Path::new(spec).exists() {
        return Backend::from_json_file(spec).map_err(|e| format!("backend load: {e}"));
    }
    Err(format!("unknown backend: {spec}"))
}

fn parse_usize(s: &str) -> Result<usize, String> {
    s.parse().map_err(|_| format!("invalid number: {s}"))
}