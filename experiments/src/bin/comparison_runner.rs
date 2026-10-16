//! Persistent JSON-lines adapter for matched, warmed compilation measurements.
//! The timed region calls the public pipeline and builds its configuration.
//! Parsing, serialization, transport, and independent verification are separate.
use q_rust::backend::Backend;
use q_rust::parser::parse_qasm;
use q_rust::transpiler::{transpile_with_report, RoutingSearchConfig, TranspilerConfig};
use serde_json::{json, Value};
use std::io::{self, BufRead, Write};
use std::time::Instant;

fn compile_one(request: &Value) -> Result<Value, String> {
    let source = request["qasm"].as_str().ok_or("missing qasm")?;
    let topology = request["topology"].as_str().ok_or("missing topology")?;
    let width = request["physical_width"]
        .as_u64()
        .ok_or("missing physical width")? as usize;
    let logical = request["n_logical"]
        .as_u64()
        .ok_or("missing logical width")? as usize;
    let edges: Vec<(usize, usize)> = request["edges"]
        .as_array()
        .ok_or("missing edges")?
        .iter()
        .map(|edge| {
            let a = edge[0].as_u64().ok_or("invalid edge")? as usize;
            let b = edge[1].as_u64().ok_or("invalid edge")? as usize;
            if a >= width || b >= width || a == b {
                return Err("invalid edge");
            }
            Ok((a, b))
        })
        .collect::<Result<_, _>>()?;
    let parse_started = Instant::now();
    let circuit = parse_qasm(source).map_err(|e| e.to_string())?;
    let parse_ns = parse_started.elapsed().as_nanos() as u64;
    if circuit.num_qubits != logical || width < logical {
        return Err("width mismatch".into());
    }

    for field in ["layout_trials", "layout_iterations", "lookahead_size"] {
        if request.get(field).is_some() && request[field].as_u64().is_none() {
            return Err(format!("{field} must be an unsigned integer"));
        }
    }
    for field in [
        "rank_progress",
        "elide_output_permutations",
        "orient_swaps",
        "profile_passes",
        "two_qubit_block_synthesis",
        "recover_phase_gadgets",
        "commute_diagonal",
        "post_block_fusion",
    ] {
        if request.get(field).is_some() && request[field].as_bool().is_none() {
            return Err(format!("{field} must be a boolean"));
        }
    }
    let mut search = RoutingSearchConfig::default();
    if let Some(value) = request["layout_trials"].as_u64() {
        search.layout_trials = Some(value as usize);
    }
    if let Some(value) = request["layout_iterations"].as_u64() {
        search.layout_iterations = Some(value as usize);
    }
    if let Some(value) = request["lookahead_size"].as_u64() {
        search.lookahead_size = value as usize;
    }
    if let Some(value) = request["rank_progress"].as_bool() {
        search.rank_progress = value;
    }
    if let Some(value) = request["elide_output_permutations"].as_bool() {
        search.elide_output_permutations = value;
    }
    if let Some(value) = request["orient_swaps"].as_bool() {
        search.orient_swaps = value;
    }
    if let Some(value) = request["commute_diagonal"].as_bool() {
        search.commute_diagonal = value;
    }
    let effective_search = json!({
        "layout_trials":search.layout_trials.unwrap_or(10),
        "layout_iterations":search.layout_iterations.unwrap_or(3),
        "lookahead_size":search.lookahead_size, "rank_progress":search.rank_progress,
        "elide_output_permutations":search.elide_output_permutations, "orient_swaps":search.orient_swaps,
        "commute_diagonal":search.commute_diagonal,
        "tie_breaking":"deterministic lexicographic", "layout_start_seeds":"0..layout_trials"
    });
    let block_synthesis_enabled = request["two_qubit_block_synthesis"]
        .as_bool()
        .unwrap_or(true);
    let phase_recovery_enabled = request["recover_phase_gadgets"].as_bool().unwrap_or(true);
    let post_block_fusion_enabled = request["post_block_fusion"].as_bool().unwrap_or(true);
    let started = Instant::now();
    let mut backend = Backend::new(topology, width);
    backend.set_coupling_map(edges.iter().flat_map(|&(a, b)| [(a, b), (b, a)]));
    let config = TranspilerConfig::builder()
        .routing_search(search)
        .optimization_level(3)
        .profile_passes(request["profile_passes"].as_bool().unwrap_or(false))
        .two_qubit_block_synthesis(block_synthesis_enabled)
        .recover_phase_gadgets(phase_recovery_enabled)
        .post_block_fusion(post_block_fusion_enabled)
        .backend(backend)
        .target_basis(["u", "cx"])
        .build();
    let (output, report) =
        transpile_with_report(&circuit, Some(config)).map_err(|e| e.to_string())?;
    let compile_ns = started.elapsed().as_nanos() as u64;
    let pass_diagnostics: Vec<Value> = report.pass_diagnostics.iter().map(|p| json!({
        "pass":p.output.stage, "elapsed_ns":p.elapsed_ns, "layout_route_evaluations":p.layout_route_evaluations, "operations":p.output.num_ops,
        "two_qubit_gates":p.output.num_2q_gates, "swaps":p.output.num_swaps, "depth":p.output.depth
    })).collect();
    let block_stats = report.two_qubit_block_stats.as_ref().map(|s| {
        json!({
            "blocks_considered":s.blocks_considered, "synthesis_attempted":s.synthesis_attempted,
            "pure_cx_fast_paths":s.pure_cx_fast_paths, "accepted":s.accepted,
            "cx_removed":s.cx_removed, "rejected_numerical":s.rejected_numerical,
            "rejected_no_cx_improvement":s.rejected_no_cx_improvement,
            "residual_sum":s.residual_sum, "max_residual":s.max_residual
        })
    });
    let fusion_stats = report.post_block_fusion_stats.as_ref().map(|s| {
        json!({
            "runs_attempted":s.runs_attempted, "accepted":s.accepted,
            "u_gates_removed":s.u_gates_removed, "rejected_numerical":s.rejected_numerical,
            "rejected_budget":s.rejected_budget, "residual_sum":s.residual_sum,
            "max_residual":s.max_residual
        })
    });
    let initial = report
        .initial_layout
        .unwrap_or_else(|| (0..logical).collect());
    let final_layout = report
        .final_layout
        .unwrap_or_else(|| (0..logical).collect());
    Ok(json!({
        "status":"compiled", "compile_ns":compile_ns, "parse_ns":parse_ns,
        "out_qasm":output.to_qasm(None), "output_qubits":output.num_qubits,
        "initial_layout":initial[..logical].to_vec(), "final_layout":final_layout[..logical].to_vec(),
        "seed_used":false, "pass_diagnostics":pass_diagnostics, "effective_search":effective_search,
        "two_qubit_block_synthesis":block_synthesis_enabled, "two_qubit_block_stats":block_stats,
        "recover_phase_gadgets":phase_recovery_enabled, "post_block_fusion":post_block_fusion_enabled,
        "phase_gadgets_recovered":report.phase_gadgets_recovered, "post_block_fusion_stats":fusion_stats,
        "pipeline":"Q-Rust public transpile_with_report O3; U/CX; optional exact phase recovery and diagonal routing dependencies; deterministic candidate ties and seeded layout starts; optional post-routing two-qubit block synthesis and checked local fusion"
    }))
}

fn main() {
    let stdin = io::stdin();
    let mut stdout = io::BufWriter::new(io::stdout().lock());
    writeln!(stdout, "{}", json!({"ready":true,"tool":"qrust"})).unwrap();
    stdout.flush().unwrap();
    for line in stdin.lock().lines() {
        let result = line.map_err(|e| e.to_string()).and_then(|line| {
            serde_json::from_str::<Value>(&line)
                .map_err(|e| e.to_string())
                .and_then(|v| compile_one(&v))
        });
        let value = match result {
            Ok(value) => value,
            Err(error) => json!({"status":"error","error":error}),
        };
        if writeln!(stdout, "{value}")
            .and_then(|_| stdout.flush())
            .is_err()
        {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request() -> Value {
        json!({
            "qasm":"OPENQASM 2.0; include \"qelib1.inc\"; qreg q[2]; cx q[0], q[1];",
            "topology":"line", "physical_width":2, "n_logical":2, "edges":[[0,1]]
        })
    }

    #[test]
    fn block_synthesis_override_has_strict_types_and_effective_metadata() {
        for bad in [json!("false"), json!(0), Value::Null] {
            let mut input = request();
            input["two_qubit_block_synthesis"] = bad;
            assert!(compile_one(&input)
                .unwrap_err()
                .contains("must be a boolean"));
        }
        for enabled in [false, true] {
            let mut input = request();
            input["two_qubit_block_synthesis"] = json!(enabled);
            let result = compile_one(&input).unwrap();
            assert_eq!(result["two_qubit_block_synthesis"], json!(enabled));
            assert_eq!(!result["two_qubit_block_stats"].is_null(), enabled);
            assert_eq!(result["pass_diagnostics"], json!([]));
        }
        assert_eq!(
            compile_one(&request()).unwrap()["two_qubit_block_synthesis"],
            json!(true)
        );
    }

    #[test]
    fn refinement_controls_reject_non_booleans_and_record_explicit_values() {
        for field in [
            "recover_phase_gadgets",
            "commute_diagonal",
            "post_block_fusion",
        ] {
            for bad in [json!("false"), json!(0), Value::Null] {
                let mut input = request();
                input[field] = bad;
                assert!(compile_one(&input)
                    .unwrap_err()
                    .contains("must be a boolean"));
            }
            for enabled in [false, true] {
                let mut input = request();
                input[field] = json!(enabled);
                let result = compile_one(&input).unwrap();
                let effective = if field == "commute_diagonal" {
                    &result["effective_search"][field]
                } else {
                    &result[field]
                };
                assert_eq!(effective, &json!(enabled));
            }
        }
    }
}
