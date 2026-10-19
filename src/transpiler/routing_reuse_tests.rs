//! Compare reuse with the same complete continuation executed independently.
use super::*;
use crate::backend::Backend;
use crate::parser::parse_qasm;

fn source(width: usize, body: &str) -> Circuit {
    parse_qasm(&format!(
        "OPENQASM 2.0; include \"qelib1.inc\"; qreg q[{width}]; creg c[2]; {body}"
    ))
    .unwrap()
}

fn config(width: usize, directed: bool) -> TranspilerConfig {
    let mut backend = Backend::new("reuse regression", width);
    let edges: Vec<_> = (0..width - 1)
        .flat_map(|q| {
            if directed {
                vec![(q, q + 1)]
            } else {
                vec![(q, q + 1), (q + 1, q)]
            }
        })
        .collect();
    backend.set_coupling_map(edges);
    TranspilerConfig::builder()
        .optimization_level(3)
        .backend(backend)
        .target_basis(["u", "cx"])
        .profile_passes(true)
        .build()
}

fn compiled(
    input: &Circuit,
    config: &TranspilerConfig,
    reuse_enabled: bool,
) -> (LoweredCandidate, Vec<PassDiagnostic>, usize) {
    validate_input(input, config).unwrap();
    let (prefix, terminal) = split_terminal_readout(input);
    let mut optimizer = build_pass_manager_for(config, Stage::Optimize).unwrap();
    let prepared = optimizer.try_run(&prefix).unwrap();
    let mut reuse = RoutingReuse::new(config);
    reuse.enabled &= reuse_enabled;
    let mut diagnostics = Vec::new();
    let output = finish_compilation_with_reuse(
        &prepared,
        &optimizer.property_set,
        config,
        terminal,
        &mut diagnostics,
        &mut reuse,
    )
    .unwrap();
    (output, diagnostics, reuse.checkpoints.len())
}

fn assert_candidate_equal(a: &LoweredCandidate, b: &LoweredCandidate) {
    assert_eq!(a.circuit, b.circuit);
    assert_eq!(a.circuit.to_qasm(None), b.circuit.to_qasm(None));
    assert_eq!(a.output_cost(), b.output_cost());
    macro_rules! property {
        ($type:ty, $name:literal) => {
            assert_eq!(
                a.properties.get::<$type>($name),
                b.properties.get::<$type>($name),
                $name
            );
        };
    }
    property!(Vec<usize>, "initial_layout");
    property!(Vec<usize>, "final_layout");
    property!(Vec<usize>, "logical_output_permutation");
    property!(usize, "swaps_inserted");
    property!(numerical_policy::NumericalBudget, "numerical_budget");
    property!(
        clifford_simplification::CliffordSimplificationStats,
        "clifford_simplification_stats"
    );
    property!(report::CliffordRetentionStats, "clifford_retention_stats");
    property!(block_synthesis::TwoQubitBlockStats, "two_qubit_block_stats");
    property!(
        local_cleanup::PostBlockFusionStats,
        "post_block_fusion_stats"
    );
}

#[test]
fn reuse_matches_independent_continuations_with_routing_readout_and_spares() {
    let cases = [
        source(6, "h q[0]; cx q[0],q[1]; cx q[0],q[2]; cx q[0],q[3]; cx q[0],q[4]; cx q[0],q[5];"),
        source(4, "cx q[0],q[1]; cx q[2],q[1]; cx q[1],q[3]; cx q[1],q[0];"),
        source(3, "swap q[0],q[1]; cx q[0],q[2]; cx q[0],q[1];"),
        source(4, "u(pi/2,0,pi) q[0]; rz(1e-9) q[2]; barrier q; cx q[0],q[1]; cx q[0],q[3]; ccx q[0],q[1],q[2]; measure q[0] -> c[0]; measure q[3] -> c[0]; measure q[2] -> c[1];"),
        parse_qasm("OPENQASM 2.0; include \"qelib1.inc\"; gate pair(a) x,y { h x; cx x,y; rz(a) y; } qreg q[4]; pair(pi/7) q[0],q[3]; cx q[0],q[1]; cx q[0],q[2];").unwrap(),
    ];
    for input in &cases {
        for directed in [false, true] {
            for extra in [0, 1] {
                let cfg = config(input.num_qubits + extra, directed);
                let (cached, diagnostics, count) = compiled(input, &cfg, true);
                let (independent, _, empty) = compiled(input, &cfg, false);
                assert_candidate_equal(&cached, &independent);
                assert!(count <= 4);
                assert_eq!(empty, 0);
                for entry in diagnostics
                    .iter()
                    .filter(|d| d.output.stage.ends_with("ReuseRoutingCheckpoint"))
                {
                    assert_eq!(entry.layout_route_evaluations, Some(0));
                }
                let mut plain = cfg.clone();
                plain.profile_passes = false;
                let (unprofiled, no_diagnostics, _) = compiled(input, &plain, true);
                assert_candidate_equal(&cached, &unprofiled);
                assert!(no_diagnostics.is_empty());
            }
        }
    }
}

#[test]
fn equal_prepared_routes_reuse_but_keep_branch_clifford_statistics() {
    let input = source(4, "h q[0]; cx q[0],q[1]; cx q[0],q[2]; cx q[0],q[3];");
    let cfg = config(5, true);
    let prefix = property_set::PropertySet::new();
    let mut reuse = RoutingReuse::new(&cfg);
    let mut first_trace = Vec::new();
    route_and_lower(&input, &prefix, &cfg, &[], "", &mut first_trace, &mut reuse).unwrap();
    assert!(!reuse.checkpoints.is_empty());
    let mut disabled = cfg.clone();
    disabled.clifford_simplification = false;
    let mut cached_trace = Vec::new();
    let cached = route_and_lower(
        &input,
        &prefix,
        &disabled,
        &[],
        "fallback",
        &mut cached_trace,
        &mut reuse,
    )
    .unwrap();
    let mut uncached_trace = Vec::new();
    let independent = route_and_lower(
        &input,
        &prefix,
        &disabled,
        &[],
        "fallback",
        &mut uncached_trace,
        &mut RoutingReuse {
            enabled: false,
            checkpoints: Vec::new(),
        },
    )
    .unwrap();
    assert_eq!(cached.len(), independent.len());
    for (a, b) in cached.iter().zip(&independent) {
        assert_candidate_equal(a, b);
        assert!(a
            .properties
            .get::<clifford_simplification::CliffordSimplificationStats>(
                "clifford_simplification_stats"
            )
            .is_none());
    }
    assert_eq!(
        cached_trace
            .iter()
            .filter(|d| d.output.stage.ends_with("ReuseRoutingCheckpoint"))
            .count(),
        1
    );
    assert!(!cached_trace
        .iter()
        .any(|d| d.output.stage.ends_with("SabreLayoutPass")
            || d.output.stage.ends_with("BeamSabrePass")));
    assert!(uncached_trace
        .iter()
        .any(|d| d.output.stage.ends_with("SabreLayoutPass")));
}

#[test]
fn route_reuse_certificate_distinguishes_bits_maps_ledgers_and_policies() {
    let positive = source(1, "u(0,0,0) q[0];");
    let mut negative = positive.clone();
    if let Operation::Gate { params, .. } = &mut negative.operations[0] {
        params[0] = -0.0;
    }
    assert_eq!(positive, negative);
    assert!(!same_circuit_bits(&positive, &negative));
    let a =
        parse_qasm("OPENQASM 2.0; gate unused(a) q { rz(a+0.0) q; } qreg q[1]; h q[0];").unwrap();
    let mut b = a.clone();
    let crate::ir::ast::ParsedStatement::Gate(_, _, params) =
        &mut b.custom_gates.defs.get_mut("unused").unwrap().body[0]
    else {
        panic!("gate body")
    };
    let crate::ir::ast::Expr::Add(_, literal) = &mut params[0] else {
        panic!("sum expression")
    };
    **literal = crate::ir::ast::Expr::Float(-0.0);
    assert_eq!(a, b);
    assert!(!same_circuit_bits(&a, &b));
    let cfg = config(3, false);
    let properties = property_set::PropertySet::new();
    let reference = RoutingContextKey::new(&cfg, &properties);
    let mut changed = property_set::PropertySet::new();
    changed.insert("logical_output_permutation", vec![1usize, 0]);
    assert!(reference != RoutingContextKey::new(&cfg, &changed));
    let mut changed = property_set::PropertySet::new();
    changed.insert(
        "numerical_budget",
        numerical_policy::NumericalBudget {
            accepted: 1,
            charged_error: 1e-13,
            ..Default::default()
        },
    );
    assert!(reference != RoutingContextKey::new(&cfg, &changed));
    for mode in 0..7 {
        let mut changed = cfg.clone();
        match mode {
            0 => changed.routing_search.elide_output_permutations = false,
            1 => changed.retain_routing_candidates = false,
            2 => changed.routing_search.layout_trials = Some(2),
            3 => changed.routing_search.commute_diagonal = false,
            4 => changed.lookahead_strategy = routing::LookaheadStrategy::DynamicV2,
            5 => changed.backend = Some(Backend::all_to_all(3)),
            6 => changed.routing_search.rank_progress = false,
            _ => unreachable!(),
        }
        assert!(reference != RoutingContextKey::new(&changed, &properties));
    }
}

#[test]
fn route_reuse_key_covers_each_numerical_and_routing_field() {
    let cfg = config(3, false);
    let empty = property_set::PropertySet::new();
    let reference = RoutingContextKey::new(&cfg, &empty);
    for mode in 0..6 {
        let mut budget = numerical_policy::NumericalBudget::default();
        match mode {
            0 => budget.accepted = 1,
            1 => budget.rejected = 1,
            2 => budget.residual_sum = 1e-13,
            3 => budget.charged_error = 1e-12,
            4 => budget.residual_sum = -0.0,
            5 => budget.charged_error = -0.0,
            _ => unreachable!(),
        }
        let mut properties = property_set::PropertySet::new();
        properties.insert("numerical_budget", budget);
        assert!(
            reference != RoutingContextKey::new(&cfg, &properties),
            "ledger field {mode}"
        );
    }
    for mode in 0..5 {
        let mut changed = cfg.clone();
        match mode {
            0 => changed.routing_search.layout_iterations = Some(2),
            1 => changed.routing_search.lookahead_size += 1,
            2 => changed.backend.as_mut().unwrap().num_qubits += 1,
            3 => {
                changed
                    .backend
                    .as_mut()
                    .unwrap()
                    .set_coupling_map([(1, 2), (2, 1), (0, 1), (1, 0)])
            }
            4 => changed
                .backend
                .as_mut()
                .unwrap()
                .set_coupling_map([(0, 1), (1, 2)]),
            _ => unreachable!(),
        }
        assert!(
            reference != RoutingContextKey::new(&changed, &empty),
            "routing field {mode}"
        );
    }
    let mut positive = cfg.clone();
    positive.lookahead_strategy = routing::LookaheadStrategy::Static { weight: 0.0 };
    let mut negative = positive.clone();
    negative.lookahead_strategy = routing::LookaheadStrategy::Static { weight: -0.0 };
    assert!(RoutingContextKey::new(&positive, &empty) != RoutingContextKey::new(&negative, &empty));
}

#[test]
fn restored_route_copies_cannot_mutate_the_cache_or_another_branch() {
    let prepared = source(2, "cx q[0],q[1];");
    let routed = source(3, "cx q[0],q[2];");
    let mut properties = property_set::PropertySet::new();
    properties.insert("initial_layout", vec![0usize, 2]);
    properties.insert("final_layout", vec![2usize, 0]);
    properties.insert("swaps_inserted", 1usize);
    properties.insert(
        "routing_output_candidates",
        vec![routing::RoutingOutputCandidate {
            circuit: source(3, "cx q[2],q[0];"),
            initial_layout: vec![2, 0],
            final_layout: vec![0, 2],
            swaps_inserted: 2,
            origin: "test alternative",
        }],
    );
    let checkpoint = RoutingCheckpoint::capture(
        &prepared,
        RoutingContextKey::new(&config(3, false), &properties),
        &routed,
        &properties,
    );
    let mut first = property_set::PropertySet::new();
    let own_budget = numerical_policy::NumericalBudget {
        accepted: 2,
        charged_error: 4e-13,
        ..Default::default()
    };
    let own_stats = clifford_simplification::CliffordSimplificationStats {
        candidates: 7,
        rewrites: 2,
        ..Default::default()
    };
    first.insert("numerical_budget", own_budget.clone());
    first.insert("clifford_simplification_stats", own_stats.clone());
    first.insert("logical_output_permutation", vec![1usize, 0]);
    let mut first_circuit = checkpoint.restore(&mut first);
    assert_eq!(numerical_policy::snapshot(&first), own_budget);
    assert_eq!(
        first.get::<clifford_simplification::CliffordSimplificationStats>(
            "clifford_simplification_stats"
        ),
        Some(&own_stats)
    );
    assert_eq!(
        first.get::<Vec<usize>>("logical_output_permutation"),
        Some(&vec![1, 0])
    );
    first_circuit.operations.clear();
    let mut first_map = first.remove::<Vec<usize>>("initial_layout").unwrap();
    first_map[0] = 99;
    first.insert("initial_layout", *first_map);
    let mut first_alternatives = first
        .remove::<Vec<routing::RoutingOutputCandidate>>("routing_output_candidates")
        .unwrap();
    first_alternatives[0].circuit.operations.clear();
    first_alternatives[0].initial_layout[0] = 99;
    first_alternatives[0].final_layout[0] = 99;
    let mut second = property_set::PropertySet::new();
    let second_circuit = checkpoint.restore(&mut second);
    assert_eq!(second_circuit, routed);
    assert_eq!(
        second.get::<Vec<usize>>("initial_layout"),
        Some(&vec![0, 2])
    );
    assert_eq!(second.get::<Vec<usize>>("final_layout"), Some(&vec![2, 0]));
    assert_eq!(second.get::<usize>("swaps_inserted"), Some(&1));
    assert!(second
        .get::<clifford_simplification::CliffordSimplificationStats>(
            "clifford_simplification_stats"
        )
        .is_none());
    let alternatives = second
        .get::<Vec<routing::RoutingOutputCandidate>>("routing_output_candidates")
        .unwrap();
    assert_eq!(alternatives[0].circuit, source(3, "cx q[2],q[0];"));
    assert_eq!(alternatives[0].initial_layout, vec![2, 0]);
    assert_eq!(alternatives[0].final_layout, vec![0, 2]);
    assert_eq!(alternatives[0].swaps_inserted, 2);
}

#[test]
fn lazy_cost_caches_unavailable_values_without_changing_selection() {
    let candidate = LoweredCandidate {
        circuit: source(3, "ccx q[0],q[1],q[2];"),
        properties: property_set::PropertySet::new(),
        cost: OnceCell::new(),
    };
    assert!(candidate.cost.get().is_none());
    assert_eq!(candidate.output_cost(), None);
    assert_eq!(candidate.cost.get(), Some(&None));
    assert!(!prefer_lowered(&candidate, &candidate, &candidate));
    let valid = LoweredCandidate {
        circuit: source(3, "cx q[0],q[1];"),
        properties: property_set::PropertySet::new(),
        cost: OnceCell::new(),
    };
    assert_eq!(
        valid.output_cost(),
        candidate_selection::output_cost(&valid.circuit)
    );
    assert!(!prefer_lowered(&valid, &candidate, &valid));
    assert!(!prefer_lowered(&valid, &valid, &candidate));
}

#[test]
fn reuse_is_not_allocated_for_unaffected_configurations() {
    let input = source(3, "h q[0]; cx q[0],q[1]; cx q[0],q[2];");
    for mode in 0..6 {
        let mut cfg = config(3, false);
        match mode {
            0 => cfg.optimization_level = 0,
            1 => cfg.optimization_level = 2,
            2 => cfg.backend = None,
            3 => cfg.retain_routing_candidates = false,
            4 => cfg.clifford_simplification = false,
            5 => cfg.routing_search.elide_output_permutations = false,
            _ => unreachable!(),
        }
        let (cached, trace, count) = compiled(&input, &cfg, true);
        let (independent, _, _) = compiled(&input, &cfg, false);
        assert_eq!(count, 0);
        assert!(!trace
            .iter()
            .any(|d| d.output.stage.ends_with("ReuseRoutingCheckpoint")));
        assert_candidate_equal(&cached, &independent);
    }
}
