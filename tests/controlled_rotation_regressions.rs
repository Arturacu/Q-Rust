use q_rust::ir::{Circuit, GateType, Operation};
use q_rust::transpiler::{
    optimization::{InverseCancellationPass, RotationMergePass},
    pass::Pass,
    property_set::PropertySet,
};
use std::f64::consts::PI;

fn cancelled_pair(gate: GateType, first: f64, second: f64, controlled: bool) -> Circuit {
    let mut circuit = Circuit::new(if controlled { 2 } else { 1 }, 0);
    for angle in [first, second] {
        circuit.add_op(Operation::Gate {
            name: gate.clone(),
            qubits: if controlled { vec![0, 1] } else { vec![0] },
            params: vec![angle],
        });
    }
    InverseCancellationPass.run(&circuit, &mut PropertySet::new())
}

#[test]
fn controlled_full_turn_is_a_relative_phase_not_an_identity() {
    // For each Pauli axis R(2π)=-I; controlled-R(2π)=diag(I,-I).
    // This is Z on the control, so deleting the pair is incorrect even up
    // to global phase. This assertion is independent of the gate simulator.
    for gate in [GateType::CRX, GateType::CRY, GateType::CRZ] {
        for sign in [-1.0, 1.0] {
            let result = cancelled_pair(gate.clone(), sign * PI, sign * PI, true);
            assert_eq!(result.operations.len(), 2, "{gate:?}: lost relative phase");
        }
    }
}

#[test]
fn controlled_true_inverses_and_double_turns_cancel() {
    for gate in [GateType::CRX, GateType::CRY, GateType::CRZ] {
        for (a, b) in [(0.73, -0.73), (2.0 * PI, 2.0 * PI), (-2.0 * PI, -2.0 * PI)] {
            assert!(cancelled_pair(gate.clone(), a, b, true)
                .operations
                .is_empty());
        }
    }
}

#[test]
fn uncontrolled_full_turns_can_cancel_up_to_global_phase() {
    for gate in [GateType::RX, GateType::RY, GateType::RZ] {
        for sign in [-1.0, 1.0] {
            assert!(cancelled_pair(gate.clone(), sign * PI, sign * PI, false)
                .operations
                .is_empty());
        }
    }
}

#[test]
fn cancellation_does_not_accept_nonfinite_or_nonperiodic_angles() {
    for gate in [GateType::RX, GateType::CRX] {
        let controlled = gate == GateType::CRX;
        for (a, b) in [
            (f64::NAN, 0.0),
            (f64::INFINITY, 0.0),
            (0.0, 1e-5),
            (PI, PI + 1e-5),
        ] {
            assert_eq!(
                cancelled_pair(gate.clone(), a, b, controlled)
                    .operations
                    .len(),
                2
            );
        }
    }
}

fn controlled_pair_with_intervening_gate(gate: GateType, intervening_qubit: usize) -> Circuit {
    let mut circuit = Circuit::new(3, 0);
    circuit.add_op(Operation::Gate {
        name: gate.clone(),
        qubits: vec![0, 1],
        params: vec![0.4],
    });
    circuit.add_op(Operation::Gate {
        name: GateType::H,
        qubits: vec![intervening_qubit],
        params: vec![],
    });
    circuit.add_op(Operation::Gate {
        name: gate,
        qubits: vec![0, 1],
        params: vec![0.7],
    });
    circuit
}

#[test]
fn controlled_rotation_merge_requires_adjacency_on_both_operands() {
    // H on either operand does not commute with these controlled rotations.
    // The other operand still gives the pair a direct DAG edge, which alone
    // must not authorize CRP(a) ... CRP(b) -> CRP(a+b).
    for gate in [GateType::CRX, GateType::CRY, GateType::CRZ] {
        for intervening_qubit in [0, 1] {
            let circuit = controlled_pair_with_intervening_gate(gate.clone(), intervening_qubit);
            let result = RotationMergePass.run(&circuit, &mut PropertySet::new());
            assert_eq!(
                result.operations, circuit.operations,
                "{gate:?} merged across H on operand {intervening_qubit}"
            );
        }
    }
}

#[test]
fn controlled_rotation_merge_keeps_adjacent_and_disjoint_wire_optimizations() {
    for gate in [GateType::CRX, GateType::CRY, GateType::CRZ] {
        for intervening_disjoint_gate in [false, true] {
            let mut circuit = controlled_pair_with_intervening_gate(gate.clone(), 2);
            if !intervening_disjoint_gate {
                circuit.operations.remove(1);
            }
            let result = RotationMergePass.run(&circuit, &mut PropertySet::new());
            assert_eq!(
                result.operations.len(),
                1 + usize::from(intervening_disjoint_gate)
            );
            let rotations: Vec<_> = result
                .operations
                .iter()
                .filter_map(|op| match op {
                    Operation::Gate {
                        name,
                        qubits,
                        params,
                    } if name == &gate => Some((qubits, params)),
                    _ => None,
                })
                .collect();
            assert_eq!(rotations.len(), 1);
            assert_eq!(rotations[0].0, &[0, 1]);
            assert!((rotations[0].1[0] - 1.1).abs() < 1e-12);
        }
    }
}
