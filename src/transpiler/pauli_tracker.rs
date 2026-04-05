use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use std::collections::HashMap;

/// A topological pass that tracks virtual Pauli frames (X, Y, Z) and propagates them
/// entirely through Clifford operations (H, S, CX, SWAP), annihilating sequential Paulis
/// globally without issuing literal physical gates unless forced by a boundary flush.
pub struct PauliTrackerPass;

impl Pass for PauliTrackerPass {
    fn name(&self) -> &str {
        "PauliTrackerPass"
    }

    fn run(
        &self,
        circuit: &Circuit,
        _property_set: &mut crate::transpiler::property_set::PropertySet,
    ) -> Circuit {
        // Tracks (has_x, has_z) per qubit wire via a virtual un-materialized frame.
        let mut frames: HashMap<usize, (bool, bool)> = HashMap::new();
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);

        // Closes the topology and natively flushes the virtual Pauli state into exactly
        // the required localized physical operators right before evaluating a non-Clifford gate.
        // It strictly drops global phase boundaries dynamically (physics-safe tracking).
        let flush = |q: usize, out_circ: &mut Circuit, f: &mut HashMap<usize, (bool, bool)>| {
            if let Some((x, z)) = f.remove(&q) {
                if x && z {
                    out_circ.add_op(Operation::Gate {
                        name: GateType::Y,
                        qubits: vec![q],
                        params: vec![],
                    });
                } else if x {
                    out_circ.add_op(Operation::Gate {
                        name: GateType::X,
                        qubits: vec![q],
                        params: vec![],
                    });
                } else if z {
                    out_circ.add_op(Operation::Gate {
                        name: GateType::Z,
                        qubits: vec![q],
                        params: vec![],
                    });
                }
            }
        };

        for op in &circuit.operations {
            match op {
                Operation::Gate { name, qubits, .. } => {
                    if qubits.len() == 1 {
                        let q = qubits[0];
                        match name {
                            GateType::X => {
                                let (x, _z) = frames.entry(q).or_insert((false, false));
                                *x = !*x;
                                continue;
                            }
                            GateType::Z => {
                                let (_x, z) = frames.entry(q).or_insert((false, false));
                                *z = !*z;
                                continue;
                            }
                            GateType::Y => {
                                let (x, z) = frames.entry(q).or_insert((false, false));
                                *x = !*x;
                                *z = !*z;
                                continue;
                            }
                            GateType::H => {
                                let (x, z) = frames.entry(q).or_insert((false, false));
                                // By Clifford definition, H exchanges X and Z frames equivalently.
                                std::mem::swap(x, z);
                                new_circuit.add_op(op.clone());
                                continue;
                            }
                            GateType::S | GateType::Sdg => {
                                let (x, z) = frames.entry(q).or_insert((false, false));
                                // S/Sdg rotates the X frame into Y (+/- XZ), while leaving Z unchanged.
                                // Global phase sign shifts are formally dropped as they track relative.
                                if *x {
                                    *z = !*z;
                                }
                                new_circuit.add_op(op.clone());
                                continue;
                            }
                            GateType::ID => {
                                continue;
                            }
                            _ => {
                                // For arbitrary parametrics (RX, U) we strictly flush exactly its wire.
                                flush(q, &mut new_circuit, &mut frames);
                                new_circuit.add_op(op.clone());
                            }
                        }
                    } else if qubits.len() == 2 && name == &GateType::CX {
                        let ctrl = qubits[0];
                        let tgt = qubits[1];

                        let (cx, cz) = frames.entry(ctrl).or_insert((false, false)).clone();
                        let (tx, tz) = frames.entry(tgt).or_insert((false, false)).clone();

                        let next_cx = cx;
                        let mut next_cz = cz;
                        let mut next_tx = tx;
                        let next_tz = tz;

                        // CX structurally forces X on the control to propagate through to the target.
                        if cx {
                            next_tx = !next_tx;
                        }
                        // CX structurally forces Z on the target to propagate backward to the control.
                        if tz {
                            next_cz = !next_cz;
                        }

                        frames.insert(ctrl, (next_cx, next_cz));
                        frames.insert(tgt, (next_tx, next_tz));

                        new_circuit.add_op(op.clone());
                    } else if qubits.len() == 2 && name == &GateType::SWAP {
                        let q1 = qubits[0];
                        let q2 = qubits[1];
                        let f1 = frames.entry(q1).or_insert((false, false)).clone();
                        let f2 = frames.entry(q2).or_insert((false, false)).clone();
                        frames.insert(q1, f2);
                        frames.insert(q2, f1);

                        new_circuit.add_op(op.clone());
                    } else {
                        // Blanket conservative flush explicitly mapping all intersecting support.
                        for &sq in qubits {
                            flush(sq, &mut new_circuit, &mut frames);
                        }
                        new_circuit.add_op(op.clone());
                    }
                }
                Operation::Measure { qubit, .. } | Operation::Reset { qubit } => {
                    flush(*qubit, &mut new_circuit, &mut frames);
                    new_circuit.add_op(op.clone());
                }
                Operation::Barrier { qubits } => {
                    for &sq in qubits {
                        flush(sq, &mut new_circuit, &mut frames);
                    }
                    new_circuit.add_op(op.clone());
                }
            }
        }

        // Final structural graph boundary limit closure.
        for q in 0..circuit.num_qubits {
            flush(q, &mut new_circuit, &mut frames);
        }

        new_circuit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pauli_tracker_annihilates_trivial_spans() {
        let mut circ = Circuit::new(1, 0);
        circ.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        circ.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circ.add_op(Operation::Gate {
            name: GateType::Z,
            qubits: vec![0],
            params: vec![],
        });

        let pass = PauliTrackerPass;
        let c2 = pass.run(
            &circ,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );
        // X passes through H becoming Z, then cancels with the trailing Z. Only H is emitted.
        assert_eq!(c2.operations.len(), 1);
        if let Operation::Gate { name, .. } = &c2.operations[0] {
            assert_eq!(*name, GateType::H);
        } else {
            panic!("Expected H gate");
        }
    }

    #[test]
    fn test_pauli_tracker_entanglement_cross_propagation() {
        let mut circ = Circuit::new(2, 0);
        circ.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        circ.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = PauliTrackerPass;
        let c2 = pass.run(
            &circ,
            &mut crate::transpiler::property_set::PropertySet::new(),
        );

        // CX followed by X on 0 and X on 1
        assert_eq!(c2.operations.len(), 3);

        if let Operation::Gate { name, .. } = &c2.operations[0] {
            assert_eq!(*name, GateType::CX);
        }
        // It flushes at the boundary
        let has_x0 = c2.operations.iter().any(
            |op| matches!(op, Operation::Gate { name: GateType::X, qubits, .. } if qubits[0] == 0),
        );
        let has_x1 = c2.operations.iter().any(
            |op| matches!(op, Operation::Gate { name: GateType::X, qubits, .. } if qubits[0] == 1),
        );
        assert!(has_x0 && has_x1);
    }
}
