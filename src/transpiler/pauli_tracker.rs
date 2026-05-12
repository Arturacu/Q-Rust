//! Pauli-frame tracking pass (DISABLED by default).

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
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
        let mut frames: HashMap<usize, (bool, bool)> = HashMap::new();
        let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        out.custom_gates = circuit.custom_gates.clone();

        fn flush(q: usize, out: &mut Circuit, f: &mut HashMap<usize, (bool, bool)>) {
            if let Some((x, z)) = f.remove(&q) {
                if x && z {
                    out.add_op(Operation::Gate {
                        name: GateType::Y,
                        qubits: vec![q],
                        params: vec![],
                    });
                } else if x {
                    out.add_op(Operation::Gate {
                        name: GateType::X,
                        qubits: vec![q],
                        params: vec![],
                    });
                } else if z {
                    out.add_op(Operation::Gate {
                        name: GateType::Z,
                        qubits: vec![q],
                        params: vec![],
                    });
                }
            }
        }

        for op in &circuit.operations {
            match op {
                Operation::Gate { name, qubits, .. } if qubits.len() == 1 => {
                    let q = qubits[0];
                    match name {
                        GateType::X => {
                            let e = frames.entry(q).or_insert((false, false));
                            e.0 = !e.0;
                        }
                        GateType::Z => {
                            let e = frames.entry(q).or_insert((false, false));
                            e.1 = !e.1;
                        }
                        GateType::Y => {
                            let e = frames.entry(q).or_insert((false, false));
                            e.0 = !e.0;
                            e.1 = !e.1;
                        }
                        GateType::H => {
                            flush(q, &mut out, &mut frames);
                            out.add_op(op.clone());
                        }
                        GateType::ID => {}
                        _ => {
                            flush(q, &mut out, &mut frames);
                            out.add_op(op.clone());
                        }
                    }
                }
                Operation::Gate { qubits, .. } => {
                    for &q in qubits {
                        flush(q, &mut out, &mut frames);
                    }
                    out.add_op(op.clone());
                }
                Operation::Measure { qubit, .. } | Operation::Reset { qubit } => {
                    flush(*qubit, &mut out, &mut frames);
                    out.add_op(op.clone());
                }
                Operation::Barrier { qubits } => {
                    for &q in qubits {
                        flush(q, &mut out, &mut frames);
                    }
                    out.add_op(op.clone());
                }
                Operation::Conditional { op: inner, .. } => {
                    for &q in inner.qubits() {
                        flush(q, &mut out, &mut frames);
                    }
                    out.add_op(op.clone());
                }
            }
        }

        for q in 0..circuit.num_qubits {
            flush(q, &mut out, &mut frames);
        }
        out
    }
}
