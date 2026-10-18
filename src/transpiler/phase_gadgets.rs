//! Recover exact parity-phase interactions before hardware routing.
//!
//! CX(c,t) RZ_t(theta) CX(c,t) = exp(-i theta Z_c Z_t / 2).
//! This is a change of representation, not an approximate rotation merge.
//! Keeping the interaction atomic exposes its diagonal commutation to routing.
//! Phase-gadget representations are established; see Cowtan et al.,
//! "Phase Gadget Synthesis for Shallow Circuits", EPTCS 318 (2020),
//! <https://arxiv.org/abs/1906.01734>. This pass implements the two-wire identity.

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::property_set::PropertySet;

#[derive(Clone, Copy, Debug)]
pub struct RecoverPhaseGadgetsPass;

impl RecoverPhaseGadgetsPass {
    fn recover(&self, circuit: &Circuit, properties: &mut PropertySet) -> Circuit {
        let n = circuit.operations.len();
        // Next operation on each operand, with every non-gate operation a
        // global fence. Disjoint gates can be skipped without repeated scans.
        let mut next = vec![n; circuit.num_qubits];
        let mut successors = vec![[n; 2]; n];
        for (i, op) in circuit.operations.iter().enumerate().rev() {
            if let Operation::Gate { qubits, .. } = op {
                for (slot, &q) in qubits.iter().take(2).enumerate() {
                    successors[i][slot] = next[q];
                }
                for &q in qubits {
                    next[q] = i;
                }
            } else {
                next.fill(i);
            }
        }
        let mut removed = vec![false; n];
        let mut replacements = vec![None; n];
        let mut recovered = 0usize;
        for (i, op) in circuit.operations.iter().enumerate() {
            if removed[i] {
                continue;
            }
            let Operation::Gate {
                name: GateType::CX,
                qubits: pair,
                ..
            } = op
            else {
                continue;
            };
            if pair.len() != 2 {
                continue;
            }
            let middle = successors[i][0].min(successors[i][1]);
            let Some(Operation::Gate {
                name: GateType::RZ,
                qubits,
                params,
            }) = circuit.operations.get(middle)
            else {
                continue;
            };
            if qubits.as_slice() != [pair[1]] || params.len() != 1 || !params[0].is_finite() {
                continue;
            }
            let end = successors[i][0].min(successors[middle][0]);
            let Some(Operation::Gate {
                name: GateType::CX,
                qubits: closing_pair,
                ..
            }) = circuit.operations.get(end)
            else {
                continue;
            };
            if closing_pair != pair {
                continue;
            }
            replacements[i] = Some(Operation::Gate {
                name: GateType::RZZ,
                qubits: pair.clone(),
                params: params.clone(),
            });
            removed[middle] = true;
            removed[end] = true;
            recovered += 1;
        }
        let mut output = circuit.clone();
        output.operations = circuit
            .operations
            .iter()
            .enumerate()
            .filter_map(|(i, op)| {
                if removed[i] {
                    None
                } else {
                    Some(replacements[i].take().unwrap_or_else(|| op.clone()))
                }
            })
            .collect();
        properties.insert("phase_gadgets_recovered", recovered);
        output
    }
}

impl Pass for RecoverPhaseGadgetsPass {
    fn name(&self) -> &str {
        "RecoverPhaseGadgetsPass"
    }

    fn run(&self, circuit: &Circuit, properties: &mut PropertySet) -> Circuit {
        match self.try_run(circuit, properties) {
            Ok(output) => output,
            Err(error) => {
                properties.insert("phase_gadgets_recovered", 0usize);
                properties.insert("phase_gadget_error", error.to_string());
                circuit.clone()
            }
        }
    }

    fn try_run(
        &self,
        circuit: &Circuit,
        properties: &mut PropertySet,
    ) -> crate::error::Result<Circuit> {
        crate::transpiler::validate_circuit_structure(circuit)?;
        properties.remove::<String>("phase_gadget_error");
        Ok(self.recover(circuit, properties))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gate(name: GateType, qubits: &[usize], params: &[f64]) -> Operation {
        Operation::Gate {
            name,
            qubits: qubits.to_vec(),
            params: params.to_vec(),
        }
    }

    #[test]
    fn malformed_direct_run_is_conservative_and_does_not_panic() {
        for operations in [
            vec![gate(GateType::CX, &[0, 9], &[])],
            vec![gate(GateType::CX, &[0, 1], &[0.3])],
            vec![
                gate(GateType::CX, &[0, 1], &[]),
                gate(GateType::RZ, &[1], &[f64::NAN]),
                gate(GateType::CX, &[0, 1], &[]),
            ],
        ] {
            let mut circuit = Circuit::new(2, 0);
            circuit.operations = operations;
            let mut properties = PropertySet::new();
            let output = RecoverPhaseGadgetsPass.run(&circuit, &mut properties);
            assert_eq!(output.operations.len(), circuit.operations.len());
            assert_eq!(properties.get::<usize>("phase_gadgets_recovered"), Some(&0));
            assert!(properties.get::<String>("phase_gadget_error").is_some());
            assert!(RecoverPhaseGadgetsPass
                .try_run(&circuit, &mut properties)
                .is_err());
        }
    }

    #[test]
    fn recovery_skips_only_disjoint_gates_and_preserves_exact_angle() {
        let mut circuit = Circuit::new(3, 0);
        circuit.operations = vec![
            gate(GateType::CX, &[1, 0], &[]),
            gate(GateType::H, &[2], &[]),
            gate(GateType::RZ, &[0], &[-9.1]),
            gate(GateType::CX, &[1, 0], &[]),
        ];
        let mut properties = PropertySet::new();
        let output = RecoverPhaseGadgetsPass.run(&circuit, &mut properties);
        assert_eq!(
            output.operations,
            vec![
                gate(GateType::RZZ, &[1, 0], &[-9.1]),
                gate(GateType::H, &[2], &[])
            ]
        );
        assert_eq!(properties.get::<usize>("phase_gadgets_recovered"), Some(&1));
    }

    #[test]
    fn control_rotation_reversed_cx_third_wire_and_barriers_are_not_recovered() {
        for middle in [
            gate(GateType::RZ, &[0], &[0.31]),
            gate(GateType::RX, &[1], &[0.31]),
            gate(GateType::CX, &[1, 2], &[]),
            Operation::Barrier { qubits: vec![2] },
        ] {
            let mut circuit = Circuit::new(3, 0);
            circuit.operations = vec![
                gate(GateType::CX, &[0, 1], &[]),
                middle,
                gate(GateType::RZ, &[1], &[0.2]),
                gate(GateType::CX, &[0, 1], &[]),
            ];
            assert_eq!(
                RecoverPhaseGadgetsPass
                    .run(&circuit, &mut PropertySet::new())
                    .operations,
                circuit.operations
            );
        }
        let mut circuit = Circuit::new(2, 0);
        circuit.operations = vec![
            gate(GateType::CX, &[0, 1], &[]),
            gate(GateType::RZ, &[1], &[0.2]),
            gate(GateType::CX, &[1, 0], &[]),
        ];
        assert_eq!(
            RecoverPhaseGadgetsPass
                .run(&circuit, &mut PropertySet::new())
                .operations,
            circuit.operations
        );
    }
}
