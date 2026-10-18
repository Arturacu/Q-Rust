//! Exact recognition of named-Clifford encodings of SWAP.
//!
//! Equality of the signed images of X0, X1, Z0 and Z1 establishes unitary
//! equality up to one global phase. No floating-point gate classification is
//! used. Dependency-local windows contain at most 64 operations, stop at any
//! third-wire interaction, and never cross barriers or dynamic operations.

use super::pass::Pass;
use super::property_set::PropertySet;
use crate::ir::{Circuit, GateType, Operation};

const WINDOW: usize = 64;

#[derive(Clone, Copy, Debug, Default)]
pub struct EncodedSwapRecognitionPass;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EncodedSwapStats {
    pub candidates: usize,
    pub swaps_recognized: usize,
    pub operations_removed: usize,
}

// Hermitian Pauli convention: (-1)^negative i^(x.z) X^x Z^z.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Pauli {
    x: u8,
    z: u8,
    negative: bool,
}

impl Pauli {
    fn h(&mut self, q: usize) {
        let bit = 1 << q;
        let x = self.x & bit != 0;
        let z = self.z & bit != 0;
        self.negative ^= x && z;
        if x != z {
            self.x ^= bit;
            self.z ^= bit;
        }
    }
    fn s(&mut self, q: usize) {
        let bit = 1 << q;
        self.negative ^= self.x & self.z & bit != 0;
        if self.x & bit != 0 {
            self.z ^= bit;
        }
    }
    fn cx(&mut self, control: usize, target: usize) {
        let xc = self.x & (1 << control) != 0;
        let xt = self.x & (1 << target) != 0;
        let zc = self.z & (1 << control) != 0;
        let zt = self.z & (1 << target) != 0;
        self.negative ^= xc && zt && (xt ^ zc ^ true);
        if xc {
            self.x ^= 1 << target;
        }
        if zt {
            self.z ^= 1 << control;
        }
    }
    fn gate(&mut self, name: &GateType, qs: &[usize]) {
        use GateType::*;
        let q = qs[0];
        match name {
            H => self.h(q),
            S => self.s(q),
            Sdg => {
                self.s(q);
                self.s(q);
                self.s(q);
            }
            SX => {
                self.h(q);
                self.s(q);
                self.h(q);
            }
            X => self.negative ^= self.z & (1 << q) != 0,
            Z => self.negative ^= self.x & (1 << q) != 0,
            Y => self.negative ^= (self.x ^ self.z) & (1 << q) != 0,
            CX => self.cx(q, qs[1]),
            CY => {
                self.s(qs[1]);
                self.s(qs[1]);
                self.s(qs[1]);
                self.cx(q, qs[1]);
                self.s(qs[1]);
            }
            CZ => {
                self.h(qs[1]);
                self.cx(q, qs[1]);
                self.h(qs[1]);
            }
            ID => {}
            _ => unreachable!("only supported named Clifford gates reach the tableau"),
        }
    }
}

fn supported(op: &Operation) -> bool {
    matches!(op, Operation::Gate { name, qubits, params }
    if params.is_empty() && match name {
        GateType::H | GateType::S | GateType::Sdg | GateType::SX |
        GateType::X | GateType::Y | GateType::Z | GateType::ID => qubits.len() == 1,
        GateType::CX | GateType::CY | GateType::CZ => qubits.len() == 2 && qubits[0] != qubits[1],
        _ => false,
    })
}

fn process_segment(
    ops: &mut [Operation],
    removed: &mut [bool],
    width: usize,
    stats: &mut EncodedSwapStats,
) {
    let mut wires = vec![Vec::new(); width];
    for (i, op) in ops.iter().enumerate() {
        for &q in op.qubits() {
            wires[q].push(i);
        }
    }
    for start in 0..ops.len() {
        if removed[start] || !supported(&ops[start]) {
            continue;
        }
        let first = ops[start].qubits()[0];
        let mut pair = None;
        let pos = wires[first].partition_point(|&i| i < start);
        for &i in wires[first][pos..]
            .iter()
            .filter(|&&i| !removed[i])
            .take(WINDOW)
        {
            if !supported(&ops[i]) {
                break;
            }
            if ops[i].qubits().len() == 2 {
                pair = Some([
                    first,
                    *ops[i].qubits().iter().find(|&&q| q != first).unwrap(),
                ]);
                break;
            }
        }
        let Some(pair) = pair else {
            continue;
        };
        stats.candidates += 1;
        let mut positions = pair.map(|q| wires[q].partition_point(|&i| i < start));
        let mut tableau = [
            Pauli {
                x: 1,
                z: 0,
                negative: false,
            },
            Pauli {
                x: 2,
                z: 0,
                negative: false,
            },
            Pauli {
                x: 0,
                z: 1,
                negative: false,
            },
            Pauli {
                x: 0,
                z: 2,
                negative: false,
            },
        ];
        let swap = [tableau[1], tableau[0], tableau[3], tableau[2]];
        let mut consumed = Vec::new();
        let mut anchor = None;
        let mut entanglers = 0;
        for _ in 0..WINDOW {
            for j in 0..2 {
                while positions[j] < wires[pair[j]].len() && removed[wires[pair[j]][positions[j]]] {
                    positions[j] += 1;
                }
            }
            let next = (0..2)
                .filter_map(|j| wires[pair[j]].get(positions[j]).copied())
                .min();
            let Some(i) = next else {
                break;
            };
            let op = &ops[i];
            if !supported(op) || op.qubits().iter().any(|q| !pair.contains(q)) {
                break;
            }
            let Operation::Gate { name, qubits, .. } = op else {
                unreachable!()
            };
            let local: Vec<_> = qubits.iter().map(|q| usize::from(*q == pair[1])).collect();
            if local.len() == 2 {
                entanglers += 1;
                anchor.get_or_insert(i);
            }
            for p in &mut tableau {
                p.gate(name, &local);
            }
            consumed.push(i);
            for j in 0..2 {
                if wires[pair[j]].get(positions[j]) == Some(&i) {
                    positions[j] += 1;
                }
            }
            if entanglers >= 3 && tableau == swap {
                let anchor = anchor.unwrap();
                // The anchor already touches both wires, preserving the wire
                // index lists. All intervening unconsumed gates are disjoint.
                for &index in &consumed {
                    removed[index] = index != anchor;
                }
                ops[anchor] = Operation::Gate {
                    name: GateType::SWAP,
                    qubits: pair.to_vec(),
                    params: vec![],
                };
                stats.swaps_recognized += 1;
                stats.operations_removed += consumed.len() - 1;
                break;
            }
        }
    }
}

impl Pass for EncodedSwapRecognitionPass {
    fn name(&self) -> &str {
        "EncodedSwapRecognitionPass"
    }
    fn run(&self, circuit: &Circuit, properties: &mut PropertySet) -> Circuit {
        if super::validate_circuit_structure(circuit).is_err() {
            return circuit.clone();
        }
        let mut out = circuit.clone();
        let mut removed = vec![false; out.operations.len()];
        let mut stats = EncodedSwapStats::default();
        let mut start = 0;
        for end in 0..=out.operations.len() {
            if end == out.operations.len() || !matches!(out.operations[end], Operation::Gate { .. })
            {
                process_segment(
                    &mut out.operations[start..end],
                    &mut removed[start..end],
                    out.num_qubits,
                    &mut stats,
                );
                start = end + 1;
            }
        }
        out.operations = out
            .operations
            .into_iter()
            .enumerate()
            .filter_map(|(i, op)| (!removed[i]).then_some(op))
            .collect();
        properties.insert("encoded_swap_stats", stats);
        out
    }
}
