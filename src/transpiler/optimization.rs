use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::pass::Pass;
use crate::transpiler::synthesis::zyz::{u_to_matrix, zyz_decomposition, Unitary2x2};

/// A pass that merges consecutive single-qubit gates on the same qubit.
pub struct GateFusionPass;

impl Pass for GateFusionPass {
    fn name(&self) -> &str {
        "GateFusionPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);

        // Pending unitary for each qubit: Option<Unitary2x2>
        // We initialize with Identity (None implies no pending gate, or Identity)
        let mut pending_gates: Vec<Option<Unitary2x2>> = vec![None; circuit.num_qubits];

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } => {
                    if qubits.len() == 1 {
                        // Single qubit gate - accumulate
                        let q = qubits[0];
                        let current_matrix = match name {
                            GateType::U => {
                                if params.len() == 3 {
                                    u_to_matrix(params[0], params[1], params[2])
                                } else {
                                    // Should not happen if parser is correct
                                    flush_qubit(&mut new_circuit, &mut pending_gates, q);
                                    new_circuit.add_op(op.clone());
                                    continue;
                                }
                            }
                            // If we encounter other single qubit gates, we should ideally decompose them first
                            // or handle them here. Assuming basis decomposition ran first, we only see U.
                            // But for robustness, let's assume we might see others if decomposition wasn't run.
                            // For now, let's only fuse U gates. If it's not U, we flush and add it.
                            _ => {
                                flush_qubit(&mut new_circuit, &mut pending_gates, q);
                                new_circuit.add_op(op.clone());
                                continue;
                            }
                        };

                        if let Some(pending) = pending_gates[q] {
                            // Multiply: New * Old (since New is applied after Old)
                            // Matrix multiplication: C = A * B
                            // [c00 c01] = [a00 a01] * [b00 b01]
                            // [c10 c11]   [a10 a11]   [b10 b11]
                            let a = current_matrix;
                            let b = pending;

                            let c00 = a[0][0] * b[0][0] + a[0][1] * b[1][0];
                            let c01 = a[0][0] * b[0][1] + a[0][1] * b[1][1];
                            let c10 = a[1][0] * b[0][0] + a[1][1] * b[1][0];
                            let c11 = a[1][0] * b[0][1] + a[1][1] * b[1][1];

                            pending_gates[q] = Some([[c00, c01], [c10, c11]]);
                        } else {
                            pending_gates[q] = Some(current_matrix);
                        }
                    } else {
                        // Multi-qubit gate - flush involved qubits
                        for &q in qubits {
                            flush_qubit(&mut new_circuit, &mut pending_gates, q);
                        }
                        new_circuit.add_op(op.clone());
                    }
                }
                _ => {
                    // Other ops (Measure, Barrier) - flush all involved or all?
                    // Measure involves specific qubits. Barrier might involve specific or all.
                    // To be safe, let's look at the operation details if possible,
                    // but Operation enum structure for Measure/Barrier is different.
                    // For simplicity/safety, flush ALL qubits on non-gate operations for now,
                    // or try to be more granular.
                    // Measure has qubits. Barrier has qubits.

                    // Let's flush all for safety in this first iteration.
                    for q in 0..circuit.num_qubits {
                        flush_qubit(&mut new_circuit, &mut pending_gates, q);
                    }
                    new_circuit.add_op(op.clone());
                }
            }
        }

        // Flush remaining gates at the end
        for q in 0..circuit.num_qubits {
            flush_qubit(&mut new_circuit, &mut pending_gates, q);
        }

        new_circuit
    }
}

fn flush_qubit(circuit: &mut Circuit, pending_gates: &mut Vec<Option<Unitary2x2>>, q: usize) {
    if let Some(matrix) = pending_gates[q] {
        // Decompose matrix back to U gate
        let (theta, phi, lambda, _gamma) = zyz_decomposition(matrix);
        // We ignore global phase _gamma for the gate parameters,
        // as U gate definition doesn't include it.
        // Note: This drops global phase, which is physically fine but mathematically lossy.

        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![q],
            params: vec![theta, phi, lambda],
        });
        pending_gates[q] = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_gate_fusion_simple() {
        let mut circuit = Circuit::new(1, 0);
        // U(pi/2, 0, pi) -> H
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        // U(pi/2, 0, pi) -> H
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });

        let pass = GateFusionPass;
        let new_circuit = pass.run(&circuit);

        // H * H = I
        // Should result in U(0, 0, 0) or similar identity-equivalent
        assert_eq!(new_circuit.operations.len(), 1);

        if let Operation::Gate { name, params, .. } = &new_circuit.operations[0] {
            if let GateType::U = name {
                assert!(params[0].abs() < 1e-6); // Identity has theta=0
            } else {
                panic!("Expected U gate");
            }
        }
    }

    #[test]
    fn test_gate_fusion_blocking() {
        let mut circuit = Circuit::new(2, 0);
        // H q0
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });
        // CX q0, q1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H q0
        circuit.add_op(Operation::Gate {
            name: GateType::U,
            qubits: vec![0],
            params: vec![PI / 2.0, 0.0, PI],
        });

        let pass = GateFusionPass;
        let new_circuit = pass.run(&circuit);

        // Should NOT merge H gates because CX blocks them
        assert_eq!(new_circuit.operations.len(), 3);
    }
}

/// A pass that removes adjacent identical CX gates (CX a,b; CX a,b -> ID).
pub struct CXCancellationPass;

impl Pass for CXCancellationPass {
    fn name(&self) -> &str {
        "CXCancellationPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let mut skip_next = false;

        for i in 0..circuit.operations.len() {
            if skip_next {
                skip_next = false;
                continue;
            }

            let op = &circuit.operations[i];

            // Check if this is a CX gate
            let is_cx = match op {
                Operation::Gate { name, .. } => *name == GateType::CX,
                _ => false,
            };

            if is_cx && i + 1 < circuit.operations.len() {
                let next_op = &circuit.operations[i + 1];
                if op == next_op {
                    // Found adjacent identical CX gates, skip both
                    skip_next = true;
                    continue;
                }
            }

            new_circuit.add_op(op.clone());
        }

        new_circuit
    }
}

/// A pass that cancels CX gates separated by commuting single-qubit gates.
/// Pattern: CX a,b; [Commuting Ops]; CX a,b -> [Commuting Ops]
pub struct CommutationCancellationPass;

impl Pass for CommutationCancellationPass {
    fn name(&self) -> &str {
        "CommutationCancellationPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let mut i = 0;

        while i < circuit.operations.len() {
            let op = &circuit.operations[i];

            // Check start of potential cancellation block: CX a,b
            if let Operation::Gate {
                name: GateType::CX,
                qubits: cx_qubits,
                ..
            } = op
            {
                let ctrl = cx_qubits[0];
                let target = cx_qubits[1];

                // Look ahead for matching CX
                let mut j = i + 1;
                let mut commute = true;
                let mut intermediate_ops = Vec::new();

                while j < circuit.operations.len() {
                    let next_op = &circuit.operations[j];

                    // If we find the matching CX, we can stop
                    if let Operation::Gate {
                        name: GateType::CX,
                        qubits: next_qubits,
                        ..
                    } = next_op
                    {
                        if next_qubits == cx_qubits {
                            // Found match!
                            break;
                        }
                    }

                    // Check commutation
                    if !commutes_with_cx(next_op, ctrl, target) {
                        commute = false;
                        break;
                    }

                    intermediate_ops.push(next_op.clone());
                    j += 1;
                }

                if commute && j < circuit.operations.len() {
                    // We found a matching CX and everything in between commutes.
                    // Skip the first CX (i), add intermediate ops, and skip the second CX (j).
                    for mid_op in intermediate_ops {
                        new_circuit.add_op(mid_op);
                    }
                    i = j + 1; // Continue after the second CX
                    continue;
                }
            }

            new_circuit.add_op(op.clone());
            i += 1;
        }

        new_circuit
    }
}

fn commutes_with_cx(op: &Operation, ctrl: usize, target: usize) -> bool {
    match op {
        Operation::Gate {
            name,
            qubits,
            params,
        } => {
            // If the gate doesn't involve the CX qubits, it commutes
            let involves_ctrl = qubits.contains(&ctrl);
            let involves_target = qubits.contains(&target);

            if !involves_ctrl && !involves_target {
                return true;
            }

            // If it involves other qubits (multi-qubit gate), assume no commutation for now
            // unless we implement more complex logic.
            if qubits.len() > 1 {
                return false;
            }

            // Single qubit gate logic
            let q = qubits[0];

            if q == ctrl {
                // Gate on control qubit: Commutes if it is Z-basis (Z, RZ, U1, S, T)
                // Anti-commutes if X-basis (X, RX)
                match name {
                    GateType::Z
                    | GateType::S
                    | GateType::Sdg
                    | GateType::T
                    | GateType::Tdg
                    | GateType::RZ => true,
                    GateType::U => {
                        // U(0, 0, lambda) is RZ(lambda) -> Commutes
                        // U(theta, phi, lambda) generally doesn't unless theta=0
                        if params.len() >= 1 && params[0].abs() < 1e-9 {
                            true
                        } else {
                            false
                        }
                    }
                    _ => false, // X, Y, H, etc.
                }
            } else if q == target {
                // Gate on target qubit: Commutes if it is X-basis (X, RX)
                // Anti-commutes if Z-basis (Z, RZ)
                match name {
                    GateType::X | GateType::RX => true,
                    GateType::U => {
                        // U(pi, -pi/2, pi/2) is X -> Commutes
                        // U(theta, -pi/2, pi/2) is RX(theta) -> Commutes
                        if params.len() >= 3 {
                            let _theta = params[0];
                            let phi = params[1];
                            let lambda = params[2];
                            // Check if it's RX-like: phi = -pi/2, lambda = pi/2
                            let pi_2 = std::f64::consts::FRAC_PI_2;
                            if (phi - -pi_2).abs() < 1e-9 && (lambda - pi_2).abs() < 1e-9 {
                                true
                            } else {
                                false
                            }
                        } else {
                            false
                        }
                    }
                    _ => false,
                }
            } else {
                true // Should be covered by first check, but safe fallback
            }
        }
        _ => false, // Measure/Barrier/etc assume no commutation
    }
}

#[cfg(test)]
mod optimization_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_cx_cancellation() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CXCancellationPass;
        let new_circuit = pass.run(&circuit);
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_commutation_cancellation_rz_control() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // RZ(pi) 0 (Commutes with Control)
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![PI],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(&circuit);

        // Should remove CXs, leaving only RZ
        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::RZ);
        } else {
            panic!("Expected RZ");
        }
    }

    #[test]
    fn test_commutation_cancellation_rx_target() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // RX(pi) 1 (Commutes with Target X-component)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![1],
            params: vec![PI],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::RX);
        } else {
            panic!("Expected RX");
        }
    }

    #[test]
    fn test_commutation_blocking() {
        let mut circuit = Circuit::new(2, 0);
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });
        // X 0 (Anti-commutes with Control Z-component)
        circuit.add_op(Operation::Gate {
            name: GateType::X,
            qubits: vec![0],
            params: vec![],
        });
        // CX 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::CX,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = CommutationCancellationPass;
        let new_circuit = pass.run(&circuit);

        // Should NOT cancel
        assert_eq!(new_circuit.operations.len(), 3);
    }
}

/// A pass that removes redundant SWAP gates (SWAP a,b; SWAP a,b -> ID).
pub struct SwapSimplificationPass;

impl Pass for SwapSimplificationPass {
    fn name(&self) -> &str {
        "SwapSimplificationPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let mut i = 0;

        while i < circuit.operations.len() {
            let op = &circuit.operations[i];

            // Check start of potential cancellation block: SWAP a,b
            if let Operation::Gate {
                name: GateType::SWAP,
                qubits: swap_qubits,
                ..
            } = op
            {
                // Look ahead for matching SWAP
                let mut j = i + 1;
                let mut commute = true;
                let mut intermediate_ops = Vec::new();

                while j < circuit.operations.len() {
                    let next_op = &circuit.operations[j];

                    // If we find the matching SWAP, we can stop
                    if let Operation::Gate {
                        name: GateType::SWAP,
                        qubits: next_qubits,
                        ..
                    } = next_op
                    {
                        // Check if qubits match (order doesn't matter for SWAP)
                        if (next_qubits[0] == swap_qubits[0] && next_qubits[1] == swap_qubits[1])
                            || (next_qubits[0] == swap_qubits[1]
                                && next_qubits[1] == swap_qubits[0])
                        {
                            // Found match!
                            break;
                        }
                    }

                    // Check commutation: SWAP(a,b) commutes with any gate that doesn't involve a or b
                    if involves_any(next_op, &swap_qubits) {
                        commute = false;
                        break;
                    }

                    intermediate_ops.push(next_op.clone());
                    j += 1;
                }

                if commute && j < circuit.operations.len() {
                    // We found a matching SWAP and everything in between commutes.
                    // Skip the first SWAP (i), add intermediate ops, and skip the second SWAP (j).
                    for mid_op in intermediate_ops {
                        new_circuit.add_op(mid_op);
                    }
                    i = j + 1; // Continue after the second SWAP
                    continue;
                }
            }

            new_circuit.add_op(op.clone());
            i += 1;
        }

        new_circuit
    }
}

fn involves_any(op: &Operation, qubits: &[usize]) -> bool {
    match op {
        Operation::Gate {
            qubits: op_qubits, ..
        } => {
            for q in op_qubits {
                if qubits.contains(q) {
                    return true;
                }
            }
            false
        }
        Operation::Measure { qubit, .. } => qubits.contains(qubit),
        Operation::Reset { qubit } => qubits.contains(qubit),
        Operation::Barrier {
            qubits: barrier_qubits,
        } => {
            for q in barrier_qubits {
                if qubits.contains(q) {
                    return true;
                }
            }
            false
        }
    }
}

#[cfg(test)]
mod swap_tests {
    use super::*;

    #[test]
    fn test_swap_cancellation_simple() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(&circuit);
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_swap_cancellation_reversed() {
        let mut circuit = Circuit::new(2, 0);
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![1, 0],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(&circuit);
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_swap_cancellation_commuting() {
        let mut circuit = Circuit::new(3, 0);
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H 2 (Commutes)
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![2],
            params: vec![],
        });
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { name, .. } = &new_circuit.operations[0] {
            assert_eq!(*name, GateType::H);
        } else {
            panic!("Expected H");
        }
    }

    #[test]
    fn test_swap_blocking() {
        let mut circuit = Circuit::new(3, 0);
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });
        // H 0 (Blocks)
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        // SWAP 0,1
        circuit.add_op(Operation::Gate {
            name: GateType::SWAP,
            qubits: vec![0, 1],
            params: vec![],
        });

        let pass = SwapSimplificationPass;
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 3);
    }
}

/// A pass that simplifies gate parameters (e.g., modulo 2pi, removing negligible rotations).
pub struct ParameterSimplificationPass {
    pub epsilon: f64,
}

impl Default for ParameterSimplificationPass {
    fn default() -> Self {
        Self { epsilon: 1e-9 }
    }
}

impl Pass for ParameterSimplificationPass {
    fn name(&self) -> &str {
        "ParameterSimplificationPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let two_pi = 2.0 * std::f64::consts::PI;

        for op in &circuit.operations {
            match op {
                Operation::Gate {
                    name,
                    qubits,
                    params,
                } => {
                    let mut new_params = params.clone();
                    let mut keep_gate = true;

                    // Normalize parameters and check for identity
                    match name {
                        GateType::RX | GateType::RY | GateType::RZ | GateType::U => {
                            for param in &mut new_params {
                                // Normalize to [0, 2pi)
                                *param = *param % two_pi;
                                if *param < 0.0 {
                                    *param += two_pi;
                                }
                            }

                            // Check if gate is effectively identity
                            match name {
                                GateType::RX | GateType::RY | GateType::RZ => {
                                    if let Some(theta) = new_params.first() {
                                        if theta.abs() < self.epsilon
                                            || (theta - two_pi).abs() < self.epsilon
                                        {
                                            keep_gate = false;
                                        }
                                    }
                                }
                                GateType::U => {
                                    // U(0, 0, 0) is Identity
                                    // U(0, 0, lambda) is RZ(lambda) - keep unless lambda ~ 0
                                    // U(theta, phi, lambda)
                                    if new_params.len() >= 3 {
                                        let theta = new_params[0];
                                        let phi = new_params[1];
                                        let lambda = new_params[2];

                                        // If theta ~ 0, it's diagonal.
                                        // If phi + lambda ~ 0 (mod 2pi), it's global phase only?
                                        // U(0, phi, lambda) = diag(1, e^(i(phi+lambda)))
                                        // If phi+lambda = 0, it's I.

                                        if theta.abs() < self.epsilon {
                                            let phase_sum = (phi + lambda) % two_pi;
                                            if phase_sum.abs() < self.epsilon
                                                || (phase_sum - two_pi).abs() < self.epsilon
                                            {
                                                keep_gate = false;
                                            }
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                        _ => {}
                    }

                    if keep_gate {
                        new_circuit.add_op(Operation::Gate {
                            name: name.clone(),
                            qubits: qubits.clone(),
                            params: new_params,
                        });
                    }
                }
                _ => new_circuit.add_op(op.clone()),
            }
        }

        new_circuit
    }
}

/// A pass that performs peephole optimization using pattern matching.
/// Currently implements simple self-inverse cancellation (H-H, X-X, etc.).
pub struct PeepholeOptimizationPass;

impl Pass for PeepholeOptimizationPass {
    fn name(&self) -> &str {
        "PeepholeOptimizationPass"
    }

    fn run(&self, circuit: &Circuit) -> Circuit {
        let mut new_circuit = Circuit::new(circuit.num_qubits, circuit.num_cbits);
        let mut skip_next = false;

        for i in 0..circuit.operations.len() {
            if skip_next {
                skip_next = false;
                continue;
            }

            let op = &circuit.operations[i];

            if i + 1 < circuit.operations.len() {
                let next_op = &circuit.operations[i + 1];

                if are_inverses(op, next_op) {
                    skip_next = true;
                    continue;
                }
            }

            new_circuit.add_op(op.clone());
        }

        new_circuit
    }
}

fn are_inverses(op1: &Operation, op2: &Operation) -> bool {
    match (op1, op2) {
        (
            Operation::Gate {
                name: n1,
                qubits: q1,
                params: p1,
            },
            Operation::Gate {
                name: n2,
                qubits: q2,
                params: p2,
            },
        ) => {
            if q1 != q2 {
                return false;
            }

            match (n1, n2) {
                // Self-inverses
                (GateType::H, GateType::H) => true,
                (GateType::X, GateType::X) => true,
                (GateType::Y, GateType::Y) => true,
                (GateType::Z, GateType::Z) => true,
                (GateType::CX, GateType::CX) => true, // Already handled by CXCancellation, but good to have
                (GateType::SWAP, GateType::SWAP) => true, // Already handled

                // Inverse pairs
                (GateType::S, GateType::Sdg) => true,
                (GateType::Sdg, GateType::S) => true,
                (GateType::T, GateType::Tdg) => true,
                (GateType::Tdg, GateType::T) => true,

                // Rotations with opposite angles
                (GateType::RX, GateType::RX)
                | (GateType::RY, GateType::RY)
                | (GateType::RZ, GateType::RZ) => {
                    if let (Some(a1), Some(a2)) = (p1.first(), p2.first()) {
                        (a1 + a2).abs() < 1e-9
                            || (a1 + a2 - 2.0 * std::f64::consts::PI).abs() < 1e-9
                    } else {
                        false
                    }
                }

                _ => false,
            }
        }
        _ => false,
    }
}

#[cfg(test)]
mod param_peephole_tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_param_simplification_identity() {
        let mut circuit = Circuit::new(1, 0);
        // RZ(0) -> Identity
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![0.0],
        });
        // RZ(2pi) -> Identity
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![2.0 * PI],
        });
        // RX(1e-10) -> Identity (negligible)
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![1e-10],
        });

        let pass = ParameterSimplificationPass::default();
        let new_circuit = pass.run(&circuit);
        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_param_simplification_normalization() {
        let mut circuit = Circuit::new(1, 0);
        // RZ(3pi) -> RZ(pi)
        circuit.add_op(Operation::Gate {
            name: GateType::RZ,
            qubits: vec![0],
            params: vec![3.0 * PI],
        });

        let pass = ParameterSimplificationPass::default();
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 1);
        if let Operation::Gate { params, .. } = &new_circuit.operations[0] {
            assert!((params[0] - PI).abs() < 1e-6);
        } else {
            panic!("Expected Gate");
        }
    }

    #[test]
    fn test_peephole_inverses() {
        let mut circuit = Circuit::new(1, 0);
        // H H -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::H,
            qubits: vec![0],
            params: vec![],
        });
        // S Sdg -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::S,
            qubits: vec![0],
            params: vec![],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::Sdg,
            qubits: vec![0],
            params: vec![],
        });

        let pass = PeepholeOptimizationPass;
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 0);
    }

    #[test]
    fn test_peephole_rotation_cancellation() {
        let mut circuit = Circuit::new(1, 0);
        // RX(pi) RX(-pi) -> ID
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![PI],
        });
        circuit.add_op(Operation::Gate {
            name: GateType::RX,
            qubits: vec![0],
            params: vec![-PI],
        });

        let pass = PeepholeOptimizationPass;
        let new_circuit = pass.run(&circuit);

        assert_eq!(new_circuit.operations.len(), 0);
    }
}
