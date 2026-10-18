OPENQASM 2.0;
include "qelib1.inc";
// Quantum Phase Estimation (toy 2-qubit version)
// Estimates eigenphase of U=T (T gate, phase = pi/4)
// Ancilla q[0] + eigenstate register q[1]
qreg q[2];
// Prepare eigenstate of T: |1> (T|1> = e^{i*pi/4}|1>)
x q[1];
// QPE step: H + controlled-T
h q[0];
// Controlled-T applied to q[1] conditioned on q[0]
crz(0.7853981634) q[0], q[1];
// Inverse QFT on 1 qubit = just Hadamard
h q[0];
