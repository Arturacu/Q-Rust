OPENQASM 2.0;
include "qelib1.inc";
// Deutsch-Jozsa: 2 query qubits + 1 ancilla
// Hidden function is balanced: f(x) = x[0] XOR x[1]
// Expected: measurement on q[0],q[1] gives |11> (balanced indicator)
qreg q[3];
// Prepare ancilla in |-> = H|1>
x q[2];
h q[0];
h q[1];
h q[2];
// Balanced oracle: CX(q0, ancilla) and CX(q1, ancilla)
cx q[0], q[2];
cx q[1], q[2];
// Hadamard on query registers
h q[0];
h q[1];
