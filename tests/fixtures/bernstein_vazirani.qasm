OPENQASM 2.0;
include "qelib1.inc";
// Bernstein-Vazirani: hidden bitstring s = "101" (3 bits)
// Query qubits: q[0], q[1], q[2]; ancilla: q[3]
// Oracle: CX for each set bit in s
qreg q[4];
x q[3];
h q[0];
h q[1];
h q[2];
h q[3];
// Oracle for s=101: CX(q0,ancilla) and CX(q2,ancilla)
cx q[0], q[3];
cx q[2], q[3];
h q[0];
h q[1];
h q[2];
// Measuring q[0..2] should yield 1,0,1 (s="101")
