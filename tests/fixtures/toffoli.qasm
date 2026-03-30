// Toffoli (CCX) gate circuit
OPENQASM 2.0;
qreg q[3];
creg c[3];
x q[0];
x q[1];
ccx q[0], q[1], q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
