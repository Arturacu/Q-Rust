// Tests the CSX (controlled sqrt-X) gate
OPENQASM 2.0;
qreg q[2];
creg c[2];
h q[0];
csx q[0], q[1];
measure q -> c;
