// Tests all newly added controlled and interaction gates
OPENQASM 2.0;
qreg q[3];
creg c[3];

// Controlled gates
cz q[0], q[1];
cy q[0], q[1];
ch q[0], q[1];

// Controlled rotations
crx(1.0) q[0], q[1];
cry(0.5) q[0], q[1];
crz(0.7) q[1], q[2];

// Ising interactions
rxx(0.3) q[0], q[1];
ryy(0.4) q[1], q[2];
rzz(0.5) q[0], q[2];

measure q -> c;
