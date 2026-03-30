// Parametric rotation gates with various angles
OPENQASM 2.0;
qreg q[2];
creg c[2];
rx(pi/4) q[0];
ry(pi/3) q[1];
rz(pi/2) q[0];
rx(1.234) q[1];
ry(0.5) q[0];
rz(2.0 * pi / 3.0) q[1];
measure q -> c;
