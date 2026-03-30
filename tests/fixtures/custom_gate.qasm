// Custom gate definition with parameter expression
OPENQASM 2.0;

gate my_rot(theta) q { U(theta, 0, 0) q; }

qreg q[2];
creg c[2];
my_rot(1.57) q[0];
h q[1];
cx q[0], q[1];
measure q -> c;
