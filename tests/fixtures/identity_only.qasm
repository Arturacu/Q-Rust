// Trivial circuit: only identity gates (optimization baseline)
OPENQASM 2.0;
qreg q[2];
creg c[2];
id q[0];
id q[1];
id q[0];
measure q -> c;
