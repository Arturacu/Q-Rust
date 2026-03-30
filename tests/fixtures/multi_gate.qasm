// Mixed standard gates covering all single-qubit gate types
OPENQASM 2.0;
qreg q[3];
creg c[3];
h q[0];
x q[1];
y q[2];
z q[0];
s q[1];
sdg q[2];
t q[0];
tdg q[1];
id q[2];
swap q[0], q[1];
cx q[1], q[2];
h q[2];
measure q -> c;
