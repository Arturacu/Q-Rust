#!/usr/bin/env python3
"""Deterministic OpenQASM 2.0 fixture generator for the Q-Rust benchmark suite.

Produces a standard suite of algorithmic benchmark circuits spanning a range of widths and depths:

  * GHZ states                 n = 3, 5, 7, 10, 14, 18, 22, 25
  * Quantum Fourier Transform  n = 3, 5, 7, 10, 14, 18, 22
  * Cuccaro ripple-carry adder n = 4, 8, 16   (total qubit count)
  * Random Clifford circuits   n = 5, 10, 15  (depth 20)

Design constraints (verified against Q-Rust's parser in `src/parser/`):
  * STRICT OpenQASM 2.0 only — never 3.0 syntax.
  * Every emitted gate name is one Q-Rust parses natively
    (`h, x, s, sdg, t, cx, ccx, swap, crz, rz, ...`).  In particular the QFT
    uses controlled-RZ (`crz`) for its phase rotations, NOT `cp`/`cu1`, which
    Q-Rust would treat as unsupported custom gates.
  * Output is fully deterministic (fixed RNG seeds) so the fixture suite — and
    therefore every downstream benchmark number — is reproducible.

The script also writes `fixtures/manifest.json`, recording for each fixture its
family, nominal `n`, true qubit count, and source-level gate tallies.  This
manifest is the single source of truth the orchestrator iterates over.

Usage:
    python3 experiments/generate_fixtures.py
"""
from __future__ import annotations

import json
import math
import os
import random
from dataclasses import dataclass, asdict, field
from typing import List

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(HERE, "fixtures")

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


@dataclass
class Fixture:
    """One generated benchmark circuit and its source-level metadata."""

    name: str
    family: str          # "ghz" | "qft" | "adder" | "clifford"
    nominal_n: int       # the `n` label used in reports and tables
    num_qubits: int      # actual qubit register width
    source_gates: int = 0
    source_2q_gates: int = 0
    source_cx_gates: int = 0
    notes: str = ""
    body: List[str] = field(default_factory=list)

    def qasm(self) -> str:
        return HEADER + f"qreg q[{self.num_qubits}];\n" + "".join(self.body)


def _tally(fx: Fixture, line: str, n_qubits_in_gate: int, is_cx: bool) -> None:
    """Record a gate into the fixture body and update source-level counters."""
    fx.body.append(line)
    fx.source_gates += 1
    if n_qubits_in_gate >= 2:
        fx.source_2q_gates += 1
    if is_cx:
        fx.source_cx_gates += 1


# --------------------------------------------------------------------------- #
# GHZ                                                                          #
# --------------------------------------------------------------------------- #
def make_ghz(n: int) -> Fixture:
    """|GHZ_n> = (|0..0> + |1..1>)/sqrt(2): one H then a CX ladder."""
    fx = Fixture(name=f"ghz_{n}", family="ghz", nominal_n=n, num_qubits=n,
                 notes="H on q[0] followed by a linear CX ladder.")
    _tally(fx, "h q[0];\n", 1, False)
    for i in range(n - 1):
        _tally(fx, f"cx q[{i}],q[{i + 1}];\n", 2, True)
    return fx


# --------------------------------------------------------------------------- #
# QFT (textbook, controlled-RZ phase rotations + final bit-reversal swaps)     #
# --------------------------------------------------------------------------- #
def make_qft(n: int) -> Fixture:
    """n-qubit QFT.

    Each qubit gets a Hadamard followed by a descending chain of controlled
    phase rotations implemented with `crz(pi / 2**k)` (Q-Rust-native), and the
    routine closes with the standard bit-reversal SWAP network.  This is a
    legitimate QFT up to single-qubit phases, which is irrelevant for the
    transpiler-equivalence harness (it compares input vs. output, both built
    from the same primitives).
    """
    fx = Fixture(name=f"qft_{n}", family="qft", nominal_n=n, num_qubits=n,
                 notes="Hadamard + crz(pi/2^k) phase ladder, then bit-reversal swaps.")
    for j in range(n):
        _tally(fx, f"h q[{j}];\n", 1, False)
        for k in range(j + 1, n):
            angle = math.pi / (2 ** (k - j))
            _tally(fx, f"crz({angle:.12g}) q[{k}],q[{j}];\n", 2, False)
    for i in range(n // 2):
        _tally(fx, f"swap q[{i}],q[{n - 1 - i}];\n", 2, False)
    return fx


# --------------------------------------------------------------------------- #
# Cuccaro ripple-carry adder (quant-ph/0410184)                               #
# --------------------------------------------------------------------------- #
def make_adder(total_qubits: int) -> Fixture:
    """Cuccaro ripple-carry adder sized to occupy exactly `total_qubits`.

    Layout (2w + 2 qubits for w-bit operands): c0 ancilla, interleaved
    a_i/b_i operand bits, and a high carry-out z.  We solve w = (n-2)/2.
    Adds register a into register b in place, with the MAJ / UMA primitives.
    A few X gates seed a non-trivial addend so the circuit is not the identity.
    """
    assert total_qubits >= 4 and total_qubits % 2 == 0, "adder needs even n >= 4"
    w = (total_qubits - 2) // 2
    fx = Fixture(name=f"adder_{total_qubits}", family="adder",
                 nominal_n=total_qubits, num_qubits=total_qubits,
                 notes=f"Cuccaro ripple-carry adder, {w}-bit operands "
                       f"({total_qubits} qubits = 2*{w}+2).")

    c0 = 0
    a = [1 + 2 * i for i in range(w)]
    b = [2 + 2 * i for i in range(w)]
    z = total_qubits - 1

    def cx(ctrl, tgt):
        _tally(fx, f"cx q[{ctrl}],q[{tgt}];\n", 2, True)

    def ccx(c1, c2, tgt):
        _tally(fx, f"ccx q[{c1}],q[{c2}],q[{tgt}];\n", 3, False)

    def x(t):
        _tally(fx, f"x q[{t}];\n", 1, False)

    def maj(cq, bq, aq):
        cx(aq, bq); cx(aq, cq); ccx(cq, bq, aq)

    def uma(cq, bq, aq):
        ccx(cq, bq, aq); cx(aq, cq); cx(cq, bq)

    # Seed an example addend (a = 0b...01, b = 0b...01) so the adder does work.
    x(a[0])
    x(b[0])

    carry = c0
    for i in range(w):
        maj(carry, b[i], a[i])
        carry = a[i]
    cx(a[w - 1], z)
    carry = a[w - 2] if w >= 2 else c0
    for i in reversed(range(w)):
        cq = c0 if i == 0 else a[i - 1]
        uma(cq, b[i], a[i])
    return fx


# --------------------------------------------------------------------------- #
# Random Clifford circuits                                                     #
# --------------------------------------------------------------------------- #
_CLIFFORD_1Q = ["h", "s", "sdg", "x", "y", "z"]


def make_clifford(n: int, depth: int = 20, seed: int = 0, tag: str = "") -> Fixture:
    """Random Clifford circuit of `n` qubits and `depth` layers.

    Each layer: a random single-qubit Clifford on every qubit, then a random
    perfect-ish matching of CX gates over a shuffled qubit list.  Seeded per-n
    for bit-for-bit reproducibility.  `tag` distinguishes independent instances
    at the same width (a second random draw).
    """
    rng = random.Random(0xC11FF0 + n * 1000 + seed)
    fx = Fixture(name=f"clifford_{n}{tag}", family="clifford", nominal_n=n,
                 num_qubits=n,
                 notes=f"Random Clifford, depth {depth}, seed {seed} "
                       f"(reproducible).")
    for _ in range(depth):
        for q in range(n):
            g = rng.choice(_CLIFFORD_1Q)
            _tally(fx, f"{g} q[{q}];\n", 1, False)
        order = list(range(n))
        rng.shuffle(order)
        for i in range(0, n - 1, 2):
            ctrl, tgt = order[i], order[i + 1]
            _tally(fx, f"cx q[{ctrl}],q[{tgt}];\n", 2, True)
    return fx


# --------------------------------------------------------------------------- #
# QAOA (MaxCut on a ring, p layers): rzz cost + rx mixer                       #
# --------------------------------------------------------------------------- #
def make_qaoa(n: int, p: int = 2) -> Fixture:
    fx = Fixture(name=f"qaoa_{n}", family="qaoa", nominal_n=n, num_qubits=n,
                 notes=f"MaxCut QAOA on a ring, p={p}: rzz(gamma) cost + rx(beta) mixer.")
    for q in range(n):
        _tally(fx, f"h q[{q}];\n", 1, False)
    for layer in range(p):
        gamma = 0.6 + 0.2 * layer
        beta = 0.4 - 0.1 * layer
        for i in range(n):
            j = (i + 1) % n
            _tally(fx, f"rzz({gamma:.6g}) q[{i}],q[{j}];\n", 2, False)
        for q in range(n):
            _tally(fx, f"rx({beta:.6g}) q[{q}];\n", 1, False)
    return fx


# --------------------------------------------------------------------------- #
# Bernstein-Vazirani (secret string = all ones): width = s+1                   #
# --------------------------------------------------------------------------- #
def make_bv(total_qubits: int) -> Fixture:
    s = total_qubits - 1  # number of input qubits; last qubit is the ancilla
    anc = s
    fx = Fixture(name=f"bv_{total_qubits}", family="bv", nominal_n=total_qubits,
                 num_qubits=total_qubits,
                 notes=f"Bernstein-Vazirani, {s}-bit secret (all ones), 1 ancilla.")
    _tally(fx, f"x q[{anc}];\n", 1, False)
    _tally(fx, f"h q[{anc}];\n", 1, False)
    for q in range(s):
        _tally(fx, f"h q[{q}];\n", 1, False)
    for q in range(s):
        _tally(fx, f"cx q[{q}],q[{anc}];\n", 2, True)
    for q in range(s):
        _tally(fx, f"h q[{q}];\n", 1, False)
    return fx


# --------------------------------------------------------------------------- #
# Ising / Trotter evolution on a linear chain                                  #
# --------------------------------------------------------------------------- #
def make_ising(n: int, steps: int = 3) -> Fixture:
    fx = Fixture(name=f"ising_{n}", family="ising", nominal_n=n, num_qubits=n,
                 notes=f"Transverse-field Ising Trotterization, {steps} steps "
                       f"(rzz couplings + rx field).")
    dt, j_coupling, h_field = 0.2, 1.0, 0.8
    for _ in range(steps):
        for i in range(n - 1):
            _tally(fx, f"rzz({2 * j_coupling * dt:.6g}) q[{i}],q[{i + 1}];\n", 2, False)
        for q in range(n):
            _tally(fx, f"rx({2 * h_field * dt:.6g}) q[{q}];\n", 1, False)
    return fx


# --------------------------------------------------------------------------- #
# W-state preparation (ry/cx cascade)                                          #
# --------------------------------------------------------------------------- #
def make_wstate(n: int) -> Fixture:
    fx = Fixture(name=f"wstate_{n}", family="wstate", nominal_n=n, num_qubits=n,
                 notes="W-state preparation via a ry/cx cascade.")
    _tally(fx, f"ry({2 * math.acos(math.sqrt(1 / n)):.9g}) q[0];\n", 1, False)
    for k in range(1, n):
        rem = n - k
        theta = 2 * math.acos(math.sqrt(1 / rem)) if rem > 0 else 0.0
        _tally(fx, f"ry({theta / 2:.9g}) q[{k}];\n", 1, False)
        _tally(fx, f"cx q[{k - 1}],q[{k}];\n", 2, True)
        _tally(fx, f"ry({-theta / 2:.9g}) q[{k}];\n", 1, False)
        _tally(fx, f"cx q[{k - 1}],q[{k}];\n", 2, True)
    return fx


# --------------------------------------------------------------------------- #
# Graph / cluster state on a ring (H + CZ mesh)                                #
# --------------------------------------------------------------------------- #
def make_graphstate(n: int) -> Fixture:
    fx = Fixture(name=f"graphstate_{n}", family="graphstate", nominal_n=n,
                 num_qubits=n, notes="Ring cluster state: H on all + CZ on ring edges.")
    for q in range(n):
        _tally(fx, f"h q[{q}];\n", 1, False)
    for i in range(n):
        j = (i + 1) % n
        if i < j:  # avoid duplicate edge on the wrap for n=2
            _tally(fx, f"cz q[{i}],q[{j}];\n", 2, False)
    if n > 2:  # closing ring edge (n-1, 0)
        _tally(fx, f"cz q[0],q[{n - 1}];\n", 2, False)
    return fx


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def build_all() -> List[Fixture]:
    # Suite design: the widths <= 18 are the *verifiable core* (every one is
    # certified by the harness on every topology it fits); the 22- and
    # 25-qubit GHZ/QFT fixtures are retained solely as compile-time *scaling*
    # probes and are flagged wherever they appear. The extra sub-14-qubit
    # widths and second Clifford seeds densify the verified core.
    fixtures: List[Fixture] = []
    for n in (3, 5, 7, 8, 10, 12, 14, 18, 22, 25):
        fixtures.append(make_ghz(n))
    for n in (3, 5, 7, 8, 10, 12, 14, 18, 22):
        fixtures.append(make_qft(n))
    for n in (4, 6, 8, 10, 12, 16):
        fixtures.append(make_adder(n))
    for n in (5, 7, 10, 12, 15):
        fixtures.append(make_clifford(n, depth=20))
    # Second independent random-Clifford instances at two widths.
    fixtures.append(make_clifford(5, depth=20, seed=1, tag="b"))
    fixtures.append(make_clifford(10, depth=20, seed=1, tag="b"))
    # Extended families for the cross-tool comparison breadth study.
    for n in (5, 7, 10, 12, 14):
        fixtures.append(make_qaoa(n))
    for n in (5, 7, 8, 10, 12, 14):
        fixtures.append(make_bv(n))
    for n in (5, 7, 10, 12, 14):
        fixtures.append(make_ising(n))
    for n in (5, 7, 10, 12, 14):
        fixtures.append(make_wstate(n))
    for n in (5, 7, 10, 12, 14):
        fixtures.append(make_graphstate(n))
    return fixtures


def main() -> None:
    os.makedirs(FIXTURE_DIR, exist_ok=True)
    fixtures = build_all()
    manifest = []
    for fx in fixtures:
        path = os.path.join(FIXTURE_DIR, f"{fx.name}.qasm")
        with open(path, "w") as fh:
            fh.write(fx.qasm())
        entry = asdict(fx)
        entry.pop("body")
        entry["file"] = f"fixtures/{fx.name}.qasm"
        manifest.append(entry)
        print(f"  wrote {fx.name:<12} n={fx.num_qubits:<3} "
              f"gates={fx.source_gates:<5} 2q={fx.source_2q_gates:<4} "
              f"cx={fx.source_cx_gates}")

    manifest_path = os.path.join(FIXTURE_DIR, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\n{len(fixtures)} fixtures written to {FIXTURE_DIR}")
    print(f"manifest -> {manifest_path}")


if __name__ == "__main__":
    main()
