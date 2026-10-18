"""Independent small-circuit checks using Qiskit, not Q-Rust matrices.

Compare full input/output isometries, including dirty-ancilla branches, and
terminal classical readout. Only numerical full-operator evidence is claimed.

Optional dependencies (outside the Rust library): qiskit and numpy. The first
validation used qiskit==2.2.3 and numpy==2.0.2; actual versions are saved in every
report. Install these in a separate Python environment, then run from the repo:

    python experiments/verify_benchmark_qiskit.py \
        --records experiments/results/correctness_20261016/raw_benchmark_data.json \
        --output experiments/results/correctness_20261016/qiskit_validation.json

Omit --records to check only corrected QFT/W fixtures and checker self-tests.
The default 8-qubit limit avoids large dense physical operators. The report
includes every skipped or failed item, hashes, numerical residuals, and versions.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform

try:
    import numpy as np
    import qiskit
    from qiskit import QuantumCircuit, qasm2
    from qiskit.quantum_info import Operator, Statevector
except ImportError as error:
    raise SystemExit("Optional checker dependencies are missing: install qiskit and numpy in a separate Python environment.") from error


def parse_qasm(text):
    return qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)


def split_readout(circuit):
    prefix = QuantumCircuit(circuit.num_qubits)
    prefix.global_phase = circuit.global_phase
    readout = {}
    measured = False
    for item in circuit.data:
        name = item.operation.name
        qubits = [circuit.find_bit(q).index for q in item.qubits]
        if name == "barrier":
            continue
        if name == "measure":
            measured = True
            cbit = circuit.find_bit(item.clbits[0]).index
            if cbit in readout:
                raise ValueError("repeated classical write")
            readout[cbit] = qubits[0]
        else:
            if measured or name in ("reset", "if_else") or item.clbits or getattr(item.operation, "condition", None) is not None:
                raise ValueError("unsupported nonunitary or mid-circuit operation")
            prefix.append(item.operation, qubits)
    return prefix, readout


def physical_indices(nlogical, nphysical, layout):
    if nphysical < nlogical or len(layout) < nlogical:
        raise ValueError("layout width mismatch")
    layout = layout[:nlogical]
    if len(set(layout)) != nlogical or any(q < 0 or q >= nphysical for q in layout):
        raise ValueError("invalid layout")
    basis = np.arange(1 << nlogical, dtype=np.int64)
    indices = np.zeros_like(basis)
    for logical, physical in enumerate(layout):
        indices |= ((basis >> logical) & 1) << physical
    return indices


def isometry_check(source, output, initial, final):
    source_unitary, source_readout = split_readout(source)
    output_unitary, output_readout = split_readout(output)
    nlogical, nphysical = source.num_qubits, output.num_qubits
    inputs = physical_indices(nlogical, nphysical, initial)
    outputs = physical_indices(nlogical, nphysical, final)
    expected_readout = {c: final[q] for c, q in source_readout.items()}
    if expected_readout != output_readout or (
        source_readout and source.num_clbits != output.num_clbits
    ):
        raise ValueError("terminal readout map or classical width mismatch")
    source_matrix = np.asarray(Operator(source_unitary).data)
    output_matrix = np.asarray(Operator(output_unitary).data)
    actual = output_matrix[:, inputs]
    # Full physical output retained: no projection, tracing, or renormalization.
    expected = np.zeros_like(actual)
    expected[outputs, :] = source_matrix
    dim = 1 << nlogical
    overlap = np.vdot(expected, actual)
    phase = overlap / abs(overlap) if abs(overlap) else 1.0
    difference = actual - phase * expected
    frobenius = float(np.linalg.norm(difference) / np.sqrt(dim))
    max_column = float(np.max(np.linalg.norm(difference, axis=0)))
    clean_probability = np.sum(np.abs(actual[outputs, :]) ** 2, axis=0)
    return {
        "passed": frobenius < 1e-8 and max_column < 1e-8,
        "isometry_fidelity": float(abs(overlap / dim) ** 2),
        "normalized_frobenius_residual": frobenius,
        "max_basis_column_residual": max_column,
        "max_ancilla_leakage_probability": float(max(0.0, 1.0 - np.min(clean_probability))),
        "logical_qubits": nlogical,
        "physical_qubits": nphysical,
    }


def self_checks():
    logical = QuantumCircuit(1)
    dirty = QuantumCircuit(2)
    dirty.h(1)
    separable = isometry_check(logical, dirty, [0], [0])
    assert not separable["passed"] and abs(separable["isometry_fidelity"] - 0.5) < 1e-12
    entangled = QuantumCircuit(2)
    entangled.cx(0, 1)
    entangled_result = isometry_check(logical, entangled, [0], [0])
    assert not entangled_result["passed"] and abs(entangled_result["isometry_fidelity"] - 0.25) < 1e-12
    logical.h(0)
    padded = QuantumCircuit(3)
    padded.h(2)
    assert isometry_check(logical, padded, [2], [2])["passed"]
    measured = QuantumCircuit(1, 2)
    measured.measure(0, 0)
    changed = QuantumCircuit(1, 2)
    changed.measure(0, 1)
    try:
        isometry_check(measured, changed, [0], [0])
    except ValueError:
        pass
    else:
        raise AssertionError("changed measurement destination accepted")
    return {"passed": True, "cases": ["separable dirty ancilla", "entangled ancilla", "clean padded permutation", "changed classical destination"]}


def fixture_checks(root, max_qubits):
    rows = []
    paths = sorted((root / "experiments/fixtures").glob("qft_*.qasm"))
    paths += sorted((root / "experiments/fixtures").glob("wstate_*.qasm"))
    paths += sorted((root / "tests/fixtures").glob("qft_*.qasm"))
    for path in paths:
        text = path.read_text()
        circuit, _ = split_readout(parse_qasm(text))
        row = {"path": str(path), "source_sha256": hashlib.sha256(text.encode()).hexdigest(), "qubits": circuit.num_qubits}
        if circuit.num_qubits > max_qubits:
            row.update(status="skipped_width", passed=None)
        elif path.name.startswith("qft_"):
            dim = 1 << circuit.num_qubits
            coordinates = np.arange(dim)
            expected = np.exp(2j * np.pi * np.outer(coordinates, coordinates) / dim) / np.sqrt(dim)
            actual = np.asarray(Operator(circuit).data)
            residual = float(np.max(np.abs(actual - expected)))
            row.update(status="checked", passed=residual < 1e-10,
                       max_amplitude_error=residual, oracle="analytic positive-sign discrete Fourier matrix, including output bit reversal")
        else:
            actual = np.asarray(Statevector.from_instruction(circuit).data)
            expected = np.zeros_like(actual)
            expected[[1 << q for q in range(circuit.num_qubits)]] = 1 / np.sqrt(circuit.num_qubits)
            residual = float(np.max(np.abs(actual - expected)))
            row.update(status="checked", passed=residual < 1e-10,
                       max_amplitude_error=residual, oracle="analytic W-state amplitudes on single-excitation basis states")
        rows.append(row)
    return rows


def trial_checks(root, records_path, max_qubits):
    records = json.loads(records_path.read_text())
    manifest = {f["name"]: f for f in json.loads((root / "experiments/fixtures/manifest.json").read_text())}
    rows = []
    for cell in records:
        for index, trial in enumerate(cell.get("trials", [cell])):
            row = {"fixture": trial.get("fixture"), "topology": trial.get("topology"),
                   "opt_level": trial.get("opt_level"), "trial_index": index, "runner_status": trial.get("status")}
            rows.append(row)
            if "out_qasm" not in trial:
                row.update(status="no_output", passed=None)
                continue
            try:
                source_path = root / "experiments" / manifest[trial["fixture"]]["file"]
                source_text = source_path.read_text()
                output_text = trial["out_qasm"]
                input_hash = hashlib.sha256(source_text.encode()).hexdigest()
                output_hash = hashlib.sha256(output_text.encode()).hexdigest()
                row.update(input_sha256=input_hash, output_sha256=output_hash)
                if trial.get("input_sha256", input_hash) != input_hash or trial.get("output_sha256", output_hash) != output_hash:
                    raise ValueError("recorded artifact hash mismatch")
                source, output = parse_qasm(source_text), parse_qasm(output_text)
                if max(source.num_qubits, output.num_qubits) > max_qubits:
                    row.update(status="skipped_width", passed=None)
                    continue
                initial = trial.get("initial_layout")
                final = trial.get("final_layout")
                if initial is None and final is None and source.num_qubits == output.num_qubits:
                    initial = final = list(range(source.num_qubits))
                if initial is None or final is None:
                    raise ValueError("missing or partial layout metadata")
                row.update(isometry_check(source, output, initial, final), status="checked")
            except Exception as error:
                row.update(status="error", passed=False, error=str(error))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="repository root (defaults to this script's repository)")
    parser.add_argument("--records", type=Path)
    parser.add_argument("--max-qubits", type=int, default=8,
                        help="maximum physical width for dense independent checks (default: 8)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_qubits < 0:
        parser.error("--max-qubits must be nonnegative")
    result = {"oracle": "Qiskit Operator / full physical isometry plus analytic source-fixture expectations",
              "qiskit_version": qiskit.__version__, "numpy_version": np.__version__,
              "python_version": platform.python_version(), "maximum_qubits": args.max_qubits,
              "tolerance": 1e-8, "self_checks": self_checks(),
              "fixtures": fixture_checks(args.root, args.max_qubits)}
    if args.records:
        result["records_path"] = str(args.records)
        result["records_sha256"] = hashlib.sha256(args.records.read_bytes()).hexdigest()
        result["trials"] = trial_checks(args.root, args.records, args.max_qubits)
    def summarize(rows):
        return {"passed": sum(r.get("passed") is True for r in rows),
                "failed": sum(r.get("passed") is False for r in rows),
                "skipped": sum(r.get("passed") is None for r in rows),
                "status_counts": dict(Counter(r["status"] for r in rows))}
    result["fixture_summary"] = summarize(result["fixtures"])
    if "trials" in result:
        result["trial_summary"] = summarize(result["trials"])
    result["summary"] = summarize(result["fixtures"] + result.get("trials", []))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"]))
    return 1 if result["summary"]["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
