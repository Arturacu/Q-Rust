"""Independent zero-input Clifford readout-distribution evidence.

Measurements are coherently copied into fresh zero-initialized record qubits.
This exactly defers Z measurements when no reset or classical feedback occurs,
including later quantum operations on a measured wire and classical overwrite.
Only the last record for each observed classical label is read; physical qubits
and old records are traced out. No input/output layout or compiler SDK is used.

The diagonal signed-Pauli subgroup on the observed records is extracted by GF(2)
elimination. Its parities characterize the entire joint distribution, which is
uniform over an affine subspace. Equality is stronger than matching individual
marginals, but establishes only this all-zero-input classical-readout contract,
not equality of the quantum state or arbitrary-input unitary/isometry.

Local Clifford recognition uses phase-aligned floating-point matrix residuals.
Their summed error budget is an ideal-real-arithmetic bound on approximation;
it is not a certified bound covering all floating-point arithmetic.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import re

from comparison_validation import Gate, parse_qasm, _angle
from scalable_clifford import State, multiply, recognize, TOTAL_TOLERANCE


@dataclass(frozen=True)
class Measurement:
    qubit: int
    classical_label: str


@dataclass(frozen=True)
class Program:
    width: int
    operations: tuple
    classical_labels: tuple[str, ...]

    @property
    def measured_labels(self):
        return tuple(sorted({op.classical_label for op in self.operations
                             if isinstance(op, Measurement)}))


def _unit_label(unit):
    if (not isinstance(unit, (list, tuple)) or len(unit) != 2
            or not isinstance(unit[0], str) or len(unit[1]) != 1
            or type(unit[1][0]) is not int or unit[1][0] < 0):
        raise ValueError('Unsupported unit identifier')
    return f'{unit[0]}[{unit[1][0]}]'


def parse_program_qasm(text):
    """Parse static gates plus measurement/barrier; unknown effects fail closed."""
    text = re.sub(r'//[^\n]*', '', text)
    statements = [s.strip() for s in text.split(';') if s.strip()]
    qregs, cregs, header, ops = {}, {}, [], []
    width = 0
    started = False
    for s in statements:
        if re.fullmatch(r'OPENQASM\s+2\.0', s) or re.fullmatch(r'include\s+"qelib1\.inc"', s):
            if started:
                raise ValueError('Late QASM header')
            header.append(s + ';')
            continue
        declaration = re.fullmatch(r'(qreg|creg)\s+(\w+)\[(\d+)\]', s)
        if declaration:
            if started:
                raise ValueError('Late register declaration')
            kind, name, count = declaration.group(1), declaration.group(2), int(declaration.group(3))
            if count < 1 or name in qregs or name in cregs:
                raise ValueError('Invalid register declaration')
            if kind == 'qreg':
                qregs[name] = (width, count)
                width += count
            else:
                cregs[name] = count
            header.append(s + ';')
            continue
        started = True
        if s.startswith('barrier '):
            for operand in s[len('barrier '):].split(','):
                operand = operand.strip()
                if operand in qregs:
                    continue
                indexed = re.fullmatch(r'(\w+)\[(\d+)\]', operand)
                if (indexed is None or indexed.group(1) not in qregs
                        or int(indexed.group(2)) >= qregs[indexed.group(1)][1]):
                    raise ValueError('Invalid barrier operand')
            continue
        measured = re.fullmatch(r'measure\s+(\w+)\[(\d+)\]\s*->\s*(\w+)\[(\d+)\]', s)
        if measured:
            qname, qi, cname, ci = measured.group(1), int(measured.group(2)), measured.group(3), int(measured.group(4))
            if qname not in qregs or cname not in cregs or qi >= qregs[qname][1] or ci >= cregs[cname]:
                raise ValueError('Out-of-range readout')
            ops.append(Measurement(qregs[qname][0] + qi, f'{cname}[{ci}]'))
        else:
            parsed_width, gates = parse_qasm('\n'.join(header) + '\n' + s + ';')
            if parsed_width != width or len(gates) != 1:
                raise ValueError('Unsupported operation statement')
            ops.append(gates[0])
    if width < 1:
        raise ValueError('No quantum register')
    return Program(width, tuple(ops), tuple(f'{name}[{i}]' for name, count in cregs.items() for i in range(count)))


def parse_tket_native(data):
    """Interpret retained Circuit.to_dict, using actual command wire identities.

    Final quantum-only implicit permutations are irrelevant because quantum
    outputs are traced out. Measurements and their order remain explicit.
    """
    units = data.get('qubits')
    if not isinstance(units, list) or not units:
        raise ValueError('Missing native quantum wires')
    labels = [_unit_label(unit) for unit in units]
    if len(set(labels)) != len(labels):
        raise ValueError('Repeated quantum wire')
    wire = {label: i for i, label in enumerate(labels)}
    bits = tuple(_unit_label(unit) for unit in data.get('bits', []))
    if len(set(bits)) != len(bits):
        raise ValueError('Repeated classical bit')
    ops = []
    for command in data.get('commands', []):
        op = command['op']
        kind = op['type']
        args = [_unit_label(unit) for unit in command['args']]
        if kind in {'Barrier', 'Phase'}:
            continue
        if kind == 'Measure':
            if len(args) != 2 or args[0] not in wire or args[1] not in bits:
                raise ValueError('Malformed native measurement')
            ops.append(Measurement(wire[args[0]], args[1]))
            continue
        if any(arg not in wire for arg in args):
            raise ValueError('Classical feedback or unsupported arguments')
        name = {'noop': 'id'}.get(kind.lower(), kind.lower())
        params = tuple(_angle(str(value)) * math.pi for value in op.get('params', []))
        # Reuse the strict, independent QASM signature/arity checker.
        statement = name + ('(' + ','.join(repr(x) for x in params) + ')' if params else '')
        statement += ' ' + ','.join(f'q[{wire[arg]}]' for arg in args) + ';'
        _, gates = parse_qasm(f'OPENQASM 2.0; include "qelib1.inc"; qreg q[{len(wire)}];' + statement)
        if len(gates) != 1:
            raise ValueError('Malformed native gate')
        ops.append(gates[0])
    return Program(len(wire), tuple(ops), bits)


def _canonical_parities(rows, count):
    rows = list(rows)
    pivot = 0
    for column in range(count):
        bit = 1 << column
        found = next((i for i in range(pivot, len(rows)) if rows[i][0] & bit), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        for i in range(len(rows)):
            if i != pivot and rows[i][0] & bit:
                rows[i] = (rows[i][0] ^ rows[pivot][0], rows[i][1] ^ rows[pivot][1])
        pivot += 1
    for mask, sign in rows[pivot:]:
        if mask or sign:
            raise ValueError('Inconsistent measured parity subgroup')
    return tuple(rows[:pivot])


def diagonal_parities(rows, width, observed_wires):
    """Kernel of X-all and Z-unobserved projections, retaining Pauli signs."""
    if len(set(observed_wires)) != len(observed_wires):
        raise ValueError('Repeated observed record wire')
    observed_mask = sum(1 << q for q in observed_wires)
    outside = ((1 << width) - 1) ^ observed_mask
    pivots, diagonal = {}, []
    for original in rows:
        row = original
        bad = row[0] | ((row[1] & outside) << width)
        while bad:
            pivot = bad.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (bad, row)
                break
            old_bad, old_row = pivots[pivot]
            row = multiply(row, old_row)
            bad ^= old_bad
        if not bad:
            if row[0] or row[1] & outside:
                raise ValueError('Projection kernel error')
            mask = sum(((row[1] >> q) & 1) << i for i, q in enumerate(observed_wires))
            diagonal.append((mask, row[2]))
    return _canonical_parities(diagonal, len(observed_wires))


def readout_signature(program, classical_order):
    if (not program.measured_labels or len(set(classical_order)) != len(classical_order)
            or set(classical_order) != set(program.measured_labels)):
        raise ValueError('Observed classical-label sets differ or are empty')
    measurements = sum(isinstance(op, Measurement) for op in program.operations)
    state = State(program.width + measurements)
    last, record, errors = {}, program.width, []
    measured_physical, later_gate_touches = set(), 0
    for op in program.operations:
        if isinstance(op, Measurement):
            state.apply('cx', op.qubit, record)
            last[op.classical_label] = record
            record += 1
            measured_physical.add(op.qubit)
        elif len(op.qubits) == 1:
            later_gate_touches += sum(q in measured_physical for q in op.qubits)
            word, residual = recognize(op.name, op.params)
            for symbol in word:
                state.apply(symbol, *op.qubits)
            errors.append(residual)
        elif op.name in {'cx', 'cy', 'cz', 'swap'} and not op.params:
            later_gate_touches += sum(q in measured_physical for q in op.qubits)
            state.apply(op.name, *op.qubits)
        else:
            raise ValueError('Non-Clifford, reset, or unsupported operation: ' + op.name)
    parities = diagonal_parities(state.rows, state.width, [last[label] for label in classical_order])
    return parities, {
        'physical_qubits_in_artifact': program.width,
        'coherent_record_qubits': measurements,
        'observed_classical_bits': len(classical_order),
        'unused_classical_declarations': sorted(set(program.classical_labels) - set(classical_order)),
        'later_quantum_gate_wire_touches_after_measurement': later_gate_touches,
        'parity_rank': len(parities),
        'uniform_support_log2_size': len(classical_order) - len(parities),
        'recognized_one_qubit_gates': len(errors),
        'max_local_recognition_residual': max(errors, default=0),
        'recognition_residual_sum': math.fsum(errors),
        'parities': [{'mask': mask, 'value': sign} for mask, sign in parities],
    }


def compare_programs(source, output):
    result = {'passed': False, 'evidence_class': 'unverified_readout_distribution',
              'contract': 'all quantum inputs zero; joint distribution of written classical labels; unmeasured quantum outputs and unused classical storage excluded',
              'arbitrary_input_equivalence': False, 'quantum_state_equivalence': False,
              'recognition_is_numerical': True, 'floating_point_bound_certified': False,
              'measurement_handling': 'coherent copy to distinct record qubits in original operation order; final signed-Z subgroup on last written records',
              'numerical_recognition_limit': TOTAL_TOLERANCE}
    try:
        labels = source.measured_labels
        expected, source_info = readout_signature(source, labels)
        actual, output_info = readout_signature(output, labels)
        error = source_info['recognition_residual_sum'] + output_info['recognition_residual_sum']
        result.update(classical_order=list(labels), source=source_info, output=output_info,
                      parity_subgroups_equal=expected == actual, recognition_residual_sum=error)
        result['passed'] = expected == actual and error <= TOTAL_TOLERANCE
        if result['passed']:
            result['evidence_class'] = 'zero_input_clifford_joint_readout_distribution'
        else:
            result['error'] = 'Readout parity subgroups differ' if expected != actual else 'Recognition residual budget exceeded'
    except (ValueError, TypeError, IndexError, SyntaxError, OverflowError, ZeroDivisionError, KeyError) as error:
        result['error'] = str(error)
    return result


def validate_native_artifact(artifact):
    try:
        source = parse_program_qasm(artifact['source_qasm'])
        if artifact['tool'] == 'qiskit':
            output = parse_program_qasm(artifact['native_qasm'])
            binding = 'retained native_qasm including operation-order measurements'
        elif artifact['tool'] == 'tket':
            output = parse_tket_native(artifact['native_dict'])
            binding = 'retained native_dict including operation-order measurements; replay maps unused'
        else:
            raise ValueError('Unsupported native artifact tool')
        result = compare_programs(source, output)
        result['artifact_interpretation'] = binding
        result['source_sha256'] = hashlib.sha256(artifact['source_qasm'].encode()).hexdigest()
        result['native_representation_sha256'] = hashlib.sha256(json.dumps(artifact.get('native_dict', artifact.get('native_qasm')), sort_keys=True).encode()).hexdigest()
        return result
    except (ValueError, TypeError, IndexError, SyntaxError, OverflowError, ZeroDivisionError, KeyError) as error:
        return {'passed': False, 'evidence_class': 'unverified_readout_distribution', 'error': str(error)}
