#!/usr/bin/env python3
"""Freeze community-sourced unitary-prefix circuits before any compiler timing.

Acquisition is a separate step: this program reads a pinned local source mirror.
It never optimizes circuits, routes them, or observes compiler output quality.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import sys

QASMBENCH_REVISION = '357b942396d5c2b7cbc1c229c585a6ef5ccaebac'
MQTBENCH_REVISION = '2cd3f3b08d5d8866595948fbce955e701e86a5df'
BENCHPRESS_REVISION = 'b695f30e83a32bac05b9b4d8e98d37ba9aae5236'
MQT_FAMILIES = ('ae', 'dj', 'grover', 'qnn', 'qpeexact', 'qpeinexact',
                'qwalk', 'vqe_real_amp', 'vqe_su2')
MQT_WIDTHS = (5, 10)
MAX_QUBITS = 14
MAX_RAW_BYTES = 250000
MAX_EXPANDED_GATES = 5000
SEED = 0x5152555354


def digest(data):
    return hashlib.sha256(data.encode() if isinstance(data, str) else data).hexdigest()


def git_blob_digest(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def qasmbench_inventory(tree):
    """Inventory *all* small/medium QASM, including declared scope exclusions."""
    rows = []
    for item in tree['tree']:
        path = item['path']
        match = re.fullmatch(r'(small|medium)/([^/]+)_n(\d+)/([^/]+\.qasm)', path)
        if not match:
            continue
        reason = None
        if '_transpiled' in match[4]:
            reason = 'upstream_compiler_transpiled_variant'
        elif int(match[3]) > MAX_QUBITS:
            reason = 'width_above_14'
        elif item['size'] > MAX_RAW_BYTES:
            reason = 'raw_source_above_250000_bytes'
        rows.append({'community_source': 'QASMBench', 'source_path': path,
                     'source_revision': QASMBENCH_REVISION,
                     'family': match[2], 'directory_width': int(match[3]),
                     'raw_source_bytes': item['size'], 'upstream_blob_sha1': item['sha'],
                     'selection_status': 'candidate' if reason is None else 'excluded',
                     'exclusion_reason': reason})
    return sorted(rows, key=lambda row: row['source_path'])


def unitary_prefix(circuit):
    """Only a global terminal measurement suffix may be removed; fail closed."""
    from qiskit import QuantumCircuit
    out = QuantumCircuit(circuit.num_qubits)
    out.global_phase = circuit.global_phase
    terminal = False
    removed = {'terminal_measurements': 0, 'barriers': 0}
    for instruction in circuit.data:
        operation = instruction.operation
        name = operation.name
        if name == 'barrier':
            removed['barriers'] += 1
            continue
        if name == 'measure':
            terminal = True
            removed['terminal_measurements'] += 1
            continue
        if name == 'reset' or instruction.clbits or getattr(operation, 'condition', None) is not None:
            raise ValueError('nonunitary_or_classically_controlled_operation:' + name)
        if hasattr(operation, 'blocks'):
            raise ValueError('dynamic_control_flow:' + name)
        if terminal:
            raise ValueError('measurement_precedes_later_quantum_operation')
        out.append(operation, [circuit.find_bit(q).index for q in instruction.qubits])
    return out, removed


def _statement(name, params, qubits):
    parameters = '(' + ','.join(repr(float(p)) for p in params) + ')' if params else ''
    return name + parameters + ' ' + ','.join(f'q[{q}]' for q in qubits) + ';'


def _phase_residual(expected, actual):
    import numpy as np
    overlap = np.vdot(expected, actual)
    phase = overlap / abs(overlap) if abs(overlap) > 0 else 1
    return float(np.max(np.linalg.norm(actual - phase * expected, axis=0)))


def canonical_source(circuit, max_gates=MAX_EXPANDED_GATES):
    """Expand definitions only; retain supported primitives after matrix check.

    No transpile(), pass manager, layout, or compiler optimization is invoked.
    Qiskit may synthesize operation definitions internally. Definitions are a trust
    boundary. For named primitives, a <=3-qubit matrix check prevents accidental
    substitution of a source-defined gate with a different standard meaning.
    """
    import numpy as np
    from qiskit.quantum_info import Operator
    from qiskit.circuit.library import PermutationGate, SwapGate
    from comparison_validation import _SIGNATURES, Gate, simulate
    lines = []
    expanded = Counter()
    aliases = {'u': 'u3', 'p': 'u1', 'cp': 'cu1'}
    cache = {}

    def emit(operation, qubits, depth=0):
        if depth > 60:
            raise ValueError('gate_definition_recursion_limit')
        name = aliases.get(operation.name, operation.name)
        signature = _SIGNATURES.get(name)
        params = tuple(float(p) for p in operation.params) if signature is not None else ()
        if not all(np.isfinite(p) for p in params):
            raise ValueError('nonfinite_parameter')
        retain = False
        if signature == (len(params), len(qubits)):
            # Operation identity is included so distinct source definitions with
            # the same name and parameters cannot accidentally share a result.
            key = (id(operation), name, params)
            cached = cache.get(key)
            retain = cached[1] if cached is not None and cached[0] is operation else None
            if retain is None:
                identity = np.eye(1 << len(qubits), dtype=complex)
                analytic = simulate(len(qubits), (Gate(name, params, tuple(range(len(qubits)))),), identity)
                actual = Operator(operation).data
                retain = _phase_residual(actual, analytic) <= 1e-12
                cache[key] = (operation, retain)
        if retain:
            lines.append(_statement(name, params, qubits))
            if len(lines) > max_gates:
                raise ValueError('expanded_source_above_5000_gates')
            return
        if isinstance(operation, PermutationGate):
            pattern = [int(value) for value in operation.pattern]
            if sorted(pattern) != list(range(len(qubits))):
                raise ValueError('invalid_permutation')
            current = list(range(len(qubits)))
            expanded['permutation_to_swaps'] += 1
            for target, source in enumerate(pattern):
                position = current.index(source)
                if target != position:
                    emit(SwapGate(), [qubits[target], qubits[position]], depth + 1)
                    current[target], current[position] = current[position], current[target]
            return
        definition = operation.definition
        if definition is None:
            raise ValueError('unsupported_gate_definition:' + operation.name)
        expanded[operation.name] += 1
        for child in definition.data:
            if child.operation.name == 'barrier':
                continue
            if child.clbits or child.operation.name in {'measure', 'reset'} or hasattr(child.operation, 'blocks'):
                raise ValueError('nonunitary_nested_definition:' + operation.name)
            emit(child.operation, [qubits[definition.find_bit(q).index] for q in child.qubits], depth + 1)

    for instruction in circuit.data:
        emit(instruction.operation, [circuit.find_bit(q).index for q in instruction.qubits])
    text = 'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[' + str(circuit.num_qubits) + '];\n' + '\n'.join(lines) + '\n'
    return text, {'definition_expansions': dict(sorted(expanded.items())),
                  'expanded_source_gate_count': len(lines),
                  'global_phase': 'discarded only as a whole-circuit equivalence convention'}


def normalization_evidence(original, source_qasm, normalized_qasm):
    """Check original Qiskit circuit against independent analytic NumPy gates."""
    import numpy as np
    from qiskit.quantum_info import Statevector
    from comparison_validation import parse_qasm, simulate, TOLERANCE
    width = original.num_qubits
    if width <= 6:
        inputs = np.eye(1 << width, dtype=complex)
        mode = 'full_unitary_numerical'
    else:
        rng = np.random.default_rng(SEED)
        inputs = np.zeros((1 << width, 8), dtype=complex)
        for index, basis in enumerate((0, (1 << width) - 1, int('10' * (width // 2) + ('1' if width % 2 else ''), 2))):
            inputs[basis, index] = 1
        inputs[:, 3] = 1 / np.sqrt(1 << width)
        for index in range(4, 8):
            values = rng.normal(size=1 << width) + 1j * rng.normal(size=1 << width)
            inputs[:, index] = values / np.linalg.norm(values)
        mode = 'sampled_8_inputs_numerical'
    expected = np.column_stack([Statevector(inputs[:, index]).evolve(original).data for index in range(inputs.shape[1])])
    residuals = {}
    for name, qasm in [('source', source_qasm), ('normalized', normalized_qasm)]:
        actual_width, gates = parse_qasm(qasm)
        if actual_width != width:
            raise ValueError('normalization_changed_width')
        residuals[name] = _phase_residual(expected, simulate(width, gates, inputs))
    passed = all(value <= TOLERANCE for value in residuals.values())
    return {'passed': passed, 'mode': mode, 'input_count': inputs.shape[1],
            'seed': SEED, 'tolerance': TOLERANCE, 'residuals': residuals,
            'comparison': 'Qiskit Statevector on retained original unitary prefix versus independent analytic NumPy interpreter; one shared global phase across inputs'}


def _record(circuit, entry, output, raw_path, raw_bytes):
    from bench_baselines import normalize_qasm
    entry.update(raw_source_sha256=digest(raw_bytes), raw_artifact=str(raw_path))
    if circuit.num_qubits != entry['directory_width']:
        raise ValueError('actual_width_disagrees_with_declared_width')
    prefix, removals = unitary_prefix(circuit)
    source, expansion = canonical_source(prefix)
    normalized = normalize_qasm(source)
    evidence = normalization_evidence(prefix, source, normalized)
    entry['normalization_evidence'] = evidence
    if not evidence['passed']:
        raise RuntimeError('normalization_equivalence_failed')
    name = entry['id']
    (output / (name + '.source.qasm')).write_text(source)
    (output / (name + '.normalized.qasm')).write_text(normalized)
    entry.update(fixture=name, n_logical=circuit.num_qubits, source_qasm=source, qasm=normalized,
                 source_sha256=digest(source), normalized_sha256=digest(normalized),
                 source_artifact=name + '.source.qasm', normalized_artifact=name + '.normalized.qasm',
                 normalization={'track': 'static unitary prefix; arbitrary logical inputs',
                                'removed': removals, **expansion,
                                'shared_gate_expansion': 'existing exact CRZ/RZZ to RZ/CX expansion',
                                'optimization': 'none; no placement or routing'},
                 selection_status='selected', exclusion_reason=None)
    return entry



def verify_source_mirror(sources):
    """Bind all generator Python and include/license inputs to pinned Git blobs."""
    receipts = {}
    for suite, revision in [('qasmbench', QASMBENCH_REVISION), ('mqtbench', MQTBENCH_REVISION)]:
        commit = json.loads((sources / (suite + '_commit.json')).read_text())
        tree = json.loads((sources / (suite + '_tree.json')).read_text())
        if commit['sha'] != revision or tree['sha'] not in {revision, commit['commit']['tree']['sha']}:
            raise ValueError('upstream_revision_mismatch:' + suite)
        for item in tree['tree']:
            path = item['path']
            relevant = path in {'LICENSE', 'NOTICE', 'qelib1.inc', 'pyproject.toml'}
            relevant |= suite == 'mqtbench' and path.startswith('src/') and path.endswith('.py')
            if not relevant or item['type'] != 'blob':
                continue
            data = (sources / suite / path).read_bytes()
            if git_blob_digest(data) != item['sha']:
                raise ValueError('upstream_support_blob_mismatch:' + suite + '/' + path)
            receipts[suite + '/' + path] = {'sha256': digest(data), 'git_blob_sha1': item['sha']}
    return receipts


def build(sources, output):
    from qiskit import qasm2
    from qiskit import qpy
    sources, output = sources.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    integrity = verify_source_mirror(sources)
    tree = json.loads((sources / 'qasmbench_tree.json').read_text())
    if tree['sha'] != QASMBENCH_REVISION or json.loads((sources / 'qasmbench_commit.json').read_text())['sha'] != QASMBENCH_REVISION:
        raise ValueError('upstream_tree_mismatch')
    ledger = qasmbench_inventory(tree)
    cases = []
    qelib = (sources / 'qasmbench/qelib1.inc').read_text()
    # Inline the exact source include to prevent the parser from silently using
    # an installed qelib with a different gate definition or phase convention.
    for entry in ledger:
        if entry['selection_status'] != 'candidate':
            continue
        path = sources / 'qasmbench' / entry['source_path']
        raw = path.read_bytes()
        if len(raw) != entry['raw_source_bytes'] or git_blob_digest(raw) != entry['upstream_blob_sha1']:
            raise ValueError('upstream_blob_mismatch:' + entry['source_path'])
        entry['id'] = 'qasmbench__' + path.stem
        entry['source_url'] = f"https://github.com/pnnl/QASMBench/blob/{QASMBENCH_REVISION}/{entry['source_path']}"
        entry['license'] = 'upstream Battelle license (LICENSE and NOTICE retained)'
        try:
            text = raw.decode()
            if text.count('include "qelib1.inc";') != 1:
                raise ValueError('unsupported_include_configuration')
            circuit = qasm2.loads(text.replace('include "qelib1.inc";', qelib), custom_instructions=())
            cases.append(_record(circuit, entry, output, path.relative_to(sources), raw))
        except Exception as error:
            entry.update(selection_status='excluded', exclusion_reason=f'{type(error).__name__}: {error}',
                         raw_source_sha256=digest(raw), raw_artifact=str(path.relative_to(sources)))
        print(json.dumps({'candidate': entry['id'], 'status': entry['selection_status'], 'reason': entry['exclusion_reason']}), flush=True)
    sys.path.insert(0, str(sources / 'mqtbench/src'))
    from mqt.bench import get_benchmark_alg
    generated = output / 'generated_raw'
    generated.mkdir()
    for family in MQT_FAMILIES:
        for width in MQT_WIDTHS:
            name = f'mqtbench__{family}_n{width}'
            source_path = f'src/mqt/bench/benchmarks/{family}.py'
            entry = {'id': name, 'community_source': 'MQT Bench', 'source_path': source_path,
                     'source_revision': MQTBENCH_REVISION, 'family': family, 'directory_width': width,
                     'source_url': f'https://github.com/munich-quantum-toolkit/bench/blob/{MQTBENCH_REVISION}/{source_path}',
                     'generator_source_sha256': digest((sources / 'mqtbench' / source_path).read_bytes()),
                     'license': 'MIT (upstream LICENSE retained)',
                     'generator': {'api': 'get_benchmark_alg', 'level': 'ALG', 'benchmark': family,
                                   'circuit_size': width, 'random_parameters': True,
                                   'parameter_seed': 10, 'other_options': 'pinned upstream defaults',
                                   'no_transpile_or_mirror': True}}
            ledger.append(entry)
            try:
                circuit = get_benchmark_alg(family, width, random_parameters=True)
                if circuit.parameters:
                    raise ValueError('unbound_symbolic_parameters')
                raw_path = generated / (name + '.qpy')
                with raw_path.open('wb') as stream:
                    qpy.dump(circuit, stream)
                raw = raw_path.read_bytes()
                entry['raw_source_format'] = 'Qiskit QPY (original generated circuit)'
                try:
                    exported = qasm2.dumps(circuit).encode()
                    (generated / (name + '.qasm')).write_bytes(exported)
                    entry['optional_raw_qasm_sha256'] = digest(exported)
                except Exception as export_error:
                    entry['optional_raw_qasm_export_error'] = f'{type(export_error).__name__}: {export_error}'
                cases.append(_record(circuit, entry, output, raw_path.relative_to(output), raw))
            except Exception as error:
                entry.update(selection_status='excluded', exclusion_reason=f'{type(error).__name__}: {error}')
            print(json.dumps({'candidate': name, 'status': entry['selection_status'], 'reason': entry['exclusion_reason']}), flush=True)
    versions = {name: importlib.metadata.version(name) for name in ('qiskit', 'numpy')}
    provenance = {
        'qasmbench': {'repository': 'https://github.com/pnnl/QASMBench', 'revision': QASMBENCH_REVISION,
                      'license_file': 'qasmbench/LICENSE', 'notice_file': 'qasmbench/NOTICE',
                      'include_sha256': digest((sources / 'qasmbench/qelib1.inc').read_bytes())},
        'mqtbench': {'repository': 'https://github.com/munich-quantum-toolkit/bench', 'revision': MQTBENCH_REVISION,
                     'release': '2.3.0', 'license_file': 'mqtbench/LICENSE',
                     'archive_sha256': digest((sources / 'mqtbench-v2.3.0.tar.gz').read_bytes())},
        'benchpress_methodology_reference': {'repository': 'https://github.com/Qiskit/benchpress',
                                             'revision': BENCHPRESS_REVISION,
                                             'scope': 'methodology reference only; not a Benchpress execution'},
        'generation_versions': versions,
        'source_integrity': integrity,
        'selection_rule': {'fixed_before_compiler_runs': True, 'qasmbench': 'all small/medium non-transpiled QASM paths',
                           'max_width': MAX_QUBITS, 'max_qasmbench_raw_bytes': MAX_RAW_BYTES,
                           'max_expanded_gates': MAX_EXPANDED_GATES,
                           'mqt_families': MQT_FAMILIES, 'mqt_widths': MQT_WIDTHS,
                           'eligibility': 'static unitary circuit after only global terminal-measurement suffix and barrier removal',
                           'selection_uses_compiler_timings_or_quality': False},
        'scope': 'bounded external-source compilation study, not a representative sample, complete suite, hardware benchmark, or independent standard',
        'upstream_artifact_root': '../sources',
        'normalization_trust_boundary': 'pinned Qiskit parser, generator definitions, and Statevector reference; independent NumPy gate matrices validate the transformation',
        'normalizer_sha256': digest(Path(__file__).read_bytes()),
        'oracle_sha256': digest((Path(__file__).parent / 'comparison_validation.py').read_bytes()),
    }
    exclusions = [entry for entry in ledger if entry['selection_status'] == 'excluded']
    # Ledger carries compact selection metadata; source text appears only in cases.
    compact = [{k: v for k, v in entry.items() if k not in {'source_qasm', 'qasm'}} for entry in ledger]
    manifest = {'schema': 1, 'cases': sorted(cases, key=lambda entry: entry['id']),
                'exclusions': exclusions, 'candidate_ledger': compact, 'provenance': provenance,
                'counts': {'selected': len(cases), 'excluded': len(exclusions), 'ledger': len(ledger),
                           'selected_by_suite': dict(Counter(case['community_source'] for case in cases)),
                           'selected_by_family': dict(Counter(case['family'] for case in cases))}}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest['counts'], indent=2))



def acquire(sources):
    """Fetch only the frozen inventories, eligible raw inputs, and source archive."""
    import io
    import tarfile
    from urllib.request import urlopen
    sources = sources.resolve()
    sources.mkdir(parents=True, exist_ok=False)

    def fetch(url, path):
        data = urlopen(url, timeout=60).read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return data

    for suite, repo, revision in [('qasmbench', 'pnnl/QASMBench', QASMBENCH_REVISION),
                                   ('mqtbench', 'munich-quantum-toolkit/bench', MQTBENCH_REVISION),
                                   ('benchpress', 'Qiskit/benchpress', BENCHPRESS_REVISION)]:
        for suffix, endpoint in [('commit', 'commits/'), ('tree', 'git/trees/')]:
            url = f'https://api.github.com/repos/{repo}/{endpoint}{revision}'
            if suffix == 'tree':
                url += '?recursive=1'
            fetch(url, sources / (suite + '_' + suffix + '.json'))
    inventory = qasmbench_inventory(json.loads((sources / 'qasmbench_tree.json').read_text()))
    paths = [row['source_path'] for row in inventory if row['selection_status'] == 'candidate']
    paths += ['LICENSE', 'NOTICE', 'README.md', 'qelib1.inc']
    for path in paths:
        fetch(f'https://raw.githubusercontent.com/pnnl/QASMBench/{QASMBENCH_REVISION}/{path}', sources / 'qasmbench' / path)
    archive = fetch(f'https://api.github.com/repos/munich-quantum-toolkit/bench/tarball/{MQTBENCH_REVISION}',
                    sources / 'mqtbench-v2.3.0.tar.gz')
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        for member in stream.getmembers():
            if not member.isfile():
                continue
            parts = Path(member.name).parts[1:]
            if not parts or '..' in parts or Path(*parts).is_absolute():
                raise ValueError('unsafe_archive_member')
            target = sources / 'mqtbench' / Path(*parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(stream.extractfile(member).read())
    for path in ['README.md', 'default.conf', 'LICENSE.txt']:
        fetch(f'https://raw.githubusercontent.com/Qiskit/benchpress/{BENCHPRESS_REVISION}/{path}', sources / 'benchpress' / path)
    verify_source_mirror(sources)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--acquire', action='store_true', help='download pinned sources to a new directory')
    args = parser.parse_args()
    if args.acquire:
        acquire(args.sources)
    if args.output:
        build(args.sources, args.output)
    elif not args.acquire:
        parser.error('--output is required unless --acquire is used')


if __name__ == '__main__':
    main()
