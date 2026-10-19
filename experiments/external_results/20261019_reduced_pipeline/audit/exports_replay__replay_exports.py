"""Replay archived reporting sources without mutating archive or original study."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

OUT = Path(__file__).resolve().parent
WORK = OUT.parent
ARCHIVE = WORK / 'qrust_reduced_archive_audit_20261019/extracted/working/qrust_reduced_ladder_20261019'
REPORTING = ARCHIVE / 'reporting'
CAMPAIGN = ARCHIVE / 'corpus_work/campaign'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_new(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def source_inventory():
    return {p.relative_to(ARCHIVE).as_posix(): {'bytes': p.stat().st_size, 'sha256': sha(p)}
            for root in (REPORTING, CAMPAIGN) for p in sorted(root.rglob('*')) if p.is_file()}


def main():
    before = source_inventory()
    write_new(OUT / 'archived_inputs_before.json', before)
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MPLCONFIGDIR'] = str(OUT / 'matplotlib_cache')
    jobs = [
        ('tables', [REPORTING/'export_tables.py', '--campaign', CAMPAIGN, '--output', OUT/'tables']),
        ('plots', [REPORTING/'plot_ladder.py', '--analysis', CAMPAIGN/'analysis.json',
                   '--protocol', CAMPAIGN/'protocol.json', '--output-dir', OUT/'exports', '--postflight-authorized']),
        ('report_v2', [REPORTING/'render_report.py', '--campaign', CAMPAIGN, '--tables', OUT/'tables',
                       '--template', REPORTING/'report_template_v2.md', '--report', OUT/'report_v2.md',
                       '--receipt', OUT/'report_render_receipt_v2.json']),
        ('stage_interpretation', [REPORTING/'interpret_stage_results.py', '--campaign', CAMPAIGN,
                                 '--output', OUT/'stage_interpretation.json']),
        ('variability', [REPORTING/'audit_variability.py', '--campaign', CAMPAIGN, '--output', OUT/'variability.json']),
    ]
    commands = []
    for name, args in jobs:
        argv = [sys.executable, *map(str, args)]
        entry = {'job': name, 'argv': argv, 'cwd': str(OUT), 'env_overrides': {k: env[k] for k in ('PYTHONDONTWRITEBYTECODE','MPLCONFIGDIR')},
                 'source_sha256': sha(Path(args[0])), 'log': f'logs/{name}.log', 'started_unix': time.time()}
        with (OUT / entry['log']).open('xb') as log:
            process = subprocess.run(argv, cwd=OUT, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        entry['returncode'] = process.returncode
        entry['elapsed_s'] = time.time() - entry['started_unix']
        entry['log_sha256'] = sha(OUT / entry['log'])
        commands.append(entry)
        with (OUT / 'commands.jsonl').open('a') as stream:
            stream.write(json.dumps(entry, sort_keys=True) + '\n')
        print(json.dumps({'job': name, 'returncode': process.returncode, 'elapsed_s': entry['elapsed_s']}), flush=True)
        if process.returncode:
            write_new(OUT / 'failure_receipt.json', {'phase':name,'commands':commands,'archived_inputs_unchanged': source_inventory() == before})
            raise SystemExit(process.returncode)
    comparisons = []
    for subdir in ('tables','exports'):
        expected = {p.name for p in (REPORTING/subdir).iterdir() if p.is_file()}
        actual = {p.name for p in (OUT/subdir).iterdir() if p.is_file()}
        if actual != expected:
            raise ValueError(f'Output inventory mismatch for {subdir}: {actual} != {expected}')
        for name in sorted(expected):
            comparisons.append((f'{subdir}/{name}', REPORTING/subdir/name, OUT/subdir/name))
    comparisons += [(name, REPORTING/name, OUT/name) for name in
                    ('report_v2.md','report_render_receipt_v2.json','stage_interpretation.json','variability.json')]
    checks = [{'path': name, 'archived_bytes': archived.stat().st_size, 'replayed_bytes': replayed.stat().st_size,
               'archived_sha256': sha(archived), 'replayed_sha256': sha(replayed),
               'byte_identical': archived.read_bytes() == replayed.read_bytes()}
              for name,archived,replayed in comparisons]
    after = source_inventory()
    write_new(OUT/'archived_inputs_after.json', after)
    receipt = {'schema':'qrust_reduced_archived_exports_replay_v1', 'archived_source_root':str(ARCHIVE),
               'output_root':str(OUT), 'python':sys.executable, 'replay_source_sha256':sha(Path(__file__)),
               'source_files_checked':len(before), 'archived_inputs_unchanged':before == after,
               'commands':commands, 'comparisons':checks, 'all_outputs_byte_identical':all(c['byte_identical'] for c in checks),
               'scope':'Saved-output reporting replay only; no compiler invocation or numerical validation repeated.',
               'template_sha256':sha(REPORTING/'report_template_v2.md')}
    write_new(OUT/'replay_receipt.json',receipt)
    print(json.dumps({'outputs_compared':len(checks),'all_outputs_byte_identical':receipt['all_outputs_byte_identical'],
                      'archived_inputs_unchanged':receipt['archived_inputs_unchanged'],'receipt_sha256':sha(OUT/'replay_receipt.json')}),flush=True)
    if not receipt['all_outputs_byte_identical'] or not receipt['archived_inputs_unchanged']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
