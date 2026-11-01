#!/usr/bin/env python3
"""Run the frozen native RICH reread with shared claims and disjoint host pools."""
import argparse
import json
import os
import socket
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
from reference_prospective_flight_schema import ROOT, OUT, REF, BINARY, check
from research_flight_conditional import write, digest
from research_conditional_mc import stamp

BATCH = REF / 'batch'
ENV = '/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'


def prepare():
    p = check()
    prep = json.loads((REF / 'preparation.json').read_text())
    BATCH.mkdir(exist_ok=True)
    if (BATCH / 'protocol.json').exists():
        raise FileExistsError('Prospective reference batch already declared')
    queued = json.loads(subprocess.check_output(['condor_q', 'hxwu', '-json']).decode() or '[]')
    if queued:
        raise ValueError('Finish or inspect the existing user jobs before allocating fifty reference workers')
    text = subprocess.check_output(['condor_status', '-constraint',
        '(OpSysAndVer=="CentOS7") && (Cpus>=32) && (Memory>=8192)', '-af', 'Machine', 'State', 'Activity']).decode()
    hosts = {}
    for line in text.splitlines():
        host, state, activity = line.split()
        hosts[host] = hosts.get(host, 0) + int(state == 'Unclaimed' and activity == 'Idle')
    ordered = sorted(hosts, key=lambda h: (-hosts[h], h))
    if len(ordered) < 50:
        raise ValueError('Fewer than fifty eligible physical hosts; inspect before reducing requested concurrency')
    groups = [ordered[i::50] for i in range(50)]
    launcher = BATCH / 'worker.sh'
    launcher.write_text('#!/usr/bin/env bash\nset -eo pipefail\ncd ' + str(ROOT) + '\nsource ' + ENV +
        '\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\nexec python3 analysis/dispatch_prospective_reference.py worker --worker "$1"\n')
    launcher.chmod(0o755)
    submit = '\n'.join(['universe = vanilla', 'executable = ' + str(launcher), 'arguments = $(worker)', 'getenv = true',
        'request_cpus = 32', 'request_memory = 8GB',
        'requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'should_transfer_files = NO', 'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = ' + str(BATCH / 'worker_$(Process).out'), 'error = ' + str(BATCH / 'worker_$(Process).err'),
        'log = ' + str(BATCH / 'condor.log'), 'queue worker,hosts from ('] +
        [str(i) + ' ' + ';'.join(groups[i]) for i in range(50)] + [')', ''])
    (BATCH / 'reference.sub').write_text(submit)
    write(BATCH / 'protocol.json', dict(created_utc=stamp(), jobs=50, processes_per_job=4,
        requested_CPUs_per_job=32, physical_host_pools=groups, sources=[r['source'] for r in prep['jobs']],
        scheduling='All workers claim from the same immutable source manifest using an atomic per-source directory. Idle workers are not required for completion when other workers can finish the shared queue.',
        failure='Never reclaim an existing source claim automatically or overwrite incomplete ROOT output. A failed source remains visible for diagnosis.',
        source_sha256=digest(Path(__file__)), launcher_sha256=digest(launcher), submit_sha256=digest(BATCH / 'reference.sub'),
        preparation_sha256=digest(REF / 'preparation.json'), reference_protocol_sha256=digest(REF / 'protocol.json'),
        schema_amendment_sha256=digest(REF / 'schema_amendment.json'),
        native_binary_sha256=p['native_binary_sha256'], no_model_or_selection_change=True))
    print('PROSPECTIVE_REFERENCE_BATCH', len(prep['jobs']), 'sources', 50, 'disjoint host pools', flush=True)


def submit():
    p = json.loads((BATCH / 'protocol.json').read_text())
    if (BATCH / 'submission.json').exists():
        raise FileExistsError('Reference batch was already submitted')
    if digest(BATCH / 'reference.sub') != p['submit_sha256']:
        raise ValueError('Reference submit file changed')
    if json.loads(subprocess.check_output(['condor_q', 'hxwu', '-json']).decode() or '[]'):
        raise ValueError('Unexpected existing user jobs before reference submission')
    answer = subprocess.check_output(['condor_submit', '-terse', str(BATCH / 'reference.sub')]).decode().strip()
    cluster = int(answer.split('.')[0])
    write(BATCH / 'submission.json', dict(created_utc=stamp(), cluster=cluster, result=answer,
        protocol_sha256=digest(BATCH / 'protocol.json')))
    print('PROSPECTIVE_REFERENCE_SUBMITTED', answer, flush=True)


def worker(number):
    p = json.loads((BATCH / 'protocol.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (BINARY, 'native_binary_sha256'),
                      (REF / 'preparation.json', 'preparation_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Native reference worker input changed: ' + key)
    prep = json.loads((REF / 'preparation.json').read_text())
    jobs = sorted(prep['jobs'], key=lambda r: (-r['events'], r['source']))

    def loop(slot):
        completed = []
        for row in jobs:
            source = row['source']
            path = REF / 'parts' / ('%05d.root' % source)
            note = path.with_suffix('.json')
            if note.exists():
                done = json.loads(note.read_text())
                if not done['complete'] or done['binary_sha256'] != p['native_binary_sha256'] or digest(path) != done['output_sha256']:
                    raise ValueError('Existing native reference completion changed')
                continue
            claim = REF / 'claims' / ('%05d' % source)
            try:
                claim.mkdir()
            except FileExistsError:
                continue
            write(claim / 'owner.json', dict(created_utc=stamp(), worker=number, slot=slot,
                host=socket.getfqdn(), pid=os.getpid(), source=source))
            try:
                if path.exists() or note.exists():
                    raise FileExistsError('Inspect unrecorded native reference output')
                for filename, key in [('targets', 'targets_sha256'), ('raw_files', 'raw_files_sha256')]:
                    if digest(row[filename]) != row[key]:
                        raise ValueError('Reference native manifest changed: ' + filename)
                with path.with_suffix('.log').open('w') as log:
                    result = subprocess.run([str(BINARY), row['raw_files'], row['targets'], str(path)], stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise RuntimeError('Native RICH reader exit %d for source %d' % (result.returncode, source))
                expected = np.loadtxt(row['targets'], dtype=np.int64, ndmin=2)
                with uproot.open(path) as f:
                    a = f['reference'].arrays(['run', 'event', 'itrtrack', 'irich', 'found'], library='np')
                if len(a['run']) != row['events'] or not a['found'].all():
                    raise ValueError('Missing native RICH event')
                for j, key in enumerate(['run', 'event', 'itrtrack', 'irich']):
                    if not np.array_equal(a[key], expected[:, j]):
                        raise ValueError('Native RICH target association mismatch')
                temporary = note.with_suffix('.pending')
                write(temporary, dict(complete=True, created_utc=stamp(), source=source, events=len(a['run']), worker=number,
                    binary_sha256=p['native_binary_sha256'], output_sha256=digest(path), targets_sha256=row['targets_sha256'],
                    protocol_sha256=digest(REF / 'protocol.json')))
                os.replace(str(temporary), str(note))
                completed.append(source)
                print('RICH_SOURCE_COMPLETE', source, len(a['run']), flush=True)
            except Exception as exc:
                write(claim / 'failure.json', dict(created_utc=stamp(), source=source, error=repr(exc)))
                raise
        return completed

    with ThreadPoolExecutor(max_workers=p['processes_per_job']) as pool:
        completed = list(pool.map(loop, range(p['processes_per_job'])))
    write(BATCH / ('worker_%02d.json' % number), dict(complete=True, created_utc=stamp(), worker=number,
        completed_sources=[s for group in completed for s in group], interpretation='Worker exhausted unclaimed sources; global completion requires every source manifest'))


def status():
    p = json.loads((REF / 'preparation.json').read_text())
    counts = dict(sources=len(p['jobs']), complete=0, claimed=0, failed=0)
    for row in p['jobs']:
        tag = '%05d' % row['source']
        counts['complete'] += int((REF / 'parts' / (tag + '.json')).exists())
        counts['claimed'] += int((REF / 'claims' / tag).exists())
        counts['failed'] += int((REF / 'claims' / tag / 'failure.json').exists())
    print('PROSPECTIVE_REFERENCE_STATUS', json.dumps(counts), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'submit', 'worker', 'status'])
    parser.add_argument('--worker', type=int)
    args = parser.parse_args()
    worker(args.worker) if args.stage == 'worker' else globals()[args.stage]()
