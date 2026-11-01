#!/usr/bin/env python3
"""Replay 128 saved flight calibrations with an immutable shared work queue."""
import argparse
import json
import os
import socket
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from bootstrap_prospective_flight import ROOT, OUT, REF, DEST, check
from research_flight_conditional import write, digest
from research_conditional_mc import stamp

BATCH = DEST / 'batch'
ENV = '/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'


def queue(cluster=None):
    return json.loads(subprocess.check_output(['condor_q', str(cluster) if cluster is not None else 'hxwu',
        '-json', '-attributes', 'ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')


def prepare():
    p, audit = check()
    endpoint = json.loads((REF / 'endpoint_sensitivity/evaluation.json').read_text())
    if not endpoint['complete'] or endpoint['primary_predictions_changed']:
        raise ValueError('Resolve and record the endpoint sensitivity without changing primary outputs')
    if queue():
        raise ValueError('Inspect existing user jobs before allocating a new fifty-host batch')
    BATCH.mkdir(exist_ok=True)
    (BATCH / 'claims').mkdir(exist_ok=True)
    if (BATCH / 'protocol.json').exists():
        raise FileExistsError('Prospective recalibration dispatch already prepared')
    text = subprocess.check_output(['condor_status', '-constraint',
        '(OpSysAndVer=="CentOS7") && (Cpus>=32) && (Memory>=8192)', '-af', 'Machine', 'State', 'Activity']).decode()
    hosts = {}
    for line in text.splitlines():
        host, state, activity = line.split()
        hosts[host] = hosts.get(host, 0) + int(state == 'Unclaimed' and activity == 'Idle')
    ordered = sorted(hosts, key=lambda host: (-hosts[host], host))
    if len(ordered) < 50:
        raise ValueError('Fewer than fifty eligible physical hosts; inspect resource allocation')
    groups = [ordered[i::50] for i in range(50)]
    worker = BATCH / 'worker.sh'
    worker.write_text('#!/usr/bin/env bash\nset -eo pipefail\ncd ' + str(ROOT) + '\nsource ' + ENV +
        '\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\nexec python3 analysis/dispatch_prospective_recalibration.py worker --worker "$1"\n')
    worker.chmod(0o755)
    submit = '\n'.join(['universe = vanilla', 'executable = ' + str(worker), 'arguments = $(worker)', 'getenv = true',
        'request_cpus = 32', 'request_memory = 8GB',
        'requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'should_transfer_files = NO', 'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = ' + str(BATCH / 'worker_$(Process).out'), 'error = ' + str(BATCH / 'worker_$(Process).err'),
        'log = ' + str(BATCH / 'condor.log'), 'queue worker,hosts from ('] +
        [str(i) + ' ' + ';'.join(groups[i]) for i in range(50)] + [')', ''])
    (BATCH / 'replay.sub').write_text(submit)
    write(BATCH / 'protocol.json', dict(created_utc=stamp(), workers=50, subprocesses_per_worker=4,
        saved_calibrations=p['training_replicates'], physical_host_pools=groups,
        source_sha256=digest(Path(__file__)), worker_sha256=digest(worker), submit_sha256=digest(BATCH / 'replay.sub'),
        inference_source_sha256=digest(ROOT / 'analysis/bootstrap_prospective_flight.py'),
        inference_protocol_sha256=digest(DEST / 'protocol.json'), reference_sha256=audit['reference_sha256'],
        endpoint_sensitivity_sha256=digest(REF / 'endpoint_sensitivity/evaluation.json'),
        scheduling='Atomic per-replicate claims in a shared queue; a worker does not read a partially written completion file of a different worker',
        no_training=True, original_frozen_predictions_changed=False))
    print('PROSPECTIVE_REPLAY_BATCH', p['training_replicates'], 'saved calibrations', 50, 'host pools', flush=True)


def submit():
    p = json.loads((BATCH / 'protocol.json').read_text())
    if (BATCH / 'submission.json').exists():
        raise FileExistsError('Prospective replay already submitted')
    if queue():
        raise ValueError('Unexpected user jobs before prospective replay submission')
    if digest(BATCH / 'replay.sub') != p['submit_sha256']:
        raise ValueError('Prospective replay submit changed')
    answer = subprocess.check_output(['condor_submit', '-terse', str(BATCH / 'replay.sub')]).decode().strip()
    cluster = int(answer.split('.')[0])
    write(BATCH / 'submission.json', dict(created_utc=stamp(), cluster=cluster, result=answer,
        protocol_sha256=digest(BATCH / 'protocol.json')))
    print('PROSPECTIVE_REPLAY_SUBMITTED', answer, flush=True)


def worker(number):
    p = json.loads((BATCH / 'protocol.json').read_text())
    for path, key in [(Path(__file__), 'source_sha256'), (ROOT / 'analysis/bootstrap_prospective_flight.py', 'inference_source_sha256'),
        (DEST / 'protocol.json', 'inference_protocol_sha256'), (REF / 'reference.npz', 'reference_sha256')]:
        if digest(path) != p[key]:
            raise ValueError('Prospective replay worker input changed: ' + key)

    def loop(slot):
        complete = []
        for rep in range(p['saved_calibrations']):
            claim = BATCH / 'claims' / ('%04d' % rep)
            try:
                claim.mkdir()
            except FileExistsError:
                continue
            write(claim / 'owner.json', dict(created_utc=stamp(), worker=number, slot=slot, host=socket.getfqdn(), pid=os.getpid()))
            try:
                path = DEST / 'parts' / ('%04d.json' % rep)
                if path.exists():
                    raise FileExistsError('Unclaimed pre-existing replay; inspect before running it')
                with (claim / 'replay.log').open('w') as log:
                    result = subprocess.run(['python3', 'analysis/bootstrap_prospective_flight.py', 'replicate', '--replicate', str(rep)],
                        stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise RuntimeError('Calibration replay %d failed with exit %d' % (rep, result.returncode))
                note = json.loads(path.read_text())
                if not note['complete'] or note['replicate'] != rep or note['protocol_sha256'] != p['inference_protocol_sha256'] or digest(path.with_suffix('.npz')) != note['output_sha256']:
                    raise ValueError('Saved calibration replay failed integrity')
                write(claim / 'complete.json', dict(complete=True, created_utc=stamp(), replicate=rep, completion_sha256=digest(path)))
                complete.append(rep)
                print('PROSPECTIVE_REPLAY_PART', rep, flush=True)
            except Exception as exc:
                write(claim / 'failure.json', dict(created_utc=stamp(), replicate=rep, error=repr(exc)))
                raise
        return complete

    with ThreadPoolExecutor(max_workers=p['subprocesses_per_worker']) as pool:
        parts = list(pool.map(loop, range(p['subprocesses_per_worker'])))
    write(BATCH / ('worker_%02d.json' % number), dict(complete=True, created_utc=stamp(), worker=number,
        replicates=[rep for group in parts for rep in group], interpretation='Shared queue has no unclaimed replica; global completion requires all 128 completion files'))


def status():
    p = json.loads((BATCH / 'protocol.json').read_text())
    counts = dict(replicates=p['saved_calibrations'], claimed=0, complete=0, failed=0)
    for rep in range(p['saved_calibrations']):
        claim = BATCH / 'claims' / ('%04d' % rep)
        counts['claimed'] += int(claim.exists())
        counts['complete'] += int((claim / 'complete.json').exists())
        counts['failed'] += int((claim / 'failure.json').exists())
    print('PROSPECTIVE_REPLAY_STATUS', json.dumps(counts), flush=True)


def retire():
    destination = BATCH / 'retirement.json'
    if destination.exists():
        raise FileExistsError('Empty replay workers already retired')
    p = json.loads((BATCH / 'protocol.json').read_text())
    evaluation = json.loads((DEST / 'evaluation.json').read_text())
    if not evaluation['complete'] or evaluation['replicates'] != p['saved_calibrations'] or evaluation['protocol_sha256'] != p['inference_protocol_sha256']:
        raise ValueError('Complete recalibration intervals required before retiring empty workers')
    cluster = json.loads((BATCH / 'submission.json').read_text())['cluster']
    before = queue(cluster)
    result = dict(created_utc=stamp(), cluster=cluster, before=before,
        inference_evaluation_sha256=digest(DEST / 'evaluation.json'), running_jobs_preserved=True, complete=False)
    write(destination, result)
    idle = [r for r in before if r['JobStatus'] == 1]
    if idle:
        result['scheduler_response'] = subprocess.check_output(['condor_rm', '-constraint',
            'ClusterId == %d && JobStatus == 1' % cluster], stderr=subprocess.STDOUT).decode()
    result['after'] = queue(cluster); result['complete'] = True
    write(destination, result)
    print('PROSPECTIVE_REPLAY_EMPTY_WORKERS_RETIRED', len(idle), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['prepare', 'submit', 'worker', 'status', 'retire'])
    parser.add_argument('--worker', type=int); args = parser.parse_args()
    worker(args.worker) if args.stage == 'worker' else globals()[args.stage]()
