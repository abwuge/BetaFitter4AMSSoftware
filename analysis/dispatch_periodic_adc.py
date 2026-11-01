#!/usr/bin/env python3
"""Run the thirty declared period-response fits on separate scheduler host pools."""
import argparse
import json
import os
import socket
import subprocess
import sys
import time
from pathlib import Path
from research_periodic_adc import ROOT, OUT, CONFIGS, SOURCES
from research_flight_conditional import write, digest
from research_conditional_mc import stamp

BATCH = OUT / 'batch'
RUNTIME = ROOT / 'results/flight_calibration_research/continuation_20261013/prospective_flight/scoring/recalibration/batch/shared_runtime'


def queue():
    return json.loads(subprocess.check_output(['condor_q', 'hxwu', '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')


def prepare():
    if BATCH.exists():
        raise FileExistsError('Inspect the existing dispatch before preparing another')
    if queue():
        raise ValueError('Inspect existing scheduler work before adding the period fits')
    tasks = [(dataset, family, variant, start) for dataset in ['flight', 'mc']
        for family, _, _ in CONFIGS[dataset] for variant in ['static', 'step', 'spline'] for start in range(2)]
    if len(tasks) > 50:
        raise ValueError('More than fifty independent fit tasks')
    for dataset in ['flight', 'mc']:
        record = json.loads((OUT / dataset / 'numerical_verification.json').read_text())
        if not record['complete']:
            raise ValueError('Numerical verification must precede dispatch')
    runtime = json.loads((RUNTIME / 'manifest.json').read_text())
    if not runtime['complete'] or not (RUNTIME / 'python/numpy/__init__.py').exists():
        raise ValueError('The previously validated shared Python runtime is absent')
    text = subprocess.check_output(['condor_status', '-constraint',
        '(OpSysAndVer=="CentOS7") && (Cpus>=4) && (Memory>=3072)', '-af', 'Machine', 'State', 'Activity']).decode()
    hosts = {}
    for line in text.splitlines():
        host, state, activity = line.split()
        hosts[host] = hosts.get(host, 0) + int(state == 'Unclaimed' and activity == 'Idle')
    ordered = sorted(hosts, key=lambda host: (-hosts[host], host))
    if len(ordered) < len(tasks):
        raise ValueError('Not enough distinct eligible hosts for the thirty fit tasks')
    pools = [ordered[i::len(tasks)] for i in range(len(tasks))]
    BATCH.mkdir(parents=True)
    worker = BATCH / 'worker.sh'
    worker.write_text('''#!/usr/bin/env bash
set -eo pipefail
cd {root}
unset PYTHONHOME
export PATH=/usr/bin:/bin:$PATH
export PYTHONPATH={runtime}/python:/home/ams/hxwu/.local/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
/usr/bin/python3 -c 'import numpy,scipy,awkward; assert numpy.__version__=="1.19.1" and scipy.__version__=="0.18.1"; assert awkward.__version__=="1.9.0"; print("PERIOD_RUNTIME",numpy.__version__,scipy.__version__,flush=True)'
/usr/bin/python3 -m unittest discover -s tests -p test_periodic_adc.py -v
exec /usr/bin/python3 analysis/dispatch_periodic_adc.py worker --index "$1"
'''.format(root=ROOT, runtime=RUNTIME))
    worker.chmod(0o755)
    submission = '\n'.join(['universe = vanilla', 'executable = ' + str(worker), 'arguments = $(number)',
        'getenv = true', 'request_cpus = 4', 'request_memory = 3GB',
        'requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'should_transfer_files = NO', 'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = ' + str(BATCH / 'worker_$(number).out'), 'error = ' + str(BATCH / 'worker_$(number).err'),
        'log = ' + str(BATCH / 'condor.log'), 'queue number,hosts from ('] +
        [str(i) + ' ' + ';'.join(pools[i]) for i in range(len(tasks))] + [')', ''])
    (BATCH / 'fits.sub').write_text(submission)
    protected = {str(path.relative_to(ROOT)): digest(path) for path in [ROOT / 'analysis' / name for name in SOURCES] +
        [Path(__file__).resolve(), ROOT / 'tests/test_periodic_adc.py', RUNTIME / 'manifest.json']}
    write(BATCH / 'protocol.json', dict(created_utc=stamp(), tasks=tasks, workers=len(tasks), maximum_workers=50,
        physical_host_pools=pools, request_cpus=4, request_memory_MB=3072, source_sha256=protected,
        worker_sha256=digest(worker), submit_sha256=digest(BATCH / 'fits.sub'),
        scope='Thirty independent fits, fewer than the fifty-host default because only thirty starts were declared; disjoint physical-host eligibility pools',
        data_scope='Original separate MC and flight measurement pools; no beta reference read',
        scientific_protocols={dataset: digest(OUT / dataset / 'protocol.json') for dataset in ['flight','mc']}))
    print('PERIOD_DISPATCH_PREPARED', len(tasks), len(ordered), flush=True)


def submit():
    if (BATCH / 'submission.json').exists():
        raise FileExistsError('Period fit batch was already submitted')
    p = json.loads((BATCH / 'protocol.json').read_text())
    if queue() or digest(BATCH / 'fits.sub') != p['submit_sha256']:
        raise ValueError('Unexpected jobs or changed submission')
    answer = subprocess.check_output(['condor_submit', '-terse', str(BATCH / 'fits.sub')]).decode().strip()
    write(BATCH / 'submission.json', dict(created_utc=stamp(), result=answer, cluster=int(answer.split('.')[0]),
        protocol_sha256=digest(BATCH / 'protocol.json')))
    print('PERIOD_SUBMITTED', answer, flush=True)


def worker(index):
    p = json.loads((BATCH / 'protocol.json').read_text())
    for name, fingerprint in p['source_sha256'].items():
        if digest(ROOT / name) != fingerprint:
            raise ValueError('Changed worker source ' + name)
    dataset, family, variant, number = p['tasks'][index]
    if digest(OUT / dataset / 'protocol.json') != p['scientific_protocols'][dataset]:
        raise ValueError('Changed scientific fit protocol')
    result_path = OUT / dataset / (family + '_' + variant) / ('start_%d.json' % number)
    marker = BATCH / ('task_%02d.json' % index)
    if marker.exists():
        raise FileExistsError('Task completion already recorded')
    if result_path.exists():
        raise FileExistsError('Inspect existing result before attempting the same task')
    begin = time.monotonic()
    command = [sys.executable, 'analysis/research_periodic_adc.py', 'start', '--dataset', dataset,
        '--family', family, '--variant', variant, '--start', str(number)]
    returncode = subprocess.call(command, cwd=str(ROOT))
    record = dict(created_utc=stamp(), index=index, task=p['tasks'][index], host=socket.getfqdn(), pid=os.getpid(),
        exit_code=returncode, elapsed_seconds=time.monotonic()-begin, protocol_sha256=digest(BATCH / 'protocol.json'))
    if result_path.exists():
        record['result_sha256'] = digest(result_path)
        record['numerically_converged'] = json.loads(result_path.read_text())['numerically_converged']
    write(marker, record)
    if returncode:
        raise RuntimeError('The period fit failed; retain its output and error logs')
    print('PERIOD_WORKER_COMPLETE', json.dumps(record), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['prepare','submit','worker'])
    parser.add_argument('--index', type=int, default=0); args = parser.parse_args()
    if args.stage == 'worker':
        worker(args.index)
    else:
        globals()[args.stage]()
