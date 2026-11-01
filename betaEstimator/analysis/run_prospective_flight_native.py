#!/usr/bin/env python3
"""Probe, distribute, and verify the native measured supplement without RICH."""
import argparse
import json
import os
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import uproot
from prepare_prospective_flight import OUT, ROOT
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp

BINARY = OUT / 'bin/native'
ENV = '/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh'
OLD = ROOT / 'results/tofcalib_validation/bin/supplement_light'


def probe_prepare():
    dest = OUT / 'probe'; dest.mkdir(exist_ok=True)
    if (dest / 'protocol.json').exists():
        raise FileExistsError('Native probe already prepared')
    a = load(OUT / 'candidates.npz'); records = json.loads((OUT / 'preparation.json').read_text())['jobs']
    sources = np.array([r['source'] for r in records]); means = np.array([np.median(a['run'][a['source'] == s]) for s in sources])
    chosen = sources[np.argsort(means)[np.linspace(0, len(sources) - 1, 4).astype(int)]]
    ix = np.concatenate([np.flatnonzero(a['source'] == s)[:8] for s in chosen])
    q = np.where(a['tk_qin'][ix, 0, 2] > 8.5, a['tk_qin'][ix, 0, 2], a['tk_q'][ix, 1])
    native = np.column_stack([a['run'][ix], a['event'][ix], a['itrtrack'][ix], a['ibetah'][ix], q])
    old = np.column_stack([native[:, :4], np.full(len(ix), -1), q])
    np.savetxt(dest / 'targets.txt', native, fmt=['%d', '%d', '%d', '%d', '%.9g'])
    np.savetxt(dest / 'old_targets.txt', old, fmt=['%d', '%d', '%d', '%d', '%d', '%.9g'])
    runs = set(int(r) for r in a['run'][ix]); raws = set()
    for s in chosen:
        for path in (OUT / 'native' / ('%05d' % s) / 'raw_files.txt').read_text().splitlines():
            if int(Path(path).name.split('.')[0]) in runs:
                raws.add(path)
    (dest / 'raw_files.txt').write_text('\n'.join(sorted(raws)) + '\n')
    np.savez_compressed(dest / 'expected.npz', **{k: v[ix] for k, v in a.items()})
    write(dest / 'protocol.json', dict(created_utc=stamp(), events=len(ix), sources=chosen,
        selection='First eight measured entries from four sources spaced across median run date; no response residual or RICH selection',
        binary_sha256=digest(BINARY), original_binary_sha256=digest(OLD), source_sha256=digest(ROOT / 'analysis/prospective_flight_native.C'),
        targets_sha256=digest(dest / 'targets.txt'), old_targets_sha256=digest(dest / 'old_targets.txt'),
        raw_files_sha256=digest(dest / 'raw_files.txt'), expected_sha256=digest(dest / 'expected.npz'), RICH_read=False))
    print('PROSPECTIVE_NATIVE_PROBE', len(ix), len(raws), flush=True)


def probe_check():
    dest = OUT / 'probe'; p = json.loads((dest / 'protocol.json').read_text()); a = load(dest / 'expected.npz')
    with uproot.open(dest / 'native.root') as f:
        n = f['measured'].arrays(library='np')
    with uproot.open(dest / 'old.root') as f:
        b = f['beta_extra'].arrays(library='np')
    identity = {}
    for key in ['run', 'event', 'itrtrack', 'ibetah']:
        identity[key] = bool(np.array_equal(a[key], n[key]) and np.array_equal(n[key], b[key]))
    if not all(identity.values()) or not n['found'].all() or not n['selected'].all() or not b['selected'].all():
        raise ValueError('Native probe identity or selection failed')
    comparisons = {}
    for label, original, new in [('R', b['tk_rigidityk'][:, :, 1], n['R']), ('sigma_r', b['tk_erigidity1'][:, 2, 1], n['sigma_r']),
        ('native_time', a['t_native'], n['time']), ('ADC', a['adc'], n['adc']), ('bar', a['bar'], n['bar']),
        ('path_cm', a['x'] * 29.9792458, n['length']), ('native_chi', a['tk_chis1'][:, 2, 2, 1, :], n['chi'])]:
        diff = np.asarray(new, float) - np.asarray(original, float)
        comparisons[label] = dict(max_abs=float(np.max(abs(diff))), exact=bool(np.array_equal(original, new)), rms=float(np.sqrt(np.mean(diff * diff))))
    required = all(comparisons[k]['exact'] for k in ['R', 'sigma_r', 'native_time', 'ADC', 'bar']) and comparisons['path_cm']['max_abs'] < 1e-5
    if not required:
        raise ValueError('Required native probe measurements differ: ' + str(comparisons))
    write(dest / 'verification.json', dict(passed=True, created_utc=stamp(), identity=identity, comparisons=comparisons,
        required_measurements_equal=required, events=len(n['run']), binary_sha256=digest(BINARY), original_binary_sha256=digest(OLD),
        native_sha256=digest(dest / 'native.root'), old_sha256=digest(dest / 'old.root'), protocol_sha256=digest(dest / 'protocol.json'), RICH_read=False))
    print('PROSPECTIVE_NATIVE_PROBE_PASSED', comparisons, flush=True)


def batch_prepare():
    proof = json.loads((OUT / 'probe/verification.json').read_text())
    if not proof['passed'] or digest(BINARY) != proof['binary_sha256']:
        raise ValueError('Native equivalence gate failed')
    directory = OUT / 'batch'; directory.mkdir(exist_ok=True)
    if (directory / 'protocol.json').exists():
        raise FileExistsError('Native batch is already declared')
    jobs = json.loads((OUT / 'preparation.json').read_text())['jobs']
    active = subprocess.check_output(['condor_q', 'hxwu', '-json']).decode()
    queued = json.loads(active) if active.strip() else []
    if queued:
        raise ValueError('Inspect existing user jobs before allocating the fifty-node batch')
    raw = subprocess.check_output(['condor_status', '-constraint', 'OpSysAndVer=="CentOS7"', '-af', 'Machine', 'State', 'Activity']).decode()
    hosts = {}; states = {}
    for line in raw.splitlines():
        h, state, activity = line.split(); states.setdefault(h, []).append([state, activity])
        hosts[h] = hosts.get(h, 0) + int(state == 'Unclaimed' and activity == 'Idle')
    ordered = sorted(hosts, key=lambda h: (-hosts[h], h))
    if len(ordered) < 50:
        raise ValueError('Fewer than fifty eligible physical hosts; inspect resource state')
    host_groups = [ordered[i::50] for i in range(50)]
    task_groups = [[] for _ in range(50)]; loads = np.zeros(50, int)
    for row in sorted(jobs, key=lambda r: -r['events']):
        i = int(np.argmin(loads)); task_groups[i].append(row['source']); loads[i] += row['events']
    for i, rows in enumerate(task_groups):
        (directory / ('tasks_%02d.txt' % i)).write_text('\n'.join(map(str, rows)) + '\n')
    launcher = directory / 'worker.sh'
    launcher.write_text('#!/usr/bin/env bash\nset -eo pipefail\ncd ' + str(ROOT) + '\nsource ' + ENV + '\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\nexec python3 analysis/run_prospective_flight_native.py worker --worker "$1"\n')
    launcher.chmod(0o755)
    submit = '\n'.join(['universe = vanilla', 'executable = ' + str(launcher), 'arguments = $(worker)', 'getenv = true',
        'request_cpus = 4', 'request_memory = 8GB',
        'requirements = (TARGET.OpSysAndVer == "CentOS7") && stringListMember(TARGET.Machine,"$(hosts)",";")',
        'should_transfer_files = NO', 'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
        'output = ' + str(directory / 'worker_$(Process).out'), 'error = ' + str(directory / 'worker_$(Process).err'),
        'log = ' + str(directory / 'condor.log'), 'queue worker,hosts from ('] +
        [str(i) + ' ' + ';'.join(host_groups[i]) for i in range(50)] + [')', ''])
    (directory / 'native.sub').write_text(submit)
    write(directory / 'protocol.json', dict(created_utc=stamp(), jobs=50, processes_per_job=4, requested_CPUs_per_job=4,
        physical_host_pools=host_groups, sources=task_groups, events_per_job=loads, eligible_host_count=len(ordered),
        guarantee='Disjoint host candidate groups prevent two of this batch jobs sharing one physical host; no claim that other users are absent from those hosts',
        preexisting_user_jobs=0, binary_sha256=digest(BINARY), wrapper_sha256=digest(Path(__file__)),
        launcher_sha256=digest(launcher), submit_sha256=digest(directory / 'native.sub'),
        probe_sha256=digest(OUT / 'probe/verification.json'), preparation_sha256=digest(OUT / 'preparation.json'), RICH_read=False))
    print('PROSPECTIVE_NATIVE_BATCH_READY', len(jobs), 50, 'host_pools', len(ordered), flush=True)


def worker(number):
    p = json.loads((OUT / 'batch/protocol.json').read_text())
    if digest(Path(__file__)) != p['wrapper_sha256'] or digest(BINARY) != p['binary_sha256']:
        raise ValueError('Native executable or wrapper changed')
    notes = json.loads((OUT / 'preparation.json').read_text())['jobs']; notes = {r['source']: r for r in notes}
    def one(source):
        directory = OUT / 'native' / ('%05d' % source); path = directory / 'measurement.root'; record = directory / 'complete.json'
        if record.exists():
            old = json.loads(record.read_text())
            if old['binary_sha256'] == p['binary_sha256'] and digest(path) == old['output_sha256']:
                return old
            raise ValueError('Native source completion fingerprint changed')
        if path.exists():
            raise FileExistsError('Inspect incomplete native output rather than overwriting it')
        row = notes[source]
        for name, key in [('targets.txt', 'targets_sha256'), ('raw_files.txt', 'raw_files_sha256')]:
            if digest(directory / name) != row[key]:
                raise ValueError('Native input manifest changed')
        with (directory / 'execution.log').open('w') as log:
            proc = subprocess.run([str(BINARY), str(directory / 'raw_files.txt'), str(directory / 'targets.txt'), str(path)], stdout=log, stderr=subprocess.STDOUT)
        if proc.returncode:
            raise RuntimeError('Native source failed ' + str(source) + ' exit ' + str(proc.returncode))
        with uproot.open(path) as f:
            a = f['measured'].arrays(['run', 'event', 'found', 'selected'], library='np')
        expected = np.loadtxt(directory / 'targets.txt', ndmin=2)
        if len(a['run']) != row['events'] or not np.array_equal(a['run'], expected[:, 0]) or not np.array_equal(a['event'], expected[:, 1]) or not a['found'].all():
            raise ValueError('Native output incomplete or misassociated')
        result = dict(complete=True, created_utc=stamp(), source=source, worker=number, events=len(a['run']), selected=int(a['selected'].sum()),
            binary_sha256=p['binary_sha256'], output_sha256=digest(path), targets_sha256=row['targets_sha256'], RICH_read=False)
        write(record, result); print('NATIVE_SOURCE_COMPLETE', source, len(a['run']), flush=True); return result
    with ThreadPoolExecutor(max_workers=p['processes_per_job']) as pool:
        results = list(pool.map(one, p['sources'][number]))
    write(OUT / 'batch' / ('worker_%02d.json' % number), dict(complete=True, worker=number, created_utc=stamp(), sources=results))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['probe_prepare', 'probe_check', 'batch_prepare', 'worker']); parser.add_argument('--worker', type=int)
    args = parser.parse_args(); worker(args.worker) if args.stage == 'worker' else globals()[args.stage]()
