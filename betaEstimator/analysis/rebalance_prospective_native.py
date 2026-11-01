#!/usr/bin/env python3
"""Reuse completed workers' host pools for idle prospective measurement jobs."""
import argparse
import json
import re
import subprocess
from pathlib import Path
from prepare_prospective_flight import OUT
from research_flight_conditional import digest, write
from research_conditional_mc import stamp

CLUSTER = 398926
ATTRS = 'ClusterId,ProcId,JobStatus,Requirements,RemoteHost,RequestCpus,RequestMemory'


def queue():
    raw = subprocess.check_output(['condor_q', str(CLUSTER), '-json', '-attributes', ATTRS]).decode()
    rows = json.loads(raw) if raw.strip() else []
    if len(rows) > 50 or any(int(r['ClusterId']) != CLUSTER for r in rows):
        raise ValueError('Unexpected prospective measurement queue')
    return rows


def main(apply=False):
    destination = OUT / 'batch/completed_pool_rebalance.json'
    if destination.exists():
        raise FileExistsError('The completed-pool adjustment is already recorded')
    protocol = json.loads((OUT / 'batch/protocol.json').read_text())
    before = queue()
    active = {int(r['ProcId']) for r in before}
    reserved = set()
    for row in before:
        reserved.update(re.findall(r'amswn\d+\.ams\.sdiat\.cn', row['Requirements']))
    completed = []
    available_pool = set()
    for worker, hosts in enumerate(protocol['physical_host_pools']):
        completion = OUT / 'batch' / ('worker_%02d.json' % worker)
        if worker not in active and completion.exists():
            done = json.loads(completion.read_text())
            if done.get('complete') is True:
                completed.append(worker)
                available_pool.update(hosts)
    available_pool.difference_update(reserved)
    text = subprocess.check_output(['condor_status', '-constraint',
        '(OpSysAndVer=="CentOS7") && (Cpus>=32) && (Memory>=8192) && (State=="Unclaimed") && (Activity=="Idle")',
        '-af', 'Machine']).decode()
    eligible = sorted(available_pool.intersection(text.split()))
    idle = sorted((r for r in before if r['JobStatus'] == 1), key=lambda r: r['ProcId'])
    count = min(len(eligible), len(idle))
    changes = []
    if count:
        groups = [eligible[i::count] for i in range(count)]
        for row, hosts in zip(idle[:count], groups):
            worker = int(row['ProcId'])
            old_hosts = ';'.join(protocol['physical_host_pools'][worker])
            old = row['Requirements']
            if old.startswith('/Expr(') and old.endswith(')/'):
                old = old[6:-2]
            if old.count(old_hosts) != 1:
                raise ValueError('Unexpected requirements for worker %d' % worker)
            new = old.replace(old_hosts, ';'.join(hosts), 1)
            changes.append(dict(worker=worker, hosts=hosts, old_requirements=old, new_requirements=new))
    flat = [host for row in changes for host in row['hosts']]
    if len(set(flat)) != len(flat) or set(flat).intersection(reserved):
        raise ValueError('New host candidates are not disjoint from active pools')
    print('COMPLETED_POOL_PLAN', json.dumps(dict(active=len(before), idle=len(idle), completed=completed,
        idle_eligible_hosts=eligible, adjusted_workers=[r['worker'] for r in changes])), flush=True)
    if not apply:
        return
    if not changes:
        raise ValueError('No idle measurement jobs can use a completed host pool')
    # The scheduler-side idle condition protects a task that starts during this read.
    for change in changes:
        constraint = 'ClusterId == %d && ProcId == %d && JobStatus == 1' % (CLUSTER, change['worker'])
        command = ['condor_qedit', '-constraint', constraint, 'Requirements', change['new_requirements']]
        change['scheduler_result'] = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
    after = queue()
    write(destination, dict(complete=True, created_utc=stamp(), cluster=CLUSTER,
        reason='Idle jobs were restricted to occupied hosts; only host pools of verified completed workers are reused.',
        changes=changes, before=before, after=after, completed_worker_ids=completed,
        original_protocol_sha256=digest(OUT / 'batch/protocol.json'), source_sha256=digest(Path(__file__)),
        protected='Running tasks, event lists, all source code and binaries used for measurements, output files, model parameters and scientific selections',
        job_count_increased=False, scientific_configuration_changed=False))
    subprocess.check_call(['condor_reschedule'])
    print('COMPLETED_POOL_ADJUSTED', len(changes), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply', action='store_true')
    main(parser.parse_args().apply)
