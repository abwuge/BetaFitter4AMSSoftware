#!/usr/bin/env python3
"""Record successive idle-only host reassignments for the fixed native batch."""
import json
import re
import subprocess
from pathlib import Path
from prepare_prospective_flight import OUT
from research_flight_conditional import write, digest
from research_conditional_mc import stamp

CLUSTER = 398926


def main():
    root = OUT / 'batch/scheduler_rounds'
    root.mkdir(exist_ok=True)
    rows = json.loads(subprocess.check_output(['condor_q', str(CLUSTER), '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,Requirements,RemoteHost']).decode() or '[]')
    protocol = json.loads((OUT / 'batch/protocol.json').read_text())
    active = {int(r['ProcId']) for r in rows}
    idle = sorted((r for r in rows if r['JobStatus'] == 1), key=lambda r: r['ProcId'])
    reserved = set()
    for r in rows:
        reserved.update(re.findall(r'amswn\d+\.ams\.sdiat\.cn', r['Requirements']))
    completed = []
    reusable = set()
    for worker, hosts in enumerate(protocol['physical_host_pools']):
        path = OUT / 'batch' / ('worker_%02d.json' % worker)
        if worker not in active and path.exists() and json.loads(path.read_text()).get('complete'):
            completed.append(worker)
            reusable.update(hosts)
    text = subprocess.check_output(['condor_status', '-constraint',
        '(OpSysAndVer=="CentOS7") && (Cpus>=32) && (Memory>=8192) && (State=="Unclaimed") && (Activity=="Idle")',
        '-af', 'Machine']).decode()
    hosts = sorted(reusable.intersection(text.split()).difference(reserved))
    count = min(len(idle), len(hosts))
    print('NATIVE_SCHEDULER', 'queued', len(rows), 'idle', len(idle), 'free_reusable_hosts', len(hosts), flush=True)
    if not count:
        return
    number = len(list(root.glob('*.json')))
    destination = root / ('%03d.json' % number)
    if destination.exists():
        raise FileExistsError('Scheduler round record already exists')
    result = dict(created_utc=stamp(), cluster=CLUSTER, before=rows, completed_workers=completed, changes=[],
        source_sha256=digest(Path(__file__)), original_protocol_sha256=digest(OUT / 'batch/protocol.json'),
        scientific_configuration_changed=False, running_jobs_modified=False, jobs_added=False, complete=False)
    write(destination, result)
    for row, group in zip(idle[:count], [hosts[i::count] for i in range(count)]):
        old = row['Requirements']
        if old.startswith('/Expr(') and old.endswith(')/'):
            old = old[6:-2]
        matches = list(re.finditer(r'"(amswn\d+\.ams\.sdiat\.cn(?:;amswn\d+\.ams\.sdiat\.cn)*)"', old))
        if len(matches) != 1:
            raise ValueError('Unexpected native host expression')
        match = matches[0]
        new = old[:match.start(1)] + ';'.join(group) + old[match.end(1):]
        constraint = 'ClusterId == %d && ProcId == %d && JobStatus == 1' % (CLUSTER, row['ProcId'])
        command = ['condor_qedit', '-constraint', constraint, 'Requirements', new]
        response = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
        result['changes'].append(dict(worker=row['ProcId'], old=old, new=new, result=response))
        write(destination, result)
    result['complete'] = True
    write(destination, result)
    subprocess.check_call(['condor_reschedule'])
    print('NATIVE_SCHEDULER_ADJUSTED', count, 'record', destination, flush=True)


if __name__ == '__main__':
    main()
