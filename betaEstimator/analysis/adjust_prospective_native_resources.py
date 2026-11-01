#!/usr/bin/env python3
"""Record and repair idle native-job resource requests without restarting work."""
import json
import subprocess
from pathlib import Path
from prepare_prospective_flight import OUT
from research_flight_conditional import write, digest
from research_conditional_mc import stamp


def queue():
    text = subprocess.check_output(['condor_q', '398926', '-json']).decode()
    return json.loads(text) if text.strip() else []


if __name__ == '__main__':
    destination = OUT / 'batch/resource_amendment.json'
    if destination.exists():
        raise FileExistsError('Resource correction already recorded')
    before = queue()
    if len(before) > 50 or any(int(row['ClusterId']) != 398926 for row in before):
        raise ValueError('Unexpected native production jobs')
    selected = [row for row in before if row['JobStatus'] == 1]
    command = ['condor_qedit', '-constraint', 'ClusterId == 398926 && JobStatus == 1', 'RequestCpus', '32']
    result = subprocess.check_output(command, stderr=subprocess.STDOUT).decode()
    after = queue()
    compact = lambda rows: [{k: r.get(k) for k in ['ClusterId', 'ProcId', 'JobStatus', 'RequestCpus', 'RequestMemory', 'RemoteHost']} for r in rows]
    write(destination, dict(created_utc=stamp(), cluster=398926, original_protocol_sha256=digest(OUT / 'batch/protocol.json'),
        reason='condor_q -better-analyze -reverse on slot7@amswn021 reported START = Target.RequestCpus == 32, rejecting the initial RequestCpus=4',
        command=command, result=result, before=compact(before), after=compact(after), eligible_idle_jobs=len(selected),
        unchanged='Running jobs, worker four-process concurrency, source assignments, disjoint physical host pools, binary, selection, model parameters, and existing output files',
        model_or_input_change=False, source_sha256=digest(Path(__file__))))
    subprocess.check_call(['condor_reschedule'])
    print('NATIVE_RESOURCE_CORRECTION', len(selected), result, flush=True)
