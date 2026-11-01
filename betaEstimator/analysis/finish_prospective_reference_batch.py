#!/usr/bin/env python3
"""Retire only unstarted shared-queue workers after all RICH sources are audited."""
import json
import subprocess
from reference_prospective_flight_schema import REF, check
from research_flight_conditional import digest, write
from research_conditional_mc import stamp


def main():
    check()
    batch = REF / 'batch'
    destination = batch / 'retirement.json'
    if destination.exists():
        raise FileExistsError('Unused reference workers were already retired')
    audit = json.loads((REF / 'audit.json').read_text())
    prep = json.loads((REF / 'preparation.json').read_text())
    if not audit['complete'] or digest(REF / 'reference.npz') != audit['reference_sha256']:
        raise ValueError('Every prospective native reference must be audited before retiring idle workers')
    if set(audit['parts']) != {str(row['source']) for row in prep['jobs']}:
        raise ValueError('Reference audit does not cover the shared source manifest')
    submission = json.loads((batch / 'submission.json').read_text())
    cluster = submission['cluster']
    before = json.loads(subprocess.check_output(['condor_q', str(cluster), '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')
    idle = [r['ProcId'] for r in before if r['JobStatus'] == 1]
    result = dict(created_utc=stamp(), cluster=cluster, before=before, idle_worker_ids=idle,
        reference_audit_sha256=digest(REF / 'audit.json'), source_manifest_complete=True,
        reason='The shared queue is fully completed and audited; unstarted workers have no remaining input to process.',
        running_jobs_preserved=True, completed_files_preserved=True, complete=False)
    write(destination, result)
    if idle:
        constraint = 'ClusterId == %d && JobStatus == 1' % cluster
        result['scheduler_response'] = subprocess.check_output(['condor_rm', '-constraint', constraint], stderr=subprocess.STDOUT).decode()
    else:
        result['scheduler_response'] = 'No idle jobs remain.'
    result['after'] = json.loads(subprocess.check_output(['condor_q', str(cluster), '-json', '-attributes',
        'ClusterId,ProcId,JobStatus,RemoteHost']).decode() or '[]')
    result['complete'] = True
    write(destination, result)
    print('PROSPECTIVE_EMPTY_WORKERS_RETIRED', len(idle), 'remaining_jobs', len(result['after']), flush=True)


if __name__ == '__main__':
    main()
