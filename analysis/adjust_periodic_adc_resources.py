#!/usr/bin/env python3
"""Correct the scheduler reservation without changing any scientific fit."""
import json
import subprocess
from pathlib import Path
from research_flight_conditional import write, digest
from research_conditional_mc import stamp
from dispatch_periodic_adc import ROOT, BATCH


def main():
    path = BATCH / 'resource_amendment.json'
    if path.exists():
        raise FileExistsError('Resource amendment already recorded')
    cluster = json.loads((BATCH / 'submission.json').read_text())['cluster']
    query = ['condor_q', str(cluster), '-json', '-attributes', 'ClusterId,ProcId,JobStatus,RequestCpus,RemoteHost']
    before = json.loads(subprocess.check_output(query).decode() or '[]')
    if len(before) != 29 or any(row['JobStatus'] != 1 or row.get('RemoteHost') or row['RequestCpus'] != 4 for row in before):
        raise ValueError('Expected the twenty-nine unstarted four-CPU jobs')
    completed = json.loads((BATCH / 'task_00.json').read_text())
    if completed['exit_code'] or not completed['numerically_converged']:
        raise ValueError('Preserve and inspect the first completed fit')
    evidence = subprocess.check_output(['condor_q', str(cluster)+'.1', '-better-analyze', '-reverse', '-machine', 'amswn001.ams.sdiat.cn']).decode()
    if 'Target.RequestCpus == 32' not in evidence:
        raise ValueError('The observed cluster CPU reservation rule changed')
    source = ROOT / 'analysis/dispatch_periodic_adc.py'
    write(path, dict(created_utc=stamp(), stage='before_scheduler_change', before=before,
        reason='Reverse matchmaking reports START = Target.RequestCpus == 32 for the ordinary worker hosts; scientific fits still use four threads.',
        original_dispatch_protocol_sha256=digest(BATCH / 'protocol.json'), original_submit_sha256=digest(BATCH / 'fits.sub'),
        dispatch_source_sha256=digest(source), completed_task_sha256=digest(BATCH / 'task_00.json'),
        scientific_thread_count_unchanged=4, request_cpus_old=4, request_cpus_new=32, matchmaking_evidence=evidence,
        no_resubmission=True, no_completed_fit_replaced=True))
    reply = subprocess.check_output(['condor_qedit', str(cluster), 'RequestCpus', '32']).decode()
    after = json.loads(subprocess.check_output(query).decode() or '[]')
    if any(row['RequestCpus'] != 32 for row in after):
        raise ValueError('Scheduler reservation update was incomplete')
    write(BATCH / 'resource_amendment_applied.json', dict(complete=True, created_utc=stamp(), scheduler_reply=reply,
        after=after, amendment_sha256=digest(path), source_sha256=digest(Path(__file__))))
    print('PERIOD_RESOURCES_CORRECTED', len(after), 'pending jobs now reserve 32 CPUs; four fit threads unchanged', flush=True)


if __name__ == '__main__':
    main()
